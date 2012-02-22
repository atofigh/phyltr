/*
 * Copyright (C) 2010, 2011 Ali Tofigh
 *
 * This file is part of PhylTr, a package for phylogenetic analysis
 * using duplications and transfers.
 *
 * PhylTr is released under the terms of the license contained in the
 * file LICENSE.
 */

/*
 * This file contains source code for a program that generates a
 * binary gene trees from a birth-death process given a species tree
 * with times on the edges. Given the intensities for duplications,
 * losses and horizontal gene transfers, the program creates a gene
 * tree that has evolved inside the species tree.
 */

#include "common.hh"

#include <NHparser/NHparser.h>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <cstdlib>
#include <errno.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>


//============================================================================
//             Namespace declarations and using directives.
//============================================================================

using namespace std;

//============================================================================
//                   Typedefs and class declarations
//============================================================================

typedef Binary_tree<double>         Tree_type;
typedef Tree_type::vid_t            vid_t;

struct Node;
typedef boost::shared_ptr<Node>     Node_ptr;
typedef boost::weak_ptr<Node>       Node_ptr_weak;

/*
 * Class Node
 *
 * Tree  structure  used  when  creating  the  gene  tree  during  the
 * birth-death process.
 */
struct Node {
    enum Event {none, duplication, transfer, speciation};
    static unsigned next_id;

    Node_ptr_weak   parent;
    Node_ptr        left, right;
    vid_t           stree_label;
    Event           event;
    double          edge_time;
    bool            edge_is_transfer;
    unsigned        id;

    Node(Node_ptr p, Node_ptr l, Node_ptr r,
         vid_t label, Event e, double t, bool transfer)
        : parent(p), left(l), right(r), stree_label(label), event(e),
          edge_time(t), edge_is_transfer(transfer), id(next_id++)
    {}
};

unsigned Node::next_id = 0;

/*
 * class Stats
 *
 * Simple data struct to keep track of some statistics during the gene
 * tree evolution process.
 */
struct Stats {
    int duplications;
    int transfers;
    int losses;
    int speciations;
    double total_time;

    Stats()
        : transfers(0), duplications(0), losses(0), speciations(0), total_time(0.0)
    {}
};

//============================================================================
//                   Global Constants and variables.
//============================================================================

const string PROG_NAME = "phyltr-gen-gtree";
const string USAGE = "Usage: " + PROG_NAME +
    " [OPTION]... TREE-FILE DUPLICATION-RATE TRANSFER-RATE DEATH-RATE NAME-PREFIX";

/* The precision with which edge times will be printed. */
const int TIME_PRECISION = 4;

/* helper constant for null shared pointer to Nodes. */
Node_ptr null;

/* g_options contains the command line options. */
boost::program_options::variables_map g_options;

/* The variable holding the maximum number of attempts of creating a
   gene tree with at least one leaf. */
unsigned max_attempts;

/* The rates of duplication, transfer, and death respectively. */
double delta, tau, mu;

/* The gene tree root node. */
Node_ptr gtree_root = null;

/* The file prefix of the output files. */
string file_prefix;

/* The "only-visible-events" flag. */
bool only_visible_events = false;

/* The "only-sane-scenarios" flag. */
bool only_sane_scenarios = false;

/* The "all-species" flag. */
bool all_species= false;

/* The seed used to initiate the random number generator. */
unsigned seed = 0;

//=============================================================================
//                        Function declarations
//=============================================================================

/*
 * advance_process()
 *
 * Advances the gene evolution process during 'time' units of time
 * assuming that the process takes place during a slice of the species
 * tree, i.e., the only edges present during this time is the edges in
 * 'slice'.
 */
void advance_process(double                 time,
                     const vector<vid_t>   &slice, 
                     vector<Node_ptr>      &cur_genes,
                     Stats                 &stats);

/*
 * speciate()
 *
 * Creates children for any gene in 'gtree_leaves' that is labeled by
 * 'slice[slice_idx]'.
 */
void speciate(const Tree_type      &stree,
              vector<vid_t>        &slice,
              unsigned              slice_idx,
              vector<Node_ptr>     &cur_genes,
              Stats &stats);

/*
 * acceptable_gene_tree()
 *
 * Returns true if the gene tree produced by the process is acceptable
 * according to the options, i.e., if the gene tree is not empty,
 * it is not the case that the --sane-scenarios-only flag has been
 * given and the scenario is not sane, and perhaps other conditions.
 */
bool acceptable_gene_tree(unsigned n_leaves,
                          const Tree_type &stree);

/*
 * print_gtree(), print_sigma(), print_gamma(), print_events()
 *
 * Functions used for outputting the result of the birth-death
 * process.
 */
void print_arguments(ostream &out, int argc, char *argv[]);
void print_gtree(ostream &out, const Tree_type &stree);
void print_sigma(ostream &out, const Tree_type &stree);
void print_gamma(ostream &out, const Tree_type &stree);
void print_events(ostream &out, const Tree_type &stree);
void print_stats(ostream &out, const Stats &stats);
//=============================================================================
//         Template and inline function and member definitions.
//=============================================================================



//=============================================================================
//                                main()
//=============================================================================

int 
main(int argc, char *argv[])
{
    namespace po = boost::program_options;

    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    string stree_filename;      // Contains the filename of the species tree.
    try
        {
            /* Declare visible options */
            visible_opts.add_options()
                ("help,h", "display this help and exit")
                ("min-size",
                 po::value<unsigned>(),
                 "The minimum number of leaves required")
                ("max-size",
                 po::value<unsigned>(),
                 "The maximum number of leaves allowed")
                ("root-length,r", po::value<double>(),
                 "If not set, the birth-death process starts according to the length of "
                 "the root of the species tree. E.g., if you want to force the process "
                 "to start at the first speciation in the species tree, then set this "
                 "parameter to zero.")
                ("sane-scenarios-only,s",
                 po::bool_switch(&only_sane_scenarios)->default_value(false),
                 "Gene trees created with strange event combinations are "
                 "discarded. A scenario is sane iff (1) the children of each "
                 "duplication are placed at or below the duplication in the "
                 "species tree and at comparable species tree vertices, "
                 "(2) the head of each "
                 "transfer edge is placed at an incomparable species tree "
                 "vertex compared to the tail, (3) each gene tree "
                 "vertex has at least one of its children placed at or "
                 "below it in the species tree, and (4) the children of "
                 "each speciation are placed at incomparable species tree "
                 "vertices and such that their lca is the speciation."
                 )
                ("all-species,a",
                 po::bool_switch(&all_species)->default_value(false),
                 "If set, gene trees that lack genes in one or more "
                 "species in the species tree are discarded."
                 )
                ("max-attempts,m",
                 po::value<unsigned>(&max_attempts)->default_value(10000),
                 "The maximum number of times the process will restart if the "
                 "gene tree goes extinct or if the scenario is not sane "
                 "and the --sane-scenarios-only flag has been given.")
                ("visible-events-only,v",
                 po::bool_switch(&only_visible_events)->default_value(false),
                 "When this switch is given, a duplication is classified "
                 "as a duplication only if both its children are placed at "
                 "comparable vertices in the species tree and are also placed "
                 "at or below the duplication in the species tree, and a "
                 "transfer edge is classified as a transfer edge only if "
                 "the head is placed at an incomparable vertex compared to "
                 "the tail. Note that this does not always transform "
                 "a scenario to a sane one since other silly things can "
                 "happen, e.g., a speciation followed by a transfer that "
                 "places both children in the same edge."
                 )
                ("seed",
                 po::value<unsigned>(),
                 "A seed to initialize the random number generator with."
                 )
                ;
            
            /* Declare positional options. */
            hidden_opts.add_options()
                ("filename", po::value<string>(&stree_filename))
                ("delta", po::value<double>(&delta))
                ("tau", po::value<double>(&tau))
                ("mu", po::value<double>(&mu))
                ("file-prefix", po::value<string>(&file_prefix));
                ;
            
            po::positional_options_description positional_options;
            positional_options.add("filename", 1);
            positional_options.add("delta", 1);
            positional_options.add("tau", 1);
            positional_options.add("mu", 1);
            positional_options.add("file-prefix", 1);

            /* Gather all options into a single options_description. */
            all_options.add(visible_opts).add(hidden_opts);
            
            /* Parse the arguments. */
            po::command_line_parser parser(argc, argv);
            parser.options(all_options);
            parser.positional(positional_options);
            po::store(parser.run(), g_options);
            po::notify(g_options);
        }
    catch (const exception &e)
        {
            cerr << PROG_NAME << ": " << e.what() << '\n'
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* Show help message if --help was given. */
    if (g_options.count("help"))
        {
            cout << USAGE << '\n'
                 << visible_opts << '\n';
            exit(EXIT_SUCCESS);
        }
    
    if (g_options.count("seed"))
        {
            seed = init_rand(g_generator, g_options["seed"].as<unsigned>());
        }
    else
        {
            seed = init_rand(g_generator);
        }



    /* Check that all required positional arguments are given. */
    if (g_options.count("filename") == 0
        || g_options.count("delta") == 0
        || g_options.count("tau") == 0
        || g_options.count("mu") == 0
        || g_options.count("file-prefix") == 0)
        {
            cerr << PROG_NAME << ": Too few arguments.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* Check that all the rates are non-negative. */
    if (delta < 0 || tau < 0 || mu < 0)
        {
            cerr << PROG_NAME << ": negative rates\n";
            exit(EXIT_FAILURE);
        }

    /* Check that root-length, if given, is non-negative. */
    if (g_options.count("root-length") > 0 && 
        g_options["root-length"].as<double>() < 0)
        {
            cerr << PROG_NAME << ": negative root-length\n";
            exit(EXIT_FAILURE);
        }

    /* Check that max-attempts is not zero. */
    if (g_options["max-attempts"].as<unsigned>() == 0)
        {
            cerr << PROG_NAME << ": "
                 << "max-attempts was set to zero.\n";
            exit(EXIT_FAILURE);
        }
    
    /* Check that min-size <= max-size if given and that max-size > 0. */
    if (g_options.count("min-size") && g_options.count("max-size"))
        {
            if (g_options["min-size"].as<unsigned>() >
                g_options["max-size"].as<unsigned>())
                {
                    cerr << PROG_NAME << ": --min-size must be less than "
                         << "or equal to --max-size\n";
                    exit(EXIT_FAILURE);
                }
        }

    if (g_options.count("max-size") &&
        g_options["max-size"].as<unsigned>() == 0)
        {
            cerr << PROG_NAME << ": --max-size must be a positive integer.\n";
            exit(EXIT_FAILURE);
        }

    /*
     * Read the species tree.
     */

    /* Try to open the file */
    errno = 0;
    FILE *stree_file = fopen(stree_filename.c_str(), "r");
    if (errno)
        {
            cerr << PROG_NAME << ": "
                 << "failed to open species tree file '"
                 << stree_filename << "'\n";
            perror(PROG_NAME.c_str());
            exit(EXIT_FAILURE);
        }
    
    /* Read the tree using the NHParser library. */
    NHNode *NH_stree_root = NH_read_tree(stree_file, stderr,
                                         stree_filename.c_str());
    if (NH_stree_root == 0)
        {
            cerr << PROG_NAME << ": "
                 << "failed to read species tree\n";
            exit(EXIT_FAILURE);
        }

    /* Check that the tree is full binary. */
    if (!NHtree_is_binary(NH_stree_root))
        {
            cerr << PROG_NAME << ": "
                 << "species tree is not binary\n";
            exit(EXIT_FAILURE);
        }

    /* transform the NH-tree to Binary_tree. */
    Tree_type stree(NH_stree_root->branchLength < 0 ? 0 : NH_stree_root->branchLength,
                    NH_stree_root->label == 0 ? "" : NH_stree_root->label);
    try
        {
            create_binary_tree(NH_stree_root, stree, 0);
        }
    catch (const exception &e)
        {
            cerr << PROG_NAME << ": " << e.what() << '\n';
            exit(EXIT_FAILURE);
        }

    /* Check that all edges have a positive length. */
    for (vid_t v = stree.preorder_begin(); v != Tree_type::NONE;
         v = stree.preorder_next(v))
        {
            if (v == 0)
                {
                    if (stree.length(v) < 0)
                        {
                            cerr << PROG_NAME << ": "
                                 << "bad length at root of species tree"
                                 << '\n';
                            exit(EXIT_FAILURE);
                        }
                }
            else if (stree.length(v) <= 0)
                {
                    cerr << PROG_NAME << ": "
                         << "bad length in species tree"
                         << '\n';
                    exit(EXIT_FAILURE);
                }
        }
    
    /* Check that all leaves have unique labels. */
    vector<string> stree_labels;
    for (vid_t v = 0; v < stree.size(); ++v)
        {
            if (stree.is_leaf(v))
                stree_labels.push_back(stree.label(v));
        }
    sort(stree_labels.begin(), stree_labels.end());
    if (unique(stree_labels.begin(), stree_labels.end()) != stree_labels.end())
        {
            cerr << PROG_NAME << ": "
                 << "species tree leaves do not have unique labels.\n";
            exit(EXIT_FAILURE);
        }

    /*
     * Now we have the species tree in stree. And the parameters are
     * in delta, tau, and mu. Go create the gene tree!
     */

    unsigned attempts = max_attempts;
    bool gene_tree_accepted = false;
    Stats stats;

    while (attempts--)
        {
            stats = Stats();

            /* Create a vector with slice times in decreasing order so that we
               can pop the times from the back efficiently). */
            vector<double> times;
            for (vid_t v = 1; v < stree.size(); ++v)
                {
                    times.push_back(stree.time(v));
                }
            sort(times.begin(), times.end(), greater<double>());

            /* Create a vector that will keep track of the leaves of the gene tree. */
            vector<Node_ptr> cur_genes;
            unsigned n_leaves = 0;

            /* Create a vector keeping track of the current edges of the slice. */
            vector<vid_t> slice;
    
            /* Start with a single gene at or before the species tree root. */
            gtree_root.reset(new Node(null, null, null, 0,
                                      Node::none, 0.0, false));
            cur_genes.push_back(gtree_root);
            slice.push_back(0);

            /* Determine the time before the first speciation. */
            double root_time = stree.time(0);
            if (g_options.count("root-length"))
                root_time = g_options["root-length"].as<double>();

            /* Run the process before the first speciation. */
            advance_process(root_time, slice, cur_genes, stats);
    
            /* The first speciation occurs. Replace the stree root edge with
               its children, and duplicate the genes. */
            if (!stree.is_leaf(0))
                {
                    speciate(stree, slice, 0, cur_genes, stats);
                }
            /* Keep advancing and taking care of speciations until the end. */
            /* Note that leaves of stree might not end at the same time. */
            double cur_time = stree.time(0);
            while (times.empty() == false)
                {
                    /* Advance the process until next speciation/leaf. */
                    advance_process(times.back() - cur_time, slice, cur_genes, stats);
                    cur_time = times.back();
                    times.pop_back();

                    /* find an stree edge that ends at the current time */
                    unsigned idx = 0;
                    for (idx = 0; idx < slice.size(); ++idx)
                        {
                            if (stree.time(slice[idx]) == cur_time)
                                break;
                        }
                    vid_t cur_snode = slice[idx];
                    /* cur_snode is either speciation or a leaf. */
                    if (stree.is_leaf(cur_snode))
                        {
                            using namespace boost::lambda;
                            vector<Node_ptr>::iterator new_end = 
                                remove_if(cur_genes.begin(),
                                          cur_genes.end(),
                                          &*_1 ->* &Node::stree_label == cur_snode);
                            n_leaves += distance(new_end, cur_genes.end());
                            cur_genes.erase(new_end, cur_genes.end());

                            swap(slice[idx], slice.back());
                            slice.pop_back();
                        }
                    else
                        {
                            speciate(stree, slice, idx, cur_genes, stats);
                        }
                }
            if (acceptable_gene_tree(n_leaves, stree))
                {
                    gene_tree_accepted = true;
                    break;
                }
        }

    if (!gene_tree_accepted)
        {
            cerr << PROG_NAME << ": "
                 << "No acceptable gene tree was produced in "
                 << max_attempts << " runs of the process\n";
            exit(EXIT_FAILURE);
        }



    /*
     * Output the result to the files prefix.args, prefix.gtree,
     * prefix.sigma, prefix.gamma, and prefix.events.
     */
    
    vector<string> filenames;
    filenames.push_back(file_prefix + ".args");
    filenames.push_back(file_prefix + ".gtree");
    filenames.push_back(file_prefix + ".sigma");
    filenames.push_back(file_prefix + ".gamma");
    filenames.push_back(file_prefix + ".events");
    filenames.push_back(file_prefix + ".stats");

    typedef boost::shared_ptr<ofstream> ofstream_ptr;
    vector<ofstream_ptr> outfiles(filenames.size());
    for (unsigned i = 0; i < outfiles.size(); ++i)
        {
            outfiles[i].reset(new ofstream(filenames[i].c_str()));
            if (outfiles[i]->fail())
                {
                    cerr << PROG_NAME << ": "
                         << "failed to open '"
                         << filenames[i] << "' for writing\n";
                    exit(EXIT_FAILURE);
                }
        }

    print_arguments(*outfiles[0], argc, argv);
    print_gtree(*outfiles[1], stree);
    print_sigma(*outfiles[2], stree);
    print_gamma(*outfiles[3], stree);
    print_events(*outfiles[4], stree);
    print_stats(*outfiles[5], stats);

    return EXIT_SUCCESS;
}

//=============================================================================
//                     Helper function declarations
//=============================================================================

void postorder_gtree(vector<Node_ptr> &gtree_postorder);
void postorder_gtree_r(Node_ptr n, vector<Node_ptr> &gtree_postorder);
void print_gtree_r(ostream &out, const Tree_type &stree,
                   Node_ptr n, vector<int> &numbers);
bool is_sane(Node_ptr n, const Tree_type &stree);
bool genes_in_all_species(const Tree_type &stree);

//=============================================================================
//                   Function and member definitions.
//=============================================================================

bool
acceptable_gene_tree(unsigned n_leaves,
                     const Tree_type &stree)
{
    if (!gtree_root)
        return false;

    if (only_sane_scenarios && !is_sane(gtree_root, stree))
        return false;

    if (all_species && !genes_in_all_species(stree))
        return false;

    if (g_options.count("min-size") &&
        n_leaves < g_options["min-size"].as<unsigned>())
        {
            return false;
        }

    if (g_options.count("max-size") &&
        n_leaves > g_options["max-size"].as<unsigned>())
        {
            return false;
        }

    return true;
}


bool
is_sane(Node_ptr n, const Tree_type &stree)
{
    if (!n)
        return false;

    vid_t x = n->stree_label;
    Node_ptr parent = n->parent.lock();

    /* The head of a transfer edge must be placed at an incomparable
       stree vertex compared to the tail. */
    if (n->edge_is_transfer && parent)
        {
            vid_t p = parent->stree_label;
            if (stree.descendant(x, p))
                return false;
        }

    if (!n->left)
        return true;

    vid_t y = n->left->stree_label;
    vid_t z = n->right->stree_label;
    vid_t lcayz = stree.lca(y, z);
    
    /* at least one of the children must be placed below x */
    if (!stree.descendant(y, x) && ! stree.descendant(z, x))
        return false;
    
    /* if n is a speciation, then the children must be strictly below
       x and placed at incomparable stree vertices. */
    if (n->event == Node::speciation)
        {
            if (lcayz != x)
                return false;
        }

    /* If n is a duplications, both children must be placed below x
       and at comparable stree vertices. */
    if (n->event == Node::duplication)
        {
            if (lcayz != y && lcayz != z)
                return false;
            if (!stree.descendant(y, x) || !stree.descendant(z, x))
                return false;
        }

    return is_sane(n->left, stree) && is_sane(n->right, stree);
}


bool
genes_in_all_species(const Tree_type &stree)
{
    // Record which stree leaves occur as stree_labels in the gene
    // tree and count them.
    int stree_leaves = (stree.size() + 1) / 2;
    vector<bool> visited_leaf(stree.size(), false);

    vector<Node_ptr> gtree_postorder;
    postorder_gtree(gtree_postorder);

    BOOST_FOREACH(Node_ptr np, gtree_postorder)
        {
            if (!np->left)
                {
                    visited_leaf[np->stree_label] = true;
                }
        }
    int visited_leaves = count(visited_leaf.begin(), visited_leaf.end(), true);

    return stree_leaves == visited_leaves;
}



void
advance_process(double time,
                const vector<vid_t> &slice, 
                vector<Node_ptr> &cur_genes,
                Stats &stats)
{
    while (!cur_genes.empty())
        {
            /* When is the next event? */
            double total_rate = cur_genes.size() * (delta + tau + mu);
            double time_used;
            if (total_rate > 0)
                {
                    boost::exponential_distribution<double> rng_exp(total_rate);
                    time_used = min(rng_exp(g_rng_d), time);
                }
            else
                {
                    time_used = time;
                }
            /* update the times of the gene tree edges. */
            vector<Node_ptr>::iterator e = cur_genes.end();
            for (vector<Node_ptr>::iterator i = cur_genes.begin();
                 i != e; ++i)
                {
                    (*i)->edge_time += time_used;
                }

            /* update the stats with TOTAL time */
            stats.total_time += time_used * cur_genes.size();

            /* if nothing happened during remaining time, then end. */
            time -= time_used;
            if (time <= 0)
                break;

            /* Choose a gene to which something will happen*/
            unsigned chosen_idx = g_rng_ui(cur_genes.size());
            Node_ptr chosen = cur_genes[chosen_idx];
            
            /* Choose an event based on rates. */
            double event = g_rng_d();
            if (event < delta / (delta + tau + mu))             // duplication
                {
                    stats.duplications++;

                    Node_ptr g1(new Node(chosen, null, null,
                                         chosen->stree_label,
                                         Node::none, 0.0, false));
                    Node_ptr g2(new Node(chosen, null, null,
                                         chosen->stree_label,
                                         Node::none, 0.0, false));
                    chosen->left = g1;
                    chosen->right = g2;
                    chosen->event = Node::duplication;
                    cur_genes[chosen_idx] = g1;
                    cur_genes.push_back(g2);
                }
            else if (event < (delta + tau) / (delta + tau + mu))// transfer
                {
                    /* Do nothing if slice has only one trunk. */
                    
                    if (slice.size() > 1)
                        {
                            stats.transfers++;

                            /* Choose which trunk to transfer to. */
                            unsigned transfer_trunk_idx
                                = g_rng_ui(slice.size() - 1);
                            if (slice[transfer_trunk_idx] == chosen->stree_label)
                                transfer_trunk_idx = slice.size() - 1;
                            
                            /* Create new children for the chosen gene. */
                            Node_ptr g1(new Node(chosen, null, null,
                                                 chosen->stree_label,
                                                 Node::none, 0.0, false));
                            Node_ptr g2(new Node(chosen, null, null,
                                                 slice[transfer_trunk_idx],
                                                 Node::none, 0.0, true));
                            chosen->left = g1;
                            chosen->right = g2;
                            chosen->event = Node::transfer;
                            cur_genes[chosen_idx] = g1;
                            cur_genes.push_back(g2);
                        }
                }
            else                                                // loss
                {
                    stats.losses++;
                    /* If only one gene remaining, process stops. */
                    if (cur_genes.size() == 1)
                        {
                            cur_genes.pop_back();
                            gtree_root.reset();
                            break;
                        }

                    /* Rewire pointers to remove chosen gene. */
                    Node_ptr parent = chosen->parent.lock();
                    Node_ptr grand_parent = parent->parent.lock();
                    Node_ptr sibling =
                        parent->left == chosen ? parent->right : parent->left;
                    sibling->parent = grand_parent;
                    sibling->edge_is_transfer =
                        sibling->edge_is_transfer || parent->edge_is_transfer;
                    sibling->edge_time += parent->edge_time;
                    if (!grand_parent)
                        {
                            gtree_root = sibling;
                        }
                    else
                        {
                            if (grand_parent->left == parent)
                                grand_parent->left = sibling;
                            else
                                grand_parent->right = sibling;
                        }
                    swap(cur_genes[chosen_idx], cur_genes.back());
                    cur_genes.pop_back();
                }
        }
}


void
speciate(const Tree_type       &stree,
         vector<vid_t>         &slice,
         unsigned               slice_idx,
         vector<Node_ptr>      &cur_genes,
         Stats                 &stats)
{
    /* Replace speciating vertex in slice with its children. */
    vid_t cur_snode = slice[slice_idx];
    slice[slice_idx] = stree.left(cur_snode);
    slice.push_back(stree.right(cur_snode));
                    
    /* Replace each gene in cur_genes labelled with cur_snode with new
       children. */
    unsigned nleaves = cur_genes.size();
    for (unsigned i = 0; i < nleaves; ++i)
        {
            if (cur_genes[i]->stree_label != cur_snode)
                continue;

            stats.speciations++;

            Node_ptr cur_gene = cur_genes[i];
            Node_ptr g1(new Node(cur_gene, null, null,
                                 stree.left(cur_snode),
                                 Node::none, 0.0, false));
            Node_ptr g2(new Node(cur_gene, null, null,
                                 stree.right(cur_snode),
                                 Node::none, 0.0, false));
            cur_gene->left = g1;
            cur_gene->right = g2;
            cur_gene->event = Node::speciation;
            cur_genes[i] = g1;
            cur_genes.push_back(g2);
        }
}


void
print_arguments(ostream &out, int argc, char *argv[])
{
    for (int i = 1; i < argc - 1; ++i)
        {
            out << argv[i] << " ";
        }
    out << argv[argc - 1] << "\n";
    out << "seed: " << seed << "\n";
}


void
print_gtree(ostream &out, const Tree_type &stree)
{
    vector<int> numbers(stree.size(), 0); // For each leaf of stree
                                          // keeps track of the last
                                          // number used.
    print_gtree_r(out, stree, gtree_root, numbers);
    out << ";\n";
}


void
print_sigma(ostream &out, const Tree_type &stree)
{
    vector<int> numbers(stree.size(), 0); // For each leaf of stree
                                          // keeps track of the last
                                          // number used.
    vector<Node_ptr> gtree_postorder;
    postorder_gtree(gtree_postorder);

    // For each leaf, print the stree label and stree label
    // concatenated with a _ and the corresponding number.
    BOOST_FOREACH(Node_ptr np, gtree_postorder)
        {
            // We only want leaves.
            if (np->left)
                continue;

            out << stree.label(np->stree_label) << "_"
                << ++numbers[np->stree_label] << "\t"
                << stree.label(np->stree_label) << "\n";
        }
}


void
print_gamma(ostream &out, const Tree_type &stree)
{
    vector<Node_ptr> gtree_postorder;
    postorder_gtree(gtree_postorder);

    // Create a mapping of vid_t -> postorder number, i.e.,
    // stree_vid_to_postorder[v] is the postorder number of v.
    vector<int> stree_vid_to_postorder(stree.size());
    int i = 0;
    for (vid_t v = stree.postorder_begin();
         v != stree.NONE;
         v = stree.postorder_next(v))
        {
            stree_vid_to_postorder[v] = i++;
        }
    
    for (unsigned i = 0; i < gtree_postorder.size(); ++i)
        {
            out << i << "\t"
                << stree_vid_to_postorder[gtree_postorder[i]->stree_label]
                << "\n";
        }
}


void
print_events(ostream &out, const Tree_type &stree)
{
    vector<Node_ptr> gtree_postorder;
    postorder_gtree(gtree_postorder);
    
    out << "Transfer edges: ";
    for (unsigned i = 0; i < gtree_postorder.size() - 1; ++i) // skip root!
        {
            if (gtree_postorder[i]->edge_is_transfer)
                {
                    if (!only_visible_events)
                        {
                            out << i << " ";
                        }
                    else 
                        {
                            Node_ptr node = gtree_postorder[i];
                            Node_ptr parent = gtree_postorder[i]->parent.lock();
                            if (!stree.descendant(node->stree_label,
                                                  parent->stree_label))
                                {
                                    out << i << " ";
                                }
                        }
                }
        }
    out << "\n";
    
    out << "Duplications:    ";
    for (unsigned i = 0; i < gtree_postorder.size(); ++i)
        {
            if (gtree_postorder[i]->event == Node::duplication)
                {
                    if (!only_visible_events)
                        {
                            out << i << " ";
                        }
                    else 
                        {
                            vid_t x = gtree_postorder[i]->stree_label;
                            vid_t y = gtree_postorder[i]->left->stree_label;
                            vid_t z = gtree_postorder[i]->right->stree_label;
                            if (stree.descendant(y, x) &&
                                stree.descendant(z, x) &&
                                (stree.descendant(y, z) || 
                                 stree.descendant(z, y)))
                                {
                                    out << i << " ";
                                }
                        }
                }
        }
    out << "\n";
}


void
print_stats(ostream &out, const Stats &stats)
{
    out << stats.speciations << "\tspeciation event"
        << (stats.speciations == 1 ? "" : "s") << "\n";
    out << stats.duplications << "\tduplication event"
        << (stats.duplications == 1 ? "" : "s") << "\n";
    out << stats.transfers << "\ttransfer event"
        << (stats.transfers == 1 ? "" : "s") << "\n";
    out << stats.losses << "\tloss event"
        << (stats.losses == 1 ? "" : "s") << "\n";
    out << stats.total_time << setprecision(TIME_PRECISION)
        << "\ttotal unpruned gene tree time\n";
}

        
void
postorder_gtree(vector<Node_ptr> &gtree_postorder)
{
    gtree_postorder.clear();
    postorder_gtree_r(gtree_root, gtree_postorder);
}


void
postorder_gtree_r(Node_ptr n, vector<Node_ptr> &gtree_postorder)
{
    if (n->left)
        {
            postorder_gtree_r(n->left, gtree_postorder);
            postorder_gtree_r(n->right, gtree_postorder);
        }
    gtree_postorder.push_back(n);
}


void
print_gtree_r(ostream &out,
              const Tree_type &stree,
              Node_ptr n,
              vector<int> &numbers)
{
    if (!n->left)
        {
            out << stree.label(n->stree_label) << "_"
                << ++numbers[n->stree_label];
        }
    else
        {
            out << "(";
            print_gtree_r(out, stree, n->left, numbers);
            out << ", ";
            print_gtree_r(out, stree, n->right, numbers);
            out << ")";
        }
    out << ":" << setprecision(TIME_PRECISION) << showpoint << n->edge_time;
}
