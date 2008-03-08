/*
 * This file contains source code for a program that generates a
 * binary gene trees from a birth-death process given a species tree
 * with times on the edges. Given the intensities for duplications,
 * losses and horizontal gene transfers, the program creates a gene
 * tree that has evolved inside the species tree.
 */

#include "common.hh"

extern "C" {
#include <NHparser.h>
}
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <cstdlib>
#include <errno.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
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

/*
 * Class Node
 *
 * Tree  structure  used  when  creating  the  gene  tree  during  the
 * birth-death process.
 */
struct Node {
    enum Event {none, duplication, transfer, speciation};
    static unsigned next_id;

    Node_ptr    parent, left, right;
    vid_t       stree_label;
    Event       event;
    bool        edge_is_transfer;
    unsigned    id;

    Node(Node_ptr p, Node_ptr l, Node_ptr r,
         vid_t label, Event e, bool transfer)
        : parent(p), left(l), right(r), stree_label(label), event(e),
          edge_is_transfer(transfer), id(next_id++)
    {}
};

unsigned Node::next_id = 0;

//============================================================================
//                   Global Constants and variables.
//============================================================================

const string PROG_NAME = "phyltr-gen-gtree";
const string USAGE = "Usage: " + PROG_NAME +
    " [OPTION]... TREE-FILE DUPLICATION-RATE TRANSFER-RATE DEATH-RATE NAME-PREFIX";

/* helper constant for null shared pointer to Nodes. */
Node_ptr null;


/* g_options contains the command line options. */
boost::program_options::variables_map g_options;

/* The rates of duplication, transfer, and death respectively. */
double delta, tau, mu;

/* The gene tree root node. */
Node_ptr gtree_root = null;

/* The file prefix of the output files. */
string file_prefix;
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
                     vector<Node_ptr>      &cur_genes);

/*
 * speciate()
 *
 * Creates children for any gene in 'gtree_leaves' that is labeled by
 * 'slice[slice_idx]'.
 */
void speciate(const Tree_type      &stree,
              vector<vid_t>        &slice,
              unsigned              slice_idx,
              vector<Node_ptr>     &cur_genes);

/*
 * print_gtree(), print_sigma(), print_gamma(), print_events()
 *
 * Functions used for outputting the result of the birth-death
 * process.
 */

void print_gtree(ostream &out, const Tree_type &stree);
void print_sigma(ostream &out, const Tree_type &stree);
void print_gamma(ostream &out, const Tree_type &stree);
void print_events(ostream &out, const Tree_type &stree);
//=============================================================================
//         Template and inline function and member definitions.
//=============================================================================



//=============================================================================
//                                main()
//=============================================================================

int 
main(int argc, char *argv[])
{
    init_rand(g_generator);

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
                ("root-length,r", po::value<double>(),
                 "If not set, the birth-death process starts according to the length of "
                 "the root of the species tree. E.g., if you want to force the process "
                 "to start at the first speciation in the species tree, then set this "
                 "paramter to zero.");
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
    catch (exception &e)
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

    /* Create a vector keeping track of the current edges of the slice. */
    vector<vid_t> slice;
    
    /* Start with a single gene at or before the species tree root. */
    gtree_root.reset(new Node(null, null, null, 0, Node::none, false));
    cur_genes.push_back(gtree_root);
    slice.push_back(0);

    /* Determine the time before the first speciation. */
    double root_time = stree.time(0);
    if (g_options.count("root-length"))
        root_time = g_options["root-length"].as<double>();

    /* Run the process before the first speciation. */
    advance_process(root_time, slice, cur_genes);
    
    /* The first speciation occurs. Replace the stree root edge with
       its children, and duplicate the genes. */
    if (!stree.is_leaf(0))
        {
            speciate(stree, slice, 0, cur_genes);
        }
    /* Keep advancing and taking care of speciations until the end. */
    /* Note that leaves of stree might not end at the same time. */
    double cur_time = stree.time(0);
    while (times.empty() == false)
        {
            /* Advance the process until next speciation/leaf. */
            advance_process(times.back() - cur_time, slice, cur_genes);
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
                    cur_genes.erase(new_end, cur_genes.end());

                    swap(slice[idx], slice.back());
                    slice.pop_back();
                }
            else
                {
                    speciate(stree, slice, idx, cur_genes);
                }
        }

    if (gtree_root == 0)
        {
            cerr << PROG_NAME << ": "
                 << "gene went extinct during the process.\n";
            exit(EXIT_FAILURE);
        }



    /*
     * Output the result to the files prefix.gtree, prefix.sigma,
     * prefix.gamma, and prefix.events.
     */
    
    vector<string> filenames;
    filenames.push_back(file_prefix + ".gtree");
    filenames.push_back(file_prefix + ".sigma");
    filenames.push_back(file_prefix + ".gamma");
    filenames.push_back(file_prefix + ".events");

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

    print_gtree(*outfiles[0], stree);
    print_sigma(*outfiles[1], stree);
    print_gamma(*outfiles[2], stree);
    print_events(*outfiles[3], stree);

    return EXIT_SUCCESS;
}

//=============================================================================
//                     Helper function declarations
//=============================================================================

void postorder_gtree(vector<Node_ptr> &gtree_postorder);
void postorder_gtree_r(Node_ptr n, vector<Node_ptr> &gtree_postorder);
void print_gtree_r(ostream &out, const Tree_type &stree,
                   Node_ptr n, vector<int> &numbers);

//=============================================================================
//                   Function and member definitions.
//=============================================================================

void
advance_process(double time,
                const vector<vid_t> &slice, 
                vector<Node_ptr> &cur_genes)
{
    while (!cur_genes.empty())
        {
            /* When is the next event? */
            double total_rate = cur_genes.size() * (delta + tau + mu);
            boost::exponential_distribution<double> rng_exp(total_rate);
            time -= rng_exp(g_rng_d);
            if (time <= 0)
                break;

            /* Choose a gene to which something will happen*/
            unsigned chosen_idx = g_rng_ui(cur_genes.size());
            Node_ptr chosen = cur_genes[chosen_idx];
            
            /* Choose an event based on rates. */
            double event = g_rng_d();
            if (event < delta / (delta + tau + mu))             // duplication
                {
                    Node_ptr g1(new Node(chosen, null, null,
                                         chosen->stree_label,
                                         Node::none, false));
                    Node_ptr g2(new Node(chosen, null, null,
                                         chosen->stree_label,
                                         Node::none, false));
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
                            /* Choose which trunk to transfer to. */
                            unsigned transfer_trunk_idx
                                = g_rng_ui(slice.size() - 1);
                            if (slice[transfer_trunk_idx] == chosen->stree_label)
                                transfer_trunk_idx = slice.size() - 1;
                            
                            /* Create new children for the chosen gene. */
                            Node_ptr g1(new Node(chosen, null, null,
                                                 chosen->stree_label,
                                                 Node::none, false));
                            Node_ptr g2(new Node(chosen, null, null,
                                                 slice[transfer_trunk_idx],
                                                 Node::none, true));
                            chosen->left = g1;
                            chosen->right = g2;
                            chosen->event = Node::transfer;
                            cur_genes[chosen_idx] = g1;
                            cur_genes.push_back(g2);
                        }
                }
            else                                                // loss
                {
                    /* If only one gene remaining, process stops. */
                    if (cur_genes.size() == 1)
                        {
                            cur_genes.pop_back();
                            gtree_root.reset();
                            break;
                        }

                    /* Rewire pointers to remove chosen gene. */
                    Node_ptr parent = chosen->parent;
                    Node_ptr grand_parent = parent->parent;
                    Node_ptr sibling =
                        parent->left == chosen ? parent->right : parent->left;
                    sibling->parent = grand_parent;
                    sibling->edge_is_transfer =
                        sibling->edge_is_transfer || parent->edge_is_transfer;
                    if (grand_parent == 0)
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
         vector<Node_ptr>      &cur_genes)
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

            Node_ptr cur_gene = cur_genes[i];
            Node_ptr g1(new Node(cur_gene, null, null,
                                 stree.left(cur_snode),
                                 Node::none, false));
            Node_ptr g2(new Node(cur_gene, null, null,
                                 stree.right(cur_snode),
                                 Node::none, false));
            cur_gene->left = g1;
            cur_gene->right = g2;
            cur_gene->event = Node::speciation;
            cur_genes[i] = g1;
            cur_genes.push_back(g2);
        }
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
    for (unsigned i = 0; i < gtree_postorder.size(); ++i)
        {
            if (gtree_postorder[i]->edge_is_transfer)
                out << i << " ";
        }
    out << "\n";
    
    out << "Duplications:    ";
    for (unsigned i = 0; i < gtree_postorder.size(); ++i)
        {
            if (gtree_postorder[i]->event == Node::duplication &&
                !gtree_postorder[i]->edge_is_transfer)
                {
                    out << i << " ";
                }
        }
    out << "\n";
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
}
