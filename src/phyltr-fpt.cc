/*
 * phyltr-ftp.cc
 *
 * Source file for the phyltr-fpt application in the PhylTr
 * package. This application implements the fixed-prameter-tractable
 * algorithm for finding most parsimonious DTL-scenarios that
 * reconcile a species tree with a gene tree. [ !!insert reference
 * here once the algorithm has been publised ] At the moment only the
 * simple version is implemented. The so called depth-improvement
 * should probably be implemented as a seperate application.
 */

#include "common.hh"
#include "utils/binary-tree.hh"

extern "C" {
#include <NHparser.h>
}
#include <boost/program_options.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/concept_check.hpp>
#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cerrno>
#include <fstream>
#include <iterator>
#include <stack>

//============================================================================
//             Namespace declarations and using directives.
//============================================================================

using namespace  std;
using namespace  boost::lambda;
using            boost::dynamic_bitset;
using            boost::shared_ptr;

//============================================================================
//                   Typedefs and class declarations
//============================================================================

class Candidate;

typedef Binary_tree<double>    Tree_type;
typedef Tree_type::vid_t       vid_t;
typedef float                  Cost_type;
typedef shared_ptr<Candidate>  Candidate_ptr;

/*
 * Class Candidate
 *
 * Used to represent candidates for DTL-scenarios. The class relies on
 * the global variables holding the input (i.e., the species/gene
 * tree, sigma, duplication- and transfer-cost) in order to do its
 * work. In all the member functions, 'v' is used to indicate a gene
 * tree vertex.
 */
struct Candidate {
    dynamic_bitset<> speciations, duplications, transfer_edges;

    Candidate();

    bool is_speciation(vid_t v) const;
    bool is_duplication(vid_t v) const;
    bool is_transfer_edge(vid_t v) const;
    bool is_transfer(vid_t v) const;
    bool is_unassigned(vid_t v) const;
    Cost_type cost() const;
    void set_speciation(vid_t v);
    void set_duplication(vid_t v);
    void set_transfer_edge(vid_t v);
    void set_final();

    
    bool operator<(const Candidate &c) const;
};

/*
 * struct Program_input
 * 
 * An object g_program_input of Program_input is used to hold input to
 * the program and should be accessible by all source files.
 */
struct Program_input {
    Tree_type      species_tree;
    Tree_type      gene_tree;
    vector<vid_t>  sigma;
    Cost_type      duplication_cost;
    Cost_type      transfer_cost;
    Cost_type      min_cost;
    Cost_type      max_cost;

    Program_input() : species_tree(), gene_tree(), sigma(),
                      duplication_cost(0), transfer_cost(0),
                      min_cost(0), max_cost(0)
    {}
};
//============================================================================
//                   Global Constants and variables.
//============================================================================

const string PROG_NAME = "phyltr-fpt";
const string USAGE = "Usage: " + PROG_NAME + 
    " [OPTION]... SPECIES_TREE_FILE GENE_TREE_FILE MAP_FILE\n"
    "                  MIN_COST MAX_COST";

Program_input g_program_input = Program_input();
/* Variable holding the command line options. */
boost::program_options::variables_map g_program_options;

//=============================================================================
//                        Function declarations
//=============================================================================

/*
 * fpt_algorithm()
 *
 * The function implementing the FPT-algorithm for reconciling gene
 * trees and species trees. The function will find and output all
 * valid final candidates with cost in the range
 * [g_program_input.min_cost, g_program_input.max_cost] to 'solution'.
 * The only operation performed on 'solutions' is
 *
 *     *solutions++ = cand_ptr,
 *
 * where cand_ptr is a Candidate_ptr pointing to a final candidate.
 */
template<typename OutputIterator>
void fpt_algorithm_simple(OutputIterator solutions);

/*
 * compute_lambda()
 *
 * Computes the least common ancestor mapping of the gene tree into
 * the species tree given a candidate, and stores the result in the
 * vector 'lambda'. The size of 'lambda' must be at least
 * g_program_input.gene_tree.size().
 */
void compute_lambda(const Candidate &c, vector<vid_t> &lambda);

/*
 * compute_highest_mapping()
 *
 * Given a candidate 'c', and precomputed 'lambda', computes the
 * highest possible mapping of gene tree vertices into the species
 * tree and stores the result in 'highest'. The sizes of both 'lambda'
 * and 'highest' must be at least g_program_input.gene_tree.size().
 */
void compute_highest_mapping(const Candidate &c,
                             const vector<vid_t> &lambda,
                             vector<vid_t> &highest);

/*
 * find_unresolved()
 *
 * Given a candidate c and precomputed lambda, finds an unassigned
 * gene tree vertex such that it and one of its children is mapped by
 * lambda to the same species tree vertex. Returns Tree_type::NONE if
 * no such vertex exists. 
 */
vid_t find_unresolved(const Candidate &c, const vector<vid_t> &lambda);

/*
 * is_elegant()
 *
 * Given a _final_ candidate c and precomputed lambda, returns true if
 * the candidate is elegant and false otherwise. A candidate is called
 * elegant if there are no unnecessary duplications or transfers,
 * i.e., no duplication can be converted to a speciation with the
 * candidate remaining valid, and no transfer can be converted to a
 * speciation with the candidate remaining valid.
 */
bool is_elegant(const Candidate &c, const vector<vid_t> &lambda);

/*
 * find_comapped_vertices()
 *
 * Finds all descendants of 'u' (including 'u') that are not seperated
 * from 'u' by a transfer edge, and are mapped by lambda to the same
 * species tree vertex. The descendants will be output via 'iter'.
 */
template<typename OutIter>
void find_comapped_vertices(vid_t u,
                            const Candidate &c,
                            const vector<vid_t> &lambda,
                            OutIter iter);

//=============================================================================
//         Template and inline function and member definitions.
//=============================================================================

void 
compute_lambda(const Candidate &c, vector<vid_t> &lambda)
{
    const Tree_type &G = g_program_input.gene_tree;
    const Tree_type &S = g_program_input.species_tree;
    const vector<vid_t> &sigma = g_program_input.sigma;

    for (vid_t u = G.postorder_begin(); 
         u != Tree_type::NONE; 
         u = G.postorder_next(u))
        {
            /* Take care of gene tree leaves and continue. */
            if (G.is_leaf(u))
                {
                    lambda[u] = sigma[u];
                    continue;
                }
            
            vid_t v = G.left(u);
            vid_t w = G.right(u);
            
            if (c.is_transfer_edge(v))
                {
                    lambda[u] = lambda[w];
                }
            else if (c.is_transfer_edge(w))
                {
                    lambda[u] = lambda[v];
                }
            else
                {
                    lambda[u] = S.lca(lambda[w], lambda[v]);
                }
        }
}

void 
compute_highest_mapping(const Candidate &c,
                        const vector<vid_t> &lambda,
                        vector<vid_t> &highest)
{
    const Tree_type &G = g_program_input.gene_tree;
    const Tree_type &S = g_program_input.species_tree;
    const vector<vid_t> &sigma = g_program_input.sigma;

    /*
     * We define a function C(x, y) : V(S) x V(S) -> V(S). y must be a
     * proper descendant of x in the species tree. The function
     * returns the unique child of x that is an ancestor of y.
     */
    struct {
        vid_t operator()(vid_t x, vid_t y) {
            vid_t left = g_program_input.species_tree.left(x);
            vid_t right = g_program_input.species_tree.right(x);
            return  g_program_input.species_tree.descendant(y, left) ? left : right;
        }
    } C;

    /*
     * First, take care of the root of G.
     */
    if (c.is_speciation(0))
        {
            highest[0] = lambda[0];
        }
    else if (c.is_duplication(0))
        {
            highest[0] = 0;
        }
    else /* If the root is a transfer vertex. */
        {
            /* Let v be the transfered child of the root */
            vid_t v = G.left(0);
            v = c.is_transfer_edge(v) ? v : G.right(0);
            highest[0] = C(S.lca(lambda[0], lambda[v]), lambda[0]);
        }
    
    /*
     * Next, take care of the rest of the vertices from the root and down.
     */
    for (vid_t u = G.preorder_next(0);
         u != G.NONE;
         u = G.preorder_next(u))
        {
            vid_t pu = G.parent(u);
            vid_t x = lambda[pu];
            vid_t y = lambda[u];

            if (G.is_leaf(u))
                {
                    highest[u] = sigma[u];
                }
            else if (c.is_speciation(u))
                {
                    highest[u] = lambda[u];
                }
            else /* If u is a duplication or a transfer vertex. */
                {
                    /*
                     * Let z be the highest possible mapping of u when
                     * considering only p(u). If p(u) is a duplication
                     * or if p(u) is a transfer but u is not the
                     * transfered vertex, z = highest[pu]. Otherwise,
                     * if p(u) is a speciation, z = C(x, y), and if
                     * u is the transfered vertex, then z = C(lca(x, y), y)
                     */
                    vid_t z = highest[pu];
                    if (c.is_speciation(pu))
                        {
                            z = C(x, y);
                        }
                    else if (c.is_transfer_edge(u))
                        {
                            z = C(S.lca(x, y), y);
                        }
                    /*
                     * Let z_prime be the highest possible mapping of
                     * u when considering its children only. z_prime
                     * is the root of x unless u is a transfer. In
                     * that case, if v is the transferred child,
                     * z_prime = C(lca(lambda[u], lambda[v]), lambda[u]).
                     */
                    vid_t z_prime = 0;
                    if (c.is_transfer(u))
                        {
                            /* Let v be the transferred child of u. */
                            vid_t v = G.left(u);
                            v = c.is_transfer_edge(v) ? v : G.right(u);
                            
                            z_prime = C(S.lca(y, lambda[v]), y);
                        }
                    /*
                     * Since z and z_prime are both ancestors of
                     * lambda[u], we know that they are
                     * comparable. The one that is minimal in S is
                     * then the highest possible mapping of u.
                     */
                    highest[u] = S.descendant(z, z_prime) ? z : z_prime;
                }
        }
}

vid_t 
find_unresolved(const Candidate &c, const vector<vid_t> &lambda)
{
    const Tree_type &G = g_program_input.gene_tree;

    for (vid_t u = G.preorder_begin();
         u != G.NONE;
         u = G.preorder_next(u))
        {
            if (G.is_leaf(u) || !c.is_unassigned(u))
                continue;
            
            vid_t v = G.left(u);
            vid_t w = G.right(u);

            if (lambda[u] == lambda[v] || lambda[u] == lambda[w])
                return u;
        }
    
    return G.NONE;
}

bool
is_elegant(const Candidate &c, const vector<vid_t> &lambda)
{
    const Tree_type &G = g_program_input.gene_tree;
    const Tree_type &S = g_program_input.species_tree;

    static vector<vid_t> highest(G.size());

    /* Note that c must be a final scenario, we are not checking this? */
    compute_highest_mapping(c, lambda, highest);
    
    /*
     * Check that the children of duplications are mapped by lambda to
     * comparable species tree vertices. Otherwise, the duplication is
     * unnecessary.
     */
    for (vid_t d = c.duplications.find_first();
         d != c.duplications.npos;
         d = c.duplications.find_next(d))
        {
            vid_t v = G.left(d);
            vid_t w = G.right(d);
            if (!S.descendant(lambda[v], lambda[w]) && 
                !S.descendant(lambda[w], lambda[v]))
                return false;
        }

    /*
     * Check if the parents of transfer vertices can be mapped high
     * enough so that the transfer can be converted to a
     * speciation. If so, the transfer is unnecessary.
     */
    for (vid_t v = c.transfer_edges.find_first();
         v != c.transfer_edges.npos;
         v = c.transfer_edges.find_next(v))
        {
            /*
             * Let (u, v) be the transfer edge we are considering, let
             * pu = p(u), and x = lca{lambda[u], lambda[v]}
             */
            vid_t u = G.parent(v);
            vid_t pu = G.parent(u);
            vid_t x = S.lca(lambda[u], lambda[v]);
            
            /* The root of G is always an unnecessary transfer vertex. */
            if (u == 0)
                return false;

            /*
             * If p(u) is a speciation and x is a proper descendant of
             * highest[p(u)] = lambda[p(u)], then the transfer is
             * unnecessary.
             */
            if (c.is_speciation(pu) && 
                S.descendant(x, lambda[pu]) && 
                x != lambda[pu])
                {
                    return false;
                }
            
            /*
             * If p(u) is not a speciation, then it is enough for x to
             * be a descendant of highest[p(u)] for the transfer to be
             * unnecessary.
             */
            if (!c.is_speciation(pu) &&
                S.descendant(x, highest[pu]))
                {
                    return false;
                }
        }

    return true;
}


template<typename OutIter>
void 
find_comapped_vertices(vid_t u,
                       const Candidate &c,
                       const vector<vid_t> &lambda,
                       OutIter iter)
{
    boost::function_requires< boost::OutputIteratorConcept<OutIter, vid_t> >();

    const Tree_type &G = g_program_input.gene_tree;

    /* u is always a comapped vertex of itself. */
    *iter++ = u;

    if (G.is_leaf(u))
        return;

    vid_t v = G.left(u);
    vid_t w = G.right(u);

    if (!c.is_transfer_edge(v) && lambda[v] == lambda[u])
        find_comapped_vertices(v, c, lambda, iter);
    if (!c.is_transfer_edge(w) && lambda[w] == lambda[u])
        find_comapped_vertices(w, c, lambda, iter);
}

template<typename OutputIterator>
void
fpt_algorithm_simple(OutputIterator solutions)
{
    boost::function_requires< boost::OutputIteratorConcept<OutputIterator, Candidate_ptr> >();

    static bool print_inelegant_warning = true;

    const Tree_type &G = g_program_input.gene_tree;
    const Tree_type &S = g_program_input.species_tree;
    const Cost_type duplication_cost = g_program_input.duplication_cost;
    const Cost_type transfer_cost = g_program_input.transfer_cost;
    const Cost_type max_cost = g_program_input.max_cost;
    const Cost_type min_cost = g_program_input.min_cost;

    stack<Candidate_ptr> Q;
    vector<vid_t> lambda(G.size());
    
    Q.push(Candidate_ptr(new Candidate()));
    while (!Q.empty())
        {
            Candidate_ptr c1 = Q.top(); Q.pop();
            
            compute_lambda(*c1, lambda);

            vid_t u = find_unresolved(*c1, lambda);
            if (u == Tree_type::NONE) /* Nothing more to do for candidate. */
                {
                    c1->set_final();
                    Cost_type cost = 
                        c1->duplications.count() * duplication_cost +
                        c1->transfer_edges.count() * transfer_cost;
                    if (cost <= max_cost &&
                        cost >= min_cost && 
                        is_elegant(*c1, lambda))
                        {
                            *solutions++ = c1;
                        }
                    else if (print_inelegant_warning)
                        {
                            clog << PROG_NAME << ": Found inelegant scenario!\n";
                            clog << c1->duplications << "\n"
                                 << c1->transfer_edges << "\n";
                            print_inelegant_warning = false;
                        }
                    continue;
                }
            
            /* Here we found an unresolved gene tree vertex u. */

            /* Do not continue if c1 is a hopeless case. */
            if (c1->cost() + min(duplication_cost, transfer_cost) > max_cost)
                continue;

            vid_t x = lambda[u];
            
            /* If u is mapped to a leaf it must be a duplication. */
            if (S.is_leaf(x))
                {
                    c1->set_duplication(u);
                    if (c1->cost() <= max_cost)
                        Q.push(c1);
                    continue;
                }

            /*
             * Find all non-transfered descendants of u that are
             * mapped to the same species tree vertex. Note that R
             * will contain at least two vertices.
             */
            vector<vid_t> R;
            find_comapped_vertices(u, *c1, lambda, back_inserter(R));

            /* If a comapped vertex is a duplication or speciation, then
               all its comapped ancestors are duplications. */

            bool forced_duplication = false;
            for (unsigned i = 0; i < R.size(); ++i)
                {
                    if (c1->is_speciation(R[i]) || c1->is_duplication(R[i]))
                        {
                            forced_duplication = true;

                            vid_t anc = G.parent(R[i]);
                            while (anc != G.NONE && lambda[anc] == x)
                                {
                                    if (c1->is_unassigned(anc))
                                        c1->set_duplication(anc);
                                    anc = G.parent(anc);
                                }
                            if (c1->cost() < max_cost)
                                Q.push(c1);
                            break;
                        }
                }

            if (forced_duplication)
                continue;

            vid_t v = 0, w = 0;

            /* Find a minimal vertex in R and resolve it in three ways. */
            for (unsigned i = 0; i < R.size(); ++i)
                {
                    v = G.left(R[i]);
                    w = G.right(R[i]);
                    if (lambda[v] != x && lambda[w] != x)
                        {
                            u = R[i];
                            break;
                        }
                }
            
            /* Here u is a minimal unresolved gene tree vertex with
               children v, w, such that lambda[v] != lambda[w] != lambda[u]. */
            Candidate_ptr c2(new Candidate(*c1));
            Candidate_ptr c3(new Candidate(*c1));

            /* Make c1 the candidate with u as duplication. */
            c1->set_speciation(u);
            vid_t anc = G.parent(u);
            while (anc != G.NONE && lambda[anc] == lambda[u])
                {
                    if (c1->is_unassigned(anc))
                        c1->set_duplication(anc);
                    anc = G.parent(anc);
                }

            /* Make c2 the candidate with (u,v) as transfer. */
            c2->set_transfer_edge(v);

            /* Make c3 the candidate with (u,w) as transfer. */
            c3->set_transfer_edge(w);

            /* Insert the candidates into Q. */
            if (c3->cost() <= max_cost)
                {
                    Q.push(c3);
                }
            if (c2->cost() <= max_cost)
                {
                    Q.push(c2);
                }
            if (c1->cost() <= max_cost)
                {
                    Q.push(c1);
                }
        }

}

//=============================================================================
//                                main()
//=============================================================================

int
main(int argc, char *argv[])
{
    /*
     * Parse the command line arguments.
     */
    namespace po = boost::program_options;

    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    try
        {
            /* Declare options that are described when --help is given. */
            visible_opts.add_options()
                ("help", "display this help and exit")
                ("transfer-cost,t",
                 po::value<Cost_type>(&g_program_input.transfer_cost)->default_value(1.0),
                 "Cost of transfer events")
                ("duplication-cost,d", 
                 po::value<Cost_type>(&g_program_input.duplication_cost)->default_value(1.0),
                 "Cost of duplication events")
                ;

            /* Declare positional options. */
            hidden_opts.add_options()
                ("species-tree-file", po::value<string>())
                ("gene-tree-file", po::value<string>())
                ("map-file", po::value<string>())
                ("min-cost", po::value<Cost_type>(&g_program_input.min_cost))
                ("max-cost", po::value<Cost_type>(&g_program_input.max_cost))
                ;
    
            po::positional_options_description positional_options;
            positional_options.add("species-tree-file", 1);
            positional_options.add("gene-tree-file", 1);
            positional_options.add("map-file", 1);
            positional_options.add("min-cost", 1);
            positional_options.add("max-cost", 1);

            /* Put all options in one options_description object. */
            all_options.add(visible_opts).add(hidden_opts);

            po::command_line_parser parser(argc, argv);
            parser.options(all_options);
            parser.positional(positional_options);
    
            /* Parse and store the results. */
            po::store(parser.run(), g_program_options);
            po::notify(g_program_options);
        }
    catch (po::error &e)
        {
            cerr << PROG_NAME << ": " << e.what() << "\n"
                 << "Try '" << PROG_NAME << " --help' "
                 << "for more information.\n";
            exit(EXIT_FAILURE);
        }


    /* Show help message if --help was given. */
    if (g_program_options.count("help"))
        {
            cout << USAGE << "\n"
                 << visible_opts << "\n";
            exit(EXIT_SUCCESS);
        }

    /* Check that all required positional arguments are given. */
    if (g_program_options.count("species-tree-file") == 0 ||
        g_program_options.count("gene-tree-file") == 0 ||
        g_program_options.count("map-file") == 0 ||
        g_program_options.count("min-cost") == 0 ||
        g_program_options.count("max-cost") == 0)
        {
            cerr << PROG_NAME << ": too few arguments.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /*
     * Check that min_cost <= max_cost.
     */
    if (g_program_input.min_cost > g_program_input.max_cost)
        {
            cerr << PROG_NAME << ": "
                 << "MIN_COST cannot be greater than MAX_COST\n";
            exit(EXIT_FAILURE);
        }

    /*
     * Read the species- and gene-tree files.
     */
    string species_tree_filename =
        g_program_options["species-tree-file"].as<string>();
    string gene_tree_filename =
        g_program_options["gene-tree-file"].as<string>();

    errno = 0;
    FILE *species_tree_file = fopen(species_tree_filename.c_str(), "r");
    if (errno)
        {
            cerr << PROG_NAME << ": "
                 << "failed to open species tree file '"
                 << species_tree_filename << "'\n";
            perror(PROG_NAME.c_str());
            exit(EXIT_FAILURE);
        }

    FILE *gene_tree_file = fopen(gene_tree_filename.c_str(), "r");
    if (errno)
        {
            cerr << PROG_NAME << ": "
                 << "failed to open gene tree file '"
                 << gene_tree_filename << "'\n";
            perror(PROG_NAME.c_str());
            exit(EXIT_FAILURE);
        }

    NHNode *NH_species_root = NH_read_tree(species_tree_file, stderr,
                                           species_tree_filename.c_str());
    if (NH_species_root == 0)
        {
            cerr << PROG_NAME << ": "
                 << "failed to read species tree\n";
            exit(EXIT_FAILURE);
        }
    NHNode *NH_gene_root = NH_read_tree(gene_tree_file, stderr,
                                        gene_tree_filename.c_str());
    if (NH_gene_root == 0)
        {
            cerr << PROG_NAME << ": "
                 << "failed to read gene tree\n";
            exit(EXIT_FAILURE);
        }

    /* check that trees are binary. */
    if (!NHtree_is_binary(NH_species_root))
        {
            cerr << PROG_NAME << ": "
                 << "species tree is not binary\n";
            exit(EXIT_FAILURE);
        }
    if (!NHtree_is_binary(NH_gene_root))
        {
            cerr << PROG_NAME << ": "
                 << "gene tree is not binary\n";
            exit(EXIT_FAILURE);
        }

    /* transform the NH-trees to Binary_trees. */
    try
        {
            create_binary_tree(NH_species_root, g_program_input.species_tree, 0);
            create_binary_tree(NH_gene_root, g_program_input.gene_tree, 0);
        }
    catch (const exception &e)
        {
            cerr << PROG_NAME << ": " << e.what() << "\n";
            exit(EXIT_FAILURE);
        }
    
    /*
     * Read map file
     */
    g_program_input.sigma.resize(g_program_input.gene_tree.size());
    try 
        {
            create_gene_species_map(g_program_input.species_tree,
                                    g_program_input.gene_tree,
                                    g_program_options["map-file"].as<string>(),
                                    g_program_input.sigma);
        }
    catch (const exception &e)
        {
            cerr << PROG_NAME << ": " << e.what() << "\n";
            exit(EXIT_FAILURE);
        }

    /*
     * Check that the costs given are not too big in the worst
     * case. No-one needs such huge costs on duplications and
     * transfers.
     */
    unsigned max_elements = 
        max(g_program_input.species_tree.size() / 2, 
            g_program_input.gene_tree.size() / 2) + 1;
    Cost_type max_allowed = numeric_limits<Cost_type>::max() / max_elements;

    if (g_program_input.transfer_cost >= max_allowed)
        {
            cerr << PROG_NAME << ": transfer cost is too big.\n";
            exit(EXIT_FAILURE);
        }
    if (g_program_input.duplication_cost >= max_allowed)
        {
            cerr << PROG_NAME << ": duplication cost is too big.\n";
            exit(EXIT_FAILURE);
        }

    vector<Candidate_ptr> solutions;
    fpt_algorithm_simple(back_inserter(solutions));

    /*
     * Switch to using the postorder numbering of vertices in the
     * scenarios.
     */
    const Tree_type &G = g_program_input.gene_tree;
    const Tree_type &S = g_program_input.species_tree;

    vector<vid_t> gene_tree_numbering(G.size());
    vector<vid_t> species_tree_numbering(S.size());
    get_postorder_numbering(G, gene_tree_numbering);
    get_postorder_numbering(S, species_tree_numbering);

    BOOST_FOREACH(Candidate_ptr cp, solutions)
        {
            dynamic_bitset<> new_duplications(cp->duplications.size());
            dynamic_bitset<> new_transfers(cp->transfer_edges.size());
            for (unsigned i = cp->duplications.find_first();
                 i != dynamic_bitset<>::npos;
                 i = cp->duplications.find_next(i))
                {
                    new_duplications.set(gene_tree_numbering[i]);
                }
            for (unsigned i = cp->transfer_edges.find_first();
                 i != dynamic_bitset<>::npos;
                 i = cp->transfer_edges.find_next(i))
                {
                    new_transfers.set(gene_tree_numbering[i]);
                }
            cp->duplications = new_duplications;
            cp->transfer_edges = new_transfers;
            // cp->set_final();
            // We won't do set_final because the indices of the bitsets
            // now do not correspond to vertex ids.
        }

    sort(solutions.begin(), solutions.end(), (*_1) < (*_2));
    cout << "Number of scenarios: " << solutions.size() << "\n\n";

    BOOST_FOREACH(Candidate_ptr cp, solutions)
        {
            cout << "Transfer edges:\t";
            for (unsigned j = cp->transfer_edges.find_first();
                 j != dynamic_bitset<>::npos;
                 j = cp->transfer_edges.find_next(j))
                {
                    cout << j << " ";
                }
            cout << '\n';
            cout << "Duplications:\t";
            for (unsigned j = cp->duplications.find_first();
                 j != dynamic_bitset<>::npos;
                 j = cp->duplications.find_next(j))
                {
                    cout << j << " ";
                }
            cout << "\n\n";
        }
    
    return EXIT_SUCCESS;
}


//=============================================================================
//                     Helper function declarations
//=============================================================================

//=============================================================================
//                   Function and member definitions.
//=============================================================================

Candidate::Candidate()
    : speciations(g_program_input.gene_tree.size()),
      duplications(g_program_input.gene_tree.size()),
      transfer_edges(g_program_input.gene_tree.size())
{
}

bool
Candidate::is_speciation(vid_t v) const
{
    return speciations[v];
}

bool
Candidate::is_duplication(vid_t v) const
{
    return duplications[v];
}

bool
Candidate::is_transfer_edge(vid_t v) const
{
    return transfer_edges[v];
}

bool
Candidate::is_transfer(vid_t v) const
{
    const Tree_type &G = g_program_input.gene_tree;

    if (G.is_leaf(v))
        return false;

    return is_transfer_edge(G.left(v)) || is_transfer_edge(G.right(v));
}

bool
Candidate::is_unassigned(vid_t v) const
{
    return !(is_speciation(v) || is_duplication(v) || is_transfer(v));
}

void
Candidate::set_speciation(vid_t v)
{
    speciations.set(v);
}

void
Candidate::set_duplication(vid_t v)
{
    duplications.set(v);
}

void
Candidate::set_transfer_edge(vid_t v)
{
    transfer_edges.set(v);
}

void
Candidate::set_final()
{
    dynamic_bitset<> transfer_vertices(duplications.size());
    for (unsigned i = transfer_edges.find_first();
         i != transfer_edges.npos;
         i = transfer_edges.find_next(i))
        {
            transfer_vertices.set(g_program_input.gene_tree.parent(i));
        }

    speciations = ~(duplications | transfer_vertices);
}

Cost_type
Candidate::cost() const
{
    return duplications.count() * g_program_input.duplication_cost
        + transfer_edges.count() * g_program_input.transfer_cost;
}

bool
Candidate::operator<(const Candidate &c) const
{
    /*
     * This function must not depend on this->speciations, since the
     * indices might not refer to vertex ids, but might rather
     * correspond to preorder numbering.
     */
    
    if (cost() < c.cost())
        return true;
    if (cost() > c.cost())
        return false;
    if (transfer_edges.count() < c.transfer_edges.count())
        return true;
    if (transfer_edges.count() > c.transfer_edges.count())
        return false;
    if (transfer_edges < c.transfer_edges)
        return true;
    if (transfer_edges > c.transfer_edges)
        return false;
    return duplications < c.duplications;
}
