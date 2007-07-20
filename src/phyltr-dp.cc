/*
 * phyltr-dp.cc
 *
 * Source file for the phyltr-dp application in the PhylTr
 * package. This application implements a dynamic programming
 * algorithm for finding the minimum cost of DTL-scenarios reconciling
 * a species tree and a gene tree. [ insert reference here once the
 * algorithm is published. ] The program also backtracks to find all
 * optimal (i.e., most parsimonious) DTL-scenarios.
 *
 * TODO:
 *
 * - Change the output to something actually meaningful. Preferably
 * using XML.
 * - Implement a true C++ tree parser for reading the phylogenetic
 * trees.
 */

#include "common.hh"
#include "utils/binary-tree.hh"

extern "C" {
#include <NHparser.h>
}
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cerrno>
#include <map>
#include <list>
#include <vector>
#include <iterator>
#include <fstream>
#include <limits>
#include <iterator>

using namespace std;
namespace po = boost::program_options;
using boost::multi_array;
using boost::dynamic_bitset;

/*
 * Typedefs and class declarations
 */

typedef Binary_tree<double> Tree_type;
typedef Tree_type::vid_t    vid_t;
typedef float               Cost_type; // Instead of double to save memory.

/*
 * Class Scenario
 *
 * Used to keep track of (possibly partial) DTL-scenarios. Since we do
 * not care about the mapping of the gene tree into the species tree,
 * the only relevant information needed is the set of duplication and
 * transfer vertices. We use bitsets in order to speed up union
 * operations, though a vector<bool> might work just as well.
 *
 * Note that a Scenario object may represent only a partial scenario,
 * in which only a subtree of the gene tree is being
 * considered.
 *
 * This class is a candidate for breaking out and making it a common
 * class to be used by programs in the PhylTr package. In that case,
 * we should probably hide the implementation and define a union
 * operation. For now, a Scenario is just a pair of bitsets.
 */
struct Scenario {
    explicit Scenario(unsigned num_vertices);
    Cost_type cost(Cost_type duplication_cost, Cost_type transfer_cost);

    dynamic_bitset<> duplications;
    dynamic_bitset<> transfer_edges;
};


/*
 * Class Below
 *
 * A two-dimensional matrix of Below-objects are used in the
 * implementation of the dynamic programming algorithm. Data member
 * 'cost' is initialized to infinity. The vector 'events' is used in
 * backtracking and contains information about what computation that
 * led to the cost. The vectors 'scenarios_below' and 'scenarios_at'
 * are used in the backtracking algorithm (which is itself a
 * DP-alg!).
 *
 * The 'cost' member of each Below-object in the DP-matrix represents
 * the minimum cost of placing a gene tree vertex u, with children
 * left_u and right_u, at a descendant of a species tree vertex x,
 * with children x_left and x_right. events[S] should be set to true
 * iff the cost of u being a speciation at x, with the left_u below
 * left_ x, and right_u placed below right_x, equals
 * 'cost'. events[S_REV] should be set to true just as for events[S]
 * except that u_left is placed below x_right and u_right is placed
 * below x_left. events[T_LEFT] should be set to true iff the cost of
 * the edge (u, left_u) being a transfer, with u placed at x, equals
 * 'cost'. The same applies to events[T_RIGHT], except that the
 * transfer edge considered is (u, right_u). events[D] should be set
 * to true iff the cost of u being a duplication at x equals
 * 'cost'. events[B_LEFT] and events[B_RIGHT] should be set iff the
 * cost of placing u below x_left and below x_right respectively equal
 * 'cost'
 *
 * During backtracking, intermediate results of scenarios are stored
 * in 'scenarios_below' and 'scenarios_at'. The vector
 * 'scenarios_below' contains all scenarios where u is mapped at or
 * below x, while 'scenarios_at' contains the scenarios in which
 * u is mapped to x only.
 */
struct Below {
    enum {S, S_REV, T_LEFT, T_RIGHT, D, B_LEFT, B_RIGHT, N_EVENTS};
    
    Below();

    Cost_type cost;
    vector<bool> events;
    vector<Scenario> scenarios_below;
    vector<Scenario> scenarios_at;
    bool scenarios_below_computed;
    bool scenarios_at_computed;
};


/*
 * Class Outside
 *
 * A two-dimensional matrix of Outside-objects are used in the
 * implementation of the dynamic programming algorithm. Data member
 * 'cost' is initialized to infinity. The vector species_vertices is
 * used in backtracking.
 *
 * The member 'cost' of an Outside-object should represents the
 * minimum cost of placing a gene tree vertex u at a species tree
 * vertex incomparable to a species tree vertex x. The vector
 * species_vertices should contain the vertex ids of every vertex y
 * incomparable to x such that the minimal cost of placing u below y
 * equals 'cost'
 */
struct Outside {
    Outside();

    Cost_type cost;
    vector<vid_t> species_vertices;
};

/*
 * Global Constants
 */

const string PROG_NAME = "phyltr";
const string USAGE = 
    "Usage: " + PROG_NAME + " [OPTION]... SPECIES_TREE_FILE GENE_TREE_FILE "
    "GENEMAP_FILE";
const Cost_type COST_INF = numeric_limits<Cost_type>::infinity();


/*
 * Function declarations
 */ 

/*
 * backtrack_scenarios_below()
 *
 * The function computes each unique (in the sense described above in
 * the description of the Scenario class) most parsimonious partial
 * DTL-scenario (a Scenario object) such that 'gene_vertex' is placed
 * at or below 'species_vertex', and inserts them into
 * below[gene_vertex][species_vertex].scenarios_below.
 */
void backtrack_scenarios_below(const Tree_type &spcies_tree,
                               const Tree_type &gene_tree,
                               const vector<vid_t> &sigma,
                               multi_array<Below, 2> &below,
                               const multi_array<Outside, 2> &outside,
                               Cost_type duplication_cost,
                               Cost_type transfer_cost,
                               vid_t gene_vertex, 
                               vid_t species_vertex);

int
main(int argc, char *argv[])
{
    /*
     * Declaration of variables holding the options from the command line.aaaaa
     */
    double transfer_cost;
    double duplication_cost;
    string species_tree_filename;
    string gene_tree_filename;
    string map_filename;

    /*
     * Parse the command line options.
     */

    /* Declare options that are described when --help is given. */
    po::options_description visible_opts("Command Line Options");
    visible_opts.add_options()
        ("help", "display this help and exit")
        ("transfer-cost,t",
         po::value<double>(&transfer_cost)->default_value(1.0),
         "Cost of transfer events")
        ("duplication-cost,d", 
         po::value<double>(&duplication_cost)->default_value(1.0),
         "Cost of duplication events")
        ;

    /* Declare positional options. */
    po::options_description hidden_opts("");
    hidden_opts.add_options()
        ("species_tree_file", po::value<string>(&species_tree_filename), "")
        ("gene_tree_file", po::value<string>(&gene_tree_filename), "")
        ("map_file", po::value<string>(&map_filename), "")
        ;
    
    po::options_description all_options("");
    all_options.add(visible_opts).add(hidden_opts);

    po::positional_options_description positional_options;
    positional_options.add("species_tree_file", 1);
    positional_options.add("gene_tree_file", 1);
    positional_options.add("map_file", 1);
 
    /* Parse the arguments. */
    po::variables_map vm;
    try
        {
            po::store(po::command_line_parser(argc, argv).options(all_options)
                      .positional(positional_options).run(), vm);
            po::notify(vm);

            /* Show help message if --help was given. */
            if (vm.count("help"))
                {
                    cout << USAGE << '\n'
                         << visible_opts << '\n';
                    exit(EXIT_SUCCESS);
                }

            /* Check that all required positional arguments are given. */
            struct Too_few_arguments : public exception {
                const char *what() const throw() {
                    return "too few arguments";
                }
            };

            if (vm.count("species_tree_file") == 0 
                || vm.count("gene_tree_file") == 0 
                || vm.count("map_file") == 0)
                {
                    throw Too_few_arguments();
                }
 
        }
    catch (exception &e)
        {
            cerr << PROG_NAME << ": " << e.what() << '\n'
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /*
     * Read the species- and gene-tree files.
     */

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
    Tree_type species_tree;
    Tree_type gene_tree;
    try
        {
            create_binary_tree(NH_species_root, species_tree, 0);
            create_binary_tree(NH_gene_root, gene_tree, 0);
        }
    catch (const exception &e)
        {
            cerr << PROG_NAME << ": " << e.what() << '\n';
            exit(EXIT_FAILURE);
        }
    
    /*
     * Read map file
     */
    
    vector<vid_t> sigma;
    try 
        {
            sigma = create_gene_species_map(species_tree, gene_tree, map_filename);
        }
    catch (const exception &e)
        {
            cerr << PROG_NAME << ": " << e.what() << '\n';
            exit(EXIT_FAILURE);
        }

    /*
     * Check that the costs given are not too big in the worst
     * case. No-one needs such huge costs on duplications and
     * transfers.
     */
    unsigned max_elements = 
        max(species_tree.size() / 2, gene_tree.size() / 2) + 1;
    if (transfer_cost >= numeric_limits<Cost_type>::max() / max_elements)
        {
            cerr << PROG_NAME << ": transfer cost is too big";
            exit(EXIT_FAILURE);
        }
    if (duplication_cost >= numeric_limits<Cost_type>::max() / max_elements)
        {
            cerr << PROG_NAME << ": duplication cost is too big";
            exit(EXIT_FAILURE);
        }
    
    
    /*
     * The dynamic programming algorithm goes here.
     *
     * The variables holding the input are:
     *   species_tree
     *   gene_tree
     *   sigma
     *   duplication_cost
     *   transfer_cost
     */
    
    typedef Tree_type::vid_t vid_t;

    /* Create the below and outside matrices. */
      
    multi_array<Below, 2> below(boost::extents[gene_tree.size() + 1]
                                              [species_tree.size() + 1]);
    multi_array<Outside, 2> outside(boost::extents[gene_tree.size() + 1]
                                                  [species_tree.size() + 1]);

    /* The main loop. */
    
    for (vid_t u = gene_tree.postorder_begin();
         u != gene_tree.NONE;
         u = gene_tree.postorder_next(u))
        {
            /* First compute below[u][*]. */
            for (vid_t x = species_tree.postorder_begin();
                 x != species_tree.NONE;
                 x = species_tree.postorder_next(x))
                {
                    /* If u is a leaf use sigma. */
                    if (gene_tree.is_leaf(u))
                        {
                            if (species_tree.descendant(sigma[u], x))
                                {
                                    below[u][x].cost = 0;
                                }
                            continue;
                        }

                    /* If u is not a leaf compute costs. */
                    vector<Cost_type> costs(Below::N_EVENTS, COST_INF);

                    vid_t v = gene_tree.left(u);
                    vid_t w = gene_tree.right(u);
                    
                    costs[Below::D] = 
                        duplication_cost + below[v][x].cost + below[w][x].cost;
                    costs[Below::T_LEFT] = 
                        transfer_cost + outside[v][x].cost + below[w][x].cost;
                    costs[Below::T_RIGHT] = 
                        transfer_cost + outside[w][x].cost + below[v][x].cost;
                    if (!species_tree.is_leaf(x))
                        {
                            vid_t y = species_tree.left(x);
                            vid_t z = species_tree.right(x);
                            costs[Below::S] = 
                                below[v][y].cost + below[w][z].cost;
                            costs[Below::S_REV]
                                = below[v][z].cost + below[w][y].cost;
                            costs[Below::B_LEFT] = below[u][y].cost;
                            costs[Below::B_RIGHT] = below[u][z].cost;
                        }
                    Cost_type min_cost = 
                        *min_element(costs.begin(), costs.end());
                    below[u][x].cost = min_cost;
                    for (unsigned i = 0; i < Below::N_EVENTS; ++i)
                        {
                            if (costs[i] == min_cost)
                                below[u][x].events[i] = true;
                        }
                }

            /* Now compute outside[u][*]. */
            
            for (vid_t x = species_tree.preorder_begin();
                 x != species_tree.NONE;
                 x = species_tree.preorder_next(x))
                {
                    if (species_tree.is_leaf(x))
                        continue;
                    
                    vid_t y = species_tree.left(x);
                    vid_t z = species_tree.right(x);
                    
                    Cost_type outside_u_x = outside[u][x].cost;
                    Cost_type below_u_y = below[u][y].cost;
                    Cost_type below_u_z = below[u][z].cost;
                    outside[u][y].cost = min(outside_u_x, below_u_z);
                    outside[u][z].cost = min(outside_u_x, below_u_y);

                    if (outside[u][y].cost == outside_u_x)
                        {
                            copy(outside[u][x].species_vertices.begin(),
                                 outside[u][x].species_vertices.end(),
                                 back_inserter(outside[u][y].species_vertices));
                        }
                    if (outside[u][y].cost == below_u_z)
                        {
                            outside[u][y].species_vertices.push_back(z);
                        }

                    if (outside[u][z].cost == outside_u_x)
                        {
                            copy(outside[u][x].species_vertices.begin(),
                                 outside[u][x].species_vertices.end(),
                                 back_inserter(outside[u][z].species_vertices));
                        }
                    if (outside[u][z].cost == below_u_y)
                        {
                            outside[u][z].species_vertices.push_back(y);
                        }
                }
        }

    /* The DP-matrix is computed. */
    
    /*
     * Output the cost.
     */
    cout << "Cost: " << below[0][0].cost << '\n';
    
    /*
     * Backtrack and output all scenarios. We should do this
     * xml-style.
     */
    
    backtrack_scenarios_below(species_tree, gene_tree, sigma,
                              below, outside, duplication_cost, transfer_cost,
                              0, 0);
    vector<Scenario> &scenarios = below[0][0].scenarios_below;
    
    cout << "Number of scenarios: " << scenarios.size() << '\n';

    cout << "S: " << species_tree << '\n';
    cout << "G: " << gene_tree << '\n' << '\n';

    for (unsigned i = 0; i < scenarios.size(); ++i)
        {
            Scenario &sc = scenarios[i];
            cout << "Duplications:\t";
            for (unsigned j = sc.duplications.find_first();
                 j != sc.duplications.npos; j = sc.duplications.find_next(j))
                {
                    cout << setw(4) << j;
                }
            cout << '\n';
            cout << "Transfer edges:\t";
            for (unsigned j = sc.transfer_edges.find_first();
                 j != sc.transfer_edges.npos; j = sc.transfer_edges.find_next(j))
                {
                    cout << setw(4) << j;
                }
            cout << "\n\n";
        }
    
//     cout << "The DP-matrices (rows = genes, columns = species)\n";
//     cout << "below:\n";
//     cout << setw(4) << " ";
//     for (unsigned i = 0; i < species_tree.size(); ++i)
//         {
//             cout << setw(4) << i;
//         }
//     cout << '\n';
    
//     for (unsigned i = 0; i < gene_tree.size(); ++i)
//         {
//             cout << setw(4) << i << ":";
//             for (unsigned j = 0; j < species_tree.size(); ++j)
//                 {
//                     cout << setw(4) << below[i][j].cost;
//                 }
//             cout << '\n';
//         }


//     cout << "\noutside:\n";
//     cout << setw(4) << " ";
//     for (unsigned i = 0; i < species_tree.size(); ++i)
//         {
//             cout << setw(4) << i;
//         }
//     cout << '\n';
//     for (unsigned i = 0; i < gene_tree.size(); ++i)
//         {
//             cout << setw(4) << i << ":";
//             for (unsigned j = 0; j < species_tree.size(); ++j)
//                 {
//                     cout << setw(4) << outside[i][j].cost;
//                 }
//             cout << '\n';
//         }

//     cout << "\nbelow[u][x].scenarios_below.size():\n";
//     cout << setw(4) << " ";
//     for (unsigned i = 0; i < species_tree.size(); ++i)
//         {
//             cout << setw(4) << i;
//         }
//     cout << '\n';
    
//     for (unsigned i = 0; i < gene_tree.size(); ++i)
//         {
//             cout << setw(4) << i << ":";
//             for (unsigned j = 0; j < species_tree.size(); ++j)
//                 {
//                     cout << setw(4) << below[i][j].scenarios_below.size();
//                 }
//             cout << '\n';
//         }
//     cout << "\nbelow[u][x].scenarios_at.size():\n";
//     cout << setw(4) << " ";
//     for (unsigned i = 0; i < species_tree.size(); ++i)
//         {
//             cout << setw(4) << i;
//         }
//     cout << '\n';
    
//     for (unsigned i = 0; i < gene_tree.size(); ++i)
//         {
//             cout << setw(4) << i << ":";
//             for (unsigned j = 0; j < species_tree.size(); ++j)
//                 {
//                     cout << setw(4) << below[i][j].scenarios_at.size();
//                 }
//             cout << '\n';
//         }
//     cout << "\noutside[u][x].species_vertices.size():\n";
//     cout << setw(4) << " ";
//     for (unsigned i = 0; i < species_tree.size(); ++i)
//         {
//             cout << setw(4) << i;
//         }
//     cout << '\n';
    
//     for (unsigned i = 0; i < gene_tree.size(); ++i)
//         {
//             cout << setw(4) << i << ":";
//             for (unsigned j = 0; j < species_tree.size(); ++j)
//                 {
//                     cout << setw(4) << outside[i][j].species_vertices.size();
//                 }
//             cout << '\n';
//         }

//     for (unsigned i = 0; i < gene_tree.size(); ++i)
//         {
//             for (unsigned j = 0; j < species_tree.size(); ++j)
//                 {
//                     cout << setw(4) << i
//                          << setw(4) << j << '\t';
//                     for (unsigned k = 0; k < below[i][j].events.size();++k)
//                         {
//                             cout << setw(6) << below[i][j].events[k];
//                         }
//                     cout << "\n";
//                 }
//         }

    return EXIT_SUCCESS;
}

/*
 * Helper function declarations.
 */


/*
 * combine_scenarios()
 *
 * Combines each scenarios in one sequence with each scenario in
 * another sequence and outputs the result to the output iterator
 * 'out'. This function is a helper function used during backtracking.
 */
template<typename Initer1, typename Initer2, typename Outiter>
void combine_scenarios(Initer1 begin1, Initer1 end1,
                       Initer2 begin2, Initer2 end2,
                       Outiter out);



/*
 * backtrack_scenarios_at()
 *
 * The function computes each unique (in the sense described above in
 * the description of the Scenario class) most parsimonious partial
 * DTL-scenario (a Scenario object) such that 'gene_vertex' is placed
 * _at_ 'species_vertex', and inserts them into
 * below[gene_vertex][species_vertex].scenarios_at.
 */
void backtrack_scenarios_at(const Tree_type &spcies_tree,
                            const Tree_type &gene_tree,
                            const vector<vid_t> &sigma,
                            multi_array<Below, 2> &below,
                            const multi_array<Outside, 2> &outside,
                            Cost_type duplication_cost,
                            Cost_type transfer_cost,
                            vid_t gene_vertex,
                            vid_t species_vertex);

/*
 * backtrack_[speciation,duplicttion,transfer]_at()
 *
 * These are helper functions called by backtrack_scenarios_at.
 */

void backtrack_speciation_at(const Tree_type &spcies_tree,
                             const Tree_type &gene_tree,
                             const vector<vid_t> &sigma,
                             multi_array<Below, 2> &below,
                             const multi_array<Outside, 2> &outside,
                             Cost_type duplication_cost,
                             Cost_type transfer_cost,
                             vid_t gene_vertex_left,
                             vid_t gene_vertex_right,
                             vid_t species_vertex);

void backtrack_duplication_at(const Tree_type &spcies_tree,
                              const Tree_type &gene_tree,
                              const vector<vid_t> &sigma,
                              multi_array<Below, 2> &below,
                              const multi_array<Outside, 2> &outside,
                              Cost_type duplication_cost,
                              Cost_type transfer_cost,
                              vid_t gene_vertex_1,
                              vid_t gene_vertex_2,
                              vid_t species_vertex);

void backtrack_transfer_at(const Tree_type &spcies_tree,
                           const Tree_type &gene_tree,
                           const vector<vid_t> &sigma,
                           multi_array<Below, 2> &below,
                           const multi_array<Outside, 2> &outside,
                           Cost_type duplication_cost,
                           Cost_type transfer_cost,
                           vid_t gene_vertex_at,
                           vid_t gene_vertex_outside,
                           vid_t species_vertex);




/*
 * Function and member definitions.
 */

// Scenario constructor
Scenario::Scenario(unsigned num_vertices)
    : duplications(num_vertices), transfer_edges(num_vertices)
{
    ; // Do nothing.
}

// Scenario cost function
Cost_type
Scenario::cost(Cost_type duplication_cost, Cost_type transfer_cost)
{
    return duplications.count() * duplication_cost +
        transfer_edges.count() * transfer_cost;
}


// Below constructor
Below::Below() : cost(COST_INF), events(N_EVENTS, false),
                 scenarios_below_computed(false), scenarios_at_computed(false)
{
    ; // Do nothing.
}


// Outside constructor
Outside::Outside() : cost(COST_INF)
{
    ; // Do nothing
}

/*
 * Template function definitions. (must be here since definition must
 * be available at point of usage)
 */

template<typename Initer1, typename Initer2, typename Outiter>
void combine_scenarios(Initer1 begin1, Initer1 end1,
                       Initer2 begin2, Initer2 end2,
                       Outiter out)
{
    if (begin1 == end1 || begin2 == end2)
        return;

    unsigned size = begin1->duplications.size();

    for (Initer1 i = begin1; i != end1; ++i)
        {
            for (Initer2 j = begin2; j != end2; ++j)
                {
                    Scenario sc(size);
                    sc.duplications = i->duplications | j->duplications;
                    sc.transfer_edges = i->transfer_edges | j->transfer_edges;
                    *out = sc; out++;
                }
        }
}



void backtrack_scenarios_below(const Tree_type &S,
                               const Tree_type &G,
                               const vector<vid_t> &sigma,
                               multi_array<Below, 2> &below,
                               const multi_array<Outside, 2> &outside,
                               Cost_type duplication_cost,
                               Cost_type transfer_cost,
                               vid_t u, 
                               vid_t x)
{
    if (below[u][x].scenarios_below_computed)
        return;

    vector<Scenario> &scenarios = below[u][x].scenarios_below;
    
    /* If u is a leaf. */
    if (G.is_leaf(u))
        {
            if (S.descendant(sigma[u], x))
                {
                    scenarios.push_back(Scenario(G.size()));
                }
            below[u][x].scenarios_below_computed = true;
            return;
        }
    
    /* If u is not a leaf. */
    backtrack_scenarios_at(S, G, sigma, below, outside,
                           duplication_cost, transfer_cost,
                           u, x);
    copy(below[u][x].scenarios_at.begin(), below[u][x].scenarios_at.end(),
         back_inserter(scenarios));

    if (!S.is_leaf(x))
        {
            vid_t y = S.left(x);
            vid_t z = S.right(x);

            if (below[u][x].events[Below::B_LEFT] == true)
                {
                    backtrack_scenarios_below(S, G, sigma, below, outside,
                                              duplication_cost, transfer_cost,
                                              u, y);
                    copy(below[u][y].scenarios_below.begin(),
                         below[u][y].scenarios_below.end(),
                         back_inserter(scenarios));
                }
            if (below[u][x].events[Below::B_RIGHT] == true)
                {
                    backtrack_scenarios_below(S, G, sigma, below, outside,
                                              duplication_cost, transfer_cost,
                                              u, z);
                    copy(below[u][z].scenarios_below.begin(),
                         below[u][z].scenarios_below.end(),
                         back_inserter(scenarios));
                }
        }
    below[u][x].scenarios_below_computed = true;
}

void backtrack_scenarios_at(const Tree_type &S,
                            const Tree_type &G,
                            const vector<vid_t> &sigma,
                            multi_array<Below, 2> &below,
                            const multi_array<Outside, 2> &outside,
                            Cost_type duplication_cost,
                            Cost_type transfer_cost,
                            vid_t u,
                            vid_t x)
{
    if (below[u][x].scenarios_at_computed)
        return;
    
    typedef vector<Scenario>::iterator Iter;
    
    vector<Scenario> &scenarios = below[u][x].scenarios_at;

    /* If u is a leaf. */
    if (G.is_leaf(u))
        {
            if (sigma[u] == x)
                {
                    scenarios.push_back(Scenario(G.size()));
                }
            below[u][x].scenarios_at_computed = true;
            return;
        }

    /* If u is not a leaf. */
    vid_t u_left = G.left(u);
    vid_t u_right = G.right(u);

    /*
     * In case of speciation being optimal
     */
    if (below[u][x].events[Below::S])
        {
            backtrack_speciation_at(S, G, sigma, below, outside,
                                    duplication_cost, transfer_cost,
                                    u_left, u_right, x);
        }

    if (below[u][x].events[Below::S_REV])
        {
            backtrack_speciation_at(S, G, sigma, below, outside,
                                    duplication_cost, transfer_cost,
                                    u_right, u_left, x);
        }
    /*
     * In case of transfer being optimal.
     */
    if (below[u][x].events[Below::T_LEFT])
        {
            backtrack_transfer_at(S, G, sigma, below, outside,
                                  duplication_cost, transfer_cost,
                                  u_right, u_left, x);
        }
    if (below[u][x].events[Below::T_RIGHT])
        {
            backtrack_transfer_at(S, G, sigma, below, outside,
                                  duplication_cost, transfer_cost,
                                  u_left, u_right, x);
        }
        
    /*
     * In case of duplication being optimal.
     */
    if (below[u][x].events[Below::D])
        {
            backtrack_duplication_at(S, G, sigma, below, outside,
                                     duplication_cost, transfer_cost,
                                     u_left, u_right, x);
        }
    below[u][x].scenarios_at_computed = true;
}

void backtrack_speciation_at(const Tree_type &S,
                             const Tree_type &G,
                             const vector<vid_t> &sigma,
                             multi_array<Below, 2> &below,
                             const multi_array<Outside, 2> &outside,
                             Cost_type duplication_cost,
                             Cost_type transfer_cost,
                             vid_t left,
                             vid_t right,
                             vid_t x)
{
    /*
     * We are sure that scenarios must be computed here. The check is
     * made in backtrack_scenarios_at() which calls this function.
     */
    if (S.is_leaf(x))
        return;
    
    vid_t x_left = S.left(x);
    vid_t x_right = S.right(x);

    backtrack_scenarios_below(S, G, sigma, below, outside,
                              duplication_cost, transfer_cost,
                              left, x_left);
    backtrack_scenarios_below(S, G, sigma, below, outside,
                              duplication_cost, transfer_cost, 
                              right, x_right);
            
    vector<Scenario> &left_scenarios = below[left][x_left].scenarios_below;
    vector<Scenario> &right_scenarios = below[right][x_right].scenarios_below;

    combine_scenarios(left_scenarios.begin(), left_scenarios.end(),
                      right_scenarios.begin(), right_scenarios.end(),
                      back_inserter(below[G.parent(left)][x].scenarios_at));
}

void backtrack_duplication_at(const Tree_type &S,
                              const Tree_type &G,
                              const vector<vid_t> &sigma,
                              multi_array<Below, 2> &below,
                              const multi_array<Outside, 2> &outside,
                              Cost_type duplication_cost,
                              Cost_type transfer_cost,
                              vid_t v,
                              vid_t w,
                              vid_t x)
{
    /*
     * We are sure that scenarios must be computed here. The check is
     * made in backtrack_scenarios_at() which calls this function.
     */
    vid_t u = G.parent(v);

    vector<Scenario> all;

    backtrack_scenarios_at(S, G, sigma, below, outside,
                           duplication_cost, transfer_cost,
                           v, x);
    backtrack_scenarios_at(S, G, sigma, below, outside,
                           duplication_cost, transfer_cost,
                           w, x);
    combine_scenarios(below[v][x].scenarios_at.begin(),
                      below[v][x].scenarios_at.end(),
                      below[w][x].scenarios_at.begin(),
                      below[w][x].scenarios_at.end(),
                      back_inserter(all));

    if (!S.is_leaf(x))
        {
            vid_t y = S.left(x);
            vid_t z = S.right(x);

            backtrack_scenarios_below(S, G, sigma, below, outside,
                                      duplication_cost, transfer_cost,
                                      w, y);
            backtrack_scenarios_below(S, G, sigma, below, outside,
                                      duplication_cost, transfer_cost,
                                      w, z);
            backtrack_scenarios_at(S, G, sigma, below, outside,
                                   duplication_cost, transfer_cost,
                                   v, x);
            combine_scenarios(below[w][y].scenarios_below.begin(),
                              below[w][y].scenarios_below.end(),
                              below[v][x].scenarios_at.begin(),
                              below[v][x].scenarios_at.end(),
                              back_inserter(all));
            combine_scenarios(below[w][z].scenarios_below.begin(),
                              below[w][z].scenarios_below.end(),
                              below[v][x].scenarios_at.begin(),
                              below[v][x].scenarios_at.end(),
                              back_inserter(all));

            backtrack_scenarios_below(S, G, sigma, below, outside,
                                      duplication_cost, transfer_cost,
                                      v, y);
            backtrack_scenarios_below(S, G, sigma, below, outside,
                                      duplication_cost, transfer_cost,
                                      v, z);
            backtrack_scenarios_at(S, G, sigma, below, outside,
                                   duplication_cost, transfer_cost,
                                   w, x);
            combine_scenarios(below[v][y].scenarios_below.begin(),
                              below[v][y].scenarios_below.end(),
                              below[w][x].scenarios_at.begin(),
                              below[w][x].scenarios_at.end(),
                              back_inserter(all));
            combine_scenarios(below[v][z].scenarios_below.begin(),
                              below[v][z].scenarios_below.end(),
                              below[w][x].scenarios_at.begin(),
                              below[w][x].scenarios_at.end(),
                              back_inserter(all));
        }
    for (unsigned i = 0; i < all.size(); ++i)
        {
            all[i].duplications.set(u);
            if (all[i].cost(duplication_cost, transfer_cost) == 
                below[u][x].cost)
                {
                    below[u][x].scenarios_at.push_back(all[i]);
                }
        }
}

void backtrack_transfer_at(const Tree_type &S,
                           const Tree_type &G,
                           const vector<vid_t> &sigma,
                           multi_array<Below, 2> &below,
                           const multi_array<Outside, 2> &outside,
                           Cost_type duplication_cost,
                           Cost_type transfer_cost,
                           vid_t u_at,
                           vid_t u_outside,
                           vid_t x)
{
    /*
     * We are sure that scenarios must be computed here. The check is
     * made in backtrack_scenarios_at() which calls this function.
     */
    vector<Scenario> all;
    
    vector<vid_t>::const_iterator i =
        outside[u_outside][x].species_vertices.begin();
    vector<vid_t>::const_iterator e =
        outside[u_outside][x].species_vertices.end();
    for ( ; i != e; ++i)
        {
            backtrack_scenarios_below(S, G, sigma, below, outside,
                                      duplication_cost, transfer_cost,
                                      u_outside, *i);
        }
    backtrack_scenarios_at(S, G, sigma, below, outside,
                           duplication_cost, transfer_cost,
                           u_at, x);

    for (i = outside[u_outside][x].species_vertices.begin(); i != e; ++i)
        {
            combine_scenarios(below[u_outside][*i].scenarios_below.begin(),
                              below[u_outside][*i].scenarios_below.end(),
                              below[u_at][x].scenarios_at.begin(),
                              below[u_at][x].scenarios_at.end(),
                              back_inserter(all));
        }
    for (unsigned j = 0; j < all.size(); ++j)
        {
            all[j].transfer_edges.set(u_outside);
            below[G.parent(u_at)][x].scenarios_at.push_back(all[j]);
        }
}

