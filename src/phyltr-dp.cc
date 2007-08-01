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
 */

#include "common.hh"
#include "utils/binary-tree.hh"

extern "C" {
#include <NHparser.h>
}
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/concept_check.hpp>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <vector>
#include <iterator>
#include <fstream>
#include <limits>

//============================================================================
//             Namespace declarations and using directives.
//============================================================================

using namespace  std;
using            boost::multi_array;
using            boost::dynamic_bitset;

//============================================================================
//                   Typedefs and class declarations
//============================================================================

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
 */
struct Scenario {
    explicit Scenario(unsigned num_vertices);
    Cost_type cost() const;

    dynamic_bitset<> duplications;
    dynamic_bitset<> transfer_edges;
    
    bool operator<(const Scenario &sc) const;
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
 * u is mapped to exactly x.
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
    bool           cost_only;

    Program_input() : species_tree(), gene_tree(), sigma(),
                      duplication_cost(0), transfer_cost(0),
                      cost_only(false)
    {}
};
//============================================================================
//                   Global Constants and variables.
//============================================================================

const string PROG_NAME = "phyltr-dp";
const string USAGE = 
    "Usage: " + PROG_NAME + " [OPTION]... SPECIES_TREE_FILE GENE_TREE_FILE "
    "GENEMAP_FILE";
const Cost_type COST_INF = numeric_limits<Cost_type>::infinity();

Program_input g_program_input = Program_input();

/* Variable holding the command line options. */
boost::program_options::variables_map g_program_options;
/* The DP-algorithm matrices. */
multi_array<Below, 2> g_below;
multi_array<Outside, 2> g_outside;

//=============================================================================
//                        Function declarations
//=============================================================================

/*
 * backtrack_scenarios_below()
 *
 * The function computes each unique (in the sense described above in
 * the description of the Scenario class) most parsimonious partial
 * DTL-scenario (a Scenario object) such that 'gene_vertex' is placed
 * at or below 'species_vertex', and inserts them into
 * g_below[gene_vertex][species_vertex].scenarios_below.
 */
void backtrack_scenarios_below(vid_t gene_vertex, 
                               vid_t species_vertex);

/*
 * dp_algorithm()
 *
 * The function performing the dynamic programming algorithm. The
 * result is stored in the global matrices 'g_below' and 'g_outside'.
 */
void dp_algorithm();

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

    /*
     * Parse the command line options.
     */

    /* Declare options that are described when --help is given. */
    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    try
        {
            visible_opts.add_options()
                ("help", "display this help and exit")
                ("transfer-cost,t",
                 po::value<Cost_type>(&g_program_input.transfer_cost)
                 ->default_value(1.0),
                 "Cost of transfer events")
                ("duplication-cost,d", 
                 po::value<Cost_type>(&g_program_input.duplication_cost)
                 ->default_value(1.0),
                 "Cost of duplication events")
                ("cost-only,c",
                 po::bool_switch(&g_program_input.cost_only)
                 ->default_value(false),
                 "If set, only the optimum cost is printed")
                ;
            
            /* Declare positional options. */
            hidden_opts.add_options()
                ("species-tree-file",po::value<string>())
                ("gene-tree-file", po::value<string>())
                ("map-file", po::value<string>())
                ;
            
            po::positional_options_description positional_options;
            positional_options.add("species-tree-file", 1);
            positional_options.add("gene-tree-file", 1);
            positional_options.add("map-file", 1);

            /* Gather all options into a single options_description. */
            all_options.add(visible_opts).add(hidden_opts);
            
            /* Parse the arguments. */
            po::command_line_parser parser(argc, argv);
            parser.options(all_options);
            parser.positional(positional_options);
            po::store(parser.run(), g_program_options);
            po::notify(g_program_options);

        }
    catch (exception &e)
        {
            cerr << PROG_NAME << ": " << e.what() << '\n'
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* Show help message if --help was given. */
    if (g_program_options.count("help"))
        {
            cout << USAGE << '\n'
                 << visible_opts << '\n';
            exit(EXIT_SUCCESS);
        }
    
    /* Check that all required positional arguments are given. */
    if (g_program_options.count("species-tree-file") == 0 
        || g_program_options.count("gene-tree-file") == 0 
        || g_program_options.count("map-file") == 0)
        {
            cerr << PROG_NAME << ": Too few arguments.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }
    
    /*
     * Read the species- and gene-tree files.
     */
    errno = 0;
    string species_tree_filename =
        g_program_options["species-tree-file"].as<string>();
    FILE *species_tree_file = fopen(species_tree_filename.c_str(), "r");
    if (errno)
        {
            cerr << PROG_NAME << ": "
                 << "failed to open species tree file '"
                 << species_tree_filename << "'\n";
            perror(PROG_NAME.c_str());
            exit(EXIT_FAILURE);
        }

    string gene_tree_filename =
        g_program_options["gene-tree-file"].as<string>();
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
            cerr << PROG_NAME << ": " << e.what() << '\n';
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
            cerr << PROG_NAME << ": " << e.what() << '\n';
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
            cerr << PROG_NAME << ": transfer cost is too big";
            exit(EXIT_FAILURE);
        }
    if (g_program_input.duplication_cost >= max_allowed)
        {
            cerr << PROG_NAME << ": duplication cost is too big";
            exit(EXIT_FAILURE);
        }
    
    
    /* Create the below and outside matrices. */
      
    g_below.resize(boost::extents[g_program_input.gene_tree.size() + 1]
                                 [g_program_input.species_tree.size() + 1]);
    g_outside.resize(boost::extents[g_program_input.gene_tree.size() + 1]
                                   [g_program_input.species_tree.size() + 1]);
    dp_algorithm();

    /*
     * If the only_cost flag was given, Output the cost.and exit.
     */
    if (g_program_options["cost-only"].as<bool>())
        {
            cout << g_below[0][0].cost << '\n';
            return EXIT_SUCCESS;
        }

    /*
     * Backtrack and output all scenarios. We should do this
     * xml-style.
     */
    backtrack_scenarios_below(0, 0);

    vector<Scenario> &scenarios = g_below[0][0].scenarios_below;

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

    BOOST_FOREACH(Scenario &sc, scenarios)
        {
            dynamic_bitset<> new_duplications(sc.duplications.size());
            dynamic_bitset<> new_transfers(sc.transfer_edges.size());
            for (unsigned i = sc.duplications.find_first();
                 i != sc.duplications.npos;
                 i = sc.duplications.find_next(i))
                {
                    new_duplications.set(gene_tree_numbering[i]);
                }
            for (unsigned i = sc.transfer_edges.find_first();
                 i != sc.transfer_edges.npos;
                 i = sc.transfer_edges.find_next(i))
                {
                    new_transfers.set(gene_tree_numbering[i]);
                }
            sc.duplications = new_duplications;
            sc.transfer_edges = new_transfers;
        }

    sort(scenarios.begin(), scenarios.end());

    cout << "Number of scenarios: " << scenarios.size() << "\n\n";
    BOOST_FOREACH(const Scenario &sc, scenarios)
        {
            cout << "Transfer edges:\t";
            for (unsigned j = sc.transfer_edges.find_first();
                 j != sc.transfer_edges.npos; j = sc.transfer_edges.find_next(j))
                {
                    cout << j << " ";
                }
            cout << '\n';
            cout << "Duplications:\t";
            for (unsigned j = sc.duplications.find_first();
                 j != sc.duplications.npos; j = sc.duplications.find_next(j))
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
 * g_below[gene_vertex][species_vertex].scenarios_at.
 */
void backtrack_scenarios_at(vid_t gene_vertex,
                            vid_t species_vertex);

/*
 * backtrack_[speciation,duplication,transfer]_at()
 *
 * These are helper functions called by backtrack_scenarios_at.
 */

void backtrack_speciation_at(vid_t gene_vertex_left,
                             vid_t gene_vertex_right,
                             vid_t species_vertex);

void backtrack_duplication_at(vid_t gene_vertex_1,
                              vid_t gene_vertex_2,
                              vid_t species_vertex);

void backtrack_transfer_at(vid_t gene_vertex_at,
                           vid_t gene_vertex_outside,
                           vid_t species_vertex);

//=============================================================================
//                   Function and member definitions.
//=============================================================================

void
dp_algorithm()
{
    const Tree_type &G = g_program_input.gene_tree;
    const Tree_type &S = g_program_input.species_tree;
    const vector<vid_t> &sigma = g_program_input.sigma;
    const Cost_type transfer_cost = g_program_input.transfer_cost;
    const Cost_type duplication_cost = g_program_input.duplication_cost;

    for (vid_t u = G.postorder_begin();
         u != G.NONE;
         u = G.postorder_next(u))
        {
            /* First compute g_below[u][*]. */
            for (vid_t x = S.postorder_begin();
                 x != S.NONE;
                 x = S.postorder_next(x))
                {
                    /* If u is a leaf use sigma. */
                    if (G.is_leaf(u))
                        {
                            if (S.descendant(sigma[u], x))
                                {
                                    g_below[u][x].cost = 0;
                                }
                            continue;
                        }

                    /* If u is not a leaf compute costs. */
                    vector<Cost_type> costs(Below::N_EVENTS, COST_INF);

                    vid_t v = G.left(u);
                    vid_t w = G.right(u);
                    
                    costs[Below::D] = duplication_cost
                        + g_below[v][x].cost + g_below[w][x].cost;
                    costs[Below::T_LEFT] = transfer_cost
                        + g_outside[v][x].cost + g_below[w][x].cost;
                    costs[Below::T_RIGHT] = transfer_cost
                        + g_outside[w][x].cost + g_below[v][x].cost;
                    if (!S.is_leaf(x))
                        {
                            vid_t y = S.left(x);
                            vid_t z = S.right(x);
                            costs[Below::S] = 
                                g_below[v][y].cost + g_below[w][z].cost;
                            costs[Below::S_REV]
                                = g_below[v][z].cost + g_below[w][y].cost;
                            costs[Below::B_LEFT] = g_below[u][y].cost;
                            costs[Below::B_RIGHT] = g_below[u][z].cost;
                        }
                    Cost_type min_cost = 
                        *min_element(costs.begin(), costs.end());
                    g_below[u][x].cost = min_cost;
                    for (unsigned i = 0; i < Below::N_EVENTS; ++i)
                        {
                            if (costs[i] == min_cost)
                                g_below[u][x].events[i] = true;
                        }
                }

            /* Now compute g_outside[u][*]. */
            
            for (vid_t x = S.preorder_begin();
                 x != S.NONE;
                 x = S.preorder_next(x))
                {
                    if (S.is_leaf(x))
                        continue;
                    
                    vid_t y = S.left(x);
                    vid_t z = S.right(x);
                    
                    Cost_type outside_u_x = g_outside[u][x].cost;
                    Cost_type below_u_y = g_below[u][y].cost;
                    Cost_type below_u_z = g_below[u][z].cost;
                    g_outside[u][y].cost = min(outside_u_x, below_u_z);
                    g_outside[u][z].cost = min(outside_u_x, below_u_y);

                    if (g_outside[u][y].cost == outside_u_x)
                        {
                            copy(g_outside[u][x].species_vertices.begin(),
                                 g_outside[u][x].species_vertices.end(),
                                 back_inserter(g_outside[u][y].species_vertices));
                        }
                    if (g_outside[u][y].cost == below_u_z)
                        {
                            g_outside[u][y].species_vertices.push_back(z);
                        }

                    if (g_outside[u][z].cost == outside_u_x)
                        {
                            copy(g_outside[u][x].species_vertices.begin(),
                                 g_outside[u][x].species_vertices.end(),
                                 back_inserter(g_outside[u][z].species_vertices));
                        }
                    if (g_outside[u][z].cost == below_u_y)
                        {
                            g_outside[u][z].species_vertices.push_back(y);
                        }
                }
        }
}


template<typename Initer1, typename Initer2, typename Outiter>
void
combine_scenarios(Initer1 begin1, Initer1 end1,
                  Initer2 begin2, Initer2 end2,
                  Outiter out)
{
    boost::function_requires< boost::InputIteratorConcept<Initer1> >();
    boost::function_requires< boost::InputIteratorConcept<Initer2> >();
    boost::function_requires< boost::OutputIteratorConcept<Outiter, Scenario> >();

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
                    *out++ = sc;
                }
        }
}

// Scenario constructor
Scenario::Scenario(unsigned num_vertices)
    : duplications(num_vertices), transfer_edges(num_vertices)
{
    ; // Do nothing.
}

// Scenario cost function
Cost_type
Scenario::cost() const
{
    return duplications.count() * g_program_input.duplication_cost +
        transfer_edges.count() * g_program_input.transfer_cost;
}

bool
Scenario::operator<(const Scenario &sc) const
{
    if (cost() < sc.cost())
        return true;
    if (cost() > sc.cost())
        return false;
    if (transfer_edges.count() < sc.transfer_edges.count())
        return true;
    if (transfer_edges.count() > sc.transfer_edges.count())
        return false;
    if (transfer_edges < sc.transfer_edges)
        return true;
    if (transfer_edges > sc.transfer_edges)
        return false;
    return duplications < sc.duplications;
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




void
backtrack_scenarios_below(vid_t u, vid_t x)
{
    const Tree_type &G = g_program_input.gene_tree;
    const Tree_type &S = g_program_input.species_tree;
    const vector<vid_t> &sigma = g_program_input.sigma;

    if (g_below[u][x].scenarios_below_computed)
        return;

    vector<Scenario> &scenarios = g_below[u][x].scenarios_below;
    
    /* If u is a leaf. */
    if (G.is_leaf(u))
        {
            if (S.descendant(sigma[u], x))
                {
                    scenarios.push_back(Scenario(G.size()));
                }
            g_below[u][x].scenarios_below_computed = true;
            return;
        }
    
    /* If u is not a leaf. */
    backtrack_scenarios_at(u, x);
    copy(g_below[u][x].scenarios_at.begin(), g_below[u][x].scenarios_at.end(),
         back_inserter(scenarios));

    if (!S.is_leaf(x))
        {
            vid_t y = S.left(x);
            vid_t z = S.right(x);

            if (g_below[u][x].events[Below::B_LEFT] == true)
                {
                    backtrack_scenarios_below(u, y);
                    copy(g_below[u][y].scenarios_below.begin(),
                         g_below[u][y].scenarios_below.end(),
                         back_inserter(scenarios));
                }
            if (g_below[u][x].events[Below::B_RIGHT] == true)
                {
                    backtrack_scenarios_below(u, z);
                    copy(g_below[u][z].scenarios_below.begin(),
                         g_below[u][z].scenarios_below.end(),
                         back_inserter(scenarios));
                }
        }
    g_below[u][x].scenarios_below_computed = true;
}

void
backtrack_scenarios_at(vid_t u, vid_t x)
{
    const Tree_type &G = g_program_input.gene_tree;
    const vector<vid_t> &sigma = g_program_input.sigma;

    if (g_below[u][x].scenarios_at_computed)
        return;
    
    typedef vector<Scenario>::iterator Iter;
    
    vector<Scenario> &scenarios = g_below[u][x].scenarios_at;

    /* If u is a leaf. */
    if (G.is_leaf(u))
        {
            if (sigma[u] == x)
                {
                    scenarios.push_back(Scenario(G.size()));
                }
            g_below[u][x].scenarios_at_computed = true;
            return;
        }

    /* If u is not a leaf. */
    vid_t u_left = G.left(u);
    vid_t u_right = G.right(u);

    /*
     * In case of speciation being optimal
     */
    if (g_below[u][x].events[Below::S])
        {
            backtrack_speciation_at(u_left, u_right, x);
        }

    if (g_below[u][x].events[Below::S_REV])
        {
            backtrack_speciation_at(u_right, u_left, x);
        }
    /*
     * In case of transfer being optimal.
     */
    if (g_below[u][x].events[Below::T_LEFT])
        {
            backtrack_transfer_at(u_right, u_left, x);
        }
    if (g_below[u][x].events[Below::T_RIGHT])
        {
            backtrack_transfer_at(u_left, u_right, x);
        }
        
    /*
     * In case of duplication being optimal.
     */
    if (g_below[u][x].events[Below::D])
        {
            backtrack_duplication_at(u_left, u_right, x);
        }
    g_below[u][x].scenarios_at_computed = true;
}

void
backtrack_speciation_at(vid_t left, vid_t right, vid_t x)
{
    const Tree_type &G = g_program_input.gene_tree;
    const Tree_type &S = g_program_input.species_tree;

    /*
     * We are sure that scenarios must be computed here. The check is
     * made in backtrack_scenarios_at() which calls this function.
     */
    if (S.is_leaf(x))
        return;
    
    vid_t x_left = S.left(x);
    vid_t x_right = S.right(x);

    backtrack_scenarios_below(left, x_left);
    backtrack_scenarios_below(right, x_right);
            
    vector<Scenario> &left_scenarios = g_below[left][x_left].scenarios_below;
    vector<Scenario> &right_scenarios = g_below[right][x_right].scenarios_below;

    combine_scenarios(left_scenarios.begin(), left_scenarios.end(),
                      right_scenarios.begin(), right_scenarios.end(),
                      back_inserter(g_below[G.parent(left)][x].scenarios_at));
}

void
backtrack_duplication_at(vid_t v, vid_t w, vid_t x)
{
    const Tree_type &G = g_program_input.gene_tree;
    const Tree_type &S = g_program_input.species_tree;

    /*
     * We are sure that scenarios must be computed here. The check is
     * made in backtrack_scenarios_at() which calls this function.
     */
    vid_t u = G.parent(v);

    vector<Scenario> all;

    backtrack_scenarios_at(v, x);
    backtrack_scenarios_at(w, x);
    combine_scenarios(g_below[v][x].scenarios_at.begin(),
                      g_below[v][x].scenarios_at.end(),
                      g_below[w][x].scenarios_at.begin(),
                      g_below[w][x].scenarios_at.end(),
                      back_inserter(all));

    if (!S.is_leaf(x))
        {
            vid_t y = S.left(x);
            vid_t z = S.right(x);

            backtrack_scenarios_below(w, y);
            backtrack_scenarios_below(w, z);
            backtrack_scenarios_at(v, x);
            combine_scenarios(g_below[w][y].scenarios_below.begin(),
                              g_below[w][y].scenarios_below.end(),
                              g_below[v][x].scenarios_at.begin(),
                              g_below[v][x].scenarios_at.end(),
                              back_inserter(all));
            combine_scenarios(g_below[w][z].scenarios_below.begin(),
                              g_below[w][z].scenarios_below.end(),
                              g_below[v][x].scenarios_at.begin(),
                              g_below[v][x].scenarios_at.end(),
                              back_inserter(all));

            backtrack_scenarios_below(v, y);
            backtrack_scenarios_below(v, z);
            backtrack_scenarios_at(w, x);
            combine_scenarios(g_below[v][y].scenarios_below.begin(),
                              g_below[v][y].scenarios_below.end(),
                              g_below[w][x].scenarios_at.begin(),
                              g_below[w][x].scenarios_at.end(),
                              back_inserter(all));
            combine_scenarios(g_below[v][z].scenarios_below.begin(),
                              g_below[v][z].scenarios_below.end(),
                              g_below[w][x].scenarios_at.begin(),
                              g_below[w][x].scenarios_at.end(),
                              back_inserter(all));
        }
    for (unsigned i = 0; i < all.size(); ++i)
        {
            all[i].duplications.set(u);
            if (all[i].cost() == g_below[u][x].cost)
                {
                    g_below[u][x].scenarios_at.push_back(all[i]);
                }
        }
}

void
backtrack_transfer_at(vid_t u_at, vid_t u_outside, vid_t x)
{
    const Tree_type &G = g_program_input.gene_tree;

    /*
     * We are sure that scenarios must be computed here. The check is
     * made in backtrack_scenarios_at() which calls this function.
     */
    vector<Scenario> all;
    
    vector<vid_t>::const_iterator i =
        g_outside[u_outside][x].species_vertices.begin();
    vector<vid_t>::const_iterator e =
        g_outside[u_outside][x].species_vertices.end();
    for ( ; i != e; ++i)
        {
            backtrack_scenarios_below(u_outside, *i);
        }
    backtrack_scenarios_at(u_at, x);

    for (i = g_outside[u_outside][x].species_vertices.begin(); i != e; ++i)
        {
            combine_scenarios(g_below[u_outside][*i].scenarios_below.begin(),
                              g_below[u_outside][*i].scenarios_below.end(),
                              g_below[u_at][x].scenarios_at.begin(),
                              g_below[u_at][x].scenarios_at.end(),
                              back_inserter(all));
        }
    for (unsigned j = 0; j < all.size(); ++j)
        {
            all[j].transfer_edges.set(u_outside);
            g_below[G.parent(u_at)][x].scenarios_at.push_back(all[j]);
        }
}
