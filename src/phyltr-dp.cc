/*
 * Copyright (C) 2010, 2011 Ali Tofigh
 *
 * This file is part of PhylTr, a package for phylogenetic analysis
 * using duplications and transfers.
 *
 * PhylTr is released under the terms of the license contained in the
 * file LICENSE.
 */

//*****************************************************************************
// phyltr-dp.cc
//
// Implementation of the dynamic programming algorithm for finding the
// minimum cost of any DTL-scenario reconciling a species tree and a
// gene tree. A backtracking algorithm is also implemented to find all
// optimal DTL-scenarios.
//*****************************************************************************

#include "common.hh"
#include "utils/binary-tree.hh"

#include <NHparser/NHparser.h>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>

#include <exception>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <vector>
#include <bitset>
#include <iterator>
#include <limits>
#include <utility>



namespace           po = boost::program_options;

using namespace     std;
using               boost::multi_array;
using               boost::dynamic_bitset;
using               boost::shared_ptr;



typedef Binary_tree<double>     tree_type;
typedef tree_type::vid_t        vid_t;
typedef float                   cost_type; // Instead of double to save memory.



const string PROGRAM_NAME = "phyltr-dp";
const string USAGE =
    "Usage: " + PROGRAM_NAME + " [OPTION]... SPECIES-TREE-FILE GENE-TREE-FILE "
    "MAP-FILE";



//*****************************************************************************
// class Scenario
//
// A Scenario object is used to keep track of the duplications and
// transfer edges of solutions during backtracking. A Scenario object
// may represent a partial solution. Throughout the program, an edge
// is represented by the vertex at its head, i.e., the vertex farthest
// away from the root.
//
// operator<() is used to sort the scenarios for printing. It sorts
//      first on the number of transfers, then on the number of losses,
//      and lastly according to lexicographic order.
//
// scenario_less_lex() is used to sort scenarios lexicographically. It
//      is used in the dated version of the algorithm when removing
//      duplicate scenarios.
//
// operator==() is needed when removing duplicated scenarios using the
//      std::unique() function.
//
// operator<<() is used to print Scenarios for the user. The syntax is
//      currently not very user friendly.
//*****************************************************************************
class Scenario {
public:
    dynamic_bitset<> duplications;
    dynamic_bitset<> transfer_edges;

    Scenario(unsigned size);
};

bool operator<(const Scenario &, const Scenario &);
bool scenario_less_lex(const Scenario&, const Scenario &);
bool operator==(const Scenario &sc1, const Scenario &sc2);
ostream &operator<<(ostream &, const Scenario &);



//*****************************************************************************
// class BacktrackElement
//
// Used to hold the information needed for backtracking after the
// dynamic programming algorithm has run. For gene tree vertex u and
// species tree vertex x, the members have the following meaning:
//
// The enums denote biological events that can be associated with a
// gene tree vertex and a species tree vertex. The below_events member
// is used to keep track of the events that led to the optimal cost of
// placing u at a descendant of x (possibly x itself):
//
// S
//      A speciation with the left child of u below left child of x,
//      and right child of u below right child of x.
//
// S_REV
//      A speciation with left child of u below right child of x, and
//      right child of u below left child of x.
//
// D
//      A duplication
//
// T_LEFT
//      The left edge of u is a transfer edge.
//
// T_RIGHT
//      The right edge of u is a transfer edge.
//
// BELOW_LEFT
//      The optimum cost can be achieved by placing u at a descendant
//      of the left child of x.
//
// BELOW_RIGHT
//      The optimum cost can be achieved by placing u at a descendant
//      of the right child of x.
//
// below_events:
//      For an Event e, below_events[e] is set iff the event
//      represented by e led to the optimal cost of placing u at a
//      descendant of x (possibly x itself). This member is set during
//      the dynamic programming algorithm.
//
// outside_sibling:
//      The vid_t of the sibling of x, if placing u below the sibling
//      gives the optimal cost of placing u outside x, and
//      tree_type::NONE otherwise. This member is set during
//      the dynamic programming algorithm.
//
// outside_ancestor:
//      The vid_t of the nearest proper ancestor of x such that
//      placing u below the sibling of the ancestor gives the optimal
//      cost of placing u outside x, and tree_type::NONE if no such
//      ancestor exists. This member is set during the dynamic
//      programming algorithm.
//
// min_transfers:
//      The minimum number of transfers required for obtaining the
//      optimum cost of placing u at a descendant of x. This is used
//      when the --minimum-transfers flag has been set.
//
// scenarios_below_needed:
//      Set to true iff we need to compute all the scenarios
//      corresponding to placing u at a descendant of x.
//
// scenarios_at_needed:
//      Set to true iff we need to compute all the scenarios
//      corresponding to placing u _at_ x.
//
// below_placements:
//      The set of descendants of x at which u can be placed to optain
//      the optimum cost. If x is itself among this set, then it must
//      always be the first element in the vector.
//
// scenarios_at:
//      The set of scenarios corresponding to placing u _at_ x.
//*****************************************************************************
class BacktrackElement {
public:
    enum Event {S, S_REV, D, T_LEFT, T_RIGHT, BELOW_LEFT, BELOW_RIGHT, N_EVENTS};

    bitset<N_EVENTS>    below_events;
    vid_t               outside_sibling;
    vid_t               outside_ancestor;
    unsigned            min_transfers;
    bool                scenarios_below_needed;
    bool                scenarios_at_needed;
    vector<vid_t>       below_placements;
    vector<Scenario>    scenarios_at;

    BacktrackElement();
};

//*****************************************************************************
// class BacktrackElementDated
//
// Used to hold the information needed for backtracking after the
// dynamic programming algorithm for the dated species trees has
// run. For gene tree vertex u and species tree vertex x, the members
// have the following meaning:
//
// The enums have the same meaning as BacktrackElement above with the
// addition of:
//
// BELOW:
//      the optimum cost can be achieved by placing u at x on a level
//      closer to the leaves.
//
// Note that BELOW_LEFT and BELOW_RIGHT are only used if x is the
// speciation event at the current level.
//
// The other fields have the following meaning:
//
// below_events
//      Same as in BacktrackElement.
//
// scenarios_needed
//      When 'true', indicates that the backtrack algorithm requires the
//      scenarios for this position (i.e., u, x, level)
//
// scenarios_at
//      Holds the scenarios at for the current position. The algorithm
//      removes duplicated scenarios at each step.
//*****************************************************************************
class BacktrackElementDated {
public:
    enum Event {S, S_REV, D, T_LEFT, T_RIGHT, BELOW_LEFT, BELOW_RIGHT, BELOW, N_EVENTS};

    bitset<N_EVENTS>    below_events;
    bool                scenarios_needed;
    vector<Scenario>    scenarios_at;

    BacktrackElementDated();
};

//*****************************************************************************
// global variables
//
// g_options
//      Contains the command line options passed to the
//      program. Should not really be needed as the options are all
//      copied to g_input.
//
// g_input
//      Holds the input to the program. This includes both flags,
//      filenames, and the data contained in the
//      files. gene_tree_numbering is used to number the gene tree
//      vertices for output, i.e., gene_tree_numbering[u] is the
//      number of vertex u when printing solutions.
//
// g_below
// g_outside
//      The DP matrices. For a gene tree vertex u and a species
//      tree vertex x, g_below[u][x] is the minimum cost of placing u
//      at a descendant of x (possibly at x itself), and
//      g_outside[u][x] is the minimum cost of placing u at a vertex
//      incomaparable to x.
//
// g_backtrack_matrix
// g_backtrack_matrix_dated
//      Holds information needed during backtracking. See description
//      of BacktrackElement and BacktrackElementDated for more details.
//
// g_level_speciations
// g_level_edges
//
//      These two variables hold information about the splitting of the
//      species tree when the dated version of the algorithm is
//      requested. A dated species tree with n internal leaves is split
//      into n + 1 levels, where level 0 is the level of the leaves and
//      level n is the level of the root. Except for level 0, each level
//      has a unique speciation vertex x; these are stored in
//      'g_level_speciations'. The set of edges at each level are stored
//      in 'g_level_edges'; these are useful to have when considering
//      transfer events in the DP algorithm as well as during the
//      backtracking. These variables are initialized and filled by the
//      'dp_algorithm_dated()' function.
//*****************************************************************************
const cost_type COST_INF = numeric_limits<cost_type>::infinity();

multi_array<cost_type, 3>               g_below;
multi_array<cost_type, 3>               g_outside;
multi_array<BacktrackElement, 3>        g_backtrack_matrix;
multi_array<BacktrackElementDated, 3>   g_backtrack_matrix_dated;
vector<vid_t>                           g_level_speciations;
vector< vector<vid_t> >                 g_level_edges;

po::variables_map                       g_options;
struct ProgramInput {
    string                  species_tree_fname;
    string                  gene_tree_fname;
    string                  sigma_fname;
    shared_ptr<tree_type>   species_tree;
    shared_ptr<tree_type>   gene_tree;
    vector<unsigned>        sigma;
    vector<unsigned>        gene_tree_numbering;
    cost_type               duplication_cost;
    cost_type               transfer_cost;
    bool                    print_only_cost;
    bool                    use_species_tree_dates;
    bool                    print_only_minimal_transfer_scenarios;
    bool                    print_only_minimal_loss_scenarios;
    bool                    unsorted;
}                                       g_input;


//*****************************************************************************
// print_error()
//
// Prints message to stderr with a suffix identifying the program.
//*****************************************************************************
void print_error(const char *);

//*****************************************************************************
// parse_options()
//
// Parses the command line options and fills in values into
// g_input. Returns an options_description suitable for printing with
// the program's help message.
//*****************************************************************************
po::options_description parse_options(int argc, char*argv[]);

//*****************************************************************************
// check_options()
//
// Should be called after parse_options. Returns true if options given
// are consistent and all arguments have been provided. If any errors
// are detected, a message is printed to stderr and false is
// returned. The caller can then simply end the program.
//*****************************************************************************
bool check_options();

//*****************************************************************************
// check_input()
//
// Should be called to check the data in g_input once everything has
// been read. Performs (consistency) checks that are not performed
// when getting/reading the data members of g_input. In this case, it
// checks that the costs of transfer and duplication are not too
// large. If any errors are detected, a message is printed to stderr
// and false is returned. Otherwise true is returned.
//*****************************************************************************
bool check_input();

//*****************************************************************************
// read_species_tree()
// read_gene_tree()
//
// The functions read the files whose name is given in
// g_input.species_tree_fname and g_input.gene_tree_fname, and
// attempt to construct binary trees. If errors are detected, an
// appropriate error message is written to stderr and false is
// returned. Otherwise, true is returned.
//*****************************************************************************
bool read_species_tree();
bool read_gene_tree();

//*****************************************************************************
// read_sigma()
//
// Constructs the mapping g_input.sigma of the gene tree leaves to the
// species tree leaves by reading the file whose filename is given in
// g_input.sigma_fname. If any errors are detected, an apropriate error
// message is written to stderr and false is returned. Otherwise, true
// is returned.
//*****************************************************************************
bool read_sigma();

//*****************************************************************************
// dp_algorithm()
// dp_algorithm_dated()
//
// dp_algorithm() runs the standard dynamic programming algorithm. For
// details see the relevant article cited in the beginning of the
// file. dp_algorithm_dated() is the modified version that uses
// timed/dated species trees.
//*****************************************************************************
void dp_algorithm();
void dp_algorithm_dated();


//*****************************************************************************
// compute_below()
// compute_below_dated()
// compute_outside()
//
// These are helper functions used by dp_algorithm() and dp_algorithm_dated().
//*****************************************************************************
void compute_below(vid_t u, vid_t x);
void compute_outside(vid_t u, vid_t x);
void compute_below_dated(vid_t u, vid_t x, int level,
                         const vector<vid_t> &edges, vid_t level_speciation);

//*****************************************************************************
// backtrack()
// backtrack_dated()
//
// Run the backtrack algorithms and finds all optimal
// DTL-scenarios. Must be called after dp_algorithm() or
// dp_algorithm_dated(), respectively. The scenarios are returned via
// the return vector that is passed as argument. The flags
// g_input.print_only_minimal_transfer_scenarios and
// g_input.print_only_minimal_loss_scenarios affect the behaviour of
// backtrack(). If both flags are set, the scenarios with minimal
// transfers are found first, and among these, the ones with minimal
// number of losses are inserted into the return vector. These flags are
// currently not implemented for the dated version of the
// algorithm.
//*****************************************************************************
void backtrack(vector<Scenario> &);
void backtrack_dated(vector<Scenario> &);

//*****************************************************************************
// backtrack_below_placements()
// backtrack_outside_placements()
// backtrack_mark_needed_scenarios_below()
// backtrack_min_transfers()
// backtrack_scenarios_at()
// combine_scenarios()
//
// These are helper functions called by backtrack(). Since no
// documentation is yet available on how the backtrack algorithm
// works, some short descriptions about what the functions do is given
// here. u denotes a gene tree vertex and x a species tree vertex.
//
// backtrack_below_placements(u, x)
//      Finds the descendants of x _at_ which u can be placed to
//      obtain the minimal cost g_below[u][x]. Inserts the vertices
//      into g_backtrack_matrix[u][x].below_placements. If x is one of
//      the descendants at which u may be placed, it will be the first
//      vertex inserted. This function assumes that the
//      below_placements of the descendants of u and descendants of x
//      have already been computed.
//
// backtrack_outside_placements()
//      Finds the vertices incomparable to x _below_ which u may be
//      placed to obtain the minimal cost g_outside[u][x]. The
//      vertices are inserted into the vector that is passed as
//      argument. This function only relies on the outside_sibling and
//      outside_ancestor members of BacktrackElement.
//
// backtrack_mark_needed_scenarios_below()
//      Determines which scenarios need to be computed if u is to be
//      placed _at_ x. This function sets the members
//      scenarios_below_needed, and in the case of duplications,
//      some scenarios_at_needed, of BacktrackElement.
//
// backtrack_min_transfers()
//      This function computes the minimal number of transfers that
//      are needed when placing u below x. It assumes that the minimum
//      transfers for descendants of u and x have already been
//      computed.
//
// backtrack_scenarios_at()
//      This function actually computes the minimal cost scenarios
//      where u is placed _at_ x, and inserts these into
//      g_backtrack_matrix[u][x].scenarios_at. This function assumes
//      that the needed scenarios for all descendants of u at all
//      species tree vertices have already been computed.
//
// combine_scenarios()
//      Combines each scenario in the first vector with each scenario
//      in the other vector (basically a union operation) and inserts
//      the resulting scenarios in the back of
//      g_backtrack_matrix[u][x].scenarios_at. The Event passed is
//      used to set the bits of transfer_edges or duplications in the
//      newly created scenarios.
//*****************************************************************************
void backtrack_below_placements(vid_t u, vid_t x);
void backtrack_outside_placements(vid_t u, vid_t x, vector<vid_t> &);
void backtrack_mark_needed_scenarios_below(vid_t u, vid_t x);
void backtrack_min_transfers(vid_t u, vid_t x);
void backtrack_scenarios_at(vid_t u, vid_t x);
void combine_scenarios(const vector<Scenario> &, const vector<Scenario> &,
                       vid_t u, vid_t x, BacktrackElement::Event);

//*****************************************************************************
// backtrack_dated_mark_needed()
// backtrack_dated_compute_scenarios()
// combine_scenarios_dated()
//
// These are helper functions used by backtrack_dated().
//
// backtrack_dated_mark_needed()
//      Recursively marks the elements of g_backtrack_matrix_dated to
//      indicate where scenarios are need to be computed. This function
//      should be called with level = level of the root, u = 0, and x = 0.
//
// backtrack_dated_compute_scenarios()
//      should be called for each element of 'g_backtrack_matrix_dated'
//      that has 'scenarios_needed' set to true (after a call to
//      'backtrack_dated_mark_needed()'). The order of the calls is
//      important and should guarantee that required scenarios have been
//      previously computed. level should start at 0 (i.e., the level of
//      the leaves) and each gene tree vertex (in postorder) should be
//      placed on all edges of the level.
//
// combine_scenarios_dated()
//      Combines each scenario in the first vector with each scenario in
//      the second according the 'Event' and places the scenarios at the
//      end of g_backtrack_matrix[level][u][x].scenarios_at. Duplicated
//      scenarios may be inserted and need to be removed by the caller.
//*****************************************************************************

void backtrack_dated_mark_needed(int level, vid_t u, vid_t x);
void backtrack_dated_compute_scenarios(int level, vid_t u, vid_t x);
void combine_scenarios_dated(const vector<Scenario> &, const vector<Scenario> &,
                             int level, vid_t u, vid_t x,
                             BacktrackElementDated::Event);


int
main(int argc, char *argv[])
{
    // Parse command line options.
    try
        {
            po::options_description visible_opts = parse_options(argc, argv);
            if (g_options.count("help"))
                {
                    cout << USAGE << "\n"
                         << visible_opts << "\n";
                    exit(EXIT_SUCCESS);
                }
        }
    catch (po::error &e)
        {
            print_error(e.what());
            cerr << "Try '" << PROGRAM_NAME << " "
                 << "--help' for more information.\n";
            exit(EXIT_FAILURE);
        }


    // Read the input and check that everything is ok.
    bool input_ok =
        check_options() &&
        read_species_tree() &&
        read_gene_tree() &&
        read_sigma() &&
        check_input();

    if (!input_ok)
        {
            exit(EXIT_FAILURE);
        }

    // We will use postorder to number the vertices of the gene tree for output.
    get_postorder_numbering(*g_input.gene_tree, g_input.gene_tree_numbering);


    // Run the main algorithm of this program.
    if (g_input.use_species_tree_dates)
        {
            dp_algorithm_dated();
        }
    else
        {
            dp_algorithm();
        }


    // Output the results.
    if (g_input.print_only_cost)
        {
            if (g_input.use_species_tree_dates)
                {
                    const tree_type &S = *g_input.species_tree;
                    cout << g_below[((S.size() + 1) / 2) - 1][0][0] << "\n";
                }
            else
                {
                    cout << g_below[0][0][0] << "\n";
                }
        }
    else
        {
            vector<Scenario> scenarios;
            if (g_input.use_species_tree_dates)
                {
                    backtrack_dated(scenarios);
                }
            else
                {
                    backtrack(scenarios);
                }
            if (!g_input.unsorted)
                {
                    sort(scenarios.begin(), scenarios.end());
                }
            BOOST_FOREACH (Scenario sc, scenarios)
                {
                    cout << sc << "\n";
                }
        }

    return EXIT_SUCCESS;
}



void
print_error(const char *msg)
{
    cerr << PROGRAM_NAME << ": " << msg << "\n";
}



po::options_description
parse_options(int argc, char*argv[])
{
    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    // Declare options that are shown when --help is given on the
    // command line.
    visible_opts.add_options()
        ("help,h", "display this help and exit")
        ("transfer-cost,t",
         po::value<cost_type>(&g_input.transfer_cost)->default_value(1.0),
         "Cost of transfer events")
        ("duplication-cost,d",
         po::value<cost_type>(&g_input.duplication_cost)->default_value(1.0),
         "Cost of duplication events")
        ("cost-only,c",
         po::bool_switch(&g_input.print_only_cost)->default_value(false),
         "If set, only the optimum cost is printed")
        ("use-species-tree-dates",
         po::bool_switch(&g_input.use_species_tree_dates)->default_value(false),
         "If set, the species tree dates are taken into account")
        ("minimum-transfers",
         po::bool_switch(&g_input.print_only_minimal_transfer_scenarios)
         ->default_value(false),
         "If set, only the scenarios with minimum number "
         "of transfers are printed. See also --minimum-losses. "
         "This is currently not implemented in combination with the "
         "--use-species-tree-dates flag.")
        ("minimum-losses",
         po::bool_switch(&g_input.print_only_minimal_loss_scenarios)
         ->default_value(false),
         "If set, only the scenarios with minimum number "
         "of losses are printed. This option may be combined with "
         "--minimum-transfers in which case only the scenarios with minimal "
         "losses among those with minimum number of transfers are written. "
         "In other words, transfers are minimized first. "
         "This is currently not implemented in combination with the "
         "--use-species-tree-dates flag.")
        ("unsorted",
         po::bool_switch(&g_input.unsorted)->default_value(false),
         "If set, the scenarios will not be sorted") ;

    // Declare the positional options
    hidden_opts.add_options()
        ("species-tree-file",po::value<string>(&g_input.species_tree_fname))
        ("gene-tree-file", po::value<string>(&g_input.gene_tree_fname))
        ("map-file", po::value<string>(&g_input.sigma_fname))
        ;

    po::positional_options_description positional_options;
    positional_options.add("species-tree-file", 1);
    positional_options.add("gene-tree-file", 1);
    positional_options.add("map-file", 1);

    // Add all options into one options_description for parsing
    all_options.add(visible_opts).add(hidden_opts);

    // Parse options.
    po::command_line_parser parser(argc, argv);
    parser.options(all_options);
    parser.positional(positional_options);
    po::store(parser.run(), g_options);
    po::notify(g_options);

    return visible_opts;
}



bool
check_options()
{
    // Check that all mandatory arguments are given.
    if (g_options.count("species-tree-file") +
        g_options.count("gene-tree-file") +
        g_options.count("map-file") != 3)
        {
            string msg = "Too few arguments\nTry '" + PROGRAM_NAME +
                " --help' for more information";
            print_error(msg.c_str());
            return false;
        }

    return true;
}



bool
check_input()
{
    // Check that the maximum possible cost fits in cost_type.
    unsigned max_events
        = g_input.gene_tree->size() / 2 + 1; // #internal vertices of G
    cost_type max_cost_allowed = numeric_limits<cost_type>::max() / max_events;
    if (g_input.duplication_cost > max_cost_allowed)
        {
            print_error("duplication cost too large");
            return false;
        }
    if (g_input.transfer_cost > max_cost_allowed)
        {
            print_error("transfer cost too large");
            return false;
        }

    // The combination of --min-transfers or --min-losses with
    // --use-species-tree-dates is currently not implemented.
    if (g_input.use_species_tree_dates &&
        (g_input.print_only_minimal_transfer_scenarios ||
         g_input.print_only_minimal_loss_scenarios))
        {
            print_error("currently, --min-transfers and --min-losses "
                        "cannot be combined with --use-species-tree-dates.");
            return false;
        }

    // If species tree dates are to be used, check that no two
    // speciations have exactly the same time and that the time of
    // leaves are after all speciations
    const tree_type &S = *g_input.species_tree;

    vector<tree_type::real_type> times;
    if (g_input.use_species_tree_dates)
        {
            for (vid_t x = 0; x < S.size(); ++x)
                {
                    if (!S.is_leaf(x))
                        {
                            times.push_back(S.time(x));
                        }
                }
            sort(times.begin(), times.end());
            if (unique(times.begin(), times.end()) != times.end())
                {
                    print_error("different speciations cannot occur at the same time");
                    return(false);
                }

            tree_type::real_type max_time = *max_element(times.begin(), times.end());
            for (vid_t x = 0; x < S.size(); ++x)
            {
                if (S.is_leaf(x) && S.time(x) <= max_time)
                {
                    print_error("leaf time is less than speciation time");
                    return(false);
                }
            }
        }


    return true;
}



bool
read_species_tree()
{
    // Open the tree file
    errno = 0;
    FILE *stree_file = fopen(g_input.species_tree_fname.c_str(), "r");
    if (errno != 0 || stree_file == 0)
        {
            string msg = "failed to open '" + g_input.species_tree_fname + "'";
            print_error(msg.c_str());
            return false;
        }

    // Read the tree
    NHNode *root = NH_read_tree(stree_file, 0, 0);
    if (root == 0)
        {
            print_error("failed to read species tree");
            return false;
        }

    // Check that the tree is binary
    if (!NHtree_is_binary(root))
        {
            print_error("species tree is not binary");
            return false;
        }

    // Create the species tree from its NH-equivalent
    g_input.species_tree.reset(new tree_type());
    create_binary_tree(root, *g_input.species_tree, 0);

    return true;
}



bool
read_gene_tree()
{
    // Open the tree file
    errno = 0;
    FILE *gtree_file = fopen(g_input.gene_tree_fname.c_str(), "r");
    if (errno != 0 || gtree_file == 0)
        {
            string msg = "failed to open '" + g_input.gene_tree_fname + "'";
            print_error(msg.c_str());
            return false;
        }

    // Read the tree
    NHNode *root = NH_read_tree(gtree_file, 0, 0);
    if (root == 0)
        {
            print_error("failed to read gene tree");
            return false;
        }

    // Check that the tree is binary
    if (!NHtree_is_binary(root))
        {
            print_error("gene tree is not binary");
            return false;
        }

    // Create the species tree from its NH-equivalent
    g_input.gene_tree.reset(new tree_type());
    create_binary_tree(root, *g_input.gene_tree, 0);

    return true;
}



bool
read_sigma()
{
    g_input.sigma.resize(g_input.gene_tree->size());
    try
        {
            create_gene_species_map(*g_input.species_tree,
                                    *g_input.gene_tree,
                                    g_input.sigma_fname,
                                    g_input.sigma);
        }
    catch (logic_error &e)
        {
            print_error(e.what());
            return false;
        }

    return true;
}



void
dp_algorithm()
{
    const tree_type        &S = *g_input.species_tree;
    const tree_type        &G = *g_input.gene_tree;
    const bool              do_backtrack = !g_input.print_only_cost;



    // Allocate memory for the global matrices and transfer destinations.
    g_below.resize(boost::extents[1][G.size()][S.size()]);
    g_outside.resize(boost::extents[1][G.size()][S.size()]);
    if (do_backtrack)
        {
            g_backtrack_matrix.resize(boost::extents[1][G.size()][S.size()]);
        }

    // Initialize below and outside to infinity.
    for (vid_t u = 0; u < G.size(); ++u)
        {
            for (vid_t x = 0; x < S.size(); ++x)
                {
                    g_below[0][u][x] = COST_INF;
                    g_outside[0][u][x] = COST_INF;
                }
        }



    // The algorithm itself is described in a published paper.
    for (vid_t u = G.postorder_begin();
         u != tree_type::NONE;
         u = G.postorder_next(u))
        {
            // First compute g_below[u][*].
            for (vid_t x = S.postorder_begin();
                 x != tree_type::NONE;
                 x = S.postorder_next(x))
                {
                    compute_below(u, x);
                }

            // Compute g_outside[u][*].
            for (vid_t x = S.preorder_begin();
                 x != S.NONE;
                 x = S.preorder_next(x))
                {
                    compute_outside(u, x);
                }
        }
}



void
compute_below(vid_t u, vid_t x)
{
    const tree_type        &S = *g_input.species_tree;
    const tree_type        &G = *g_input.gene_tree;
    const vector<vid_t>    &sigma = g_input.sigma;
    const cost_type         tcost = g_input.transfer_cost;
    const cost_type         dcost = g_input.duplication_cost;
    const bool              do_backtrack = !g_input.print_only_cost;

    if (G.is_leaf(u))
        {
            if (S.descendant(sigma[u], x))
                {
                    g_below[0][u][x] = 0;
                }
        }
    else
        {
            vector<cost_type> costs(BacktrackElement::N_EVENTS, COST_INF);

            costs[BacktrackElement::D] =
                dcost + g_below[0][G.left(u)][x] + g_below[0][G.right(u)][x];
            costs[BacktrackElement::T_LEFT] =
                tcost + g_outside[0][G.left(u)][x] + g_below[0][G.right(u)][x];
            costs[BacktrackElement::T_RIGHT] =
                tcost + g_outside[0][G.right(u)][x] + g_below[0][G.left(u)][x];

            if (!S.is_leaf(x))
                {
                    costs[BacktrackElement::S] =
                        g_below[0][G.left(u)][S.left(x)] +
                        g_below[0][G.right(u)][S.right(x)];
                    costs[BacktrackElement::S_REV] =
                        g_below[0][G.left(u)][S.right(x)] +
                        g_below[0][G.right(u)][S.left(x)];

                    costs[BacktrackElement::BELOW_LEFT] = g_below[0][u][S.left(x)];
                    costs[BacktrackElement::BELOW_RIGHT] = g_below[0][u][S.right(x)];
                }

            cost_type min_cost = *min_element(costs.begin(), costs.end());
            g_below[0][u][x] = min_cost;

            // Save the optimal events for backtracking.
            if (do_backtrack)
                {
                    bitset<BacktrackElement::N_EVENTS> &events =
                        g_backtrack_matrix[0][u][x].below_events;
                    for (unsigned e = 0; e < BacktrackElement::N_EVENTS; ++e)
                        {
                            if (min_cost != COST_INF && costs[e] == min_cost)
                                {
                                    events.set(BacktrackElement::Event(e));
                                }
                        }
                }
        }
}



void
compute_outside(vid_t u, vid_t x)
{
    const tree_type        &S = *g_input.species_tree;
    const bool              do_backtrack = !g_input.print_only_cost;

    // Cannot place u outside the root of S.
    if (x == 0)
        return;

    vid_t x_parent = S.parent(x);
    vid_t x_sibling =
        S.left(x_parent) == x ? S.right(x_parent) : S.left(x_parent);

    cost_type min_cost = min(g_below[0][u][x_sibling], g_outside[0][u][x_parent]);

    g_outside[0][u][x] = min_cost;

    // Save info for backtracking.
    if (do_backtrack)
        {
            if (g_below[0][u][x_sibling] == min_cost)
                {
                    g_backtrack_matrix[0][u][x].outside_sibling = x_sibling;
                }

            if (g_outside[0][u][x_parent] == min_cost)
                {
                    if (g_backtrack_matrix[0][u][x_parent].outside_sibling !=
                        tree_type::NONE)
                        {
                            g_backtrack_matrix[0][u][x].outside_ancestor =
                                x_parent;
                        }
                    else
                        {
                            g_backtrack_matrix[0][u][x].outside_ancestor =
                                g_backtrack_matrix[0][u][x_parent].outside_ancestor;
                        }
                }
        }
}



void
dp_algorithm_dated()
{
    const tree_type        &S = *g_input.species_tree;
    const tree_type        &G = *g_input.gene_tree;
    const vector<vid_t>    &sigma = g_input.sigma;
    const bool              do_backtrack = !g_input.print_only_cost;
    const int               n_levels = (S.size() + 1) / 2;


    // Allocate memory for the global matrices and transfer destinations.
    g_below.resize(boost::extents[n_levels][G.size()][S.size()]);
    if (do_backtrack)
        {
            g_backtrack_matrix_dated.resize(boost::extents[n_levels][G.size()][S.size()]);
        }


    // initialize below to infinity
    fill(g_below.data(), g_below.data() + g_below.num_elements(), COST_INF);

    // get the unique speciation for each level in the species tree. The
    // levels go from 0 to n_levels, where 0 represents the level of the leaves.
    g_level_speciations.resize(n_levels, tree_type::NONE);

    vector< pair<tree_type::real_type, vid_t> > pairs;
    for (vid_t x = 0; x < S.size(); ++x)
        {
            if (!S.is_leaf(x))
                {
                    pairs.push_back(make_pair(S.time(x), x));
                }
        }
    sort(pairs.begin(), pairs.end());
    for (int i = 0; i < pairs.size(); ++i)
        {
            g_level_speciations[n_levels - 1 - i] = pairs[i].second;
        }

    // for each level, get the set of vertices that represent the
    // contemporary edges.
    g_level_edges.resize(n_levels);
    g_level_edges[n_levels - 1].push_back(0);
    for (int i = n_levels - 2; i >= 0; --i)
        {
            vid_t prev_speciation = g_level_speciations[i + 1];
            BOOST_FOREACH(vid_t x, g_level_edges[i + 1])
                {
                    if (x != prev_speciation)
                        {
                            g_level_edges[i].push_back(x);
                        }
                }
            g_level_edges[i].push_back(S.left(prev_speciation));
            g_level_edges[i].push_back(S.right(prev_speciation));
        }

    // run the DP algorithm.
    for (int level = 0; level < n_levels; ++level)
        {
            for (vid_t u = G.postorder_begin();
                 u != tree_type::NONE;
                 u = G.postorder_next(u))
                {
                    BOOST_FOREACH(vid_t x, g_level_edges[level])
                        {
                            if (G.is_leaf(u))
                                {
                                    if (S.descendant(sigma[u], x))
                                        {
                                            g_below[level][u][x] = 0;
                                        }
                                }
                            else
                                {
                                    compute_below_dated(u, x, level,
                                                        g_level_edges[level],
                                                        g_level_speciations[level]);
                                }
                        }
                }
        }
}



void
compute_below_dated(vid_t u, vid_t x,
                    int level,
                    const vector<vid_t> &edges,
                    vid_t level_speciation)
{
    typedef BacktrackElementDated BED;

    const tree_type        &S = *g_input.species_tree;
    const tree_type        &G = *g_input.gene_tree;
    const vector<vid_t>    &sigma = g_input.sigma;
    const cost_type         tcost = g_input.transfer_cost;
    const cost_type         dcost = g_input.duplication_cost;
    const bool              do_backtrack = !g_input.print_only_cost;

    // We are assuming that the caller takes care of gene tree leaves,
    // so u is assumed never to be a leaf.
    if (G.is_leaf(u))
        {
            print_error("internal logical error: compute_below_dated() called "
                        "with a gene tree leaf");
            exit(EXIT_FAILURE);
        }

    // compute the cost of each possible event and store in the vector 'costs'.
    vector<cost_type> costs(BED::N_EVENTS, COST_INF);
    costs[BED::D] =
        dcost + g_below[level][G.left(u)][x] + g_below[level][G.right(u)][x];
    BOOST_FOREACH(vid_t y, edges)
        {
            if (y == x)
                continue;

            costs[BED::T_LEFT] = min(costs[BED::T_LEFT],
                                     g_below[level][G.left(u)][y] +
                                     g_below[level][G.right(u)][x]);
            costs[BED::T_RIGHT] = min(costs[BED::T_RIGHT],
                                      g_below[level][G.right(u)][y] +
                                      g_below[level][G.left(u)][x]);
        }
    costs[BED::T_LEFT] += tcost;
    costs[BED::T_RIGHT] += tcost;

    if (level != 0 && x != level_speciation)
        {
            costs[BED::BELOW] = g_below[level - 1][u][x];
        }
    else if (level != 0 && x == level_speciation)
        {
            costs[BED::S] =
                g_below[level - 1][G.left(u)][S.left(x)] +
                g_below[level - 1][G.right(u)][S.right(x)];
            costs[BED::S_REV] =
                g_below[level - 1][G.right(u)][S.left(x)] +
                g_below[level - 1][G.left(u)][S.right(x)];
            costs[BED::BELOW_LEFT] =
                g_below[level - 1][u][S.left(x)];
            costs[BED::BELOW_RIGHT] =
                g_below[level - 1][u][S.right(x)];
        }

    cost_type min_cost = *min_element(costs.begin(), costs.end());
    g_below[level][u][x] = min_cost;

    // store information for the backtracking.
    if (do_backtrack)
        {
            bitset<BED::N_EVENTS> &events =
                g_backtrack_matrix_dated[level][u][x].below_events;
            for (unsigned e = 0; e < BED::N_EVENTS; ++e)
                {
                    if (min_cost != COST_INF && costs[e] == min_cost)
                        {
                            events.set(BED::Event(e));
                        }
                }
        }
}



void
backtrack(vector<Scenario> &scenarios)
{
    const tree_type                    &S = *g_input.species_tree;
    const tree_type                    &G = *g_input.gene_tree;
    const vector<vid_t>                &sigma = g_input.sigma;
    multi_array<BacktrackElement, 3>   &matrix = g_backtrack_matrix;

    // Backtrack the placements for each u and x.
    for (vid_t u = G.postorder_begin(); u != tree_type::NONE; u = G.postorder_next(u))
        {
            for (vid_t x = S.postorder_begin(); x != tree_type::NONE; x = S.postorder_next(x))
                {
                    backtrack_below_placements(u, x);
                }
        }

    // Mark the sets of scenarios that we need to compute.
    matrix[0][0][0].scenarios_below_needed = true;
    for (vid_t u = G.preorder_begin(); u != tree_type::NONE; u = G.preorder_next(u))
        {
            for (vid_t x = 0; x < S.size(); ++x)
                {
                    if (matrix[0][u][x].scenarios_below_needed)
                        {
                            BOOST_FOREACH (vid_t y, matrix[0][u][x].below_placements)
                                {
                                    matrix[0][u][y].scenarios_at_needed = true;
                                }
                        }
                }

            if (G.is_leaf(u))
                continue;

            for (vid_t x = 0; x < S.size(); ++x)
                {
                    if (matrix[0][u][x].scenarios_at_needed)
                        {
                            backtrack_mark_needed_scenarios_below(u, x);
                        }
                }
        }

    // Compute the minimum number of transfer events for each u and x.
    for (vid_t u = G.postorder_begin(); u != tree_type::NONE; u = G.postorder_next(u))
        {
            for (vid_t x = S.postorder_begin(); x != tree_type::NONE; x = S.postorder_next(x))
                {
                    backtrack_min_transfers(u, x);
                }
        }

    // Backtrack the needed scenarios.
    for (vid_t u = G.postorder_begin(); u < G.size(); u = G.postorder_next(u))
        {
            for (vid_t x = 0; x < S.size(); ++x)
                {
                    if (matrix[0][u][x].scenarios_at_needed)
                        {
                            backtrack_scenarios_at(u, x);
                        }
                }

            // Remove the unneeded sets of scenarios to conserve memory.
            if (!G.is_leaf(u))
                {
                    for (vid_t x = 0; x < S.size(); ++x)
                        {
                            vector<Scenario>().swap(matrix[0][G.left(u)][x].scenarios_at);
                            vector<Scenario>().swap(matrix[0][G.right(u)][x].scenarios_at);
                        }
                }
        }

    // Find the minimum number of losses of placing root of G below
    // root of S.
    int max_losses = numeric_limits<int>::max();
    if (g_input.print_only_minimal_loss_scenarios)
        {
            BOOST_FOREACH(vid_t x, matrix[0][0][0].below_placements)
                {
                    BOOST_FOREACH(Scenario &sc, matrix[0][0][x].scenarios_at)
                        {
                            max_losses = min(max_losses,
                                             count_losses(S, G, sigma, sc.transfer_edges));
                        }
                }
        }

    // Find the final sets of scenarios.
    BOOST_FOREACH(vid_t x, matrix[0][0][0].below_placements)
        {
            // Take only scenarios with minimal transfers if the flag is set.
            if (g_input.print_only_minimal_transfer_scenarios &&
                !matrix[0][0][x].scenarios_at.empty() &&
                matrix[0][0][x].scenarios_at[0].transfer_edges.count() >
                matrix[0][0][0].min_transfers)
                {
                    vector<Scenario>().swap(matrix[0][0][x].scenarios_at);
                    continue;
                }
            // Take only scenarios with minimal losses if the flag is set.
            BOOST_FOREACH(Scenario &sc, matrix[0][0][x].scenarios_at)
                {
                    if (count_losses(S, G, sigma, sc.transfer_edges) <= max_losses)
                        {
                            scenarios.push_back(sc);
                        }
                }
            vector<Scenario>().swap(matrix[0][0][x].scenarios_at);
        }

    return;
}



void
backtrack_dated(vector<Scenario> &scenarios)
{
    const tree_type        &S = *g_input.species_tree;
    const tree_type        &G = *g_input.gene_tree;
    const vector<vid_t>    &sigma = g_input.sigma;
    const int               n_levels = (S.size() + 1) / 2;

    backtrack_dated_mark_needed(n_levels - 1, 0, 0);

    for (int level = 0; level < n_levels; ++level)
        {
            for (vid_t u = G.postorder_begin();
                 u != tree_type::NONE;
                 u = G.postorder_next(u))
                {
                    BOOST_FOREACH(vid_t x, g_level_edges[level])
                        {
                            if (g_backtrack_matrix_dated[level][u][x].scenarios_needed)
                                {
                                    backtrack_dated_compute_scenarios(level, u, x);
                                }
                        }
                }

            // once scenarios at level i have been created, we will
            // never need scenarios from level i - 1. So we deallocate
            // all those scenarios from previous levels to conserve
            // memory.
            if (level > 0)
                {
                    for (vid_t x = 0; x < S.size(); ++x)
                        {
                            for (vid_t u = 0; u < G.size(); ++u)
                                {
                                    vector<Scenario>().swap(g_backtrack_matrix_dated[level - 1][u][x].scenarios_at);
                                }
                        }
                }
        }
    scenarios = g_backtrack_matrix_dated[n_levels - 1][0][0].scenarios_at;
}



void
backtrack_dated_mark_needed(int level, vid_t u, vid_t x)
{
    typedef BacktrackElementDated BED;

    const tree_type        &S = *g_input.species_tree;
    const tree_type        &G = *g_input.gene_tree;

    // return if we have already been here
    if (g_backtrack_matrix_dated[level][u][x].scenarios_needed)
        return;

    g_backtrack_matrix_dated[level][u][x].scenarios_needed = true;
    bitset<BED::N_EVENTS> below_events =
        g_backtrack_matrix_dated[level][u][x].below_events;
    if (below_events[BED::S])
        {
            backtrack_dated_mark_needed(level - 1, G.left(u), S.left(x));
            backtrack_dated_mark_needed(level - 1, G.right(u), S.right(x));
        }
    if (below_events[BED::S_REV])
        {
            backtrack_dated_mark_needed(level - 1, G.left(u), S.right(x));
            backtrack_dated_mark_needed(level - 1, G.right(u), S.left(x));
        }
    if (below_events[BED::D])
        {
            backtrack_dated_mark_needed(level, G.left(u), x);
            backtrack_dated_mark_needed(level, G.right(u), x);
        }
    if (below_events[BED::T_LEFT])
        {
            BOOST_FOREACH(vid_t y, g_level_edges[level])
                {
                    if (y == x)
                        continue;

                    cost_type cost = g_below[level][G.left(u)][y] +
                        g_below[level][G.right(u)][x] + g_input.transfer_cost;
                    if (cost == g_below[level][u][x])
                        {
                            backtrack_dated_mark_needed(level, G.left(u), y);
                            backtrack_dated_mark_needed(level, G.right(u), x);
                        }
                }
        }
    if (below_events[BED::T_RIGHT])
        {
            BOOST_FOREACH(vid_t y, g_level_edges[level])
                {
                    if (y == x)
                        continue;

                    cost_type cost = g_below[level][G.right(u)][y] +
                        g_below[level][G.left(u)][x] + g_input.transfer_cost;
                    if (cost == g_below[level][u][x])
                        {
                            backtrack_dated_mark_needed(level, G.right(u), y);
                            backtrack_dated_mark_needed(level, G.left(u), x);
                        }
                }
        }
    if (below_events[BED::BELOW_LEFT])
        {
            backtrack_dated_mark_needed(level - 1, u, S.left(x));
        }
    if (below_events[BED::BELOW_RIGHT])
        {
            backtrack_dated_mark_needed(level - 1, u, S.right(x));
        }
    if (below_events[BED::BELOW])
        {
            backtrack_dated_mark_needed(level - 1, u, x);
        }
}



void
backtrack_dated_compute_scenarios(int level, vid_t u, vid_t x)
{
    typedef BacktrackElementDated BED;

    multi_array<BacktrackElementDated, 3>  &matrix = g_backtrack_matrix_dated;
    const tree_type        &S = *g_input.species_tree;
    const tree_type        &G = *g_input.gene_tree;
    vector<Scenario>       &scenarios_at = matrix[level][u][x].scenarios_at;
    bitset<BED::N_EVENTS>   below_events = matrix[level][u][x].below_events;

    if (G.is_leaf(u))
        {
            scenarios_at.push_back(Scenario(G.size()));
            return;
        }

    if (below_events[BED::S])
        {
            combine_scenarios_dated(matrix[level - 1][G.left(u)][S.left(x)].scenarios_at,
                                    matrix[level - 1][G.right(u)][S.right(x)].scenarios_at,
                                    level, u, x, BED::S);
        }
    if (below_events[BED::S_REV])
        {
            combine_scenarios_dated(matrix[level - 1][G.left(u)][S.right(x)].scenarios_at,
                                    matrix[level - 1][G.right(u)][S.left(x)].scenarios_at,
                                    level, u, x, BED::S_REV);
        }
    if (below_events[BED::D])
        {
            combine_scenarios_dated(matrix[level][G.left(u)][x].scenarios_at,
                                    matrix[level][G.right(u)][x].scenarios_at,
                                    level, u, x, BED::D);
        }
    if (below_events[BED::T_LEFT])
        {
            BOOST_FOREACH(vid_t y, g_level_edges[level])
                {
                    if (y == x)
                        continue;

                    cost_type cost = g_below[level][G.left(u)][y] +
                        g_below[level][G.right(u)][x] + g_input.transfer_cost;
                    if (cost == g_below[level][u][x])
                        {
                            combine_scenarios_dated(matrix[level][G.left(u)][y].scenarios_at,
                                                    matrix[level][G.right(u)][x].scenarios_at,
                                                    level, u, x, BED::T_LEFT);
                        }
                }
        }
    if (below_events[BED::T_RIGHT])
        {
            BOOST_FOREACH(vid_t y, g_level_edges[level])
                {
                    if (y == x)
                        continue;

                    cost_type cost = g_below[level][G.right(u)][y] +
                        g_below[level][G.left(u)][x] + g_input.transfer_cost;
                    if (cost == g_below[level][u][x])
                        {
                            combine_scenarios_dated(matrix[level][G.right(u)][y].scenarios_at,
                                                    matrix[level][G.left(u)][x].scenarios_at,
                                                    level, u, x, BED::T_RIGHT);
                        }
                }
        }
    if (below_events[BED::BELOW_LEFT])
        {
            const vector<Scenario> &to_add = matrix[level - 1][u][S.left(x)].scenarios_at;
            scenarios_at.insert(scenarios_at.end(), to_add.begin(), to_add.end());
        }
    if (below_events[BED::BELOW_RIGHT])
        {
            const vector<Scenario> &to_add = matrix[level - 1][u][S.right(x)].scenarios_at;
            scenarios_at.insert(scenarios_at.end(), to_add.begin(), to_add.end());
        }
    if (below_events[BED::BELOW])
        {
            const vector<Scenario> &to_add = matrix[level - 1][u][x].scenarios_at;
            scenarios_at.insert(scenarios_at.end(), to_add.begin(), to_add.end());
        }

    // remove duplicated scenarios
    sort(scenarios_at.begin(), scenarios_at.end(), scenario_less_lex);
    scenarios_at.erase(unique(scenarios_at.begin(), scenarios_at.end()),
                       scenarios_at.end());
}



void
backtrack_below_placements(vid_t u, vid_t x)
{
    const tree_type &S = *g_input.species_tree;
    const tree_type &G = *g_input.gene_tree;

    BacktrackElement &elem = g_backtrack_matrix[0][u][x];

    if (g_below[0][u][x] == COST_INF) // If no solutions exist
        return;

    if (G.is_leaf(u))
        {
            elem.below_placements.push_back(g_input.sigma[u]);
        }
    else if (S.is_leaf(x))
        {
            elem.below_placements.push_back(x);
        }
    else
        {
            BacktrackElement &left_elem = g_backtrack_matrix[0][G.left(u)][x];
            BacktrackElement &right_elem = g_backtrack_matrix[0][G.right(u)][x];
            const bitset<BacktrackElement::N_EVENTS> &e = elem.below_events;

            // First, determine if u is placed _at_ x.
            if (e[BacktrackElement::S] || e[BacktrackElement::S_REV] ||
                (e[BacktrackElement::T_LEFT] && right_elem.below_placements[0] == x) ||
                (e[BacktrackElement::T_RIGHT] && left_elem.below_placements[0] == x) ||
                (e[BacktrackElement::D] && right_elem.below_placements[0] == x) ||
                (e[BacktrackElement::D] && left_elem.below_placements[0] == x))
                {
                    elem.below_placements.push_back(x);
                }


            // Then, see where u is placed below x.
            if (elem.below_events[BacktrackElement::BELOW_LEFT])
                {
                    vector<vid_t> &left_placements =
                        g_backtrack_matrix[0][u][S.left(x)].below_placements;
                    copy(left_placements.begin(), left_placements.end(),
                         back_inserter(elem.below_placements));
                }
            if (elem.below_events[BacktrackElement::BELOW_RIGHT])
                {
                    vector<vid_t> &right_placements =
                        g_backtrack_matrix[0][u][S.right(x)].below_placements;
                    copy(right_placements.begin(), right_placements.end(),
                         back_inserter(elem.below_placements));
                }
        }
}



void
backtrack_mark_needed_scenarios_below(vid_t u, vid_t x)
{
    const tree_type                     &S = *g_input.species_tree;
    const tree_type                     &G = *g_input.gene_tree;
    multi_array<BacktrackElement, 3>    &matrix = g_backtrack_matrix;

    const bitset<BacktrackElement::N_EVENTS> &e = matrix[0][u][x].below_events;
    if (e[BacktrackElement::S])
        {
            matrix[0][G.left(u)][S.left(x)].scenarios_below_needed = true;
            matrix[0][G.right(u)][S.right(x)].scenarios_below_needed = true;
        }
    if (e[BacktrackElement::S_REV])
        {
            matrix[0][G.left(u)][S.right(x)].scenarios_below_needed = true;
            matrix[0][G.right(u)][S.left(x)].scenarios_below_needed = true;
        }
    if (e[BacktrackElement::D])
        {
            // This is the only time we need to set a scenarios_at_needed.
            if (matrix[0][G.right(u)][x].below_placements[0] == x)
                {
                    matrix[0][G.right(u)][x].scenarios_at_needed = true;
                    matrix[0][G.left(u)][x].scenarios_below_needed = true;
                }

            if (matrix[0][G.left(u)][x].below_placements[0] == x)
                {
                    matrix[0][G.left(u)][x].scenarios_at_needed = true;
                    matrix[0][G.right(u)][x].scenarios_below_needed = true;
                }
        }
    if (e[BacktrackElement::T_LEFT])
        {
            matrix[0][G.right(u)][x].scenarios_below_needed = true;

            vector<vid_t> outside_placements;
            backtrack_outside_placements(G.left(u), x, outside_placements);
            BOOST_FOREACH (vid_t y, outside_placements)
                {
                    matrix[0][G.left(u)][y].scenarios_below_needed = true;
                }
        }
    if (e[BacktrackElement::T_RIGHT])
        {
            matrix[0][G.left(u)][x].scenarios_below_needed = true;

            vector<vid_t> outside_placements;
            backtrack_outside_placements(G.right(u), x, outside_placements);
            BOOST_FOREACH (vid_t y, outside_placements)
                {
                    matrix[0][G.right(u)][y].scenarios_below_needed = true;
                }
        }
}



void
backtrack_scenarios_at(vid_t u, vid_t x)
{
    const tree_type                     &G = *g_input.gene_tree;
    const tree_type                     &S = *g_input.species_tree;
    multi_array<BacktrackElement, 3>    &matrix = g_backtrack_matrix;

    if (G.is_leaf(u))
        {
            // It must be the case that sigma(u) = x, otherwise the
            // algorithm is corrupt.
            matrix[0][u][x].scenarios_at.push_back(Scenario(G.size()));
            return;
        }

    const bitset<BacktrackElement::N_EVENTS> &events = matrix[0][u][x].below_events;

    if (events[BacktrackElement::S])
        {
            BOOST_FOREACH (vid_t y1, matrix[0][G.left(u)][S.left(x)].below_placements) {
                BOOST_FOREACH (vid_t y2, matrix[0][G.right(u)][S.right(x)].below_placements) {
                    combine_scenarios(matrix[0][G.left(u)][y1].scenarios_at,
                                      matrix[0][G.right(u)][y2].scenarios_at,
                                      u, x, BacktrackElement::S);
                }
            }
        }
    if (events[BacktrackElement::S_REV])
        {
            BOOST_FOREACH (vid_t y1, matrix[0][G.left(u)][S.right(x)].below_placements) {
                BOOST_FOREACH (vid_t y2, matrix[0][G.right(u)][S.left(x)].below_placements) {
                    combine_scenarios(matrix[0][G.left(u)][y1].scenarios_at,
                                      matrix[0][G.right(u)][y2].scenarios_at,
                                      u, x, BacktrackElement::S_REV);
                }
            }
        }
    if (events[BacktrackElement::D])
        {
            // Here we have to perform more work to ensure we do not
            // get duplicate scenarios. The only way that u is mapped
            // _at_ x is if at least one of the children of u is also
            // placed _at_ x.
            if (matrix[0][G.left(u)][x].below_placements[0] == x &&
                matrix[0][G.right(u)][x].below_placements[0] == x)
                {
                    combine_scenarios(matrix[0][G.left(u)][x].scenarios_at,
                                      matrix[0][G.right(u)][x].scenarios_at,
                                      u, x, BacktrackElement::D);
                }

            if (matrix[0][G.left(u)][x].below_placements[0] == x)
                {
                    BOOST_FOREACH (vid_t y, matrix[0][G.right(u)][x].below_placements) {
                        if (y == x)
                            continue;

                        combine_scenarios(matrix[0][G.right(u)][y].scenarios_at,
                                          matrix[0][G.left(u)][x].scenarios_at,
                                          u, x, BacktrackElement::D);
                    }
                }
            if (matrix[0][G.right(u)][x].below_placements[0] == x)
                {
                    BOOST_FOREACH (vid_t y, matrix[0][G.left(u)][x].below_placements) {
                        if (y == x)
                            continue;

                        combine_scenarios(matrix[0][G.left(u)][y].scenarios_at,
                                          matrix[0][G.right(u)][x].scenarios_at,
                                          u, x, BacktrackElement::D);
                    }
                }
        }
    if (events[BacktrackElement::T_LEFT])
        {
            vector<vid_t> outside_placements;
            backtrack_outside_placements(G.left(u), x, outside_placements);

            BOOST_FOREACH (vid_t y, outside_placements) {
                BOOST_FOREACH (vid_t y1, matrix[0][G.left(u)][y].below_placements) {
                    combine_scenarios(matrix[0][G.left(u)][y1].scenarios_at,
                                      matrix[0][G.right(u)][x].scenarios_at,
                                      u, x, BacktrackElement::T_LEFT);
                }
            }
        }
    if (events[BacktrackElement::T_RIGHT])
        {
            vector<vid_t> placements;
            backtrack_outside_placements(G.right(u), x, placements);

            BOOST_FOREACH (vid_t y, placements) {
                BOOST_FOREACH (vid_t y1, matrix[0][G.right(u)][y].below_placements) {
                    combine_scenarios(matrix[0][G.right(u)][y1].scenarios_at,
                                      matrix[0][G.left(u)][x].scenarios_at,
                                      u, x, BacktrackElement::T_RIGHT);
                }
            }
        }
}



void
backtrack_min_transfers(vid_t u, vid_t x)
{
    const tree_type &G = *g_input.gene_tree;
    const tree_type &S = *g_input.species_tree;
    const bitset<BacktrackElement::N_EVENTS> &events = g_backtrack_matrix[0][u][x].below_events;

    // The base case when u is a leaf.
    if (G.is_leaf(u) && g_below[0][u][x] != COST_INF)
        {
            g_backtrack_matrix[0][u][x].min_transfers = 0;
            return;
        }

    unsigned min_transfers = G.size() + 1; // Max possible number of transfers.
    if (events[BacktrackElement::D])
        {
            min_transfers = min(min_transfers,
                                g_backtrack_matrix[0][G.left(u)][x].min_transfers +
                                g_backtrack_matrix[0][G.right(u)][x].min_transfers);
        }
    if (events[BacktrackElement::T_LEFT])
        {
            vector<vid_t> placements;
            backtrack_outside_placements(G.left(u), x, placements);
            unsigned transfers = G.size() + 1;
            BOOST_FOREACH (vid_t y, placements)
                {
                    transfers = min(transfers, g_backtrack_matrix[0][G.left(u)][y].min_transfers);
                }
            transfers += 1 + g_backtrack_matrix[0][G.right(u)][x].min_transfers;
            min_transfers = min(min_transfers, transfers);
        }
    if (events[BacktrackElement::T_RIGHT])
        {
            vector<vid_t> placements;
            backtrack_outside_placements(G.right(u), x, placements);
            unsigned transfers = G.size() + 1;
            BOOST_FOREACH (vid_t y, placements)
                {
                    transfers = min(transfers, g_backtrack_matrix[0][G.right(u)][y].min_transfers);
                }
            transfers += 1 + g_backtrack_matrix[0][G.left(u)][x].min_transfers;
            min_transfers = min(min_transfers, transfers);
        }
    if (events[BacktrackElement::S])
        {
            min_transfers = min(min_transfers,
                                g_backtrack_matrix[0][G.left(u)][S.left(x)].min_transfers +
                                g_backtrack_matrix[0][G.right(u)][S.right(x)].min_transfers);
        }
    if (events[BacktrackElement::S_REV])
        {
            min_transfers = min(min_transfers,
                                g_backtrack_matrix[0][G.left(u)][S.right(x)].min_transfers +
                                g_backtrack_matrix[0][G.right(u)][S.left(x)].min_transfers);
        }
    if (events[BacktrackElement::BELOW_LEFT])
        {
            min_transfers = min(min_transfers,
                                g_backtrack_matrix[0][u][S.left(x)].min_transfers);
        }
    if (events[BacktrackElement::BELOW_RIGHT])
        {
            min_transfers = min(min_transfers,
                                g_backtrack_matrix[0][u][S.right(x)].min_transfers);
        }
    g_backtrack_matrix[0][u][x].min_transfers = min_transfers;

}



void
backtrack_outside_placements(vid_t u, vid_t x, vector<vid_t> &placements)
{
    vid_t outside_sibling = g_backtrack_matrix[0][u][x].outside_sibling;
    if (outside_sibling != tree_type::NONE)
        {
            placements.push_back(outside_sibling);
        }

    for(vid_t cur = g_backtrack_matrix[0][u][x].outside_ancestor;
        cur != tree_type::NONE;
        cur = g_backtrack_matrix[0][u][cur].outside_ancestor)
        {
            placements.push_back(g_backtrack_matrix[0][u][cur].outside_sibling);
        }
}



void
combine_scenarios(const vector<Scenario> &vec1,
                  const vector<Scenario> &vec2,
                  vid_t u, vid_t x, BacktrackElement::Event e)
{
    const tree_type &G = *g_input.gene_tree;

    if (vec1.empty() || vec2.empty())
        return;

    if (g_input.print_only_minimal_transfer_scenarios)
        {
            unsigned transfers =
                vec1[0].transfer_edges.count() +
                vec2[0].transfer_edges.count();
            if (e == BacktrackElement::T_LEFT || e == BacktrackElement::T_RIGHT)
                {
                    transfers += 1;
                }
            if (transfers > g_backtrack_matrix[0][u][x].min_transfers)
                {
                    return;
                }
        }


    BOOST_FOREACH (const Scenario &sc1, vec1)
        {
            BOOST_FOREACH (const Scenario &sc2, vec2)
                {
                    Scenario new_sc(g_input.gene_tree->size());
                    new_sc.transfer_edges =
                        sc1.transfer_edges | sc2.transfer_edges;
                    new_sc.duplications =
                        sc1.duplications | sc2.duplications;
                    if (e == BacktrackElement::D)
                        {
                            new_sc.duplications.set(u);
                        }
                    if (e == BacktrackElement::T_LEFT)
                        {
                            new_sc.transfer_edges.set(G.left(u));
                        }
                    if (e == BacktrackElement::T_RIGHT)
                        {
                            new_sc.transfer_edges.set(G.right(u));
                        }
                    g_backtrack_matrix[0][u][x].scenarios_at.push_back(new_sc);
                }
        }
}



void
combine_scenarios_dated(const vector<Scenario> &vec1,
                        const vector<Scenario> &vec2,
                        int level, vid_t u, vid_t x,
                        BacktrackElementDated::Event e)
{
    typedef BacktrackElementDated BED;

    const tree_type &G = *g_input.gene_tree;

    if (vec1.empty() || vec2.empty())
        return;

    vector<Scenario> &scenarios_at = g_backtrack_matrix_dated[level][u][x].scenarios_at;
    BOOST_FOREACH (const Scenario &sc1, vec1)
        {
            BOOST_FOREACH (const Scenario &sc2, vec2)
                {
                    Scenario new_sc(G.size());
                    new_sc.transfer_edges =
                        sc1.transfer_edges | sc2.transfer_edges;
                    new_sc.duplications =
                        sc1.duplications | sc2.duplications;
                    if (e == BED::D)
                        {
                            new_sc.duplications.set(u);
                        }
                    if (e == BED::T_LEFT)
                        {
                            new_sc.transfer_edges.set(G.left(u));
                        }
                    if (e == BED::T_RIGHT)
                        {
                            new_sc.transfer_edges.set(G.right(u));
                        }
                    scenarios_at.push_back(new_sc);
                }
        }
}


ostream &
operator<<(ostream &out, const Scenario &sc)
{
    vector<unsigned> transfer_edges;
    for (vid_t u = 0; u < g_input.gene_tree->size(); ++u)
        {
            if (sc.transfer_edges[u])
                {
                    transfer_edges.push_back(g_input.gene_tree_numbering[u]);
                }
        }
    sort(transfer_edges.begin(), transfer_edges.end());

    vector<unsigned> duplications;
    for (vid_t u = 0; u < g_input.gene_tree->size(); ++u)
        {
            if (sc.duplications[u])
                {
                    duplications.push_back(g_input.gene_tree_numbering[u]);
                }
        }
    sort(duplications.begin(), duplications.end());

    out << "Transfer edges:\t";
    copy(transfer_edges.begin(), transfer_edges.end(),
         ostream_iterator<unsigned>(out, " "));
    out << "\nDuplications:\t";
    copy(duplications.begin(), duplications.end(),
         ostream_iterator<unsigned>(out, " "));
    out << "\nNumber of losses: " << count_losses(*g_input.species_tree,
                                                  *g_input.gene_tree,
                                                  g_input.sigma,
                                                  sc.transfer_edges);
    out << "\n";

    return out;
}



bool
operator<(const Scenario &sc1, const Scenario &sc2)
{
    if (sc1.transfer_edges.count() < sc2.transfer_edges.count())
        {
            return true;
        }
    if (sc1.transfer_edges.count() > sc2.transfer_edges.count())
        {
            return false;
        }

    int losses1 = count_losses(*g_input.species_tree, *g_input.gene_tree,
                               g_input.sigma, sc1.transfer_edges);
    int losses2 = count_losses(*g_input.species_tree, *g_input.gene_tree,
                               g_input.sigma, sc2.transfer_edges);
    if (losses1 < losses2)
        {
            return true;
        }
    if (losses1 > losses2)
        {
            return false;
        }

    dynamic_bitset<> transfer_edges1(sc1.transfer_edges.size());
    dynamic_bitset<> transfer_edges2(sc1.transfer_edges.size());
    for (unsigned i = 0; i < sc1.transfer_edges.size(); ++i)
        {
            if (sc1.transfer_edges[i])
                {
                    transfer_edges1.set(g_input.gene_tree_numbering[i]);
                }
        }
    for (unsigned i = 0; i < sc2.transfer_edges.size(); ++i)
        {
            if (sc2.transfer_edges[i])
                {
                    transfer_edges2.set(g_input.gene_tree_numbering[i]);
                }
        }
    if (transfer_edges1 < transfer_edges2)
        {
            return true;
        }
    if (transfer_edges1 > transfer_edges2)
        {
            return false;
        }

    dynamic_bitset<> duplications1(sc1.duplications.size());
    dynamic_bitset<> duplications2(sc1.duplications.size());
    for (unsigned i = 0; i < sc1.duplications.size(); ++i)
        {
            if (sc1.duplications[i])
                {
                    duplications1.set(g_input.gene_tree_numbering[i]);
                }
        }
    for (unsigned i = 0; i < sc2.duplications.size(); ++i)
        {
            if (sc2.duplications[i])
                {
                    duplications2.set(g_input.gene_tree_numbering[i]);
                }
        }
    if (duplications1 < duplications2)
        {
            return true;
        }
    if (duplications1 > duplications2)
        {
            return false;
        }

    return false;
}



bool
scenario_less_lex(const Scenario& sc1, const Scenario &sc2)
{
    if (sc1.duplications < sc2.duplications)
        {
            return true;
        }
    if (sc1.duplications == sc2.duplications && sc1.transfer_edges < sc2.transfer_edges)
        {
            return true;
        }
    return false;
}


bool
operator==(const Scenario &sc1, const Scenario &sc2)
{
    return sc1.duplications == sc2.duplications &&
        sc1.transfer_edges == sc2.transfer_edges;
}


Scenario::Scenario(unsigned size) :
    duplications(size),
    transfer_edges(size)
{
}



BacktrackElement::BacktrackElement() :
    outside_sibling(tree_type::NONE),
    outside_ancestor(tree_type::NONE),
    min_transfers(0),
    scenarios_below_needed(false),
    scenarios_at_needed(false)
{
}



BacktrackElementDated::BacktrackElementDated() :
    scenarios_needed(false)
{
}
