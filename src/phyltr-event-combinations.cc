//*****************************************************************************
// phyltr-event-combinations.cc
//
// Implementation of the dynamic programming algorithm for finding all
// optimal combinations of the number of duplications and transfers that are
// possible under different cost settings.
//*****************************************************************************

#include <iostream>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <utility>
#include <vector>
#include <iterator>
#include <limits>
#include <iomanip>

#include <boost/rational.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <boost/foreach.hpp>

#include "common.hh"
#include "utils/binary-tree.hh"
#include <NHparser/NHparser.h>



namespace           po = boost::program_options;

using namespace     std;
using               boost::rational;
using               boost::multi_array;
using               boost::shared_ptr;



typedef Binary_tree<double>             tree_type;
typedef tree_type::vid_t                vid_t;
typedef pair<int, int>                  combo_t; // first=#dups, second=#transfers
typedef vector<combo_t>                 combo_seq_t;
typedef rational<long>                  rational_t;
typedef pair<rational_t, rational_t>    interval_t;

const string PROGRAM_NAME = "phyltr-event-combinations";
const string USAGE = "Usage: " + PROGRAM_NAME + 
    " [OPTION]... SPECIES-TREE-FILE GENE-TREE-FILE MAP-FILE\n"
    "\n"
    "Output consists of one line for each possible event count pair.\n"
    "Each line consists of five columns: The first column contains\n"
    "a '+', '.', or '-' depending on whether the event count pair is\n"
    "optimal in an interval, a single point, or is never optimal,\n"
    "respectively. The other columns comtain the number of duplications,\n"
    "the number of transfers, and the lower and upper limits of the\n" 
    "duplication cost in which the event count pair is optimal.\n";


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
//      files.
//
// g_below
// g_outside
// 
//      The DP matrices. For a gene tree vertex u and a species tree
//      vertex x, g_below[u][x] is the set of all minimal event
//      combinations of reconciling G_u with S_x, and g_outside[u][x]
//      is the set of all minimal event combinations when reconciling
//      G_u with S_y for any incomparable species tree vertex y.
//*****************************************************************************
multi_array<combo_seq_t, 2>           g_below;
multi_array<combo_seq_t, 2>           g_outside;

po::variables_map                   g_options;
struct ProgramInput {
    string                  species_tree_fname;
    string                  gene_tree_fname;
    string                  sigma_fname;
    shared_ptr<tree_type>   species_tree;
    shared_ptr<tree_type>   gene_tree;
    vector<unsigned>        sigma;
} g_input;


//*****************************************************************************
// print_error()
//
// Prints message to stderr with a suffix identifying the program.
//*****************************************************************************
void print_error(const char *);

//*****************************************************************************
// parse_options()
//
// Parses the command-line options and fills in values from the
// command-line into g_input. Returns an options_description suitable
// for printing with the program's help message. Throws po::error if
// an error is detected
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
// when getting/reading the data members of g_input. If any errors are
// detected, a message is printed to stderr and false is
// returned. Otherwise, true is returned.
//*****************************************************************************
bool check_input();

//*****************************************************************************
// read_species_tree()
// read_gene_tree()
//
// The functions read the files whose name is given in
// g_input.species_tree_fname and g_input.gene_tree_fname, and
// attempt to construct binary trees. If errors are detected, an
// apropriate error message is written to stderr and false is
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
//
// Runs the dynamic programming algorithm.
//*****************************************************************************
void dp_algorithm();

//*****************************************************************************
// compute_below()
// compute_outside()
//
// These are helper functions used by dp_algorithm.
//*****************************************************************************
void compute_below(vid_t u, vid_t x);
void compute_outside(vid_t u, vid_t x);

//*****************************************************************************
// compute_intervals()
//
// Computes, for each event pair in 'combos', the interval for the
// duplication cost in which the cost of the event pair is less than
// or equal to that of any other event pair in 'combos'. It assumes
// that the sum of the duplication cost and transfer cost is one. The
// results are returned in the vector 'results'. The elements of
// combos must be sorted by the number of duplications (i.e., the
// first coordinate of the pairs).
//*****************************************************************************
void compute_intervals(const combo_seq_t &combos,
                       vector<interval_t> &result);


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

    // Run the dynamic programming algorithm.
    dp_algorithm();



    // Compute the cost intervals in which each minimal
    // event-combination is optimal.
    vector<interval_t> intervals;
    compute_intervals(g_below[0][0], intervals);



    // Output the results.
    vector<interval_t>::iterator iter = intervals.begin();
    BOOST_FOREACH (combo_t &combo, g_below[0][0])
    {
        int duplications = combo.first;
        int transfers = combo.second;
        rational_t lower_r = iter->first;
        rational_t upper_r = iter->second;
        double lower_d = boost::rational_cast<double>(lower_r);
        double upper_d = boost::rational_cast<double>(upper_r);
        
        if (lower_r < upper_r)
        {
            cout << "+";
        }
        else if (lower_r == upper_r)
        {
            cout << ".";
        }
        else
        {
            cout << "-";
        }

        rational_t::int_type common_denominator =
            boost::lcm(lower_r.denominator(), upper_r.denominator());
        rational_t::int_type lower_numerator = lower_r.numerator() *
            common_denominator / lower_r.denominator();
        rational_t::int_type upper_numerator = upper_r.numerator() *
            common_denominator / upper_r.denominator();

        cout << setw(6) << duplications
             << setw(6) << transfers

             << setw(6) << lower_numerator
             << "/"
             << setw(6) << left << common_denominator << right

             << setw(6) << upper_numerator
             << "/"
             << setw(6) << left << common_denominator << right

             << setw(8) << setprecision(3) << lower_d
             << setw(8) << setprecision(3) << upper_d;
        
        cout << "\n";
        ++iter;
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
        ;

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
    const tree_type &S = *g_input.species_tree;
    const tree_type &G = *g_input.gene_tree;

    // Allocate memory for the matrices.
    g_below.resize(boost::extents[G.size()][S.size()]);
    g_outside.resize(boost::extents[G.size()][S.size()]);

    // The algorithm itself is described in an article
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

    if (G.is_leaf(u))
    {
        if (S.descendant(sigma[u], x))
        {
            g_below[u][x].push_back(combo_t(0, 0));
        }
    }
    else
    {
        // min_transfers[d] is the minimum number of transfers
        // required for any possible combination of events with d
        // duplications.
        //
        // We use G.size() as a placeholder for 'none' since there can
        // never be G.size() or more events in any scenario.
        static vector<int> min_transfers(G.size());

        vid_t v = G.left(u);
        vid_t w = G.right(u);

        min_transfers.clear();
        min_transfers.resize(G.size(), G.size());
        
        // Event pairs when u is a duplication at x.
        BOOST_FOREACH (combo_t cv, g_below[v][x])
        {
            BOOST_FOREACH (combo_t cw, g_below[w][x])
            {
                int dups = cv.first + cw.first + 1;
                int transfers = cv.second + cw.second;
                min_transfers[dups] = min(transfers, min_transfers[dups]);
            }
        }
        
        // Event pairs when u is a transfer with w the transferred vertex
        BOOST_FOREACH (combo_t cv, g_below[v][x])
        {
            BOOST_FOREACH (combo_t cw, g_outside[w][x])
            {
                int dups = cv.first + cw.first;
                int transfers = cv.second + cw.second + 1;
                min_transfers[dups] = min(transfers, min_transfers[dups]);
            }
        }
        
        // Event pairs when u is a transfer with v the transferred vertex
        BOOST_FOREACH (combo_t cv, g_outside[v][x])
        {
            BOOST_FOREACH (combo_t cw, g_below[w][x])
            {
                int dups = cv.first + cw.first;
                int transfers = cv.second + cw.second + 1;
                min_transfers[dups] = min(transfers, min_transfers[dups]);
            }
        }

        if (!S.is_leaf(x))
        {
            vid_t y = S.left(x);
            vid_t z = S.right(x);
            
            // Event pairs when u is a speciation with v below y and w below z.
            BOOST_FOREACH (combo_t cv, g_below[v][y])
            {
                BOOST_FOREACH (combo_t cw, g_below[w][z])
                {
                    int dups = cv.first + cw.first;
                    int transfers = cv.second + cw.second;
                    min_transfers[dups] = min(transfers, min_transfers[dups]);
                }
            }

            // Event pairs when u is a speciation with v below z and w below y.
            BOOST_FOREACH (combo_t cv, g_below[v][z])
            {
                BOOST_FOREACH (combo_t cw, g_below[w][y])
                {
                    int dups = cv.first + cw.first;
                    int transfers = cv.second + cw.second;
                    min_transfers[dups] = min(transfers, min_transfers[dups]);
                }
            }

            // Event pairs when u is placed strictly below x.
            BOOST_FOREACH (combo_t c, g_below[u][y])
            {
                min_transfers[c.first] = min(c.second, min_transfers[c.first]);
            }
            BOOST_FOREACH (combo_t c, g_below[u][z])
            {
                min_transfers[c.first] = min(c.second, min_transfers[c.first]);
            }
        }

        // Insert the resulting minimal event pairs into g_below.
        // Note how this procedures ensures that the event pairs are
        // sorted by the number of duplications.
        int prev_transfers = G.size();
        for (unsigned i = 0; i < min_transfers.size(); ++i)
        {
            // Inserting only event pairs that are incomparable.
            if (min_transfers[i] < prev_transfers)
            {
                g_below[u][x].push_back(combo_t(i, min_transfers[i]));
                prev_transfers = min_transfers[i];
            }
        }
    }
}



void
compute_outside(vid_t u, vid_t x)
{
    const tree_type &S = *g_input.species_tree;
    const tree_type &G = *g_input.gene_tree;

    // Cannot place u outside the root of S.
    if (x == 0)
        return;

    vid_t x_parent = S.parent(x);
    vid_t x_sibling =
        S.left(x_parent) == x ? S.right(x_parent) : S.left(x_parent);

    static vector<int> min_transfers(G.size());
    min_transfers.clear();
    min_transfers.resize(G.size(), G.size());
    
    BOOST_FOREACH (combo_t c, g_outside[u][x_parent])
    {
        min_transfers[c.first] = c.second;
    }

    BOOST_FOREACH (combo_t c, g_below[u][x_sibling])
    {
        min_transfers[c.first] = min(c.second, min_transfers[c.first]);
    }
    
    int prev_transfers = G.size();
    for (unsigned i = 0; i < min_transfers.size(); ++i)
    {
        if (min_transfers[i] < prev_transfers)
        {
            g_outside[u][x].push_back(combo_t(i, min_transfers[i]));
            prev_transfers = min_transfers[i];
        }
    }
}


void
compute_intervals(const combo_seq_t &combos,
                  vector<interval_t> &result)
{
    if (combos.size() == 0)
    {
        return;
    }
    else if (combos.size() == 1)
    {
        interval_t interval(rational_t(0), rational_t(1));
        result.push_back(interval);
        return;
    }
    else if (combos.size() == 2)
    {
        int duplication_diff = combos[1].first - combos[0].first;
        int transfer_diff = combos[0].second - combos[1].second;
        rational_t delta(transfer_diff, duplication_diff + transfer_diff);
        result.push_back(interval_t(delta, rational_t(1)));
        result.push_back(interval_t(rational_t(0), delta));
    }
    else
    {
        vector< vector<rational_t> > delta(combos.size());
        for (unsigned i = 1; i < delta.size(); ++i)
        {
            delta[i].resize(i);
            for (unsigned j = 0; j < i; ++j)
            {
                int duplication_diff = combos[i].first - combos[j].first;
                int transfer_diff = combos[j].second - combos[i].second;
                delta[i][j].assign(transfer_diff,
                                   duplication_diff + transfer_diff);
            }
        }

        for (unsigned i = 0; i < combos.size(); ++i)
        {
            rational_t lower(0);
            rational_t upper(1);

            for (unsigned j = i + 1; j < combos.size(); ++j)
            {
                lower = max(delta[j][i], lower);
            }
            for (unsigned j = 0; j < i; ++j)
            {
                upper = min(delta[i][j], upper);
            }
            result.push_back(interval_t(lower, upper));
        }
    }
}
