/*
 * This file contains source code for a program that generates binary
 * species trees from a birth-death process. Given the two paramters,
 * birth- and death-rate, the program creates a binary tree with times
 * on the edges. The output is in newick format.
 */

#include "common.hh"

#include <sys/time.h>
#include <boost/shared_ptr.hpp>
#include <boost/random.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <cstdlib>               // At least for EXIT_SUCCESS and EXIT_FAILURE
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>

//============================================================================
//             Namespace declarations and using directives.
//============================================================================

using namespace std;

//============================================================================
//                   Typedefs and class declarations
//============================================================================

struct Node;

typedef boost::shared_ptr<Node> Node_ptr;

struct Node {
    double time;
    Node_ptr parent, left, right;
    
    Node(double t, Node_ptr p, Node_ptr l, Node_ptr r)
        : time(t), parent(p), left(l), right(r)
    {
    }
};


//============================================================================
//                   Global Constants and variables.
//============================================================================

const string PROG_NAME = "phyltr-gen-stree";
const string USAGE = "Usage: " + PROG_NAME +
    " [OPTION]... TIME BIRTH-RATE DEATH-RATE";
const unsigned LENGTH_PRECISION = 4;    // The precision with which
                                        // times on the trees will be
                                        // printed.

boost::program_options::variables_map g_options;

unsigned max_attempts = 1000000;             // The maximum number of times the
                                        // process will be run if the
                                        // species goes extinct.
//=============================================================================
//                        Function declarations
//=============================================================================

/*
 * print_tree()
 *
 * Prints the tree rooted at 'root' to cout. If 'output_lengths' is
 * true, then the lenghts printed for each edge is the time of the
 * node minus the time of the parent's time. Otherwise, if
 * 'output_lengths' is false, then the time of each node is printed.
 */
void print_tree(Node_ptr root, bool output_lengths = true);

//=============================================================================
//         Template and inline function and member definitions.
//=============================================================================


//=============================================================================
//                                main()
//=============================================================================

int 
main(int argc, char *argv[])
{
    /* Seed the random number generator(s). */
    init_rand(g_generator);

    namespace po = boost::program_options;

    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    double lambda, mu;  // The birth- and death-rates respectively
    double T;           // The length of time of the process

    try
        {
            /* Declare visible options */
            visible_opts.add_options()
                ("help,h", "display this help and exit")
                ("start-with-speciation,s",
                 "If set, the birth-death process starts with a speciation")
                ("max-attempts,a",
                 po::value<unsigned>(&max_attempts)->default_value(10000),
                 "The maximum number of times the process will run if the "
                 "species goes extinct during.")
                ;
            
            /* Declare positional options. */
            hidden_opts.add_options()
                ("time", po::value<double>(&T))
                ("birth-rate", po::value<double>(&lambda))
                ("death-rate", po::value<double>(&mu))
                ;
            
            po::positional_options_description positional_options;
            positional_options.add("time", 1);
            positional_options.add("birth-rate", 1);
            positional_options.add("death-rate", 1);

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
            cout
                << USAGE << '\n'
                << visible_opts << '\n';
            exit(EXIT_SUCCESS);
        }
    
    /* Check that all required positional arguments are given. */
    if (g_options.count("birth-rate") == 0
        || g_options.count("death-rate") == 0
        || g_options.count("time") == 0)
        {
            cerr << PROG_NAME << ": Too few arguments.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* Check that max_attempt is not zero. */
    if (g_options["max-attempts"].as<unsigned>() == 0)
        {
            cerr << PROG_NAME << ": "
                 << "max-attempts was set to zero.\n";
            exit(EXIT_FAILURE);
        }

    /* Check that the rates and time are non-negative. */
    if (lambda < 0 || mu < 0)
        {
            cerr << PROG_NAME << ": negative rate\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    if (T < 0)
        {
            cerr << PROG_NAME << ": negative time\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /*
     * We now have everthing we need. Create the tree!
     */


    Node_ptr null((Node *)0);   // A constant to ease passing of zero-pointers.
    vector<Node_ptr> leaves;    // The current set of leaves at cur_time.
    Node_ptr root(new Node(T, null, null, null));

    int attempts = max_attempts;
    while (leaves.size() == 0 && attempts--)
        {
            /*
             * We will keep track of the root and the set of leaves at
             * cur_time during the process. Note that we keep the leaves at
             * time T and therefore we don't need to set their times at the
             * end of the process.
             */
            leaves.push_back(root);
            double cur_time = 0;        // Keeps track of the current time.
            // cur_time <= T at all times.

            /* If start-with-speciation is set, we want to start with a
               speciation at time 0. */
            if (g_options.count("start-with-speciation") > 0)
                {
                    Node_ptr np1(new Node(T, root, null, null));
                    Node_ptr np2(new Node(T, root, null, null));
                    root->time = 0;
                }


            /*
             * In each iteration the loop advances cur_time according to the
             * distribution of time-intervals (i.e., an exponential
             * distribution with the sum of the intensities as paramter) and
             * decide the event according to the intensities. We prune the
             * tree immediately upon a death event. Keep doing this until
             * cur_time surpasses T.
             */
            while (cur_time < T && !leaves.empty())
                {
                    double birth = lambda * leaves.size();
                    double death = mu * leaves.size();

                    if (birth + death > 0)
                        {
                            boost::exponential_distribution<double> rng_exp(birth + death);
                            double elapsed_time = rng_exp(g_rng_d);
                            cur_time = min(cur_time + elapsed_time, T);
                        }
                    else
                        {
                            cur_time = T;
                        }

                    /* If elapsed time continued passed T, we are done. */
                    if (cur_time == T)
                        break;

                    if (g_rng_d() < birth / (birth + death))    // birth
                        {
                            unsigned choice = g_rng_ui(leaves.size());
                            Node_ptr np1(new Node(T, leaves[choice], null, null));
                            Node_ptr np2(new Node(T, leaves[choice], null, null));
                            leaves[choice]->time = cur_time;
                            leaves[choice]->left = np1;
                            leaves[choice]->right = np2;
                            leaves[choice] = np1;
                            leaves.push_back(np2);
                        }
                    else                                        //death
                        {
                            if (leaves.size() == 1)
                                {
                                    leaves.pop_back();
                                    break;
                                }
                            unsigned choice = g_rng_ui(leaves.size());
                            Node_ptr np = leaves[choice];
                            Node_ptr parent = np->parent;
                            Node_ptr sibling = parent->left == np ? parent->right : parent->left;
                            sibling->parent = sibling->parent->parent;
                            if (sibling->parent)
                                {
                                    if (sibling->parent->left == parent)
                                        sibling->parent->left = sibling;
                                    else
                                        sibling->parent->right = sibling;
                                }
                            else
                                root = sibling;
                            swap(leaves[choice], leaves.back());
                            leaves.pop_back();
                        }
                }
        }

    /* Print error message if all species went extinct. (what else?) */
    if (leaves.empty())
        {
            cerr << PROG_NAME
                 << ": The species became extinct during all "
                 << max_attempts << " runs of the process\n";
            exit(EXIT_FAILURE);
        }
    
    /* Print the tree that is reachable from the root. */
    print_tree(root);
 
    return EXIT_SUCCESS;
}

//=============================================================================
//                     Helper function declarations
//=============================================================================

unsigned print_tree_r(Node_ptr node, unsigned next_leaf_id, bool output_lengths);

//=============================================================================
//                   Function and member definitions.
//=============================================================================

void
print_tree(Node_ptr root, bool output_lengths)
{
    unsigned next_leaf_id = 1;
    /* if node is a leaf. */
    if (!root->left)
        {
            cout << "s" << next_leaf_id;
            cout << ":" << setprecision(LENGTH_PRECISION) << root->time;
            cout << ";\n";
            
            return;
        }

    /* if node is not a leaf. */
    cout << "(";
    next_leaf_id = print_tree_r(root->left, next_leaf_id, output_lengths);
    cout << ", ";
    next_leaf_id = print_tree_r(root->right, next_leaf_id, output_lengths);
    cout << ")";
    cout << ":" << setprecision(LENGTH_PRECISION) << root->time;
    cout << ";\n";
}

unsigned
print_tree_r(Node_ptr node, unsigned next_leaf_id, bool output_lengths)
{
    /* if node is a leaf. */
    if (!node->left)
        {
            cout << "s" << next_leaf_id;
            if (output_lengths)
                cout << ":" << setprecision(LENGTH_PRECISION) << (node->time - node->parent->time);
            else
                cout << ":" << setprecision(LENGTH_PRECISION) << node->time;
            return next_leaf_id + 1;
        }

    /* if node is not a leaf. */
    cout << "(";
    next_leaf_id = print_tree_r(node->left, next_leaf_id, output_lengths);
    cout << ", ";
    next_leaf_id = print_tree_r(node->right, next_leaf_id, output_lengths);
    cout << ")";

    ostringstream os;
    if (output_lengths)
        os << ":" << setprecision(LENGTH_PRECISION) << showpoint << (node->time - node->parent->time);
    else
        os << ":" << setprecision(LENGTH_PRECISION) << showpoint << node->time;
    string s = os.str();

    cout << s;
    if (s[s.size() - 1] == '.')
    {
        cout << '0';
    }
    
    return next_leaf_id;
}
