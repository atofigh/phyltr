//*****************************************************************************
// phyltr-fpt.cc
//
// Implementation of the fixed-parameter-tractable (fpt) algorithm for
// enumerating DTL-scenarios reconciling a species tree and a gene tree.
//*****************************************************************************

#include <cstdlib>
#include <cerrno>
#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <stack>
#include <set>

#include <boost/program_options.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>

#include "common.hh"
#include "utils/binary-tree.hh"
extern "C" {
#include <NHparser.h>
}



namespace           po = boost::program_options;

using namespace     std;
using               boost::shared_ptr;
using               boost::dynamic_bitset;



typedef Binary_tree<double>     tree_type;
typedef tree_type::vid_t        vid_t;



const string PROGRAM_NAME = "phyltr-fpt";
const string USAGE =
    "Usage: " + PROGRAM_NAME + " [OPTION]... SPECIES-TREE-FILE GENE-TREE-FILE "
    "MAP-FILE MIN-COST MAX-COST";



//*****************************************************************************
// class Scenario
//
// A Scenario object is used to keep track of a final solution. Since
// the number of solutions may be large, the Scenario class has been
// optimized on memory. This is also why we do not keep final
// solutions in Candidate objects which are much more memory hungry.
//
// Throughout the program, an edge is represented by the vertex at its
// head, i.e., the vertex farthest away from the root.
// 
// The operator< is used to sort the scenarios for printing. It sorts
// first on cost, then the number of transfers, then on the number of
// losses, and lastly according to lexicographic order.
//*****************************************************************************
class Scenario {
public:
    dynamic_bitset<> duplications;
    dynamic_bitset<> transfer_edges;

    Scenario(unsigned size);
};

bool operator<(const Scenario &, const Scenario &);
ostream &operator<<(ostream &, const Scenario &);



//*****************************************************************************
// class Candidate
//
// A Candidate object represents a partial solution of a
// reconciliation. It keeps track of duplications, transfers, cost,
// s-moves, and the lca-mapping of the gene tree into the species tree
// given the events. It also keeps track of the forest that is
// obtained from the gene tree by removing transfer edges and
// contracting any vertices with only one child. Finally, it is also
// able to compute whether or not the candidate is elegant, i.e.,
// whether or not a duplication or transfer is unnecessary.
//
// Most member functions are self explanatory. Here are some notes
// that should be mentioned.
//
// It is important to know that the candidate is designed to behave
// intelligently with respect to duplications and transfers. Already
// during construction of a candidate, the object computes the forced
// duplications (i.e., the internal gene tree vertices that are mapped
// to species tree leaves). When marking an edge as a transfer edge
// and when marking a vertex as a duplication, any resulting forced
// duplications are also set.
//
// Note that is_elegant() assumes that no moves remain, i.e., it
// assumes that the candidate is final.
//
// The member functions parent(), left(), and right() are used to gain
// information about the gene tree forest described above. parent(u)
// is defined for all gene tree vertices u and returns the least
// proper ancestor of u that is not a transfer vertex (i.e., no
// outgoing edge is a transfer edge). left(u) and right(u) are only
// defined for non-transfer vertices, and return tree_type::NONE when
// u is a transfer vertex.
//*****************************************************************************
class Candidate {
public:
    class bad_transfer_exception : public exception {};
    class bad_duplication_exception : public exception {};

    Candidate();

    void    set_transfer_edge(vid_t);
    void    set_duplication(vid_t);
    bool    is_transfer_edge(vid_t) const;
    bool    is_duplication(vid_t)   const;
    double  cost()                  const;
    bool    is_elegant()            const;
    vid_t   get_s_move()            const;
    vid_t   lambda(vid_t)           const;
    vid_t   parent(vid_t)           const;
    vid_t   left(vid_t)             const;
    vid_t   right(vid_t)            const;

private:
    dynamic_bitset<>        duplications_;
    dynamic_bitset<>        transfer_edges_;
    double                  cost_;
    vector<vid_t>           lambda_;
    mutable vector<vid_t>   s_moves_;
    vector<vid_t>           P_;
    vector<vid_t>           left_;
    vector<vid_t>           right_;

    bool is_s_move_(vid_t) const;
    void compute_highest_mapping_(vector<vid_t> &) const;

    friend ostream &operator<<(ostream &, const Candidate &);
};
ostream &operator<<(ostream &, const Candidate &);



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
//*****************************************************************************
po::variables_map                   g_options;
struct ProgramInput {
    string                  species_tree_fname;
    string                  gene_tree_fname;
    string                  sigma_fname;
    double                  min_cost;
    double                  max_cost;
    shared_ptr<tree_type>   species_tree;
    shared_ptr<tree_type>   gene_tree;
    vector<unsigned>        sigma;
    vector<unsigned>        gene_tree_numbering;
    double                  duplication_cost;
    double                  transfer_cost;
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
// when getting/reading the individual data members of g_input. 
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
// g_input.sigma_fname. If any erros are detected, an apropriate error
// message is written to stderr and fase is returned. Otherwise, true
// is returned.
//*****************************************************************************
bool read_sigma();


//*****************************************************************************
// fpt_algorithm()
//
// Runs the fixed-parameter-tractable algorithm. For details see the
// relevant article cited in the beginning of the file. The vector
// passed to the function is filled with Scenario options.
// //*****************************************************************************
void fpt_algorithm(vector<Scenario> &);


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

    

    // Run the algorithm and get all final, elegant, candidates
    vector<Scenario> scenarios;
    fpt_algorithm(scenarios);

    // Sort and print the scenarios
    sort(scenarios.begin(), scenarios.end());
    BOOST_FOREACH (Scenario &sc, scenarios)
        {
            cout << sc << "\n";
        }



    return EXIT_SUCCESS;
}



void
fpt_algorithm(vector<Scenario> &scenarios)
{
    typedef shared_ptr<Candidate> cand_ptr;
    
    scenarios.clear();

    // We will do a depth first search, so we need a stack.
    stack<cand_ptr> cand_stack;

    // Push the initial candidate onto the stack. Note that the
    // initial candidate may have some duplications set already!
    cand_ptr initial_candidate(new Candidate());
    if (initial_candidate->cost() <= g_input.max_cost)
        {
            cand_stack.push(initial_candidate);
        }

    // Do the depth-first search.
    unsigned gene_tree_size = g_input.gene_tree->size();
    while (!cand_stack.empty())
        {
            cand_ptr cp = cand_stack.top(); cand_stack.pop();


            vid_t s_move = cp->get_s_move();
            double cp_cost = cp->cost();

            if (s_move != tree_type::NONE)
                {
                    // Resvole the s-move in three ways.
                    if (cp_cost + g_input.duplication_cost <= g_input.max_cost)
                        {
                            cand_ptr cp1(new Candidate(*cp));
                            cp1->set_duplication(cp->parent(s_move));
                            if (cp1->cost() <= g_input.max_cost)
                                {
                                    cand_stack.push(cp1);
                                }
                        }

                    if (cp_cost + g_input.transfer_cost <= g_input.max_cost)
                        {
                            cand_ptr cp2(new Candidate(*cp));
                            cand_ptr cp3(new Candidate(*cp));

                            cp2->set_transfer_edge(g_input.gene_tree->left(s_move));
                            cp3->set_transfer_edge(g_input.gene_tree->right(s_move));
                            
                            if (cp2->cost() <= g_input.max_cost)
                                {
                                    cand_stack.push(cp2);
                                }
                            if (cp3->cost() <= g_input.max_cost)
                                {
                                    cand_stack.push(cp3);
                                }
                        }
                }
            else
                {
                    // Insert elegant final candidates with cost in
                    // the given range into the return-vector.
                    if (cp->cost() >= g_input.min_cost &&
                        cp->cost() <= g_input.max_cost &&
                        cp->is_elegant())
                        {
                            Scenario sc(gene_tree_size);
                            for (vid_t u = 0; u < gene_tree_size; ++u)
                                {
                                    sc.duplications[u] = cp->is_duplication(u);
                                    sc.transfer_edges[u] = cp->is_transfer_edge(u);
                                }
                            scenarios.push_back(sc);
                        }
                }
        }
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
         po::value<double>(&g_input.transfer_cost)->default_value(1.0),
         "Cost of transfer events")
        ("duplication-cost,d",
         po::value<double>(&g_input.duplication_cost)->default_value(1.0),
         "Cost of duplication events")
        ;

    // Declare the positional options
    hidden_opts.add_options()
        ("species-tree-file",po::value<string>(&g_input.species_tree_fname))
        ("gene-tree-file", po::value<string>(&g_input.gene_tree_fname))
        ("map-file", po::value<string>(&g_input.sigma_fname))
        ("min-cost", po::value<double>(&g_input.min_cost))
        ("max-cost", po::value<double>(&g_input.max_cost))
        ;

    po::positional_options_description positional_options;
    positional_options.add("species-tree-file", 1);
    positional_options.add("gene-tree-file", 1);
    positional_options.add("map-file", 1);
    positional_options.add("min-cost", 1);
    positional_options.add("max-cost", 1);

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
        g_options.count("map-file") +
        g_options.count("min-cost") +
        g_options.count("max-cost") != 5)
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
    double cost1, cost2;
    cost1 =
        sc1.transfer_edges.count() * g_input.transfer_cost +
        sc1.duplications.count() * g_input.duplication_cost;
    cost2 =
        sc2.transfer_edges.count() * g_input.transfer_cost +
        sc2.duplications.count() * g_input.duplication_cost;
    if (cost1 < cost2)
        {
            return true;
        }
    if (cost1 > cost2)
        {
            return false;
        }

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



Scenario::Scenario(unsigned size) :
    duplications(size),
    transfer_edges(size)
{
}



Candidate::Candidate() :
    duplications_(g_input.gene_tree->size()),
    transfer_edges_(g_input.gene_tree->size()),
    cost_(0.0),
    lambda_(g_input.gene_tree->size()),
    P_(g_input.gene_tree->size()),
    left_(g_input.gene_tree->size()),
    right_(g_input.gene_tree->size())
{
    const tree_type &G = *g_input.gene_tree;
    const tree_type &S = *g_input.species_tree;
    
    compute_lambda(S, G, g_input.sigma, transfer_edges_, lambda_);

    for (vid_t u = 0; u < G.size(); ++u)
        {
            P_[u] = G.parent(u);
            left_[u] = G.left(u);
            right_[u] = G.right(u);
        }
    
    // Find forced duplications, i.e., internal gene tree vertices
    // that are mapped to leaves of S.
    for (vid_t u = 0; u < G.size(); ++u)
        {
            if (!G.is_leaf(u) && S.is_leaf(lambda_[u]))
                {
                    duplications_.set(u);
                    cost_ += g_input.duplication_cost;
                }
        }

    // Find all s-moves.
    for (vid_t u = 0; u < G.size(); ++u)
        {
            if (is_s_move_(u))
                {
                    s_moves_.push_back(u);
                }
        }
}



void
Candidate::set_transfer_edge(vid_t u)
{
    const tree_type &G = *g_input.gene_tree;
    const tree_type &S = *g_input.species_tree;

    if (u == 0)
        {
            throw bad_transfer_exception();
        }

    vid_t parent_u = G.parent(u);
    vid_t sibling_u = G.left(parent_u) == u ? G.right(parent_u) : G.left(parent_u);

    // parent_u must be an anchor
    if (lambda_[parent_u] == lambda_[u] ||
        lambda_[parent_u] == lambda_[sibling_u])
        {
            throw bad_transfer_exception();
        }

    // Set the transfer and update the cost.
    transfer_edges_.set(u);
    cost_ += g_input.transfer_cost;

    // update P_, left_, and right_
    vid_t v = u == G.left(parent_u) ? left_[parent_u] : right_[parent_u];
    vid_t w = u == G.left(parent_u) ? right_[parent_u] : left_[parent_u];
    for (vid_t a = v; a != parent_u; a = G.parent(a))
        {
            P_[a] = tree_type::NONE;
        }
    for (vid_t a = w; a != parent_u; a = G.parent(a))
        {
            P_[a] = P_[parent_u];
        }
    if (P_[parent_u] != tree_type::NONE)
        {
            if (left_[P_[parent_u]] == parent_u)
                {
                    left_[P_[parent_u]] = w;
                }
            else
                {
                    right_[P_[parent_u]] = w;
                }
        }
    left_[parent_u] = tree_type::NONE;
    right_[parent_u] = tree_type::NONE;


    // update lambda and s_moves, and find the vertices with
    // new positions that are forced duplications.
    lambda_[parent_u] = lambda_[sibling_u];
    vid_t last_updated_vertex = parent_u;
    for (vid_t a = G.parent(parent_u); a != tree_type::NONE; a = G.parent(a))
        {
            // Compute the new placement of a
            vid_t old_lambda = lambda_[a];
            vid_t new_lambda = S.lca(lambda_[G.left(a)], lambda_[G.right(a)]);
            if (is_transfer_edge(G.left(a)))
                {
                    new_lambda = lambda_[G.right(a)];
                }
            else if (is_transfer_edge(G.right(a)))
                {
                    new_lambda = lambda_[G.left(a)];
                }
            lambda_[a] = new_lambda;

            if (old_lambda == new_lambda)
                break;

            last_updated_vertex = a;

            // If a is not a transfer vertex and not a duplication
            if (left_[a] != tree_type::NONE && !duplications_[a])
                {
                    // Is 'a' a forced duplication?
                    if (S.is_leaf(new_lambda) ||
                        (lambda_[left_[a]] == new_lambda &&
                         duplications_[left_[a]]) ||
                        (lambda_[right_[a]] == new_lambda &&
                         duplications_[right_[a]]))
                        {
                            duplications_.set(a);
                            cost_ += g_input.duplication_cost;
                        }
                    // Otherwise its children are potential s-moves
                    else
                        {
                            s_moves_.push_back(left_[a]);
                            s_moves_.push_back(right_[a]);
                        }
                }
        }

    if (P_[last_updated_vertex] != tree_type::NONE)
        {
            s_moves_.push_back(P_[last_updated_vertex]);
        }
}



void
Candidate::set_duplication(vid_t u)
{
    // u must not be a duplication or transfer vertex
    if (is_duplication(u) || left_[u] == tree_type::NONE)
        {
            throw bad_duplication_exception();
        }

    duplications_.set(u);
    cost_ += g_input.duplication_cost;

    // Find any forced duplications as a result of u becoming a duplication.
    for (vid_t v = P_[u]; v != tree_type::NONE; v = P_[v])
        {
            if (!duplications_[v] && lambda_[v] == lambda_[u])
                {
                    duplications_.set(v);
                    cost_ += g_input.duplication_cost;
                }
            else
                {
                    break;
                }
        }
}



bool
Candidate::is_transfer_edge(vid_t u) const
{
    return transfer_edges_.test(u);
}



bool
Candidate::is_duplication(vid_t u) const
{
    return duplications_.test(u);
}



double
Candidate::cost() const
{
    return cost_;
}



bool
Candidate::is_elegant() const
{
    const tree_type &G = *g_input.gene_tree;
    const tree_type &S = *g_input.species_tree;

    vector<vid_t> highest(G.size());

    compute_highest_mapping_(highest);
    
    // Check that the children of duplications are mapped by lambda to
    // comparable species tree vertices. Otherwise, the duplication is
    // unnecessary.
    for (vid_t d = duplications_.find_first();
         d != duplications_.npos;
         d = duplications_.find_next(d))
        {
            vid_t v = G.left(d);
            vid_t w = G.right(d);
            if (!S.descendant(lambda_[v], lambda_[w]) && 
                !S.descendant(lambda_[w], lambda_[v]))
                return false;
        }

    // Check if the parents of transfer vertices can be mapped high
    // enough so that the transfer can be converted to a
    // speciation. If so, the transfer is unnecessary.
    for (vid_t v = transfer_edges_.find_first();
         v != transfer_edges_.npos;
         v = transfer_edges_.find_next(v))
        {
            // Let (u, v) be the transfer edge we are considering, let
            // pu = p(u), and x = lca{lambda_[u], lambda_[v]}
            vid_t u = G.parent(v);
            vid_t pu = G.parent(u);
            vid_t x = S.lca(lambda_[u], lambda_[v]);
            
            // The root of G is always an unnecessary transfer vertex.
            if (u == 0)
                return false;

            // If p(u) is a speciation and x is a proper descendant of
            // highest[p(u)] = lambda_[p(u)], then the transfer is
            // unnecessary.
            if (!is_duplication(pu) &&
                !is_transfer_edge(G.left(pu)) &&
                !is_transfer_edge(G.right(pu)) &&
                S.descendant(x, lambda_[pu]) && 
                x != lambda_[pu])
                {
                    return false;
                }
            
            // If p(u) is not a speciation, then it is enough for x to
            // be a descendant of highest[p(u)] for the transfer to be
            // unnecessary.
            if ((is_duplication(pu) ||
                 is_transfer_edge(G.left(pu)) ||
                 is_transfer_edge(G.right(pu))) &&
                S.descendant(x, highest[pu]))
                {
                    return false;
                }
        }

    return true;
}



vid_t
Candidate::get_s_move() const
{
    // Check vertices in s_move and find one that really is an s-move.
    while (!s_moves_.empty() && !is_s_move_(s_moves_.back()))
        {
            s_moves_.pop_back();
        }
    
    return s_moves_.empty() ? tree_type::NONE : s_moves_.back();
}


vid_t
Candidate::lambda(vid_t u) const
{
    return lambda_[u];
}



vid_t
Candidate::parent(vid_t u) const
{
    return P_[u];
}



vid_t
Candidate::left(vid_t u) const
{
    return left_[u];
}



vid_t
Candidate::right(vid_t u) const
{
    return right_[u];
}



bool
Candidate::is_s_move_(vid_t u) const
{
    return
        !g_input.gene_tree->is_leaf(u) &&
        !is_duplication(u) &&
        P_[u] != tree_type::NONE &&
        lambda_[g_input.gene_tree->left(u)] != lambda_[u] && 
        lambda_[g_input.gene_tree->right(u)] != lambda_[u] &&
        lambda_[P_[u]] == lambda_[u] &&
        !is_duplication(P_[u]);
}



void 
Candidate::compute_highest_mapping_(vector<vid_t> &highest) const
{
    const tree_type &G = *g_input.gene_tree;
    const tree_type &S = *g_input.species_tree;
    const vector<vid_t> &sigma = g_input.sigma;

    highest.resize(G.size());

    // We define a function C(x, y) : V(S) x V(S) -> V(S). y must be a
    // proper descendant of x in the species tree. The function
    // returns the unique child of x that is an ancestor of y.
    struct {
        vid_t operator()(vid_t x, vid_t y) {
            vid_t left = g_input.species_tree->left(x);
            vid_t right = g_input.species_tree->right(x);
            return  g_input.species_tree->descendant(y, left) ? left : right;
        }
    } C;

    // First, take care of the root of G.
    if (!is_duplication(0) &&
        !is_transfer_edge(G.left(0)) &&
        !is_transfer_edge(G.right(0))) // if root is a speciation
        {
            highest[0] = lambda_[0];
        }
    else if (is_duplication(0))
        {
            highest[0] = 0;
        }
    else // If the root is a transfer vertex. 
        {
            // Let v be the transfered child of the root.
            vid_t v = is_transfer_edge(G.left(0)) ? G.left(0) : G.right(0);
            highest[0] = C(S.lca(lambda_[0], lambda_[v]), lambda_[0]);
        }
        
    // Next, take care of the rest of the vertices from the root and down.
    for (vid_t u = G.preorder_next(0);
         u != G.NONE;
         u = G.preorder_next(u))
        {
            if (u == 0)
                continue;

            vid_t pu = G.parent(u);
            vid_t x = lambda_[pu];
            vid_t y = lambda_[u];

            if (G.is_leaf(u))
                {
                    highest[u] = sigma[u];
                }
            else if (!is_duplication(u) &&
                     !is_transfer_edge(G.left(u)) &&
                     !is_transfer_edge(G.right(u)))
                {
                    highest[u] = lambda_[u];
                }
            else // If u is a duplication or a transfer vertex.
                {
                    // Let z be the highest possible mapping of u when
                    // considering only p(u). If p(u) is a duplication
                    // or if p(u) is a transfer but u is not the
                    // transfered vertex, z = highest[pu]. Otherwise,
                    // if p(u) is a speciation, z = C(x, y), and if
                    // u is the transfered vertex, then z = C(lca(x, y), y)
                    vid_t z = highest[pu];
                    if (!is_duplication(pu) &&
                        !is_transfer_edge(G.left(pu)) &&
                        !is_transfer_edge(G.right(pu))) // If pu is speciation.
                        {
                            z = C(x, y);
                        }
                    else if (is_transfer_edge(u))
                        {
                            z = C(S.lca(x, y), y);
                        }

                    // Let z_prime be the highest possible mapping of
                    // u when considering its children only. z_prime
                    // is the root of x unless u is a transfer. In
                    // that case, if v is the transferred child,
                    // z_prime = C(lca(lambda_[u], lambda_[v]), lambda_[u]).
                    vid_t z_prime = 0;
                    if (is_transfer_edge(G.left(u)) ||
                        is_transfer_edge(G.right(u)))
                        {
                            // Let v be the transferred child of u.
                            vid_t v = is_transfer_edge(G.left(u)) ? G.left(u) : G.right(u);
                            
                            z_prime = C(S.lca(y, lambda_[v]), y);
                        }
                    // Since z and z_prime are both ancestors of
                    // lambda_[u], we know that they are
                    // comparable. The one that is minimal in S is
                    // then the highest possible mapping of u.
                    highest[u] = S.descendant(z, z_prime) ? z : z_prime;
                }
        }
}


ostream &
operator<<(ostream &out, const Candidate &c)
{
    out << "duplications:\t" << c.duplications_ << "\n";
    out << "transfers:\t" << c.transfer_edges_ << "\n";
    out << "potential smoves:\t";
    copy(c.s_moves_.begin(), c.s_moves_.end(),
         ostream_iterator<vid_t>(out, " "));
    out << "\n";
    out << "real smoves:\t\t";
    for (vid_t u = 0; u < g_input.gene_tree->size(); ++u)
        {
            if (c.is_s_move_(u))
                out << u << " ";
        }
    out << "\n";

    return out;
}
