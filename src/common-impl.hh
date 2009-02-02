#ifndef COMMON_IMPL_HH
#define COMMON_IMPL_HH

#include <iterator>
#include <algorithm>
#include <ostream>
#include <fstream>
#include <stdexcept>
#include <sys/time.h>

template<typename T>
void 
create_binary_tree(NHNode *root, Binary_tree<T> &tree,
                   typename Binary_tree<T>::vid_t ext)
{
    if (root->children == 0)
        return;

    char *left_label = root->children->node->label;
    char *right_label = root->children->next->node->label;
    double left_length = root->children->node->branchLength;
    double right_length = root->children->next->node->branchLength;
    tree.add_children(ext,
                      left_length <= 0 ? T(1) : T(left_length),
                      right_length <= 0 ? T(1) : T(right_length),
                      left_label == 0 ? "" : left_label,
                      right_label == 0 ? "" : right_label);

    create_binary_tree(root->children->node, tree, tree.left(ext));
    create_binary_tree(root->children->next->node, tree, tree.right(ext));
}

template<typename T>
void
create_gene_species_map(const Binary_tree<T> &species_tree, 
                        const Binary_tree<T> &gene_tree, 
                        std::string map_filename,
                        std::vector<typename Binary_tree<T>::vid_t> &sigma)
{
    using namespace std;
    typedef typename Binary_tree<T>::vid_t vid_t;

    vector<string> map_file_content;
    
    ifstream map_file(map_filename.c_str());
    copy(istream_iterator<string>(map_file), istream_iterator<string>(),
         back_inserter(map_file_content));

    /* Make sure there are even number of strings in map file. */
    if (map_file_content.size() == 0 || map_file_content.size() % 2 != 0)
        {
            throw logic_error("error reading map file.");
        } 

    /* Create a map from gene name to species name. */
    map<string, string> str_sigma;
    for (unsigned i = 0; i < map_file_content.size(); i += 2)
        str_sigma[map_file_content[i]] = map_file_content[i+1];
    
    /* Create the final map called sigma mapping vid_t to vid_t */
    for (vid_t v = 0; v <= gene_tree.last(); ++v)
        {
            if (!gene_tree.is_leaf(v))
                continue;
            
            string gene_label = gene_tree.label(v);
            string species_label = str_sigma[gene_label];
            if (species_label == "")
                {
                    string message = 
                        "gene label '" + gene_label + "' "
                        "is missing in map file.";
                    throw logic_error(message);
                }
            if (species_tree[species_label] == Binary_tree<T>::NONE)
                {
                    string message = 
                        "species label '" + species_label + "' "
                        "which occurs in map file "
                        "does not exist in species tree.";
                    throw logic_error(message);
                }
            sigma[v] = species_tree[species_label];
        }
}

template<typename T>
void
output_subtree(std::ostream &out, const Binary_tree<T> &tree,
               typename Binary_tree<T>::vid_t subtree)
{
    if (tree.is_leaf(subtree))
        {
            out << subtree << "_" << tree.label(subtree);
            return;
        }
    
    out << "(";
    output_subtree(out, tree, tree.left(subtree));
    out << ",";
    output_subtree(out, tree, tree.right(subtree));
    out << ")" << subtree << "_" << tree.label(subtree);
}

template<typename T>
void
get_postorder_numbering(const Binary_tree<T> &tree,
                        std::vector<typename Binary_tree<T>::vid_t> &numbering)
{
    typedef typename Binary_tree<T>::vid_t vid_t;
    
    numbering.resize(tree.size());

    unsigned i = 0;
    for (vid_t u = tree.postorder_begin();
         u != tree.NONE;
         u = tree.postorder_next(u))
        {
            numbering[u] = i;

            ++i;
        }
}


template<typename T>
std::ostream &
operator<<(std::ostream &out, const Binary_tree<T> &tree)
{
    output_subtree(out, tree, 0);
    return out;
}


template <class RandomNumberGenerator>
void
init_rand(RandomNumberGenerator &engine)
{
    /*
     * Seed the random number generator.
     */
    timeval tv;
    gettimeofday(&tv, 0);
    unsigned seed = tv.tv_sec * 1000000 + tv.tv_usec;
    engine.seed(seed);
}

template <class RandomNumberGenerator>
void
init_rand(RandomNumberGenerator &engine, unsigned seed)
{
    /*
     * Seed the random number generator.
     */
    engine.seed(seed);
}


template<typename T>
void compute_lambda(const Binary_tree<T> &S,
                    const Binary_tree<T> &G,
                    const std::vector<typename Binary_tree<T>::vid_t> &sigma,
                    const boost::dynamic_bitset<> &transfer_edges,
                    std::vector<typename Binary_tree<T>::vid_t> &lambda)
{
    typedef Binary_tree<T> Tree_type;
    typedef typename Tree_type::vid_t vid_t;

    lambda.resize(G.size());
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
            
            if (transfer_edges[v])
                {
                    lambda[u] = lambda[w];
                }
            else if (transfer_edges[w])
                {
                    lambda[u] = lambda[v];
                }
            else
                {
                    lambda[u] = S.lca(lambda[w], lambda[v]);
                }
        }
}

template<typename T>
int count_losses(const Binary_tree<T> &S,
                 const Binary_tree<T> &G,
                 const std::vector<typename Binary_tree<T>::vid_t> &sigma,
                 const boost::dynamic_bitset<> &transfer_edges)
{
    typedef Binary_tree<T> Tree_type;
    typedef typename Binary_tree<T>::vid_t vid_t;

    // Compute lambda
    std::vector<vid_t> lambda;
    compute_lambda(S, G, sigma, transfer_edges, lambda);

    // For each non-transfer edge (u, v) in G, count the number of
    // speciations that we pass from lambda(u) to lambda(v).  A loss
    // is also incurred if u is a duplication and lambda(u) !=
    // lambda(v).
    int losses = 0;
    for (vid_t u = 1; u < G.size(); ++u)
        {
            vid_t p = G.parent(u);
            if (transfer_edges[u] || lambda[p] == lambda[u])
                continue;

            vid_t x = S.parent(lambda[u]);

            vid_t u_sibling = G.left(p) == u ? G.right(p) : G.left(p);
            if (lambda[u_sibling] == lambda[p]) // we know that lampda(u) != lambda(p)!
                {
                    losses += 1;
                }

            while (x != lambda[p])
                {
                    losses += 1;
                    x = S.parent(x);
                }
        }
    
    return losses;
}


#endif /* COMMON_IMPL_HH */
