/*
 * common.hh
 *
 * Header file containing declarations of functions and classes that
 * are commonly used by programs in the PhylTr package.
 */

#ifndef COMMON_HH
#define COMMON_HH

#include "utils/binary-tree.hh"
extern "C" {
#include <NHtypes.h>
}
#include <map>
#include <string>
#include <ostream>

/*
 * Some_exception
 *
 * An exception class used for easily throwing new exceptions with a
 * meaningful string as the 'what' argument. This is instead of
 * defining a new exception type for each kind of exception we want to
 * throw. This might change in the future when the program is somewhat
 * more mature.
 */
struct Some_exception : public std::exception {
    const char *w_;
    Some_exception(const char *w) : w_(w) {}
    const char *what() const throw() {return w_;}
};


/*
 * NHtree_is_binary()
 *
 * Returns true iff the subtree rooted at 'root' is a full binary
 * tree.
 */
bool NHtree_is_binary(NHNode *root);


/*
 * create_binary_tree()
 *
 * The function will create the NH-subtree rooted at 'root' at the
 * vertex 'ext' in the Binary_tree 'tree'. Note that this means that
 * the 'root' vertex in the NHtree corresponds to the 'ext' vertex in
 * the Binary_tree. This solution is not as elegant a Binary_tree, but
 * we don't really want to do that. Also, the Binary_tree class does
 * not support creation of empty trees. Given these facts, this design
 * was the best we could do.
 */
template<typename T>
void create_binary_tree(NHNode *root, Binary_tree<T> &tree, 
                        typename Binary_tree<T>::vid_t ext);


/*
 * create_gene_species_map()
 *
 * Given the name of the gene-to-species map-file (i.e., sigma), a
 * gene-to-species std::map is returned. It might be bettor to
 * implement this as taking a map and filling it rather than returning
 * the entire map.
 */
template<typename T>
std::vector<typename Binary_tree<T>::vid_t>
create_gene_species_map(const Binary_tree<T> &species_tree, 
                        const Binary_tree<T> &gene_tree, 
                        std::string map_filename);


/*
 * ostream << Binary_tree<T>
 *
 * This is a helper function useful during debugging that prints the
 * Binary_tree 'tree' to the ostream specified in newick format.
 */
template<typename T>
std::ostream &operator<<(std::ostream &out, const Binary_tree<T> &tree);



#include "common-impl.hh"


#endif /* COMMON_HH */
