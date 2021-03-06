/*
 * Copyright (C) 2010, 2011 Ali Tofigh
 *
 * This file is part of PhylTr, a package for phylogenetic analysis
 * using duplications and transfers.
 *
 * PhylTr is released under the terms of the license contained in the
 * file LICENSE.
 */

/*
 * common.hh
 *
 * Header file containing declarations of functions and classes that
 * are commonly used by programs in the PhylTr package.
 */

#ifndef COMMON_HH
#define COMMON_HH

#include "utils/binary-tree.hh"
#include <NHparser/NHtypes.h>
#include <string>
#include <iosfwd>
#include <vector>
#include <boost/random.hpp>
#include <boost/dynamic_bitset.hpp>


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
 * This function computes the mapping (sigma) between gene tree leaves
 * and species tree leaves, and stores the result in 'sigma'. The size
 * of 'sigma' must be at least gene_tree.size(). 
 *
 * May throw std::logic_error in case problems occur while reading the
 * map file. The what() function will then return a meaningful
 * description that can be output by the program to the user. This is
 * more convenient at the moment than providing different error
 * classes for each possible error.
 */
template<typename T>
void create_gene_species_map(const Binary_tree<T> &species_tree, 
                             const Binary_tree<T> &gene_tree, 
                             std::string map_filename,
                             std::vector<typename Binary_tree<T>::vid_t> &sigma);

/*
 * compute_lambda()
 *
 * Computes the least common ancestor mapping of the gene tree into
 * the species tree given a transfer edges, and stores the result in
 * the vector 'lambda'. The size of 'lambda' is adjusted by the
 * function to the size of the gene tree.
 */
template<typename T>
void compute_lambda(const Binary_tree<T> &species_tree,
                    const Binary_tree<T> &gene_tree,
                    const std::vector<typename Binary_tree<T>::vid_t> &sigma,
                    const boost::dynamic_bitset<> &transfer_edges,
                    std::vector<typename Binary_tree<T>::vid_t> &lambda,
                    bool dated_stree = false);


/*
 * count_losses()
 *
 * computes the number of losses given a species tree, gene tree, and
 * a set of transfer edges.
 */
template<typename T>
int count_losses(const Binary_tree<T> &species_tree,
                 const Binary_tree<T> &gene_tree,
                 const std::vector<typename Binary_tree<T>::vid_t> &sigma,
                 const boost::dynamic_bitset<> &transfer_edges,
                 bool dated_stree = false);

/*
 * get_postorder_numbering()
 *
 * Provides a postorder numbering of the vertices in 'tree'. The size
 * of the vector 'numbering' will be 'tree.size()' after the
 * call. After a call to this function, numbering[v] is the postorder
 * number of vertex v.
 */
template<typename T>
void get_postorder_numbering(const Binary_tree<T> &tree,
                             std::vector<typename Binary_tree<T>::vid_t> &numbering);

/*
 * ostream << Binary_tree<T>
 *
 * This is a helper function useful during debugging that prints the
 * Binary_tree 'tree' to the ostream specified in newick format.
 */
template<typename T>
std::ostream &operator<<(std::ostream &out, const Binary_tree<T> &tree);


/*
 * Global random number generators:
 *
 * 'g_rng_ui': used to get a uniformly distributed random
 * unsigned. g_rng_ui(n) gives a random number in the range [0,n).
 *
 * 'g_rng_d': used to get a uniformly distributed random double in the
 * range [0,1).
 *
 * The random number engine 'g_generator' should be seeded in main.
 */
extern boost::mt19937    g_generator;
extern boost::random_number_generator<boost::mt19937, unsigned> g_rng_ui;
extern boost::variate_generator<boost::mt19937&, boost::uniform_real<> > g_rng_d;


/*
 * init_rand()
 *
 * initializes the random number generator using gettimeofday() or a
 * supplied seed. 'engine' must have a seed-method taking an unsigned
 * int as argument. The seed used is returned.
 */
template <class RandomNumberGenerator>
unsigned init_rand(RandomNumberGenerator &engine);

template <class RandomNumberGenerator>
unsigned init_rand(RandomNumberGenerator &engine, unsigned seed);





#include "common-impl.hh"


#endif /* COMMON_HH */
