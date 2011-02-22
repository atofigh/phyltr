/*
 * Copyright (C) 2010, 2011 Ali Tofigh
 *
 * This file is part of PhylTr, a package for phylogenetic analysis
 * using duplications and transfers.
 *
 * PhylTr is released under the terms of the license contained in the
 * file LICENSE.
 */

#ifndef BINARY_TREE_HH
#define BINARY_TREE_HH

#include "assert.hh"

#include <exception>
#include <vector>
#include <string>
#include <limits>
#include <map>
#include <boost/multi_array.hpp>
#include <boost/static_assert.hpp>

/*
 * [ type Binary_tree ]: The class Binary_tree implements the full
 * rooted binary tree concept. Furthermore, a unique string label may
 * be attached to any vertex. An empty string implies no label. Every
 * edge has a non-negative length. In addition, the root vertex is
 * considered as having an incoming edge.
 *
 * The template parameter must be a type behaving like a
 * floating-point type. For now, the code checks that
 *
 * std::numeric_limits<Real>::is_specialized == true, and
 * std::numeric_limits<Real>::is_integer == false.
 *
 * This might change in the future when TR1 becomes part of the
 * standard, in which case the check would be:
 *
 * std::is_floating_point<Real>::value == true.
 *
 * Invariants:
 * 1. size() >= 1 (There is always a root.)
 * 2. All vertices have unique non-negative ids.
 * 3. The root's id is always zero.
 * 4. Vertex labels are unique (except for the empty string).
 * 5. if v is a valid vertex, either left(v)==right(v)==NONE, or
 *    left(v) and right(v) are both valid vertices.
 * 6. label_map[str] == v  <==>  label(v) == str.
 * 7. All edge lengths are valid non-negative real numbers (except the
 *    root edge which may be zero). Valid means not +-inf or Nan.
 */

template <typename Real>
class Binary_tree {
public:
    /*
     * [ type vid_t ]: An unsigned integral type used for vertex
     * ids. A constant value 'NONE' is reserved to signal a
     * non-existant vertex id and is guaranteed not to be allocated
     * for a vertex.
     */
    typedef unsigned                            vid_t;

    typedef std::vector<vid_t>::size_type       size_type;
    typedef Real                                real_type;

    /*
     * We cannot use numeric_limits<vid_t>::max() to initialize NONE.
     * It is illegal to call a function here, and when initializing it
     * outside, we run into a bug in gcc.
     */
    static const vid_t                          NONE = -1;

    /*
     * [ Exception classes ]: Exceptions that may be thrown by the
     * functions in this class. All derived from std::exception.
     */

    struct Invalid_length 
        : public std::exception {const char *what() const throw();};
    struct Invalid_id 
        : public std::exception {const char *what() const throw();};
    struct Invalid_label 
        : public std::exception {const char *what() const throw();};
    struct Invalid_parent : 
        public std::exception {const char *what() const throw();};
    
    /*
     * [ Constructors ]: The optional paramters set the the length of
     * the incoming edge of the root, and the root's label.
     *
     * Preconditions:
     *   1. 'length' is non-negative and valid.
     *
     * Default Copying and assignment are used.
     */
    explicit Binary_tree(real_type length = 0, std::string label = "");

    /*
     * [ add_children ]: Creates children for existing vertex.
     *
     * Preconditions:
     * 1. 'parent' exists.
     * 2. 'left_len' and 'right_len' are valid non-negative reals.
     * 3. Labels assigned to vertices must be unique (except for "").
     * 4. 'parent' must be a leaf.
     */
    void add_children(vid_t parent,
                      real_type left_len,
                      real_type right_len,
                      std::string left_label = "",
                      std::string right_label = "");

    /*
     * [ size ]: Returns the number of vertices in the tree.
     */
    size_type size() const;

    /*
     * [ last ]: Returns the largest valid vertex id in the tree.
     */
    vid_t last() const;

    /*
     * [ operator[] ]: Returns the id for the vertex with label
     * 'label'. If no such vertex exists, it will return 'NONE'
     */
    vid_t operator[](std::string label) const;

    /*
     * [ left/right/parent ]: These functions return the id of the left
     * child, right child, and parent of 'vertex', respectively. If such
     * a vertex does not exist, NONE is returned.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    vid_t left(vid_t vertex) const;
    vid_t right(vid_t vertex) const;
    vid_t parent(vid_t vertex) const;

    /*
     * [ label ]: Returns the unique label of 'vertex'. If no label
     * was set, the empty string is returned.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    std::string label(vid_t vertex) const;

    /*
     * [ length ]: Returns the length of the incoming edge of
     * 'vertex'.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    real_type length(vid_t vertex) const;

    /*
     * [ time ]: Returns the sum of edge lengths from the root to
     * 'vertex'. NB!!! This does not include the length of the edge
     * coming in to the root.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    real_type time(vid_t vertex) const;

    /*
     * [ is_leaf ]: Returns true if vertex has no children.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    bool is_leaf(vid_t vertex) const;

    /*
     * [ preorder_begin/preorder_next]: preorder_begin returns the
     * first vertex of a preorder traversal of the tree (this is
     * always the root, but was included for symmetri with the
     * postorder functions). preorder_next returns the next vertex in
     * a preorder traversal or NONE if there are no more vertices to visit.
     *
     * A preorder traversal of the tree
     *        ---4
     *       |
     *   ----2
     *   |   |
     *---0    ---3
     *   |
     *   --------1
     *
     * would be 0, 1, 2, 3, 4.
     *
     * The time it takes for these functions to return is not
     * constant, but the total contribution from these two functions
     * in a loop that visits all vertices is O(n), where n is the
     * number of vertices in the tree.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    vid_t preorder_begin() const;
    vid_t preorder_next(vid_t vertex) const;

    /*
     * [ postorder_begin/postorder_next]: postorder_begin returns the
     * first vertex of a postorder traversal of the
     * tree. postorder_next returns the next vertex in a postorder
     * traversal or NONE if there are no more vertices to visit.
     *
     * A postorder traversal of the tree
     *        ---4
     *       |
     *   ----2
     *   |   |
     *---0    ---3
     *   |
     *   --------1
     *
     * would be 1, 3, 4, 2, 0.
     *
     * The comments on time complexity for the preorder-functions
     * above apply for the postorder-functions as well.
     * 
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    vid_t postorder_begin() const;
    vid_t postorder_next(vid_t vertex) const;

    /*
     * [ lca ]: Returns the the least common ancestor of vertices v1
     * and v2. The first time this function is called, some
     * precomputation is done that takes time O(n log n). Thereafter
     * the lca is computed in constant time. Note, however, that after
     * each change to the tree (e.g. after a call to add_children),
     * the precomputation has to be redone!
     *
     * The algorithm is from the article 'Lowest common ancestors in
     * trees and directed acyclic graphs' by M Bender et al., 2005.
     *
     * Precondition:
     * 1. 'v1' and 'v2' are valid vertex ids.
     */
    vid_t lca(vid_t v1, vid_t v2) const;

    /*
     * [ descendant ]: Returns true if 'v1' is a descendant of
     * 'v2'. This function is simply a wrapper around lca, since 'v1'
     * is a descendant of 'v2' if and only if lca(v1,v2) = v2.
     */
    bool descendant(vid_t v1, vid_t v2) const;


private:
    BOOST_STATIC_ASSERT(std::numeric_limits<real_type>::is_specialized == true &&
                        std::numeric_limits<real_type>::is_integer == false);
    struct vertex_t {
        real_type length;
        vid_t left, right, parent;
        std::string label;
        vertex_t(vid_t l, vid_t r, vid_t p, real_type len, std::string label);
    };

    /* Helper function needed for computing lca. */
    template<typename T> static unsigned most_significant_bit(T v);
    /* Function for building lca. */
    void build_lca() const;
    
    std::vector<vertex_t> vertices;
    std::map<std::string, vid_t> label_map;

    mutable std::vector<vid_t> E; /* Euler-path for lca-compuatation. */
    mutable std::vector<unsigned> L; /* Level array corresponding to E. */
    mutable std::vector<std::vector<vid_t>::size_type> Ref;
                             /* Representative array for lca-computation */
    mutable boost::multi_array<std::vector<vid_t>::size_type, 2> M; 
                                     /* M-matrix for the RMQ-algorithm. */
    mutable bool lca_is_valid;
};


#include "binary-tree-impl.hh"



#endif /* not BINARY_TREE_HH. */
