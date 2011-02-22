/*
 * Copyright (C) 2010, 2011 Ali Tofigh
 *
 * This file is part of PhylTr, a package for phylogenetic analysis
 * using duplications and transfers.
 *
 * PhylTr is released under the terms of the license contained in the
 * file LICENSE.
 */

#include "assert.hh"
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <numeric>
#include <functional>
#include <limits>

/*
 * Static variable definitions.
 */

template<typename R>
const typename Binary_tree<R>::vid_t Binary_tree<R>::NONE;

/*
 * Definitions of exception classes used in Binary_tree
 */

template<typename R>
const char *
Binary_tree<R>::Invalid_length::what() const throw()
{
    return "Invalid branch length passed to Binary_tree class.";
}

template<typename R>
const char *
Binary_tree<R>::Invalid_id::what() const throw()
{
    return "Invalid vertex id passed to Binary_tree class.";
}

template<typename R>
const char *
Binary_tree<R>::Invalid_label::what() const throw()
{
    return "Invalid vertex label passed to Binary_tree class. "
        "All labels in a tree must be unique.";
}

template<typename R>
const char *
Binary_tree<R>::Invalid_parent::what() const throw()
{
    return "Invalid parent passed to Binary_tree::add_children. "
        "It is not possible to add children to vertices who "
        "already have children.";
}

/*
 * Constructor of internal vertex_t class.
 */

template<typename R>
Binary_tree<R>::vertex_t::vertex_t(vid_t l, vid_t r, vid_t p,
                                real_type len, std::string label)
    : length(len), left(l), right(r), parent(p), label(label)
{}


/*
 * Constructor and functions definitions.
 */

template<typename R>
Binary_tree<R>::Binary_tree(real_type length, std::string label)
    : lca_is_valid(false)
{
    Assert<Invalid_length>(length >= 0);
    
    vertex_t v(NONE, NONE, NONE, length, label);
    vertices.push_back(v);
    if (label != "") label_map[label] = 0;
}

template<typename R>
void
Binary_tree<R>::add_children(vid_t parent,
                          real_type l_len, real_type r_len,
                          std::string l_label, std::string r_label)
{
    Assert<Invalid_id>(parent >= 0 && parent < size());
    Assert<Invalid_length>(l_len >= 0 && r_len >= 0);
    Assert<Invalid_label>(l_label == ""
                          || label_map.find(l_label) == label_map.end());
    Assert<Invalid_label>(r_label == "" ||
                          label_map.find(r_label) == label_map.end());
    Assert<Invalid_label>(l_label == "" || l_label != r_label);
    Assert<Invalid_parent>(is_leaf(parent));
    
    vertex_t v1(NONE, NONE, parent, l_len, l_label);
    vertex_t v2(NONE, NONE, parent, r_len, r_label);
    vertices[parent].left = vertices.size();
    vertices.push_back(v1);
    vertices[parent].right = vertices.size();
    vertices.push_back(v2);
    label_map[l_label] = left(parent);
    label_map[r_label] = right(parent);

    lca_is_valid = false;
}
    
template<typename R>
typename Binary_tree<R>::size_type
Binary_tree<R>::size() const
{
    return vertices.size();
}

template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::last() const
{
    return static_cast<vid_t>(size() - 1);
}

template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::operator[](std::string label) const
{
    std::map<std::string, vid_t>::const_iterator iter = label_map.find(label);
    if (iter == label_map.end()) return NONE;

    return iter->second;
}

template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::left(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return vertices[v].left;
}

template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::right(vid_t v) const 
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return vertices[v].right;
}

template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::parent(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return vertices[v].parent;
}

template<typename R>
std::string
Binary_tree<R>::label(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());
    
    return vertices[v].label;
}

template<typename R>
typename Binary_tree<R>::real_type
Binary_tree<R>::length(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return vertices[v].length;
}

template<typename R>
typename Binary_tree<R>::real_type
Binary_tree<R>::time(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    std::vector<real_type> lengths;
    while (v != NONE) {
        lengths.push_back(vertices[v].length);
        v = vertices[v].parent;
    }

    std::sort(lengths.begin(), lengths.end());
    return std::accumulate(lengths.begin(), lengths.end(), 0.0);
}


template<typename R>
bool
Binary_tree<R>::is_leaf(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    /* We have a full binary tree, so we only need to check one child. */
    return vertices[v].left == NONE;
}

template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::preorder_begin() const
{
    return 0;
}

template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::preorder_next(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    if (!is_leaf(v))
        return left(v);

    /* special case when the tree consists of one vertex only. */
    if (v == 0)
        return NONE;
    
    while (v != 0 && right(parent(v)) == v)
        v = parent(v);
    
    if (v == 0)
        return NONE;
    
    return right(parent(v));
}
 
template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::postorder_begin() const
{
    vid_t v = 0;
    while (!is_leaf(v))
        v = left(v);

    return v;
}

template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::postorder_next(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());
    
    if (v == 0)
        return NONE;

    vid_t p = parent(v);
    vid_t r = right(p);

    if (v == r)
        return p;

    while (!is_leaf(r))
        r = left(r);
    return r;
}

/*
 * The algorithm is from the article 'Lowest common ancestors in trees
 * and directed acyclic graphs' by M Bender et al, 2005. The tables
 * and arrays have the same names as in the article, except for the
 * R-array which here is renamed to Ref.
 */
template<typename R>
typename Binary_tree<R>::vid_t
Binary_tree<R>::lca(vid_t v1, vid_t v2) const
{
    Assert<Invalid_id>(v1 >= 0 && v1 < size());
    Assert<Invalid_id>(v2 >= 0 && v2 < size());

    if (!lca_is_valid)
        build_lca();
    
    /* Get the representatives (i.e, indexes into L) of v1 and v2. */
    unsigned r1 = Ref[v1];
    unsigned r2 = Ref[v2];

    /* Make sure that r2 is the bigger one, so that the range
       of indices is [r1, r2]. */
    if (r1 > r2)
        std::swap(r1, r2);

    unsigned k = most_significant_bit(r2 - r1 + 1);
    unsigned idx1 = M[r1][k];
    unsigned idx2 = M[r2-(1u<<k)+1][k]; /* 1u<<k == 2^k */
    if (L[idx2] < L[idx1])
        idx1 = idx2;

    return E[idx1];
}
    
template<typename R>
void
Binary_tree<R>::build_lca() const
{
    /* First, we build the E and L arrays, i.e., the Euler tour and
       the level array. */
    
    /* We define a recursive function to compute E and L. */
    struct Create_EL {
        void operator()(const Binary_tree<R> &tree,
                        vid_t cur_vertex, unsigned cur_level,
                        std::vector<vid_t> &E, std::vector<vid_t> &L)
        {
            E.push_back(cur_vertex);
            L.push_back(cur_level);
            
            if (tree.is_leaf(cur_vertex))
                return;
            
            this->operator()(tree, tree.left(cur_vertex), cur_level + 1, E, L);
            E.push_back(cur_vertex);
            L.push_back(cur_level);
            this->operator()(tree,tree.right(cur_vertex), cur_level + 1, E, L);
            E.push_back(cur_vertex);
            L.push_back(cur_level);
        }
    };

    E.clear(); E.reserve(this->size()); 
    L.clear(); L.reserve(this->size());
    
    Create_EL create_EL;
    create_EL(*this, 0, 0, E, L);
    
    /* Create the R-vector. */
    Ref.clear(); Ref.resize(this->size(), NONE);
    for (unsigned i = 0; i < E.size(); ++i)
        {
            if (Ref[E[i]] == NONE)
                Ref[E[i]] = i;
        }


    /* Initialize and build the M-matrix. */
    unsigned rows = L.size();
    unsigned cols = most_significant_bit(L.size()) + 1;
    M.resize(boost::extents[rows][cols]);

    for (unsigned i = 0; i < rows; ++i)
        M[i][0] = i;
    
    for (unsigned j = 1; j < cols; ++j)
        {
            unsigned stride = 1u << j;
            for (unsigned i = 0; i < rows; ++i)
                {
                    if (i + stride/2 >= rows)
                        {
                            M[i][j] = M[i][j-1];
                            continue;
                        }
                    unsigned idx1 = M[i][j-1];
                    unsigned idx2 = M[i + stride/2][j-1];
                    
                    M[i][j] = L[idx2] < L[idx1] ? idx2 : idx1;
                }
        }

    lca_is_valid = true;
}


template<typename R>
template<typename T>
unsigned
Binary_tree<R>::most_significant_bit(T v)
{
    static const char LogTable256[] = 
        {
            0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
        };
    
    unsigned bits = std::numeric_limits<T>::digits;
    unsigned c = 0; // c will be lg(v) truncated
    while (bits > 8)
        {
            bits /= 2;
            if (v >> bits)
                {
                    c += bits;
                    v >>= bits;
                }
        }
    c += LogTable256[v];
    return c;
}


template<typename R>
bool
Binary_tree<R>::descendant(vid_t v1, vid_t v2) const
{
    /* Is v1 descendant of v2? */

    return lca(v1, v2) == v2;
}
