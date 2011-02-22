/*
 * Copyright (C) 2010, 2011 Ali Tofigh
 *
 * This file is part of PhylTr, a package for phylogenetic analysis
 * using duplications and transfers.
 *
 * PhylTr is released under the terms of the license contained in the
 * file LICENSE.
 */

#include "common.hh"
#include <boost/random.hpp>

bool
NHtree_is_binary(NHNode *root)
{
    /* true if root is a leaf. */
    if (root->children == 0)
        return true;

    /* false if root does not have exactly two children. */
    if (root->children->next == 0
        || root->children->next->next != 0)
        return false;

    /* otherwise it depends on the children of root. */
    return NHtree_is_binary(root->children->node) 
        && NHtree_is_binary(root->children->next->node);
}


boost::mt19937           g_generator;
boost::random_number_generator<boost::mt19937, unsigned>
                         g_rng_ui(g_generator);
boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
                         g_rng_d(g_generator, boost::uniform_real<>());
