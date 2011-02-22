/*
 * Copyright (C) 2010, 2011 Ali Tofigh
 *
 * This file is part of PhylTr, a package for phylogenetic analysis
 * using duplications and transfers.
 *
 * PhylTr is released under the terms of the license contained in the
 * file LICENSE.
 */

#ifndef ASSERT_HH
#define ASSERT_HH


template<class E, class A>
inline void Assert(A assertion)
{
    if (!assertion) throw E();
}

#endif /* not ASSERT_HH. */
