#ifndef ASSERT_HH
#define ASSERT_HH


template<class E, class A>
inline void Assert(A assertion)
{
    if (!assertion) throw E();
}

#endif /* not ASSERT_HH. */
