#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([phyltr], [0.9], [alix.tofigh@gmail.com])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config/m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([src/phyltr-dp.cc])
AC_CONFIG_HEADERS([src/config.h])
AC_LANG([C++])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.

## we test the result of checking of boost outside of AX_BOOST_BASE
## becuase if the user sets --with-boost=no, then niether the
## ACTION-IF-FOUND nor ACTION-IF-NOT-FOUND is executed.
AX_BOOST_BASE([1.34], [ax_have_boost=yes], [])
if test "x$ax_have_boost" != xyes; then
    AC_MSG_FAILURE([cannot build without boost >= 1.34])
fi
AX_BOOST_PROGRAM_OPTIONS
CPPFLAGS_SAVED="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
AC_CHECK_HEADERS([boost/dynamic_bitset.hpp boost/foreach.hpp \
                  boost/lambda/lambda.hpp boost/multi_array.hpp \
                  boost/program_options.hpp boost/random.hpp \
                  boost/rational.hpp boost/shared_ptr.hpp boost/smart_ptr.hpp \
                  boost/static_assert.hpp boost/weak_ptr.hpp],
                 [], [AC_MSG_FAILURE([Required header not found])])
CPPFLAGS="$CPPFLAGS_SAVED" 

# Check for NHparser library
NHPARSER_CPPFLAGS=
NHPARSER_LDFLAGS=
CPPFLAGS_SAVED="$CPPFLAGS"
LDFLAGS_SAVED="$LDFLAGS"
AC_ARG_WITH([NHparser], 
            [AS_HELP_STRING([--with-NHparser[=PATH]],
                            [NHparser installation directory])],
            [AS_IF([test "x$with_NHparser" = xno],
                   [AC_MSG_ERROR([NHparser library is required])])],
            [])
AS_IF([test "x$with_NHparser" = xno],
      [
        AC_MSG_ERROR([NHparser library is required])
      ],
      [test "x$with_NHparser" != yes],
      [
        with_NHparser=`echo "$with_NHparser" | sed 's:/$::'`
        NHPARSER_CPPFLAGS="-I$with_NHparser/include"
        NHPARSER_LDFLAGS="-L$with_NHparser/lib"
      ])
CPPFLAGS="$CPPFLAGS_SAVED $NHPARSER_CPPFLAGS"
AC_CHECK_HEADERS([NHparser/NHparser.h NHparser/NHtypes.h], [],
                 [AC_MSG_FAILURE([NHparser headers not found])])
LDFLAGS="$LDFLAGS_SAVED $NHPARSER_LDFLAGS"
AC_CHECK_LIB([NHparser], [NH_read_tree], [:],
             [
               NHPARSER_LDFLAGS="-L$with_NHparser/lib64"
               LDFLAGS="$LDFLAGS_SAVED $NHPARSER_LDFLAGS"
               AC_CHECK_LIB([NHparser], [NH_read_tree], [:],
                            [AC_MSG_FAILURE([NHparser library not found])])
             ])
NHPARSER_LIB=-lNHparser
AC_DEFINE([HAVE_LIBNHPARSER],,[Define if the NHparser library is available])
AC_SUBST(NHPARSER_CPPFLAGS)
AC_SUBST(NHPARSER_LDFLAGS)
AC_SUBST(NHPARSER_LIB)
CPPFLAGS="$CPPFLAGS_SAVED"
LDFLAGS="$LDFLAGS_SAVED"

# Checks for misc header files.
AC_CHECK_HEADERS_ONCE([algorithm bitset cerrno cstdio cstdlib errno.h \
                       exception fstream functional iomanip iosfwd iostream \
                       iterator limits map numeric ostream set sstream stack \
                       stdexcept string sys/time.h utility vector], [],
                      [AC_MSG_FAILURE([Required header not found])])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday])


AC_CONFIG_FILES([
    Makefile
    src/Makefile
])
AC_OUTPUT
