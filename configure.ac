#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(BTL_BLOOMFILTER, 1.0.1, cjustin@bcgsc.ca)
AC_CANONICAL_SYSTEM

AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([BloomFilter.hpp])
AC_CONFIG_HEADER([config.h])
AC_PROG_RANLIB

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
AC_CHECK_LIB([dl], [dlopen])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Set compiler flags.
# make sure -O3 by default for speed
AC_SUBST(CPPFLAGS,"-O3 -std=c++11 $CPPFLAGS")
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Werror $AM_CXXFLAGS")

# Options to configure.
AC_ARG_WITH(sdsl, AS_HELP_STRING([--with-sdsl=PATH],
			[specify directory for the sdsl header files]))
if test "$with_sdsl" -a -d "$with_sdsl"; then
	sdsl_cppflags="-isystem$with_sdsl"
fi

# Find the absolute paths.
srcdir=$(cd $srcdir; pwd)

# Set compiler flags.
AC_SUBST(CPPFLAGS,"-I$srcdir $sdsl_cppflags $CPPFLAGS")

# Checks for header files.
AC_LANG([C++])	#sets the auto conf to look for c++ headers


AC_CHECK_HEADERS([sdsl/bit_vector_il.hpp])
if test $ac_cv_header_sdsl_bit_vector_il_hpp != yes; then
	AC_MSG_ERROR([Requires sdsl])
fi

AC_CONFIG_FILES([
		Makefile
		Tests/Unit/Makefile
		Tests/AdHoc/Makefile
])
AC_OUTPUT

