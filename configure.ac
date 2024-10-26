#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(EnumerationIndependentSet, 1.0, mdsikir@irb.hr)
AM_INIT_AUTOMAKE(EnumerationIndependentSet, 1.0)
AC_OUTPUT(Makefile src/Makefile)
AC_CONFIG_SRCDIR([src/CanonicIndependFamily.c])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgmp':
AC_CHECK_LIB([gmp], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_OUTPUT
