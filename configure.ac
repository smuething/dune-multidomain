# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dune/pdelab/multidomain/multidomaingridfunctionspace.hh])
AC_CONFIG_HEADERS([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-istl dune-localfunctions dune-pdelab dune-multidomaingrid
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-localfunctions], [dune-pdelab], [dune-multidomaingrid]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, '$(DUNE_LDFLAGS)')
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  test/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/pdelab/Makefile
  dune/pdelab/multidomain/Makefile
  m4/Makefile
  dune-multidomain.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
