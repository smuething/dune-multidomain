# we need the module file to be able to build via dunecontrol
EXTRA_DIST =					\
	CMakeLists.txt				\
	COPYING.md				\
	dune.module				\
	GPL-2					\
	GPL-3					\
	LGPL-3					\
	README					\
	README.md				\
	README.GIT.md				\
	RELEASE_NOTES.md

SUBDIRS = m4 dune doc test

# don't follow the full GNU-standard
# we need automake 1.9
AUTOMAKE_OPTIONS = foreign 1.9

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT) --with-dune-multidomaingrid=$(DUNE_MULTIDOMAINGRID_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules
