dunedir = $(includedir)/dune

include $(top_srcdir)/am/global-rules
AM_CPPFLAGS += $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS += $(ALL_PKG_LDFLAGS) $(DUNEMPILDFLAGS)
LDADD = $(ALL_PKG_LIBS) $(DUNEMPILIBS)
AM_LIBS = $(ALL_PKG_LIBS) $(DUNEMPILIBS)

check_PROGRAMS = testmultidomaingridfunctionspace \
	reactiondiffusion \
	explicitlycoupledpoisson \
	testmatrixconstruction \
	testcouplinggfs \
	stokesdarcy \
	stokesdarcy2 \
	testlocalfunctionspace

testmultidomaingridfunctionspace_SOURCES = testmultidomaingridfunctionspace.cc
#testmultidomaingridfunctionspace_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
#testmultidomaingridfunctionspace_LDFLAGS = $(AM_LDFLAGS) $(ALL_PKG_LDFLAGS)
#testmultidomaingridfunctionspace_LIBS = $(ALL_PKG_LIBS)

reactiondiffusion_SOURCES = reactiondiffusion.cc

explicitlycoupledpoisson_SOURCES = explicitlycoupledpoisson.cc

testmatrixconstruction_SOURCES = testmatrixconstruction.cc

testlocalfunctionspace_SOURCES = testlocalfunctionspace.cc

stokesdarcy_SOURCES = stokesdarcy.cc
stokesdarcy_CPPFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
stokesdarcy_LDFLAGS = $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
stokesdarcy_LDADD = $(AM_LIBS) $(SUPERLU_LIBS)

stokesdarcy2_SOURCES = stokesdarcy2.cc
stokesdarcy2_CPPFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
stokesdarcy2_LDFLAGS = $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
stokesdarcy2_LDADD = $(AM_LIBS) $(SUPERLU_LIBS)

testcouplinggfs_SOURCES = testcouplinggfs.cc

check_PROGRAMS += testinstationarypoisson
testinstationarypoisson_SOURCES = testinstationarypoisson.cc

check_PROGRAMS += testoverlappinginstationarypoisson
testoverlappinginstationarypoisson_SOURCES = testoverlappinginstationarypoisson.cc

check_PROGRAMS += testmortarpoisson
testmortarpoisson_SOURCES = testmortarpoisson.cc

check_PROGRAMS += testpowermortarpoisson
testpowermortarpoisson_SOURCES = testpowermortarpoisson.cc

check_PROGRAMS += testpoisson
testpoisson_SOURCES = testpoisson.cc
testpoisson_LDADD = $(AM_LIBS) $(DUNEMPILIBS)

check_PROGRAMS += testpoisson-dirichlet-neumann
testpoisson_dirichlet_neumann_SOURCES = testpoisson-dirichlet-neumann.cc
testpoisson_dirichlet_neumann_LDADD = $(AM_LIBS) $(DUNEMPILIBS)

check_PROGRAMS += testpoisson-assembly
testpoisson_assembly_SOURCES = testpoisson-assembly.cc
testpoisson_assembly_LDADD = $(AM_LIBS) $(DUNEMPILIBS)

MOSTLYCLEANFILES = *.vtu *.pvtu *.pvd
