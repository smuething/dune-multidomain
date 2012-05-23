dunedir = $(includedir)/dune

include $(top_srcdir)/am/global-rules
AM_CPPFLAGS += $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS) $(DUNE_CPPFLAGS) $(DUNEMPIFLAGS)
AM_LDFLAGS += $(ALL_PKG_LDFLAGS) $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS)
LDADD = $(ALL_PKG_LIBS) $(DUNE_LIBS) $(DUNEMPILIBS)

check_PROGRAMS = testmultidomaingridfunctionspace \
	poisson \
	reactiondiffusion \
	explicitlycoupledpoisson \
	instationarypoisson \
	instationaryexplicitlycoupledpoisson \
	testmatrixconstruction \
	testcouplinggfs \
	stokesdarcy \
	stokesdarcy2 \
	navierstokesdarcy-analytic \
	testlocalfunctionspace

testmultidomaingridfunctionspace_SOURCES = testmultidomaingridfunctionspace.cc
#testmultidomaingridfunctionspace_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
#testmultidomaingridfunctionspace_LDFLAGS = $(AM_LDFLAGS) $(ALL_PKG_LDFLAGS)
#testmultidomaingridfunctionspace_LIBS = $(ALL_PKG_LIBS)

poisson_SOURCES = poisson.cc


reactiondiffusion_SOURCES = reactiondiffusion.cc

explicitlycoupledpoisson_SOURCES = explicitlycoupledpoisson.cc

instationarypoisson_SOURCES = instationarypoisson.cc

instationaryexplicitlycoupledpoisson_SOURCES = instationaryexplicitlycoupledpoisson.cc

testmatrixconstruction_SOURCES = testmatrixconstruction.cc

testlocalfunctionspace_SOURCES = testlocalfunctionspace.cc

stokesdarcy_SOURCES = stokesdarcy.cc
stokesdarcy_CPPFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
stokesdarcy_LDFLAGS = $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
stokesdarcy_LDADD = $(LDADD) $(AM_LIBS) $(SUPERLU_LIBS)

stokesdarcy2_SOURCES = stokesdarcy2.cc
stokesdarcy2_CPPFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
stokesdarcy2_LDFLAGS = $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
stokesdarcy2_LDADD = $(LDADD) $(AM_LIBS) $(SUPERLU_LIBS)

navierstokesdarcy_analytic_SOURCES = navierstokesdarcy-analytic.cc
navierstokesdarcy_analytic_CPPFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
navierstokesdarcy_analytic_LDFLAGS = $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
navierstokesdarcy_analytic_LDADD = $(LDADD) $(AM_LIBS) $(SUPERLU_LIBS)

testcouplinggfs_SOURCES = testcouplinggfs.cc

check_PROGRAMS += testvisitor
testvisitor_SOURCES = testvisitor.cc

check_PROGRAMS += testassembler
testassembler_SOURCES = testassembler.cc

check_PROGRAMS += testinstationarypoisson
testinstationarypoisson_SOURCES = testinstationarypoisson.cc

check_PROGRAMS += testoverlappinginstationarypoisson
testoverlappinginstationarypoisson_SOURCES = testoverlappinginstationarypoisson.cc

check_PROGRAMS += testxfempattern
testxfempattern_SOURCES = testxfempattern.cc

check_PROGRAMS += testmortarpoisson
testmortarpoisson_SOURCES = testmortarpoisson.cc

check_PROGRAMS += testpowermortarpoisson
testpowermortarpoisson_SOURCES = testpowermortarpoisson.cc