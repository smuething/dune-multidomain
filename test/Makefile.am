dunedir = $(includedir)/dune

include $(top_srcdir)/am/global-rules
AM_CPPFLAGS += $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS += $(ALL_PKG_LDFLAGS)
AM_LIBS = $(ALL_PKG_LIBS)
CXXFLAGS += -g3

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
stokesdarcy_LDADD = $(AM_LIBS) $(SUPERLU_LIBS)

stokesdarcy2_SOURCES = stokesdarcy2.cc
stokesdarcy2_CPPFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
stokesdarcy2_LDFLAGS = $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
stokesdarcy2_LDADD = $(AM_LIBS) $(SUPERLU_LIBS)

testcouplinggfs_SOURCES = testcouplinggfs.cc

check_PROGRAMS += testvisitor
testvisitor_SOURCES = testvisitor.cc

check_PROGRAMS += testassembler
testassembler_SOURCES = testassembler.cc

check_PROGRAMS += testinstationarypoisson
testinstationarypoisson_SOURCES = testinstationarypoisson.cc
