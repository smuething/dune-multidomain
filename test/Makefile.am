dunedir = $(includedir)/dune

include $(top_srcdir)/am/global-rules
AM_CPPFLAGS += $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS += $(ALL_PKG_LDFLAGS)
AM_LIBS = $(ALL_PKG_LIBS)

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

MOSTLYCLEANFILES = *.vtu *.pvtu *.pvd
