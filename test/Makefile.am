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

testcouplinggfs_SOURCES = testcouplinggfs.cc
