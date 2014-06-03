set(TESTS
  reactiondiffusion.cc
  stokesdarcy.cc
  stokesdarcy2.cc
  testcouplinggfs.cc
  testinstationarypoisson.cc
  testlocalfunctionspace.cc
  testmatrixconstruction.cc
  testmortarpoisson.cc
  testmultidomaingridfunctionspace.cc
  testoverlappinginstationarypoisson.cc
  testpoisson.cc
  testpowermortarpoisson.cc
  testvisitor.cc
  testxfempattern.cc
)

add_executable(reactiondiffusion reactiondiffusion.cc)
add_executable(stokesdarcy stokesdarcy.cc)
add_dune_superlu_flags(stokesdarcy)
add_executable(stokesdarcy2 stokesdarcy2.cc)
add_dune_superlu_flags(stokesdarcy2)
add_executable(testcouplinggfs testcouplinggfs.cc)
add_executable(testinstationarypoisson testinstationarypoisson.cc)
add_executable(testlocalfunctionspace testlocalfunctionspace.cc)
add_executable(testmatrixconstruction testmatrixconstruction.cc)
add_executable(testmortarpoisson testmortarpoisson.cc)
add_executable(testmultidomaingridfunctionspace testmultidomaingridfunctionspace.cc)
add_executable(testoverlappinginstationarypoisson testoverlappinginstationarypoisson.cc)
add_executable(testpoisson testpoisson.cc)
add_executable(testpowermortarpoisson testpowermortarpoisson.cc)
add_executable(testvisitor testvisitor.cc)
add_executable(testxfempattern testxfempattern.cc)

# install th header as done for the auto-tools
install(FILES checkgeometry.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/test)

foreach(_test ${TESTS})
  add_test(${_test} ${_test})
#   target_link_libraries(${_test} "dunecommon")
endforeach(_test ${TESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})
