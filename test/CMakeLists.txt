set(TESTS
  reactiondiffusion
  stokesdarcy
  stokesdarcy2
  testcouplinggfs
  testinstationarypoisson
  testlocalfunctionspace
  testmatrixconstruction
  testmortarpoisson
  testmultidomaingridfunctionspace
  testoverlappinginstationarypoisson
  testpoisson
  testpowermortarpoisson
)

add_executable(reactiondiffusion EXCLUDE_FROM_ALL reactiondiffusion.cc)
add_executable(stokesdarcy EXCLUDE_FROM_ALL stokesdarcy.cc)
add_dune_superlu_flags(stokesdarcy)
add_executable(stokesdarcy2 EXCLUDE_FROM_ALL stokesdarcy2.cc)
add_dune_superlu_flags(stokesdarcy2)
add_executable(testcouplinggfs EXCLUDE_FROM_ALL testcouplinggfs.cc)
add_executable(testinstationarypoisson EXCLUDE_FROM_ALL testinstationarypoisson.cc)
add_executable(testlocalfunctionspace EXCLUDE_FROM_ALL testlocalfunctionspace.cc)
add_executable(testmatrixconstruction EXCLUDE_FROM_ALL testmatrixconstruction.cc)
add_executable(testmortarpoisson EXCLUDE_FROM_ALL testmortarpoisson.cc)
add_executable(testmultidomaingridfunctionspace EXCLUDE_FROM_ALL testmultidomaingridfunctionspace.cc)
add_executable(testoverlappinginstationarypoisson EXCLUDE_FROM_ALL testoverlappinginstationarypoisson.cc)
add_executable(testpoisson EXCLUDE_FROM_ALL testpoisson.cc)
add_executable(testpowermortarpoisson EXCLUDE_FROM_ALL testpowermortarpoisson.cc)

foreach(_test ${TESTS})
  add_test(${_test} ${_test})
  target_link_libraries(${_test} ${DUNE_LIBS})
endforeach(_test ${TESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})
