set(tests
  stokesdarcy2
  testcouplinggfs
  testlocalfunctionspace
  testmatrixconstruction
  testmortarpoisson
  testmultidomaingridfunctionspace
  testoverlappinginstationarypoisson
  testpowermortarpoisson
)

set(tests_with_arguments
  testinstationarypoisson
  testpoisson
)

set(all_tests "${tests};${tests_with_arguments}")

add_executable(stokesdarcy2 EXCLUDE_FROM_ALL stokesdarcy2.cc)
add_dune_superlu_flags(stokesdarcy2)
add_executable(testcouplinggfs EXCLUDE_FROM_ALL testcouplinggfs.cc)
add_executable(testinstationarypoisson EXCLUDE_FROM_ALL testinstationarypoisson.cc)
add_test(NAME testinstationarypoisson
  COMMAND testinstationarypoisson "3" "2.0" "0.04" "0.8")
add_executable(testlocalfunctionspace EXCLUDE_FROM_ALL testlocalfunctionspace.cc)
add_executable(testmatrixconstruction EXCLUDE_FROM_ALL testmatrixconstruction.cc)
add_executable(testmortarpoisson EXCLUDE_FROM_ALL testmortarpoisson.cc)
add_executable(testmultidomaingridfunctionspace EXCLUDE_FROM_ALL testmultidomaingridfunctionspace.cc)
add_executable(testoverlappinginstationarypoisson EXCLUDE_FROM_ALL testoverlappinginstationarypoisson.cc)
add_executable(testpoisson EXCLUDE_FROM_ALL testpoisson.cc)
add_test(NAME testpoisson
  COMMAND testpoisson "4" "2.0")
add_executable(testpowermortarpoisson EXCLUDE_FROM_ALL testpowermortarpoisson.cc)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${all_tests})

foreach(_test ${tests})
  add_test(${_test} ${_test})
endforeach()

foreach(_test ${all_tests})
  target_link_libraries(${_test} ${DUNE_LIBS})
endforeach()
