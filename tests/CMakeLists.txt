# Implement test logic
function(add_tpcc_test ccfile)
  get_filename_component(file ${ccfile} NAME_WE)
  add_executable(${file} ${ccfile})
  target_link_libraries(${file} PUBLIC tpcc)
  add_test(${file} ${file})
endfunction()

# Explicitly add the test
add_tpcc_test(binomial.cc)
add_tpcc_test(combination_01.cc)
add_tpcc_test(combination_02.cc)
add_tpcc_test(element_01.cc)
add_tpcc_test(lexicographic_01.cc)
add_tpcc_test(lexicographic_02.cc)
add_tpcc_test(lexicographic_03.cc)

# These tests currently do not pass
#add_tpcc_test(slab_01.cc)
#add_tpcc_test(slab_02.cc)
#add_tpcc_test(slab_03.cc)
