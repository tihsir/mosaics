if(EXISTS "/workspaces/CS225/cs225git/mp_mosaics/build/test_tests-b12d07c.cmake")
  include("/workspaces/CS225/cs225git/mp_mosaics/build/test_tests-b12d07c.cmake")
else()
  add_test(test_NOT_BUILT-b12d07c test_NOT_BUILT-b12d07c)
endif()
