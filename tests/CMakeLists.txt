
set(TEST_FILES catch.hpp CatchRequire.h MoleculeTest.cpp SimulationTest.cpp CellMatrixTest.cpp)
add_executable(all_tests ${TEST_FILES})
target_link_libraries(all_tests molly)
