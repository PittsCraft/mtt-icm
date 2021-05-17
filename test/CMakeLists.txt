set (Boost_USE_STATIC_LIBS OFF)
find_package (Boost REQUIRED COMPONENTS unit_test_framework)
include_directories (${Boost_INCLUDE_DIRS})

add_executable (Boost_Tests_run test-main.cpp
        test-math.cpp ../src/math.hpp ../src/math.cpp
        test-icm.cpp ../include/mtt-icm/icm.hpp ../src/icm.cpp
        test-monte-carlo-icm.cpp ../include/mtt-icm/monte-carlo-icm.hpp ../src/monte-carlo-icm.cpp ../src/monte-carlo-icm-util.hpp)
target_link_libraries (Boost_Tests_run ${Boost_LIBRARIES})
target_link_libraries(Boost_Tests_run wtl::sfmt)