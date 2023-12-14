cmake_minimum_required( VERSION 3.1 FATAL_ERROR)
project(containers LANGUAGES CXX)
set(CMAKE_VERBOSE_MAKEFILE true)

# Enable ASan and LSan sanitizers
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

add_executable(test_shapes_2d app/test_shapes_2d.cpp)
target_include_directories(test_shapes_2d PUBLIC include)
target_link_libraries(test_shapes_2d gmp)

add_executable(test_shapes_3d app/test_shapes_3d.cpp)
target_include_directories(test_shapes_3d PUBLIC include)
target_link_libraries(test_shapes_3d gmp)

add_executable(detect_collision app/detect_collision.cpp)
target_include_directories(detect_collision PUBLIC include)
target_link_libraries(detect_collision gmp)

install(TARGETS detect_collision test_shapes_2d test_shapes_3d DESTINATION bin)
install(PROGRAMS demo PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION bin)
install(DIRECTORY shapes_2d shapes_3d DESTINATION bin)