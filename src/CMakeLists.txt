add_library(SE3 SHARED SE3.cpp)

add_executable(test_SE3 main.cpp)
target_link_libraries(test_SE3 SE3)
