
cmake_minimum_required( VERSION 3.14 )
project( helens )

file(GLOB_RECURSE SRCS src/*.cpp)

add_subdirectory(ThirdParty/cxxopts)

add_executable(Helens ${SRCS})
target_compile_features(Helens PUBLIC cxx_std_17)

target_link_libraries(Helens PUBLIC cxxopts::cxxopts)
