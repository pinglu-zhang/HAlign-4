cmake_minimum_required(VERSION 3.28)
project(HAlign4)

set(CMAKE_CXX_STANDARD 20)

set(SEQAN3_DIR "${CMAKE_SOURCE_DIR}/extern/seqan3/share/cmake/seqan3")
# set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/extern/seqan3/share/cmake/seqan3")
set(SEQAN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/extern/seqan3/include")

list(PREPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/extern/seqan3/share/cmake")

include_directories(
        ${PROJECT_SOURCE_DIR}/extern/seqan3/include
)

# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)

add_executable(halign4 main.cpp)

target_link_libraries (halign4 seqan3::seqan3)
