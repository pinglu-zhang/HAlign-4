cmake_minimum_required(VERSION 3.28)
project(HAlign4)

set(CMAKE_CXX_STANDARD 20)

set(SEQAN3_DIR "${CMAKE_SOURCE_DIR}/extern/seqan3/share/cmake/seqan3")
# set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/extern/seqan3/share/cmake/seqan3")
set(SEQAN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/extern/seqan3/include")

list(PREPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/extern/seqan3/share/cmake")



include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/extern/seqan3/include
        ${PROJECT_SOURCE_DIR}/extern/CLI11/include
        ${PROJECT_SOURCE_DIR}/extern/spdlog/include
)

# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)

add_executable(halign4
        src/halign4.cpp
        src/utils/preprocess.cpp
        src/utils/seq_io.cpp
        src/utils/cmd.cpp
        src/consensus/selector.cpp
        src/utils/file_io.cpp
)

find_package(CURL)
if (CURL_FOUND)
    target_link_libraries(halign4 PRIVATE CURL::libcurl)
    target_compile_definitions(halign4 PRIVATE HALIGN4_HAVE_LIBCURL=1)
else()
    message(STATUS "libcurl not found; downloadFile will fall back to calling curl/wget at runtime")
endif()


find_package(ZLIB)
if (ZLIB_FOUND)
    target_link_libraries(halign4 PRIVATE ZLIB::ZLIB)
    target_compile_definitions(halign4 PRIVATE HALIGN4_HAVE_ZLIB=1)
else()
    message(STATUS "zlib not found; compression support will be disabled or limited")
endif()

target_link_libraries(halign4 PRIVATE seqan3::seqan3)



target_compile_definitions(halign4 PRIVATE
        $<$<CONFIG:Debug>:_DEBUG>
)
