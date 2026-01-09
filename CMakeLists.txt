cmake_minimum_required(VERSION 3.28)
project(HAlign4)

set(CMAKE_CXX_STANDARD 20)

# 如果用户没有指定 CMAKE_BUILD_TYPE，则默认设置为 Release，以便默认开启优化构建
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (default: Release)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not set -> defaulting to: ${CMAKE_BUILD_TYPE}")
else()
    message(STATUS "CMAKE_BUILD_TYPE is: ${CMAKE_BUILD_TYPE}")
endif()


include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/extern/CLI11/include
        ${PROJECT_SOURCE_DIR}/extern/spdlog/include
        ${PROJECT_SOURCE_DIR}/extern/cereal/include
)

add_executable(halign4
        src/halign4.cpp
        src/utils/preprocess.cpp
        src/utils/seq_io.cpp
        src/utils/cmd.cpp
        src/utils/file_io.cpp
        src/mash/hash.cpp
        src/mash/xxhash.c
        src/consensus/selector.cpp
        src/consensus/consensus.cpp
        src/cluster/cluster.cpp
        src/seed/minimizer.cpp
        src/mash/mash.cpp
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

find_package(OpenMP REQUIRED)

target_link_libraries(halign4 PRIVATE OpenMP::OpenMP_CXX)



# 在 Release 时添加激进的优化标志：GCC/Clang 与 MSVC 分别处理
if (MSVC)
    # MSVC 优化：/O2 开启优化，/arch:AVX2 请求 AVX2 指令集（可根据需要改为 AVX 或 AVX512），/GL 全局优化
    target_compile_options(halign4 PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/GL>
    )
    # 链接时启用链接时期优化（LTCG）
    target_link_options(halign4 PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
else()
    # GCC / Clang: -O3, -march=native, -flto, -pipe
    target_compile_options(halign4 PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
    )
endif()


target_compile_definitions(halign4 PRIVATE
        $<$<CONFIG:Debug>:_DEBUG>
)
