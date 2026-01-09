cmake_minimum_required(VERSION 3.28)
project(HAlign4Tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include(CTest)  # 提供 ctest 基础设施与 add_test :contentReference[oaicite:1]{index=1}

# 如果用户没有指定 CMAKE_BUILD_TYPE，则默认设置为 Release，以便默认开启优化构建
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (default: Release)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not set -> defaulting to: ${CMAKE_BUILD_TYPE}")
else()
    message(STATUS "CMAKE_BUILD_TYPE is: ${CMAKE_BUILD_TYPE}")
endif()
get_filename_component(HALIGN4_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)


find_package(OpenMP REQUIRED)
find_package(CURL)
find_package(ZLIB)

# --- 复用你的实现：做一个 support 静态库 ---
add_library(halign4_test_support STATIC
        "${HALIGN4_ROOT}/src/utils/preprocess.cpp"
        "${HALIGN4_ROOT}/src/utils/seq_io.cpp"
        "${HALIGN4_ROOT}/src/utils/cmd.cpp"
        "${HALIGN4_ROOT}/src/utils/file_io.cpp"
        "${HALIGN4_ROOT}/src/utils/xxhash.c"
        "${HALIGN4_ROOT}/src/utils/hash.cpp"
        "${HALIGN4_ROOT}/src/consensus/selector.cpp"
        "${HALIGN4_ROOT}/src/consensus/consensus.cpp"
        "${HALIGN4_ROOT}/src/seed/minimizer.cpp"
        "${HALIGN4_ROOT}/src/mash/mash.cpp"
)

target_include_directories(halign4_test_support PUBLIC
        "${HALIGN4_ROOT}/include"
        "${HALIGN4_ROOT}/extern/CLI11/include"
        "${HALIGN4_ROOT}/extern/spdlog/include"
        "${HALIGN4_ROOT}/extern/cereal/include"
        "${HALIGN4_ROOT}/extern/doctest/include"
)

target_link_libraries(halign4_test_support PUBLIC
        OpenMP::OpenMP_CXX
)

# 在 Release 时对测试支持库也启用与主项目一致的激进优化
if (MSVC)
    target_compile_options(halign4_test_support PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/GL>
    )
    target_link_options(halign4_test_support PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
else()
    target_compile_options(halign4_test_support PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
    )
endif()

if (CURL_FOUND)
    target_link_libraries(halign4_test_support PUBLIC CURL::libcurl)
    target_compile_definitions(halign4_test_support PUBLIC HALIGN4_HAVE_LIBCURL=1)
endif()

if (ZLIB_FOUND)
    target_link_libraries(halign4_test_support PUBLIC ZLIB::ZLIB)
    target_compile_definitions(halign4_test_support PUBLIC HALIGN4_HAVE_ZLIB=1)
endif()

target_compile_definitions(halign4_test_support PUBLIC
        $<$<CONFIG:Debug>:_DEBUG>
)

# 收集当前目录下所有 test_*.cpp


add_executable(halign4_tests
        "${CMAKE_CURRENT_LIST_DIR}/doctest_main.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_consensus.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_read_fasta.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_minimizer.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_mash.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_jaccard.cpp"
)

target_link_libraries(halign4_tests PRIVATE halign4_test_support)

# 对测试可执行文件也应用 Release 优化选项
if (MSVC)
    target_compile_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/GL>
    )
    target_link_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
else()
    target_compile_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
        $<$<CONFIG:Release>:-flto>
        $<$<CONFIG:Release>:-pipe>
    )
    target_link_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:-flto>
    )
endif()

# 建议显式启用（即使 include(CTest) 已经做了也不冲突）
enable_testing()

# 注册一个测试：运行这个可执行文件
# add_test(NAME halign4_tests COMMAND halign4_tests)

add_test(NAME consensus
        COMMAND halign4_tests --test-suite=consensus)

add_test(NAME read_fasta
        COMMAND halign4_tests --test-suite=read_fasta)
set_tests_properties(consensus PROPERTIES TIMEOUT 0)

add_test(NAME minimizer
        COMMAND halign4_tests --test-suite=minimizer)

add_test(NAME mash
        COMMAND halign4_tests --test-suite=mash)

add_test(NAME jaccard
        COMMAND halign4_tests --test-suite=jaccard)
