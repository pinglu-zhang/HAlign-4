cmake_minimum_required(VERSION 3.28)
project(HAlign4Tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include(CTest)  # 提供 ctest 基础设施与 add_test :contentReference[oaicite:1]{index=1}


get_filename_component(HALIGN4_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

# 关键：让 version check 能找到 version.hpp
set(SEQAN3_INCLUDE_DIR "${HALIGN4_ROOT}/extern/seqan3/include" CACHE PATH "" FORCE)

# 关键：显式指向 seqan3-config.cmake 所在目录（避免找错）
set(seqan3_DIR "${HALIGN4_ROOT}/extern/seqan3/share/cmake/seqan3" CACHE PATH "" FORCE)
# 你也可以同时保留 prefix path（可选）
list(PREPEND CMAKE_PREFIX_PATH "${HALIGN4_ROOT}/extern/seqan3/share/cmake")

find_package(seqan3 3.0 REQUIRED)
find_package(seqan3 3.0 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(CURL)
find_package(ZLIB)

# --- 复用你的实现：做一个 support 静态库 ---
add_library(halign4_test_support STATIC
        "${HALIGN4_ROOT}/src/utils/preprocess.cpp"
        "${HALIGN4_ROOT}/src/utils/seq_io.cpp"
        "${HALIGN4_ROOT}/src/utils/cmd.cpp"
        "${HALIGN4_ROOT}/src/utils/file_io.cpp"
        "${HALIGN4_ROOT}/src/consensus/selector.cpp"
        "${HALIGN4_ROOT}/src/consensus/consensus.cpp"
)

target_include_directories(halign4_test_support PUBLIC
        "${HALIGN4_ROOT}/include"
        "${HALIGN4_ROOT}/extern/seqan3/include"
        "${HALIGN4_ROOT}/extern/CLI11/include"
        "${HALIGN4_ROOT}/extern/spdlog/include"
        "${HALIGN4_ROOT}/extern/cereal/include"
        "${HALIGN4_ROOT}/extern/doctest/include"
)

target_link_libraries(halign4_test_support PUBLIC
        seqan3::seqan3
        OpenMP::OpenMP_CXX
)

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
)

target_link_libraries(halign4_tests PRIVATE halign4_test_support)

# 建议显式启用（即使 include(CTest) 已经做了也不冲突）
enable_testing()

# 注册一个测试：运行这个可执行文件
# add_test(NAME halign4_tests COMMAND halign4_tests)

add_test(NAME consensus
        COMMAND halign4_tests --test-suite=consensus)
set_tests_properties(consensus PROPERTIES TIMEOUT 0)