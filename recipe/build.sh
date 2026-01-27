#!/usr/bin/env bash
set -euxo pipefail

# conda-build 环境下不要用 -march=native（可移植性灾难）
# 若你按我建议在 CMakeLists 增加了 HALIGN4_NATIVE_ARCH 选项，这里保持 OFF
# 同时也避免子项目被注入 native flags
export CFLAGS="${CFLAGS} -O3"
export CXXFLAGS="${CXXFLAGS} -O3"

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DBUILD_TESTING=OFF \
  -DWFA2LIB_BUILD_BENCHMARK=OFF \
  -DWFA2LIB_BUILD_TESTS=OFF \
  -DHALIGN4_NATIVE_ARCH=OFF

cmake --build build -j "${CPU_COUNT}"

# 你的工程目前没有 install() 规则，所以这里手动安装可执行文件
install -d "${PREFIX}/bin"
install -m 0755 build/halign4 "${PREFIX}/bin/halign4"
