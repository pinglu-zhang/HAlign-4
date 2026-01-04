#!/usr/bin/env bash
set -euo pipefail

# 运行位置无关：总是以脚本所在目录为 test 源码目录
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

BUILD_DIR_DEFAULT="${ROOT_DIR}/build-test"
BUILD_DIR="${BUILD_DIR_DEFAULT}"
BUILD_TYPE="Debug"
GENERATOR=""
JOBS=""
CLEAN=0
VERBOSE=1
RUN_PERF=0
CTEST_ARGS=()

usage() {
  cat <<EOF
Usage: $(basename "$0") [options] [-- <ctest args...>]

Options:
  -B, --build-dir <dir>     Build directory (default: ${BUILD_DIR_DEFAULT})
  -t, --type <Debug|Release>Build type (default: Debug)
  -G, --generator <name>    CMake generator (e.g. "Ninja", "Unix Makefiles")
  -j, --jobs <n>            Parallel build jobs (passed to cmake --build -j)
  --clean                   Remove build directory before configuring
  --no-verbose              Run ctest without -V
  --perf                    Enable perf tests (export HALIGN4_RUN_PERF=1)
  --suite) SUITE="$2"; shift 2;;
  --file)  SOURCE_FILE="$2"; shift 2;;
  -h, --help                Show help

Examples:
  ./run_tests.sh
  ./run_tests.sh -t Release -j 16
  ./run_tests.sh --perf -- -R consensus   (pass extra args to ctest)
  ./run_tests.sh -- --output-on-failure
EOF
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -B|--build-dir) BUILD_DIR="$2"; shift 2;;
    -t|--type) BUILD_TYPE="$2"; shift 2;;
    -G|--generator) GENERATOR="$2"; shift 2;;
    -j|--jobs) JOBS="$2"; shift 2;;
    --clean) CLEAN=1; shift;;
    --no-verbose) VERBOSE=0; shift;;
    --perf) RUN_PERF=1; shift;;
    -h|--help) usage; exit 0;;
    --) shift; CTEST_ARGS+=("$@"); break;;
    *) echo "Unknown option: $1"; usage; exit 1;;
  esac
done

if [[ "${CLEAN}" -eq 1 ]]; then
  echo "[run_tests] cleaning build dir: ${BUILD_DIR}"
  rm -rf -- "${BUILD_DIR}"
fi

mkdir -p -- "${BUILD_DIR}"

# Configure (always run configure to pick up new/changed test sources reliably)
echo "[run_tests] configuring..."
CMAKE_CONFIG_ARGS=(
  -S "${SCRIPT_DIR}"
  -B "${BUILD_DIR}"
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
  -DBUILD_TESTING=ON
)

if [[ -n "${GENERATOR}" ]]; then
  CMAKE_CONFIG_ARGS+=(-G "${GENERATOR}")
fi

cmake "${CMAKE_CONFIG_ARGS[@]}"

# Build
echo "[run_tests] building..."
CMAKE_BUILD_ARGS=(--build "${BUILD_DIR}")
if [[ -n "${JOBS}" ]]; then
  CMAKE_BUILD_ARGS+=(-j "${JOBS}")
fi
cmake "${CMAKE_BUILD_ARGS[@]}"

# Run doctest directly (bypass ctest) when filtering by suite and/or source file.
if [[ -n "${SUITE:-}" || -n "${SOURCE_FILE:-}" ]]; then
  ...
  echo "[run_tests] running doctest (filters: suite='${SUITE:-}', file='${SOURCE_FILE:-}')"
  ...
  if [[ -n "${SUITE:-}" ]]; then
    DOCTEST_ARGS+=(-ts="${SUITE}")
  fi
  if [[ -n "${SOURCE_FILE:-}" ]]; then
    DOCTEST_ARGS+=(--source-file="${SOURCE_FILE}")   # doctest 支持 --source-file 过滤 :contentReference[oaicite:1]{index=1}
  fi
  ...
fi



# Run tests
echo "[run_tests] running ctest..."
if [[ "${RUN_PERF}" -eq 1 ]]; then
  export HALIGN4_RUN_PERF=1
fi

CTEST_CMD=(ctest --test-dir "${BUILD_DIR}")
if [[ "${VERBOSE}" -eq 1 ]]; then
  CTEST_CMD+=(-V)
fi
# Forward extra ctest args
CTEST_CMD+=("${CTEST_ARGS[@]}")

"${CTEST_CMD[@]}"
