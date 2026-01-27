# HAlign-4
## Requirements
### Build tools

- CMake >= 3.28
- A C++20 compiler (GCC/Clang/MSVC)

### Optional dependencies

- OpenMP (currently REQUIRED by the CMake configuration)
- libcurl (optional: used for downloading URL inputs; if missing, runtime may fall back to shell tools like `curl`/`wget` depending on your environment)
- zlib (optional: used for reading `.gz` FASTA)

---

## Build

This project uses CMake. A typical out-of-source build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The produced binary is:

- `build/halign4`

> Tip
>
> The root `CMakeLists.txt` defaults `CMAKE_BUILD_TYPE` to `Release` if you don’t set it.

---

## Run

`halign4` requires three mandatory arguments:

- `-i, --input` : input FASTA path (**must exist**; validated by CLI11 and runtime checks)
- `-o, --output`: output FASTA path
- `-w, --workdir`: working directory for intermediate files

Example (with the bundled example data):

```bash
./build/halign4 \
  -i example/data/mt1x.fasta.gz \
  -o out.fasta \
  -w work
```

### Important workdir behavior

In the current implementation, whether a non-empty `workdir` is allowed depends on the build mode (see `checkOption` in `src/halign4.cpp`):

- **Release**: `workdir` must be empty (to avoid overwriting prior results)
- **Debug**: `workdir` may be non-empty (more convenient for iterative development)

### Important CLI options

Defined in `include/config.hpp`:

- `-t, --thread <int>`: number of threads (default: hardware concurrency)
- `--cons-n <int>`: number of sequences selected in preprocess Top-K (by length) (default: 1000)
- `--kmer-size <int>`: minimizer k-mer size (default: 15; range: 4..31)
- `--kmer-window <int>`: minimizer window size (default: 10)
- `--sketch-size <int>`: Mash/MinHash sketch size (default: 2000)
- `-c, --center-path <path>`: optional; a center/reference FASTA file path (must exist)
- `-p, --msa-cmd <path>`: **MSA command template file path** (must exist; note: this is not a raw template string)
- `--keep-first-length`: keep only the first/center sequence length unchanged
- `--keep-all-length`: keep all center sequences’ lengths unchanged

### MSA command template (`-p/--msa-cmd`)

#### 1) Default template

The code ships with a built-in default template (see `DEFALT_MSA_CMD` in `include/config.hpp`):

```text
minipoa {input} -S -t {thread} -r1 > {output}
```

#### 2) Custom template file

You can pass a **text file path** via `-p/--msa-cmd`. The program will:

1. Validate the file exists (CLI11 + `file_io::requireRegularFile`)
2. Run a template self-check during argument validation (`cmd::testCommandTemplate`)

Template substitution rules:

- `{input}` and `{output}` are required
- `{thread}` (if present) will be replaced by the thread count
- Execution uses the system shell (`std::system`), so input should be considered trusted

> Note: Because `-p/--msa-cmd` is validated with `ExistingFile` in the current implementation,
> you can’t pass a raw string like `"minipoa ..."` directly; you must pass a template **file path**.

---

## Tests

Unit tests are implemented with **doctest** and built via `test/CMakeLists.txt`.

### Recommended: use the helper script

```bash
cd test
./run_tests.sh -t Release -j 8
```

The script will:

- configure and build `build-test/`
- run CTest
- when `--perf` is provided, export `HALIGN4_RUN_PERF=1` (to enable longer perf tests)

Example:

```bash
cd test
./run_tests.sh --perf -- -R align
```

### Run CTest directly

```bash
ctest --test-dir build-test -V
```

Registered test suites currently include (see `test/CMakeLists.txt`):

- `consensus`
- `read_fasta`
- `minimizer`
- `mash`
- `jaccard`
- `write_fasta`
- `align`
- `anchor`

---

## Citation

If you use HAlign-4 in academic work, please cite:

https://doi.org/10.1093/bioinformatics/btae718
