# Pangenome-Index

A toolkit for building and querying pangenome indexes using an r-index and tag arrays. It is the indexing and MEM-finding backend for [MEM-Projection](https://github.com/hlakshmidevi10/MEM-Projection), forked from [pangenome-index](https://github.com/parsaeskandar/pangenome-index).

It finds MEMs between reads and the pangenome using a lightweight tag index (`.ltags`), and writes them in a binary format that [path-walker](https://github.com/hlakshmidevi10/path-walker) reads directly to compute node coverage.

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Building an Index](#building-an-index)
- [Finding MEMs](#finding-mems)
- [Known Limitations](#known-limitations)

## Dependencies

### Required Dependencies

- **[SDSL](https://github.com/vgteam/sdsl-lite)** (vgteam fork) - For low-level data structures. Use the vgteam fork, the original simongog repo no longer builds with current CMake.
- **[HandleGraph](https://github.com/vgteam/libhandlegraph)** - For the handle graph interface
- **[GBWT](https://github.com/jltsiren/gbwt)** - For the backend
- **[GBWTGraph](https://github.com/jltsiren/gbwtgraph)** - For the graph backend; also provides `gbz_extract` and `gbz_stats`
- **[grlBWT](https://github.com/ddiazdom/grlBWT)** - For building the run-length BWT (included as a submodule and built automatically, including `grlbwt-cli`)

### System Dependencies

- **C++ compiler** with OpenMP (GCC, or Homebrew LLVM on macOS)
- **CMake**, **Make** and **pkg-config**
- **OpenSSL**, **zstd**, **zlib**

**macOS**:
```bash
brew install llvm@17 libomp cmake zstd zlib openssl@3 pkg-config
```

**Linux (Debian/Ubuntu)**:
```bash
sudo apt install build-essential cmake pkg-config libssl-dev libzstd-dev zlib1g-dev
```

## Installation

**Using MEM-Projection?** Its `install.sh` builds this repository and all its dependencies at pinned versions; you don't need the steps below.

All dependencies install into `$HOME/include` and `$HOME/lib`. Keep that prefix: the grlBWT submodule only looks for SDSL there. The makefiles also expect `sdsl-lite` to sit next to the other repos, so clone everything into one directory.

### Step 1: Set Up the Compiler (macOS only)

Apple Clang doesn't ship OpenMP, so point the build at Homebrew LLVM before building anything. SDSL records the compiler at install time, and the other repos inherit it.

```bash
export PATH="$(brew --prefix llvm@17)/bin:$PATH"
export CC=$(brew --prefix llvm@17)/bin/clang
export CXX=$(brew --prefix llvm@17)/bin/clang++
export LDFLAGS="-L$(brew --prefix libomp)/lib -L$(brew --prefix llvm@17)/lib -L$(brew --prefix zstd)/lib -L$(brew --prefix zlib)/lib"
export CPPFLAGS="-I$(brew --prefix libomp)/include -I$(brew --prefix llvm@17)/include -I$(brew --prefix zstd)/include -I$(brew --prefix zlib)/include"
export CFLAGS="$CPPFLAGS"
```

### Step 2: Clone the Repositories

```bash
mkdir -p ~/src && cd ~/src
git clone https://github.com/vgteam/sdsl-lite.git
git clone https://github.com/vgteam/libhandlegraph.git
git clone https://github.com/jltsiren/gbwt.git
git clone https://github.com/jltsiren/gbwtgraph.git
git clone --recursive https://github.com/hlakshmidevi10/pangenome-index-pvt.git
```

### Step 3: Build the Dependencies

Order matters: SDSL first, then HandleGraph and GBWT, then GBWTGraph.

```bash
# SDSL (also writes sdsl-lite/Make.helper, which the other makefiles include)
cd ~/src/sdsl-lite
./install.sh $HOME

# HandleGraph
cd ~/src/libhandlegraph
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME ..
make -j 8 && make install

# GBWT (no install target, so copy headers and library by hand)
cd ~/src/gbwt
make -j 8
cp -r include/gbwt $HOME/include/
cp lib/libgbwt.a $HOME/lib/

# GBWTGraph (same, no install target)
cd ~/src/gbwtgraph
make -j 8
cp -r include/gbwtgraph $HOME/include/
cp lib/libgbwtgraph.a $HOME/lib/
```

### Step 4: Build Pangenome-Index

```bash
cd ~/src/pangenome-index-pvt
make -j 8
```

The build will:
1. Build the grlBWT submodule's library and `grlbwt-cli` (`deps/grlBWT/build/`)
2. Create the `lib/` directory with the library
3. Create the `bin/` directory with all executables

Running `./bin/find_mems` with no arguments should print its usage line.

### Step 5: Add the Tools to Your PATH

```bash
export PATH="$HOME/src/pangenome-index-pvt/bin:$PATH"
export PATH="$HOME/src/pangenome-index-pvt/deps/grlBWT/build:$PATH"
export PATH="$HOME/src/gbwtgraph/bin:$PATH"

# Linux only: libhandlegraph is a shared library in $HOME/lib
export LD_LIBRARY_PATH="$HOME/lib:$LD_LIBRARY_PATH"
```

## Building an Index

You need a graph in `.gbz` format. MEM-Projection's `build_index.sh` runs these steps end to end; the commands below are what it runs, for a graph `graph.gbz` with output prefix `graph`.

```bash
# Number of sequences in the graph (needed by convert_tags)
gbz_stats -i graph.gbz          # read the "Sequences:" line -> NUM_SEQ

# Extract the bidirectional sequence text
gbz_extract -b -t 16 -p graph.gbz > graph.seq

# Build the run-length BWT
grlbwt-cli graph.seq -t 16 -T ./tmp -o graph.rl_bwt

# Build the r-index
build_rindex graph.rl_bwt > graph.ri

# Build the tag arrays, then convert them to the compressed format
build_tags -k 31 graph.gbz graph.rl_bwt graph.tags
convert_tags graph.tags graph_compressed.tags --num-seq $NUM_SEQ

# Build the lightweight tag index used by find_mems
build_lightweight_tags graph_compressed.tags graph.ltags

# Path names, needed by path-walker
path_extract graph.gbz graph.paths

# GFA with the GBZ's node IDs, needed by path-walker (requires vg)
vg convert -fW --no-translation graph.gbz > graph.gfa

# Index summary
print_stats graph.ri graph_compressed.tags
```

**Outputs used by queries**: `graph.ri`, `graph.ltags`, `graph.paths`, `graph.gfa`. Keep `graph_compressed.tags` too if you want to run `find_mems` without `--lightweight-tags`.

**Optional** tag-head samples (`--tag-head-samples` in `find_mems`):
```bash
build_tag_head_samples graph.ri graph.ltags graph --s 8    # writes graph.tag_samples.s8
```

## Finding MEMs

**Purpose**: Finds Maximal Exact Matches (MEMs) between reads and the pangenome and reports where each MEM lands in the graph.

**Input Files Required**:
- R-index file (`.ri`)
- Lightweight tag index (`.ltags`), or `_compressed.tags` without `--lightweight-tags`
- Reads file (plain text, one sequence per line)

**Usage**:
```bash
find_mems <r_index.ri> <tag_index> <reads.txt> <min_mem_length> <min_occurrences> [output_prefix] [options]
```

**Example** (what MEM-Projection runs by default):
```bash
find_mems graph.ri graph.ltags reads.txt 25 1 mems --lightweight-tags --use-flipped-mems
```

**Options**:
- `--lightweight-tags`: Read a `.ltags` file. MEM-Projection always uses this.
- `--use-flipped-mems`: Right-to-left MEM finder with toehold recovery, so each MEM comes with SA[sp]. Requires `--lightweight-tags`.
- `--tag-head-samples=<file>`: Use SA samples at tag-run heads (`.tag_samples.s<N>` from `build_tag_head_samples`). Requires `--lightweight-tags`.
- `--tsv`: Also write a human-readable `_path_pos.tsv`.
- `--all-positions`: Verification only. Skips the tag array and emits every BWT position of every MEM; output gets very large.
- `--debug-stats`, `--verbose`: Extra logging.

**Output Files**:
- `<output_prefix>_path_pos_v2.bin`: MEM positions as 16-byte records, bucketed by path and orientation and sorted by path position (read by path-walker)
- `<output_prefix>_seq_id_starts.out`: Start offset of each path's bucket in the `.bin` file

The run summary (MEM counts, tag runs, timings, memory) goes to stdout.

## Known Limitations

- **Node length**: `build_tags` aborts on any node longer than 1024 bp. When building the GBZ from a GFA, use `gfa2gbwt -m 1024`. Minigraph-Cactus graphs are already chopped at 1024 bp.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
