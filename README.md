# Pangenome-Index

A comprehensive toolkit for building and querying pangenome indices using tag arrays and r-index structures. This project provides efficient algorithms for indexing pangenome graphs and performing various queries including MEM (Maximal Exact Match) finding and tag-based queries.

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Large Graph Processing](#large-graph-processing)
- [Executables](#executables)
    - [build_tags](#build_tags)
    - [build_rindex](#build_rindex)
    - [merge_tags](#merge_tags)
    - [convert_tags](#convert_tags)
    - [compress_tags](#compress_tags)
    - [find_mems](#find_mems)
    - [query_tags](#query_tags)
    - [path_extract](#path_extract)
- [File Formats](#file-formats)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

## Dependencies

### Required Dependencies

- **[SDSL](https://github.com/vgteam/sdsl-lite)** (vgteam fork) - For low-level data structures
- **[GBWT](https://github.com/jltsiren/gbwt)** - For the backend
- **[GBWTGraph](https://github.com/jltsiren/gbwtgraph)** - For the graph-based backend
- **[HandleGraph](https://github.com/vgteam/libhandlegraph)** - For the handle graph interface
- **[grlBWT](https://github.com/ddiazdom/grlBWT)** - For the run-length BWT interface (automatically built during compilation)

Please refer to each dependency's repository for installation instructions.

### System Dependencies

- **C++ Compiler**: GCC 7+ or Clang 6+
- **OpenMP**: For parallel processing
- **CMake**: For building grlBWT
- **Make**: For building the project

**macOS Users**: Install OpenMP via Homebrew:
```bash
brew install libomp
```

## Installation

### Step 1: Install Dependencies

First, install all required dependencies in accessible locations. The makefile looks for dependencies in `$(LIB_DIR)` which should be set to the location where you installed the dependencies. Follow the installation instructions provided in the [Dependencies](#dependencies) section above.

### Step 2: Clone and Build Pangenome-Index

```bash
# Clone the repository with submodules
git clone --recursive https://github.com/parsaeskandar/pangenome-index.git
cd pangenome-index

# Edit the Makefile to set the correct path to your SDSL installation
# Change this line in the Makefile:
# SDSL_DIR ?= /Users/seeskand/Documents/sdsl-lite
# To point to your SDSL installation location

# Build the project
make -j 8
```

The build process will:
1. Automatically build grlBWT dependency
2. Create the `bin/` directory with all executables
3. Create the `lib/` directory with the library
4. Create the `obj/` directory with object files

## Quick Start

For small to medium-sized graphs, you can run the pipeline with `build_tags`. If you plan to use `find_mems`, you must convert the tag arrays to the compressed format first using `convert_tags` or `compress_tags`.

### Step 1: Prepare Your Graph

You need a graph file in `.gbz` format. If you don't have one, you can create it from your graph using GBWTGraph tools.

### Step 2: Extract Sequences and Create RL-BWT

```bash
# Extract sequences from the graph (requires gbwtgraph installation)
gbz_extract -t 8 -b your_graph.gbz > graph_info

# Create the run-length BWT file using grlbwt
grlbwt-cli -t 8 graph_info
# This creates graph_info.rl_bwt
```

### Step 3: Build Tag Arrays

```bash
# Build tag arrays using the graph and RL-BWT files (algorithm format)
./bin/build_tags your_graph.gbz graph_info.rl_bwt output.tags

# Convert to compressed format required by find_mems (choose one of the following):
# Option 1: Convert existing tags
./bin/convert_tags output.tags output_compressed.tags

# Option 2: Compress tags with r-index integration
./bin/build_rindex graph_info.rl_bwt > output.ri
./bin/compress_tags output.ri output.tags output_compressed
```

## Large Graph Processing

For very large graphs (e.g., whole human genome), the pipeline should be run per-chromosome and then merged. This approach reduces memory usage and enables parallel processing.

### Step 1: Process Each Chromosome

Create a script to process each chromosome independently:

```bash
#!/bin/bash
# Process a single chromosome
CHR=$1
THREADS=16

# Create chromosome-specific directory
mkdir -p "${CHR}"
cd "${CHR}"

# Extract sequences from chromosome graph
gbz_extract -t ${THREADS} -b -p ${CHR}.gbz > ${CHR}_info

# Create run-length BWT for chromosome
grlbwt-cli -t ${THREADS} -T ${PWD}/../tmp -f 0 ${CHR}_info

# Build tag arrays for chromosome
./bin/build_tags ${CHR}.gbz ${CHR}_info.rl_bwt ${CHR}.tags

cd ..
```

Run this for each chromosome (e.g., chr1, chr2, ..., chrX, chrY).

### Step 2: Create Whole-Genome R-Index

```bash
# Extract sequences from whole-genome graph
gbz_extract -t 32 -b -p whole_genome.gbz > whole_genome_info

# Create whole-genome run-length BWT
grlbwt-cli -t 32 -T ${PWD}/tmp whole_genome_info

# Build whole-genome r-index
./bin/build_rindex whole_genome_info.rl_bwt > whole.ri
```

### Step 3: Merge Chromosome Tag Arrays

```bash
# Merge all chromosome tag arrays into whole-genome index
./bin/merge_tags whole_genome.gbz whole.ri tags/
```

**Directory Structure After Processing:**
```
project/
├── whole_genome.gbz
├── whole.ri
├── tags/
│   ├── chr1.tags
│   ├── chr2.tags
│   ├── ...
│   └── chrY.tags
└── whole_genome_tag_array_compressed.tags
```

## Executables

### build_tags

**Purpose**: Constructs tag arrays over the pangenome graph.

**Input Files Required**:
- `.gbz` file: The pangenome graph in GBZ format
- `.rl_bwt` file: Run-length BWT file created by grlbwt

**How to Create Input Files**:
```bash
# Create .gbz file (if you don't have one)
# This depends on your graph format - convert to GBZ using GBWTGraph tools

# Create .rl_bwt file from .gbz
gbz_extract -t 8 -b your_graph.gbz > graph_info
grlbwt-cli -t 8 graph_info
```

**Usage**:
```bash
./bin/build_tags <graph.gbz> <graph_info.rl_bwt> <output.tags>
```

**Example**:
```bash
./bin/build_tags test_data/bidirectional_test/xy.gbz test_data/bidirectional_test/contigs_xy.rl_bwt xy_bidirectional.tags
```

**Output**:
- `output.tags`: Binary file containing the tag arrays index (algorithm format)
- To use with `find_mems`, convert to compressed format: `./bin/convert_tags output.tags output_compressed.tags`

**What it does**:
1. Loads the pangenome graph from the GBZ file
2. Computes unique k-mers in the graph (default k=31)
3. Builds a B+ tree structure for efficient k-mer lookup
4. Creates tag arrays that map k-mers to their positions in the graph
5. Serializes the tag arrays to the output file

---

### build_rindex

**Purpose**: Builds the r-index over the pangenome from an RL-BWT file.

**Input Files Required**:
- `.rl_bwt` file: Run-length BWT file created by grlbwt

**Usage**:
```bash
./bin/build_rindex <graph_info.rl_bwt>
```

**Example**:
```bash
./bin/build_rindex test_data/bidirectional_test/contigs_xy.rl_bwt
```

**Output**:
- R-index data printed to stdout (can be redirected to a file)

**What it does**:
1. Loads the run-length BWT file
2. Constructs the r-index data structure
3. Serializes the r-index to stdout

---

### merge_tags

**Purpose**: Merges multiple tag arrays from different chromosomes into a whole-genome tag array index.

**Input Files Required**:
- Whole-genome r-index file
- Directory containing chromosome-specific `.rl_bwt` files
- Multiple `.tags` files (one per chromosome)

**Usage**:
```bash
./bin/merge_tags <whole_genome_rindex> <rl_bwt_directory> <output_merged.tags>
```

**Output**:
- `merged_output.tags`: Merged tag arrays file

**What it does**:
1. Loads the whole-genome r-index
2. Reads multiple chromosome-specific RL-BWT files
3. Merges tag arrays from different chromosomes
4. Creates a unified tag array index for the whole genome

---

### convert_tags

**Purpose**: Converts tag arrays written in the algorithm format to the compressed format required by `find_mems`.

**Input Files Required**:
- Tag arrays file (`.tags`) produced by `build_tags`

**Usage**:
```bash
./bin/convert_tags <input.tags> <output_compressed.tags>
```

**Example**:
```bash
./bin/convert_tags test_data/bidirectional_test/xy_bidirectional.tags xy_bidirectional_compressed.tags
```

**Output**:
- `output_compressed.tags`: Compressed tag arrays index

**Notes**:
- This step is required before running `find_mems` if you built the tags directly using `build_tags`.

---

### compress_tags

**Purpose**: Compresses tag arrays using r-index integration for optimized storage and query performance.

**Input Files Required**:
- R-index file (`.ri`) - Used for endmarker information and sequence count
- Tag arrays file (`.tags`) produced by `build_tags`

**Usage**:
```bash
./bin/compress_tags <r_index_file> <input_tags_file> <output_prefix>
```

**Example**:
```bash
# Build r-index first
./bin/build_rindex graph_info.rl_bwt > output.ri

# Compress tags with r-index integration
./bin/compress_tags output.ri input.tags compressed_output 2> compression.log
```

**Output Files**:
- `<output_prefix>_compressed.tags`: Main compressed tag arrays file
- Temporary files (automatically cleaned up):
    - `<output_prefix>_encoded_starts.bin`: Encoded starting positions
    - `<output_prefix>_bwt_intervals.bin`: BWT interval mappings

**What it does**:
1. Loads the r-index to get sequence and endmarker information
2. Processes tag arrays in batches to manage memory usage
3. Applies advanced compression using r-index structure
4. Merges compressed data into final optimized format
5. Produces highly compressed output suitable for `find_mems`

---

### find_mems

**Purpose**: Finds Maximal Exact Matches (MEMs) between query sequences and the pangenome, reporting unique tags and positions of those MEMs with advanced filtering and output formatting.

**Input Files Required**:
- R-index file (`.ri`)
- Compressed tag arrays index file (`.tags` created using `compress_tags` or `convert_tags`)
- Reads file (text file with one sequence per line)

**Usage**:
```bash
./bin/find_mems <r_index.ri> <compressed_tags.tags> <reads.txt> <min_mem_length> <min_occurrences> [output_file] [--debug-stats] [--verbose]
```

**Parameters**:
- `min_mem_length`: Minimum length of MEMs to report
- `min_occurrences`: Minimum number of occurrences in the pangenome
- `output_file` (optional): File to write results to (otherwise prints to stdout)
- `--debug-stats` (optional): Enable detailed duplicate statistics reporting
- `--verbose` or `--debug` (optional): Enable verbose debug output (read sequences, MEM details)

**Example**:
```bash
# Basic usage with output to stdout
./bin/find_mems output.ri compressed.tags reads.txt 30 1

# With output file and debug statistics
./bin/find_mems output.ri compressed.tags reads.txt 30 1 results --debug-stats

# With verbose debug output for profiling
./bin/find_mems output.ri compressed.tags reads.txt 30 1 results --verbose

# With both debug flags
./bin/find_mems output.ri compressed.tags reads.txt 30 1 results --debug-stats --verbose

# Real example from your usage
./bin/find_mems s28cc_flo1_chr1_186000_214000.ri s28cc_flo1_chr1_186000_214000_compressed.tags s28cc_flo1_chr1_186000_214000_N500_R1_200_reads.txt 30 1 s28cc_flo1_chr1_186000_214000_N500_R1_200
```

**Output Files** (when output_file specified):
- `<output_file>_path_pos.tsv`: Sorted MEM positions by sequence ID and node ID
- `<output_file>_seq_id_starts.out`: Starting positions for each sequence ID in the sorted file

**Output Format**:
Each line contains: `seq_id \t node_id \t offset \t is_reverse \t mem_length \t mem_start \t read_id`


**What it does**:
1. Loads the r-index and compressed tag arrays
2. Processes each read to find all MEMs above the minimum length
3. Uses r-index to locate MEM positions in the pangenome
4. Queries tag arrays to get graph positions (node, offset, strand)
5. Filters duplicates and applies occurrence thresholds
6. Sorts results by sequence ID and node ID for efficient access
7. Generates summary files for quick sequence ID lookups

---

### query_tags

**Purpose**: Queries the pangenome using tag arrays for full-read lookups and prints per-read results.

**Input Files Required**:
- R-index file (`.ri`)
- Compressed tag arrays index file (`.tags` created either using merge_tags or convert_tags)
- Reads file (text file with one sequence per line)

**Usage**:
```bash
./bin/query_tags <r_index.ri> <compressed_tags.tags> <reads.txt>
```

**Example**:
```bash
./bin/query_tags test_data/bidirectional_test/xy.ri test_data/bidirectional_test/xy_bidirectional_compressed.tags test_data/bidirectional_test/test_reads.txt
```

**Output**:
- Per-read summary lines printed to stdout (read index, length, BWT interval, number of tag runs)

**What it does**:
1. Loads the r-index and compressed tag arrays
2. Reads query sequences (one per line)
3. Performs lookups of each read in the r-index
4. Queries tag arrays for the BWT interval of each read
5. Prints a summary per read

Note: query_tags requires the compressed tag arrays format. If you built tags directly with build_tags, run:
```bash
./bin/convert_tags output.tags output_compressed.tags
```

---

### path_extract

**Purpose**: Extracts all path names from a GBZ pangenome graph file and writes them to an output file. This is useful for understanding the structure of your pangenome and identifying available paths/sequences.

**Input Files Required**:
- `.gbz` file: The pangenome graph in GBZ format

**Usage**:
```bash
./bin/path_extract <gbz_file> <output_file>
```

**Example**:
```bash
# Extract paths from a chromosome graph
./bin/path_extract s28cc_flo1_chr1_186000_214000.gbz s28cc_flo1_chr1_186000_214000.paths
```

**Output**:
- `<output_file>`: Text file containing one path name per line
- Console output: Statistics about the graph (sequence count, path count, node count, edge count)

**What it does**:
1. Loads the GBZ graph file
2. Reports graph statistics to stderr:
    - Total number of sequences
    - Number of paths in the graph
    - Number of nodes and edges
3. Extracts all path names using GBWT metadata
4. Handles different path types (reference samples, generic paths)
5. Writes path names to the output file (one per line)
6. Also prints path names to stdout for immediate viewing

---

## File Formats

### Input Files

1. **`.gbz` files**: Pangenome graphs in GBZ format (GBWTGraph format)
2. **`.rl_bwt` files**: Run-length BWT files created by grlbwt
3. **`.ri` files**: R-index files (serialized r-index data)
4. **`.tags` files**: Tag arrays index files (binary format, both algorithm and compressed formats)
5. **Text files**: Plain text files with one sequence per line

### Output Files

1. **`.tags` files**: Binary tag arrays index files (algorithm format from build_tags)
2. **`_compressed.tags` files**: Compressed tag arrays (from compress_tags or convert_tags)
3. **`.ri` ufiles**: R-index data files
4. **`.tsv` files**: Tab-separated MEM results from find_mems
5. **`.paths` files**: Path name lists from path_extract
6. **Console output**: MEM results, query results, validation reports, graph statistics

## Examples

### Complete Pipeline Example

```bash
# 1. Prepare your graph (assuming you have a graph in GBZ format)
GRAPH_FILE="your_graph.gbz"

# 2. Extract path information (optional, for understanding your graph)
./bin/path_extract $GRAPH_FILE graph_paths.txt

# 3. Extract sequences and create RL-BWT
gbz_extract -t 8 -b $GRAPH_FILE > graph_info
grlbwt-cli -t 8 graph_info

# 4. Build tag arrays and r-index
./bin/build_tags $GRAPH_FILE graph_info.rl_bwt output.tags
./bin/build_rindex graph_info.rl_bwt > output.ri

# 5. Compress tags for optimal performance
./bin/compress_tags output.ri output.tags output_compressed 2> compression.log

# 6. Find MEMs with detailed output
./bin/find_mems output.ri output_compressed_compressed.tags your_reads.txt 30 1 results --debug-stats
```

### Working with Test Data

```bash
# Use the provided test data
./bin/build_tags test_data/x.giraffe.gbz test_data/x.rl_bwt test_output.tags
./bin/build_rindex test_data/x.rl_bwt > test_output.ri
./bin/compress_tags test_output.ri test_output.tags test_output_compressed
./bin/find_mems test_output.ri test_output_compressed_compressed.tags test_data/small_test_nl.txt 30 1 test_results
```

### Processing Real Genomic Data

```bash
# Example with chromosome data (based on your usage)
PREFIX="s28cc_flo1_chr1_186000_214000"

# Extract paths to understand the graph structure
./bin/path_extract ${PREFIX}.gbz ${PREFIX}.paths

# Build indices
./bin/build_rindex ${PREFIX}.rl_bwt > ${PREFIX}.ri
./bin/build_tags ${PREFIX}.gbz ${PREFIX}.rl_bwt ${PREFIX}.tags

# Compress tags for optimal performance
./bin/compress_tags ${PREFIX}.ri ${PREFIX}.tags ${PREFIX} 2> compress-tags-${PREFIX}.log

# Find MEMs with specific parameters
./bin/find_mems ${PREFIX}.ri ${PREFIX}_compressed.tags ${PREFIX}_N500_R1_200_reads.txt 30 1 ${PREFIX}_N500_R1_200_results
```

## Troubleshooting

### Common Issues

1. **Dependency not found errors**:
    - Ensure all dependencies are installed and built
    - Check that `SDSL_DIR` in the Makefile points to the correct location
    - Verify that library files are in the expected locations

2. **OpenMP errors on macOS**:
    - Install libomp: `brew install libomp`
    - The Makefile should automatically detect and configure OpenMP

3. **Memory issues with large graphs**:
    - For very large graphs, consider processing per-chromosome
    - Use the `merge_tags` tool to combine chromosome-specific results
    - Use `compress_tags` instead of `convert_tags` for better memory efficiency

4. **File format errors**:
    - Ensure `.gbz` files are valid GBWTGraph format
    - Verify `.rl_bwt` files are created correctly with grlbwt
    - Check that input text files have proper line endings

### Getting Help

If you encounter issues:

1. Check that all dependencies are properly installed
2. Verify file formats and paths
3. Try with the provided test data first
4. Check the console output for specific error messages

## Citation

If you use this software in your research, please cite:

Parsa Eskandar, Benedict Paten, and Jouni Sirén: Lossless Pangenome Indexing Using Tag Arrays. bioRxiv 2025.05.12.653561, 2025. DOI: 10.1101/2025.05.12.653561

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.