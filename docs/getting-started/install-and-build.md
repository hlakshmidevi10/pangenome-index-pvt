# Install and Build

## Dependencies

Required libraries:
- [SDSL (vgteam fork)](https://github.com/vgteam/sdsl-lite)
- [GBWT](https://github.com/jltsiren/gbwt)
- [GBWTGraph](https://github.com/jltsiren/gbwtgraph)
- [HandleGraph](https://github.com/vgteam/libhandlegraph)
- [grlBWT](https://github.com/ddiazdom/grlBWT) (built automatically as part of this project)

System requirements:
- C++ compiler (GCC 7+ or Clang 6+)
- OpenMP
- CMake
- Make

macOS:
```bash
brew install libomp
```

## Build

```bash
git clone --recursive https://github.com/parsaeskandar/pangenome-index.git
cd pangenome-index
make -j 8
```

If needed, update dependency paths in `makefile` (for example `SDSL_DIR`).

Build outputs:
- `bin/` executables
- `lib/` library artifacts
- `obj/` object files

