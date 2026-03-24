# Troubleshooting

## Build issues
- Verify dependency installs (SDSL/GBWT/GBWTGraph/HandleGraph).
- Check dependency paths in `makefile` (for example `SDSL_DIR`).
- On macOS, install OpenMP (`brew install libomp`).

## Runtime issues
- Confirm all input files are from compatible builds (graph/index/tags/tables).
- Ensure the tool gets the right tag format (algorithm vs compressed).
- Start with a small interval/test dataset to validate wiring.

## Large graph runs
- Prefer per-contig processing and merge workflow:
  - [Whole-Genome Pipeline](../pipelines/whole-genome-pipeline.md)

