# File Formats

## Inputs
- `.gbz`: pangenome graph (GBWTGraph format)
- `.rl_bwt`: run-length BWT produced by `grlbwt-cli`
- `.ri`: serialized r-index / fast locate index
- `.tags`: tag array files (algorithm format or compressed format depending on tool)
- text reads: one sequence per line
- translation tables:
  - Table 1 (name + global interval -> path_id + local interval)
  - Table 2 (src_path_id + target haplotype -> target path IDs by source interval)

## Outputs
- `.ri`: built r-index
- `.tags`: built or converted tag arrays
- tool-specific console outputs (query summaries, MEMs, coordinate translation mappings)

