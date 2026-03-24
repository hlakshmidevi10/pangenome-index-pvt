# find_mems

Finds maximal exact matches (MEMs) between reads and the indexed pangenome.

## Usage

```bash
./bin/find_mems <r_index.ri> <compressed_tags.tags> <reads.txt> <min_mem_length> <min_occurrences>
```

## Inputs
- `.ri` from `build_rindex`
- compressed `.tags` from `convert_tags` or merged whole-genome tags
- reads file (one sequence per line)

