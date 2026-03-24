# query_tags

Queries compressed tags for read intervals and prints per-read summaries.

## Usage

```bash
./bin/query_tags <r_index.ri> <compressed_tags.tags> <reads.txt>
```

## Inputs
- `.ri` from `build_rindex`
- compressed `.tags` from `convert_tags` or merged whole-genome tags
- reads file (one sequence per line)

