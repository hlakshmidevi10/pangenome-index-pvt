# merge_tags

Merges per-contig tags into a whole-genome tag array.

## Usage

```bash
./bin/merge_tags <whole_genome.gbz> <whole_genome.ri> <tags_dir>
```

## Inputs
- Whole-genome graph (`.gbz`)
- Whole-genome r-index (`.ri`)
- Directory containing per-contig tag outputs

## Output
- Merged whole-genome tags file

See full prep steps:
- [Whole-Genome Pipeline](../pipelines/whole-genome-pipeline.md)

