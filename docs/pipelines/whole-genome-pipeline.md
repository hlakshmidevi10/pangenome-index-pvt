# Whole-Genome Pipeline (Per-Contig + Merge)

Use this for large graphs where per-contig processing is more practical.

## 1) Process each contig separately

For each contig (for example `chr1`, `chr2`, ...):

```bash
CHR=chr10
THREADS=16

gbz_extract -t ${THREADS} -b -p ${CHR}.gbz > ${CHR}_info
grlbwt-cli -t ${THREADS} -T ${PWD}/tmp -f 0 ${CHR}_info
./bin/build_tags ${CHR}.gbz ${CHR}_info.rl_bwt ${CHR}.tags
```

## 2) Build whole-genome r-index

```bash
gbz_extract -t 32 -b -p whole_genome.gbz > whole_genome_info
grlbwt-cli -t 32 -T ${PWD}/tmp whole_genome_info
./bin/build_rindex whole_genome_info.rl_bwt > whole.ri
```

## 3) Merge contig tags

```bash
./bin/merge_tags whole_genome.gbz whole.ri tags/
```

The merged output is used by downstream query tools.

Next:
- [coordinate_translation Guide](../tools/coordinate_translation.md)

