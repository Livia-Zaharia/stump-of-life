# Suogang: two ancient individuals, no MT reads in the public BAMs

`SuogangB1` (`HRR051945`) and `SuogangB3` (`HRR051946`) are the two Penghu
individuals in Yang et al., Science 2020, DOI
[10.1126/science.aba0909](https://doi.org/10.1126/science.aba0909).
Their identities and accession relationships are verified from saved official
[B1](https://ngdc.cncb.ac.cn/gsa-human/browse/runDetail/HRR051945) and
[B3](https://ngdc.cncb.ac.cn/gsa-human/browse/runDetail/HRR051946) run pages.

Both complete public BAMs were downloaded (5,864,714 bytes combined), and their
MD5 checksums match the archive's `md5sum.txt`. **Neither BAM contains any MT
alignments.** Their 62,825 and 38,651 alignment records cover autosomes 1–22 only,
although their inherited headers list other contigs, including MT. The extracted
`*.MT.bam` files are consequently empty of reads. No complete mitochondrial
FASTA or mitochondrial read file was identified in the publicly listed deposit.
No mitochondrial consensus is generated or reference-filled.

For the artwork these can be observed ancient individuals on a site/date diagram,
with mitochondrial sequence unavailable in the retrieved deposit. They cannot
yet supply measured ancient DNA mutation rings.

## Dates and location

The primary paper's Table 1 gives both individuals as approximately **4,800–4,300
cal BP**, at rounded site coordinates **23.5° N, 119.6° E**. These coordinates are
appropriate for an approximate regional map, not an excavation plan.

The footnote and supplementary archaeological account state that the bones were
**not directly dated**. Associated shells were dated to **4,633–4,287 cal BP** and
**4,793–4,407 cal BP**. We retain both contextual ranges without assigning either
range to a particular individual. BP is relative to 1950. These are archaeological
context dates, not precise birth/death dates or ages of mitochondrial haplogroups.

The supplement reports contamination in both samples; published population
analyses used damage-restricted fragments. The fetched BAMs are not represented
as uncontaminated mitochondrial sequences.

Primary sources saved here:

- [Main paper](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/YangSciecne2020.pdf), Table 1 on page 12.
- [Supplement](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/YangSciecne2020Supplement.pdf), archaeological account on S10 and contamination methods on S2–S3.
- [Public archive directory](https://download.cncb.ac.cn/gsa-human/HRA000123/).

`manifest.json` records checksums, coverage results, exact source URLs and
limitations. `individuals.csv` provides the two map-ready metadata rows.

Reproduce from the repository root, without changing project dependencies:

```sh
uv run --no-project --with pysam --with typer python scripts/fetch_suogang.py
```
