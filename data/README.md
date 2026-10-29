# Research data

Raw downloads are local files, excluded from Git. JSON/CSV/TSV/Markdown manifests
and the retrieval scripts are versioned. Git LFS is reserved for large artwork
assets, rather than rehosting public sequencing archives.

The acquisition completed on 30 October 2026 (local date): **550 modern complete
mitogenomes; one published complete ancient mitogenome; 47 additional ancient
mitochondrial read subsets**. Both Suogang individuals' public BAMs were also
downloaded, but contain no mitochondrial alignments. The local research data,
source metadata and map cache occupy about **293 MB**. See [inventory.json](inventory.json).

- `modern/taiwan-panel/`: Ko et al. 2014 panel, **one ancient Liangdao sequence
  plus 550 modern complete mitochondrial sequences**. The downloaded FASTA is
  about 9 MB. `uv run stump-of-life fetch-taiwan` retrieves it again.
- `ancient/liangdao/`: Liangdao individuals, Matsu; availability and data type are
  recorded in the site manifest.
- `ancient/suogang/`: Suogang individuals, Penghu; availability and data type are
  recorded in the site manifest.
- `ancient/hanben_gongguan/`: Hanben and Gongguan individuals; availability and
  data type are recorded in the site manifest.

A complete mitochondrial FASTA, an extracted mitochondrial BAM, and a complete
genome BAM are different data products. A BAM may have sparse mitochondrial
coverage. Unknown or unsupported nucleotides must never be filled using the
reference merely to create a complete-looking art input.
