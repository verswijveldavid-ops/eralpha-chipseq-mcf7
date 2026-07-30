# Data provenance

| File | Source and purpose | Reference |
|---|---|---|
| `GSM365925_ER_minus_ligand_peaks.bed.gz` | Author-processed mock ERα ChIP peaks | hg18 |
| `GSM365926_ER_E2_peaks.bed.gz` | Author-processed E2 ERα ChIP peaks | hg18 |
| `Homo_sapiens.NCBI36.54.gtf.gz` | Ensembl release 54 gene annotation | NCBI36/hg18 |
| `MSigDB_Hallmark_2020.gmt` | Hallmark gene-set definitions | gene symbols |

SHA-256 checksums:

```text
a4c716754e9b313da8829e7a2dbc4ae1f90586e20cc1bcfdfae9cfab55a979b5  GSM365925_ER_minus_ligand_peaks.bed.gz
dffb9c5cf2cf8a6c7e0b010554764def9186b18d290a562f478a4bd33a000856  GSM365926_ER_E2_peaks.bed.gz
e8c8429741fe23f453e476f2d698d9eae3ae747a3d36d3cb372aa21f5ac07066  Homo_sapiens.NCBI36.54.gtf.gz
4275592957a1587652092bb398cf77216fde5b8daa2aedaa0e016f7d10bbdb81  MSigDB_Hallmark_2020.gmt
```

The peak files were deposited by the study authors in GEO. The annotation is
from Ensembl release 54. Keeping peak and gene coordinates on the same hg18
reference avoids cross-build overlap errors.

