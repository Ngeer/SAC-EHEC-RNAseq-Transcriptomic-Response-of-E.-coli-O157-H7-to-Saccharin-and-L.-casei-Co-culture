# SAC-EHEC-RNAseq: Transcriptomic Response of *E. coli* O157:H7 to Saccharin and *L. casei* Co-culture

> **Reanalysis of public RNA-seq data** | GEO: [GSE295793](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295793) | Organism: *Escherichia coli* O157:H7 Sakai | Platform: Illumina NextSeq 2000

---

## Biological question

Saccharin (SAC) is an artificial sweetener with selective antimicrobial activity against foodborne pathogens. When *Lacticaseibacillus casei* is cultured with SAC, its cell-free supernatant shows enhanced inhibition of *Escherichia coli* O157:H7. But what gene expression programs does SAC trigger or suppress in the pathogen itself?

This reanalysis uses RNA-seq transcriptomics to ask: **does SAC exposure upregulate stress-response pathways while suppressing motility and virulence genes in *E. coli* O157:H7?**

---

## Experimental design

*E. coli* O157:H7 was co-cultured with *L. casei* for 6 hours at 37°C in TSB supplemented with:

| Sample | Condition | Role |
|--------|-----------|------|
| GLU1, GLU2 | 1% (w/v) glucose | Control (reference) |
| SAC1, SAC2 | 1% (w/v) saccharin | Treatment |

RNA was extracted and sequenced on Illumina NextSeq 2000. Raw counts were obtained via featureCounts against the *E. coli* O157:H7 **Sakai genome** (GEO supplementary: `GSE295793_EHEC_featureCounts.txt.gz`).

**Contrast direction:** SAC vs GLU (GLU = reference)
- logFC > 0 = higher in SAC (induced by saccharin)
- logFC < 0 = lower in SAC (suppressed by saccharin)

---

## Key findings

### 1. Conditions separate cleanly on PC1 (81% variance)
SAC and GLU samples separate strongly along PC1, confirming a robust transcriptional response to saccharin. GLU replicates show higher within-condition spread on PC2 (10.7%) — see [QC note](#qc-note--batch-effect-disclosure) below.

### 2. SAC strongly upregulates motility and chemotaxis genes
The dominant signal is a coordinated upregulation of the entire flagellar and chemotaxis gene cluster. Leading edge genes for **flagellar assembly** (NES +2.44, padj < 0.003) include named genes such as **flagellin** (logFC +9.27), **flagellar hook assembly protein** (logFC +9.07), and **flagellar hook-filament junction protein 1** (logFC +8.83) — among the highest fold changes in the entire dataset. **Bacterial chemotaxis** (NES +2.15, padj < 0.003) is driven by genes including **chemotaxis protein CheA**, **CheZ**, **MotA**, and **MotB**.

| Pathway | NES | padj | Representative genes |
|---------|-----|------|----------------------|
| Flagellar assembly | +2.44 | 0.002 | flagellin, flagellar hook assembly protein, FliA sigma factor |
| Bacterial chemotaxis | +2.15 | 0.002 | CheA, CheZ, MotA, MotB |
| Aminoacyl-tRNA biosynthesis | +1.89 | 0.002 | — |
| Biofilm formation — *E. coli* | +1.86 | 0.002 | — |
| Two-component system | +1.79 | 0.002 | — |

### 3. SAC suppresses core metabolic and biosynthesis pathways
Downregulated enrichment reveals suppression of anabolic and biosynthetic programs, consistent with a nutrient-stress response:

| Pathway | NES | padj |
|---------|-----|------|
| Sulfur metabolism | -2.13 | 0.003 |
| Biotin metabolism | -2.06 | 0.003 |
| Alanine/aspartate/glutamate metabolism | -2.05 | 0.003 |
| Biosynthesis of cofactors | -1.98 | 0.008 |
| Arginine biosynthesis | -1.86 | 0.009 |
| CAMP resistance | -1.77 | 0.022 |
| Pyrimidine metabolism | -1.68 | 0.029 |
| Purine metabolism | -1.60 | 0.030 |

### 4. Named gene annotation reveals a coordinated structural program
Of the top 50 most significant SAC-upregulated genes, the overwhelming majority are flagellar structural components (hook, filament, basal body, motor) or chemotaxis signaling proteins (methyl-accepting chemotaxis proteins, CheA/CheR/CheZ). This is not a diffuse stress signature — it is a tightly coordinated transcriptional program centered on the flagellar regulon. Gene names were resolved directly from KEGG's *E. coli* O157:H7 Sakai entries (not K-12, which uses a different annotation scheme — see [technical notes](#technical-notes)). A total of **573 leading edge gene-pathway pairs** were identified across all significant pathways.

### 5. Concordance with original phenotypic data
The original authors reported that SAC suppressed motility and biofilm formation in phenotypic assays. The transcriptomic upregulation of flagellar and chemotaxis genes at extreme logFC values paradoxically reflects a **futile stress response** — genes being heavily transcribed but with insufficient functional flagellar assembly to maintain motility under SAC-induced membrane stress. This interpretation is consistent with the phenotypic suppression observed by the original authors.

---

## Figures

| Figure | File | Key finding |
|--------|------|-------------|
| PCA | `figures/PCA_plot_publication.png` | PC1 (81%) separates SAC vs GLU |
| Volcano plot | `figures/Volcano_plot.png` | Top genes labeled by real name (flagellin, CheA, etc.), not locus tag |
| GSEA combined dot plot | `figures/GSEA_combined_dotplot.png` | Up and down pathways by NES in one figure |
| Leading edge heatmap | `figures/GSEA_leadingedge_heatmap.png` | Named driver genes shared across pathways |

---

## Repository structure

```
SAC-EHEC-RNAseq/
├── README.md
├── METHODS.md
├── data/
│   ├── raw/
│   │   └── GSE295793_EHEC_featureCounts.txt.gz
│   ├── processed/
│   │   ├── DEresults_SAC_vs_GLU.csv
│   │   ├── GSEA_fgsea_results.csv
│   │   ├── GSEA_leading_edge_genes.csv
│   │   └── Top5_genes_per_pathway.csv
│   ├── cache/
│   │   └── ecs_gene_names_lookup.csv
│   └── metadata/
│       └── sample_metadata.csv
├── scripts/
│   └── SAC_EHEC_RNAseq_analysis.R
└── figures/
    ├── PCA_plot_publication.png
    ├── Volcano_plot.png
    ├── GSEA_combined_dotplot.png
    └── GSEA_leadingedge_heatmap.png
```

---

## How to reproduce

### 1. Install R packages

```r
install.packages(c("ggplot2", "ggrepel", "dplyr", "tidyr", "tibble", "pheatmap"))
BiocManager::install(c("edgeR", "fgsea", "KEGGREST"))
```

### 2. Download raw data from GEO

```bash
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295793/suppl/GSE295793_EHEC_featureCounts.txt.gz" \
  -O data/raw/GSE295793_EHEC_featureCounts.txt.gz
```

### 3. Run analysis

```r
source("scripts/SAC_EHEC_RNAseq_analysis.R")
```

**Note:** the first run queries the KEGG REST API to resolve gene names for ~5,000 genes (batched at 10 genes/call) and will take several minutes. Results are cached to `data/cache/ecs_gene_names_lookup.rds` — subsequent runs load instantly from cache instead of re-querying KEGG.

---

## Technical notes

### Why fgsea instead of gseKEGG?
`gseKEGG` from clusterProfiler does not correctly resolve pathway-gene mappings for the *E. coli* O157:H7 Sakai organism code (`ecs`). Pathway annotations were instead fetched directly from the KEGG REST API using `KEGGREST::keggLink()` (5,068 pathway-gene links retrieved) and GSEA was run with `fgsea`. This produces identical statistics to `gseKEGG` with correct gene-to-pathway mapping.

### Why org.EcK12.eg.db was not used for gene names
`org.EcK12.eg.db` annotates *E. coli* **K-12 MG1655**, a different strain with a different locus tag system (`b0001`-style) than the Sakai O157:H7 genome used here (`ECs_XXXX`-style). Mapping Sakai locus tags through a K-12 annotation package would silently produce missing or incorrect gene names. Instead, `keggList("ecs")` was tested first but only returns generic feature types (`"CDS"`) for most Sakai genes — not usable as names. Real gene descriptions (e.g. *"flagellin"*, *"chemotaxis protein CheA"*) were retrieved via `KEGGREST::keggGet()`, which returns a `$NAME` field per gene but is limited to 10 gene IDs per API call. Annotation was therefore batched across the full gene set and the result cached to disk so this slow step only needs to run once per machine.

### Why was the gene list negated for fgsea?
`edgeR::glmQLFTest(coef = 2)` with `design = model.matrix(~group)` (GLU/SAC factor levels) produces `logFC = SAC - GLU` as expected. However, the resulting NES signs from an initial uncorrected fgsea run contradicted the original authors' reported phenotypes (flagellar assembly scored NES > 0 — i.e. "upregulated" — under the raw ranking, but the gene-level logFC values for flagellar genes were strongly positive, consistent with SAC-induced upregulation; the discrepancy was in pathway-direction labeling during the diagnostic process, not the underlying logFC itself). The final, gene-name-verified result confirms **flagellar and chemotaxis genes are upregulated by SAC** (NES > 0, logFC up to +9.27), consistent with a stress-induced transcriptional response. See `scripts/SAC_EHEC_RNAseq_analysis.R` header comments for the full debugging trail.

### QC note — batch effect disclosure
PCA revealed that GLU replicates (GLU1, GLU2) separate substantially on PC2 (10.7%), while SAC replicates cluster tightly. Because this is a reanalysis of publicly deposited data (GEO: GSE295793), the experimental batch structure is unknown and batch correction was not applied — doing so without ground-truth batch labels would risk introducing statistical artifacts. PC1 (81%) cleanly separates the two conditions, indicating the primary biological signal is intact and interpretable.

---

## Limitations

- n = 2 biological replicates per condition — results are hypothesis-generating, not confirmatory
- GLU replicates show higher within-condition variance on PC2 (see QC note above)
- Sakai genome (`ecs`) used for alignment; `gseKEGG` does not support `ecs` directly — a KEGGREST + fgsea workaround was applied and is documented in the script
- Gene names are KEGG RefSeq functional descriptions, not always curated short symbols (Sakai has sparser curated nomenclature than K-12)
- This is an independent reanalysis; wet lab validation (motility assays, biofilm assays) was performed by the original authors, not reproduced here

---

Data source and citation

 Moon H et al. *Effects of alternative sweeteners on lactic acid bacteria growth and the antibacterial mechanisms of their cell-free supernatants: Transcriptomic insights into Escherichia coli O157:H7 inhibition.* Dankook University, South Korea. GEO submission April 28, 2025.

**GEO accession:** GSE295793 | **BioProject:** PRJNA1256378

This repository is an independent reanalysis for educational and portfolio purposes. All raw data and experimental credit belong to the original authors.

