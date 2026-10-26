# APOE Genotype and Neural Gene Expression in Cerebral Organoids

## Project Overview

This project investigates transcriptomic differences between APOE4 and APOE3 isogenic cerebral organoids using publicly available RNA-sequencing data. The analysis focuses on neuronal differentiation, synaptic gene expression, immune-associated pathways, extracellular-structure organization, and the potential relevance of REST-associated neural gene regulation to Alzheimer's disease research.

The project demonstrates a reproducible bioinformatics workflow from processed gene-count data to differential expression, pathway enrichment, visualization, and biological interpretation.

## Research Question

**How does APOE genotype associate with neuronal, synaptic, and other biological gene-expression programs in cerebral organoids?**

A secondary objective is to explore whether the observed expression patterns are relevant to neural gene-network dysregulation involving REST.

## Dataset

- **GEO accession:** GSE117588
- **Dataset title:** REST and Neural Gene Network Dysregulation in iPS Cell Models of Alzheimer's Disease (RNA-seq data set)
- **Organism:** *Homo sapiens*
- **Data type:** RNA-sequencing gene counts
- **Samples analysed:** 6 cerebral-organoid samples
- **Group sizes:** 3 APOE4 samples and 3 APOE3 isogenic samples

Dataset source: [NCBI GEO — GSE117588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117588)

### Sample Groups

| Sample | Group |
|---|---|
| ADR_11 | APOE4 (E4E4_CTR) |
| ADR_12 | APOE4 (E4E4_CTR) |
| ADR_5 | APOE4 (E4E4_CTR) |
| ADR_13 | APOE3 isogenic (E3E3_ISO) |
| ADR_14 | APOE3 isogenic (E3E3_ISO) |
| ADR_6 | APOE3 isogenic (E3E3_ISO) |

The original study investigated APOE-associated gene regulation and neural differentiation in cerebral organoids. Its findings provide the biological context for this independent computational analysis.

## Workflow

### 1. Data preparation

- Imported the processed gene-count matrix.
- Verified sample identities and genotype assignments against GEO metadata.
- Converted the gene-count matrix to a format suitable for DESeq2.
- Filtered low-expression genes before differential expression analysis.

### 2. Sample-level quality assessment

- Calculated library sizes and detected-gene counts.
- Examined principal component analysis (PCA).
- Evaluated sample-to-sample expression correlations.
- Inspected selected neuronal, immune-associated, and REST-related markers.

### 3. Differential expression analysis

Differential expression was analysed with DESeq2 using a group-based design.

The principal contrast was **APOE4 versus APOE3 isogenic organoids**.

- Positive log₂ fold change indicates higher estimated expression in APOE4.
- Negative log₂ fold change indicates higher estimated expression in APOE3.
- Differentially expressed genes were identified using an adjusted p-value threshold of 0.05.

Ensembl gene identifiers were annotated with human gene symbols and gene descriptions.

### 4. Functional enrichment

Gene set enrichment analysis (GSEA) was applied to a ranked gene list based on the DESeq2 Wald statistic.

The analysis included:

- Gene Ontology Biological Process (GO BP)
- Kyoto Encyclopedia of Genes and Genomes (KEGG)

Selected leading-edge genes were examined to identify genes contributing to representative enrichment signals.

### 5. Visualization

The project includes visualizations for:

- Sample-level PCA
- Sample-to-sample expression correlation
- Differential-expression volcano plot
- Selected marker-gene expression
- Synaptic vesicle cycle enrichment
- Myeloid leukocyte migration enrichment
- Extracellular-structure organization enrichment

## Principal Results

### 1. Differential expression

After expression filtering, 16,984 genes were retained for analysis.

At an adjusted p-value threshold of 0.05, 4,573 genes were identified as differentially expressed:

- 2,120 had positive log₂ fold changes in APOE4.
- 2,453 had negative log₂ fold changes in APOE4.

Selected results included higher expression estimates for H19 and NLRP2 in APOE4, whereas C1QA, C1QB, C1QC, CSF1R, CYBB and MRC1 showed higher expression in APOE3.

### 2. Neuronal and synaptic gene expression

The synaptic vesicle cycle (GO:0099504) showed positive enrichment toward APOE4:

- NES: +2.353
- Adjusted p-value: 2.91 × 10⁻⁹

Other enriched neuronal processes included synaptic vesicle exocytosis, synapse assembly, axon development, and vesicle-mediated transport in synapses.

Several neuronal markers, including DCX, MAPT, SYP, DLG4 and RBFOX3, also showed higher expression estimates in APOE4 organoids.

These findings suggest a difference in neuronal differentiation or synaptic gene-expression programs. They do not demonstrate that neuronal function is improved or that synaptic activity is directly increased.

### 3. Immune-associated gene expression

Myeloid leukocyte migration (GO:0097529) showed negative enrichment toward APOE3:

- NES: −2.789
- Adjusted p-value: 2.91 × 10⁻⁹

Complement-associated genes and other myeloid-associated markers also showed higher expression estimates in the APOE3 group.

These results indicate a pronounced difference in immune-associated transcription. The underlying mechanism remains uncertain and may involve cellular composition, cellular state, or regulation within particular cell populations.

### 4. Extracellular-structure organization

Extracellular structure organization (GO:0043062) was negatively enriched toward APOE3:

- NES: −2.638
- Adjusted p-value: 2.91 × 10⁻⁹

This suggests a difference in extracellular-structure-related gene expression between the genotype groups.

### 5. REST-related findings

REST was examined because of its relevance to neural gene-network regulation in the original study.

In this analysis, REST had an estimated log₂ fold change of approximately −0.523 and an adjusted p-value of approximately 0.0528.

This result did not meet the conventional adjusted significance threshold of 0.05. Moreover, REST regulatory activity, nuclear localization, and chromatin binding were not directly measured.

## Conclusion

This exploratory analysis identified differences in neuronal and synaptic, immune-associated, and extracellular-structure-related gene-expression programs between APOE4 and APOE3 isogenic cerebral organoids.

The findings provide a foundation for investigating APOE-associated transcriptomic variation in a model relevant to Alzheimer's disease research.

Further validation is needed to determine the mechanisms underlying these patterns and their functional consequences.

## Limitations

- Only three biological replicates were available per genotype.
- Differences in cellular composition or maturation may contribute to the observed expression patterns.
- RNA-sequencing and GSEA do not directly measure protein activity, neuronal function, or cellular behaviour.
- The REST expression result did not meet the adjusted significance threshold.
- The findings require validation in additional samples or independent datasets before strong mechanistic conclusions can be made.

## Reproducibility and Software

The analysis uses R and Bioconductor tools, including DESeq2, AnnotationDbi, org.Hs.eg.db, clusterProfiler and enrichplot, along with ggplot2 and pheatmap for visualization.

The project repository should contain the input-data reference, sample metadata, R analysis scripts, result tables, figures, saved R objects, and session information.

Software versions should be recorded in `session_info.txt` to support reproducibility.

## References

1. Meyer K, et al. *REST and Neural Gene Network Dysregulation in iPSC Models of Alzheimer's Disease.* Cell Reports. 2019;26(5):1112–1127.e9. [PubMed](https://pubmed.ncbi.nlm.nih.gov/30699343/)
2. NCBI Gene Expression Omnibus. GSE117588: *REST and Neural Gene Network Dysregulation in iPS Cell Models of Alzheimer's Disease (RNA-seq data set).* [Dataset record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117588)

---

**Project status:** Core exploratory differential-expression and pathway-enrichment analyses have been completed. Final reproducibility checks, script organization, and repository documentation remain part of project finalization.