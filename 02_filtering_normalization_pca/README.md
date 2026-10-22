# Stage 2A: Filtering normalization
## 1. Code Breakdown

- Step 1: Input Validation (`stopifnot`)
The script ensures the `count_mat` matrix consists only of whole numbers ($\(\ge 0\)$) without missing values (NA), and verifies that the `gene_symbols` vector is clean. This prevents down-stream math errors in DESeq2.
- Step 2 & 3: Metadata Setup & Filtering
It groups your samples into three biological factor levels: Old, Young, and AD. It then applies an expression filter: it removes any gene that fails to reach at least 10 raw reads in at least 8 individual samples. This aggressively strips away low-count noise (e.g., genes that aren't expressed in brain tissue), boosting our statistical power for subsequent differential testing.
- Step 5: Library Size Normalization (`estimateSizeFactors`)
Raw count data cannot be compared directly because different samples are sequenced at different depths (total read counts). The estimateSizeFactors function calculates a scale factor for each sample relative to a geometric pseudo-reference sample. Our provided console output shows these factors:
  - Most samples hover closely around 1.000 (e.g., 26-3A-AD = 1.006). This means sequencing depth is highly consistent across our library.
  - Sample 2-12A-Young has a high size factor of 1.420, meaning it was sequenced deeper than average, so its raw counts will be scaled down. Sample 3-17T-Young is at 0.746 (shallower sequencing) and its counts will be scaled up.
- Step 6: Variance-Stabilizing Transformation (vst)
In raw RNA-seq data, **high-expression genes naturally have much higher variance than low-expression genes** (heteroscedasticity). The vst() algorithm decouples variance from the mean, bringing the counts onto a log2-like scale. This transformation is required before applying distance-based metrics like PCA.

## Deconstructing the Sample Pattern: 3-17T-Young

```text
   3   -   17      T   -   Young
  ───     ────    ───     ───────
   │       │       │         │
   │       │       │         └─► Clinical Condition (Young, Old, or AD)
   │       │       └───────────► Brain Region Tissue (T = Temporal Cortex, A = Anterior Cingulate)
   │       └───────────────────► Unique Patient/Subject ID (Patient #17)
   └───────────────────────────► Technical Sequencing/Library Run Number
```

- 3 (Sequencing Run/Library ID): This tracking number is assigned by the sequencing facility. It represents the physical well or index slot on the sequencing machine.
- 17 (Patient ID): This identifies the individual human brain donor. Notice how you have a 3-17T-Young sample, which means this data comes from Patient #17.
- T or A (Brain Region): In Alzheimer's transcriptomics, it is highly common to sequence multiple areas from the exact same brain to see where pathology strikes first.
	- T stands for the Temporal Cortex (or Superior Temporal Gyrus), a region heavily damaged early in Alzheimer's disease that controls memory and language processing.
	- A stands for the Anterior Cingulate Cortex (or Amygdala), a separate functional zone.
	- Proof in our data: Look at Patient #1. You have sample 21-1A-AD and sample 20-1T-AD. This means you are comparing the Cingulate vs. Temporal tissue from the exact same Alzheimer's patient.
- Young (Clinical Cohort): The diagnostic bracket determined by pathology and age at death. our dataset classifies individuals into Young (healthy young controls), Old (age-matched healthy controls), and AD (confirmed Alzheimer's cases).


## Output:

```bash
 print(round(sizeFactors(dds), 3))
   21-1A-AD    20-1T-AD    23-2A-AD    22-2T-AD    26-3A-AD    24-3T-AD
      0.896       0.995       0.994       1.090       1.006       0.897
   27-5A-AD    25-5T-AD    29-6T-AD    31-7T-AD    28-8T-AD    30-9T-AD
      0.756       0.928       0.924       0.953       0.892       1.131
 18-10A-Old  11-10T-Old  19-11A-Old  13-11T-Old  15-13T-Old  16-14T-Old
      1.056       1.061       0.898       1.240       1.037       1.121
  12-6A-Old   14-7A-Old   10-8A-Old   17-9A-Old 2-12A-Young 4-13A-Young
      1.080       0.907       1.105       1.178       1.420       1.115
6-14A-Young 8-15A-Young 9-16A-Young 3-17T-Young 5-18T-Young 7-19T-Young
      1.042       0.930       0.939       0.746       1.158       1.054
```

`DESeq2` divides our raw counts by the size factor to normalize them:
> Normalized Count = Raw Count/ Size Factor


1. The Normal Baseline (Size Factor $\(\approx \)$ 1.000)
- Sample: 26-3A-AD (Size Factor = 1.006)
- Meaning: This sample was sequenced perfectly to the average depth of our experiment.
- Calculation: If a gene has 100 raw counts, its normalized count is $\(100 / 1.006 = 99.4\)$. The numbers barely change because the library size is already standard.

2. The High Size Factor (Scale Down)
- Sample: 2-12A-Young (Size Factor = 1.420)
- Meaning: This sample got 42% more reads than the average sample just by chance during sequencing. Because the machine read deeper, every gene in this sample will artificially look like it has higher expression than it should.
- Calculation: To fix this artificial boost, we divide by 1.420. If a gene has 142 raw counts, its normalized count becomes $\(142 / 1.420 = 100\)$. We scaled it down so it can be fairly compared to the other samples.

3. The Low Size Factor (Scale Up)
- Sample: 3-17T-Young (Size Factor = 0.746)
- Meaning: This sample was sequence-starved. It only gathered roughly 74.6% of the average read volume. Every gene here looks artificially low because the machine didn't read it for long enough.
- Calculation: To fix this, we divide by 0.746 (which is mathematically the same as multiplying by ~1.34). If a gene has 74 raw counts, its normalized count becomes $\(74 / 0.746 = 99.2\)$. We scaled it up to reflect what the count would have been if it had been sequenced at a normal depth.



# 2. Interpreting the PCA Visual Output
A PCA reduces thousands of dimensions (individual gene expression profiles) down to principal components that explain the largest variations in the data.

```text
PCA Interpretation Guide
  ┌────────────────────────────────────────────────────────┐
  │  ◄── PC1: 62% Variance (Dominant Biological Split) ──► │
  │                                                        │
  │  [Cluster A: Mixed/Aged]         [Cluster B: Young]    │
  │     ● AD samples (Blue)             ● Young (Green)    │
  │     ● Old samples (Red)                                │
  └────────────────────────────────────────────────────────┘
```

- **PC1 Explains 62% of Variance:**
This is a remarkably dominant first component. It dictates the horizontal separation across the plot.
- The Main Biological Story (Horizontal Split):
	- Looking at the far right of the chart (PC1 score > 30), we will find a distinct, isolated cluster of Young samples (Green dots).
	- Looking at the left side of the chart (PC1 score < 10). Here, the AD samples (Blue dots) and the Old healthy controls (Red dots) heavily overlap.
- What this tells us biologically:
**The absolute greatest driver of global gene expression differences in this brain dataset is chronological age/senescence, not the Alzheimer's pathology itself**. The molecular profiles of aging healthy brains and AD brains are fundamentally much closer to each other than either is to a young brain.

- **PC2 Explains 9% of Variance:**
This dictates vertical separation. While there is minor separation here, the AD and Old groups are still fundamentally intermingled, indicating high sample heterogeneity or potential confounding biological variables (like sex, post-mortem interval, or brain region batch effects).




# Stage 2B: Investigate PCA structure by sample-name suffix

## Code Breakdown




## Output explanations
```bash
sample group sample_code
21-1A-AD       21-1A-AD    AD           A
20-1T-AD       20-1T-AD    AD           T
23-2A-AD       23-2A-AD    AD           A
22-2T-AD       22-2T-AD    AD           T
26-3A-AD       26-3A-AD    AD           A
24-3T-AD       24-3T-AD    AD           T
27-5A-AD       27-5A-AD    AD           A
25-5T-AD       25-5T-AD    AD           T
29-6T-AD       29-6T-AD    AD           T
31-7T-AD       31-7T-AD    AD           T
28-8T-AD       28-8T-AD    AD           T
30-9T-AD       30-9T-AD    AD           T
18-10A-Old   18-10A-Old   Old           A
11-10T-Old   11-10T-Old   Old           T
19-11A-Old   19-11A-Old   Old           A
13-11T-Old   13-11T-Old   Old           T
15-13T-Old   15-13T-Old   Old           T
16-14T-Old   16-14T-Old   Old           T
12-6A-Old     12-6A-Old   Old           A
14-7A-Old     14-7A-Old   Old           A
10-8A-Old     10-8A-Old   Old           A
17-9A-Old     17-9A-Old   Old           A
2-12A-Young 2-12A-Young Young           A
4-13A-Young 4-13A-Young Young           A
6-14A-Young 6-14A-Young Young           A
8-15A-Young 8-15A-Young Young           A
9-16A-Young 9-16A-Young Young           A
3-17T-Young 3-17T-Young Young           T
5-18T-Young 5-18T-Young Young           T
7-19T-Young 7-19T-Young Young           T
```
















