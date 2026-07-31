# CoMIEA

**CoMIEA: Continuous Phenotype-Coherent Molecular Interplay Enrichment Analysis**

---

## Overview

This repository provides the R implementation and toy dataset for **Continuous Phenotype-Coherent Molecular Interplay Enrichment Analysis (CoMIEA)**.

CoMIEA identifies biological pathways enriched with molecular interplays that change coherently along a continuous phenotypic spectrum. The method integrates sample-specific regulatory coefficients, gene expression, and multiple network-topological characteristics to calculate topology-weighted regulatory activities. It then evaluates whether genes belonging to a biological pathway exhibit stronger molecular-interplay changes between phenotypically distant samples than between phenotypically similar samples.

Statistical significance is assessed using a gene permutation strategy that compares the observed pathway score with scores obtained from randomly selected gene sets of the same size.

---

## Repository Structure

```text
CoMIEA/
├── README.md
├── sampleR_CoMIEA.r
└── ToyDATA_CoMIEA/
    ├── EXP.csv
    ├── M.csv
    ├── PWgenes.csv
    ├── BETA_Sample1.csv
    ├── BETA_Sample2.csv
    ├── ...
    └── BETA_Sample50.csv
```

---

## CoMIEA Workflow

```text
Gene expression data and continuous phenotype values
                        ↓
Sample-specific regulatory coefficient matrices
                        ↓
Regulatory Effect (RE)
                        ↓
Outgoing Eigenvector Centrality (EC)
                        ↓
Louvain module detection and Participation Coefficient (PC)
                        ↓
Outgoing Communicability (CM)
                        ↓
Topology weight
                        ↓
Topology-weighted regulatory activity (THETA)
                        ↓
Pairwise phenotype-similarity weighting
                        ↓
Gene-level phenotype-coherence score
                        ↓
Pathway-level mean coherence score
                        ↓
Gene-set permutation test
                        ↓
Empirical permutation p-value
```

---

## Software Requirements

The analysis is implemented in **R** and requires the following packages:

```r
install.packages("igraph")
install.packages("expm")
install.packages("data.table")
```

Load the packages using:

```r
library(igraph)
library(expm)
library(data.table)
```

---

## Input Files

| File | Description |
|---|---|
| `EXP.csv` | Gene expression matrix used to calculate sample-specific regulatory effects |
| `M.csv` | Continuous phenotype or modulator values for the analyzed samples |
| `PWgenes.csv` | Gene symbols belonging to the biological pathway being tested |
| `BETA_Sample1.csv`–`BETA_Sample50.csv` | Sample-specific directed regulatory coefficient matrices |

The row names of `EXP.csv` and `M.csv` should represent the same samples and should appear in the same order. The row and column names of each `BETA` matrix should correspond to the analyzed genes.

---

## Running the Example

Run the analysis from R:

```r
source("sampleR_CoMIEA.r")
```

or from a terminal:

```bash
Rscript sampleR_CoMIEA.r
```

The current script uses Windows-style paths:

```r
ToyDATA_CoMIEA\\EXP.csv
```

For macOS or Linux, replace `\\` with `/`:

```r
ToyDATA_CoMIEA/EXP.csv
```

---

# Complete CoMIEA Analysis Code

## Step 0: Define Network Topology Functions

This step defines three utility functions used to quantify the topology of each sample-specific regulatory network.

- `directed_outgoing_EC()` calculates regulator-oriented outgoing eigenvector centrality from the absolute coefficient matrix.
- `calculate_PC_with_louvain()` detects network modules using the Louvain algorithm and calculates each gene's participation coefficient.
- `calculate_outgoing_communicability()` uses the matrix exponential to quantify the total outgoing communicability of each gene.

Each topology measure is min-max normalized to make the different network characteristics comparable.

```r
########################################################################################################################
##############################
## Eigenvector centrality 
##############################
directed_outgoing_EC <- function(A) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- abs(A)
  diag(A) <- 0

  if (sum(A) == 0) {
    EC <- rep(0, nrow(A))
    names(EC) <- rownames(A)
    return(EC)
  }

  eig <- eigen(A)
  idx <- which.max(Re(eig$values))
  EC <- abs(Re(eig$vectors[, idx]))
  EC_min <- min(EC)
  EC_max <- max(EC)
  if (EC_max > EC_min) {
    EC <- (EC - EC_min)/(EC_max - EC_min)
  } else {
    EC[] <- 0
  }
  names(EC) <- rownames(A)
  EC
}
##############################
## Participation Coefficient 
##############################
calculate_PC_with_louvain <- function(A) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- abs(A)
  diag(A) <- 0
  ##############################
  ## Step 1: Louvain module detection
  ##############################
  A_undirected <- (A + t(A)) / 2
  if (sum(A_undirected) == 0) {
    module_membership <- seq_len(nrow(A))
    names(module_membership) <- rownames(A)
  } else {

    g <- graph_from_adjacency_matrix(
      A_undirected,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )

    lv <- cluster_louvain(
      g,
      weights = E(g)$weight
    )
    module_membership <- membership(lv)
    module_membership <- module_membership[rownames(A)]
  }
  ##############################
  ## Step 2: Participation Coefficient
  ##############################
  modules <- sort(unique(module_membership))
  PC_raw <- numeric(nrow(A))
  names(PC_raw) <- rownames(A)
  for (j in seq_len(nrow(A))) {
    k_j <- sum(A[j, ])
    if (k_j == 0) {
      PC_raw[j] <- 0
      next
    }
    PC_raw[j] <- 1 - sum(
      sapply(
        modules,
        function(m) {
          genes_in_module <- which(module_membership == m)
          k_jm <- sum(A[j, genes_in_module])
          (k_jm / k_j)^2
        }
      )
    )
  }
  ##############################
  ## Step 3: Min-max normalization
  ##############################
  if (max(PC_raw) > min(PC_raw)) {
    PC_normalized <-
      (PC_raw - min(PC_raw)) /
      (max(PC_raw) - min(PC_raw))
  } else {
    PC_normalized <- rep(0, length(PC_raw))
    names(PC_normalized) <- names(PC_raw)
  }
  list(
    module = module_membership,
    PC_raw = PC_raw,
    PC_normalized = PC_normalized
  )
}
##############################
## Communicability 
##############################
calculate_outgoing_communicability <- function(A) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- abs(A)
  diag(A) <- 0
  if (is.null(rownames(A))) {
    rownames(A) <- paste0("G", seq_len(nrow(A)))
  }
  if (is.null(colnames(A))) {
    colnames(A) <- rownames(A)
  }
  G <- expm::expm(A)
  diag(G) <- 0
  CM_raw <- rowSums(G)
  names(CM_raw) <- rownames(A)
  ##############################
  ## Min-max normalization
  ##############################
  CM_min <- min(CM_raw, na.rm = TRUE)
  CM_max <- max(CM_raw, na.rm = TRUE)
  if (is.finite(CM_min) &&
      is.finite(CM_max) &&
      CM_max > CM_min) {
    CM_normalized <- (CM_raw - CM_min) /
                     (CM_max - CM_min)
  } else {
    CM_normalized <- rep(0, length(CM_raw))
    names(CM_normalized) <- names(CM_raw)
  }
  list(
    communicability_matrix = G,
    CM_raw = CM_raw,
    CM_normalized = CM_normalized
  )
}
```

---

## Step 1: Load Input Data and Required Packages

This step loads the toy dataset and the R packages required for the CoMIEA analysis.

- `EXP.csv` contains gene expression values.
- `M.csv` contains the continuous phenotype or modulator values.
- `PWgenes.csv` contains the genes assigned to the tested biological pathway.
- `igraph` is used for network analysis and Louvain module detection.
- `expm` is used to calculate matrix-exponential-based communicability.

```r
##################################################
### Infile Data
##################################################
EXP<-read.table(paste("ToyDATA_CoMIEA\\EXP.csv",sep=""),sep=",")
Modulator<-read.table(paste("ToyDATA_CoMIEA\\M.csv",sep=""),sep=",")
PWgenes<-read.table(paste("ToyDATA_CoMIEA\\PWgenes.csv",sep=""),sep=",")[,1]
library(igraph)
library(expm)
```

---

## Step 2: Set Analysis Parameters and Initialize Output

This step specifies the number of genes, samples, and permutation scores used in the example. It also initializes the `SCORE` matrix that stores the observed pathway score, permutation scores, and the final p-value column.

```r
## Parameters
p <- 250
n <- 50
NoPM<-201
eps <- 1e-8

SCORE<-matrix(seq(0,0,length=(NoPM+1)),ncol=(NoPM+1))
colnames(SCORE)<-c(paste("pm",c(1:NoPM),sep=""),"Pvalue")
```

---

## Step 3: Load Sample-Specific Regulatory Coefficient Matrices

For each sample, this step loads a sample-specific `BETA` matrix. Rows represent regulators, columns represent target genes, and nonzero coefficients represent estimated directed regulatory effects. The `THETA` matrix is initialized to store the topology-weighted regulatory activity of every gene across all samples.

```r
for (i in 1:n) { # Open loof i
BETA<-as.matrix(read.table(paste("ToyDATA_CoMIEA\\BETA_Sample",i,".csv",sep=""),sep=","))
gene_names <-colnames(BETA)

THETA<-matrix(numeric(p*n),ncol=n)
rownames(THETA)<-gene_names
colnames(THETA)<-paste0("Sample",1:n)
```

---

## Step 4: Calculate Regulatory Effect

This step extracts the nonzero regulator-target coefficients from the sample-specific `BETA` matrix. Each coefficient is multiplied by the expression level of its regulator in the corresponding sample. The resulting edge-level regulatory effects are summed by regulator to obtain a gene-level regulatory effect (`RE`). Both raw and min-max-normalized RE values are retained.

```r
##############################
## Regulatory Effect
##############################
library(data.table)
estNP<-c(0,0,0)
for (c in 1:p){
if (sum(BETA[,c]!=0)>0){
estNP<-rbind(estNP,cbind(rownames(BETA)[BETA[,c]!=0],colnames(BETA)[c],BETA[,c][BETA[,c]!=0]))
}
}
estNP<-estNP[-1,]
estNP<-cbind(estNP,abs(as.numeric(estNP[,3])))
colnames(estNP)<-c("RG","TG","COEF","absCOEF")
rownames(estNP)<-paste(estNP[,1],"_",estNP[,2],sep="")
estNP<-cbind(estNP,as.numeric(estNP[,"absCOEF"])*EXP[rownames(Modulator)[i],estNP[,1]])
colnames(estNP)[5]<-"RE"
library(data.table)
dt <- as.data.table(estNP)
result <- dt[, .(RE_sum = sum(as.numeric(RE), na.rm = TRUE)), by = RG]
mat_result <- as.matrix(result)
RE<-matrix(numeric(p*1),ncol=1)
rownames(RE)<-colnames(BETA)
RE[mat_result[,1],1]<-as.matrix(as.numeric(mat_result[,2]),ncol=1)

## Min-max normalization of RE
RE <- cbind(RE_raw = RE[,1],
RE_norm = (RE[,1] - min(RE[,1])) /(max(RE[,1]) - min(RE[,1])))
```

---

## Step 5: Calculate Eigenvector Centrality

The absolute sample-specific coefficient matrix is converted into a weighted adjacency matrix. Outgoing eigenvector centrality is then calculated to quantify whether a regulator connects to other topologically influential genes. The resulting values are normalized to the interval from 0 to 1.

```r
##############################
## Eigenvector centrality 
##############################
directed_outgoing_EC <- function(A) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- abs(A)
  diag(A) <- 0

  if (sum(A) == 0) {
    EC <- rep(0, nrow(A))
    names(EC) <- rownames(A)
    return(EC)
  }

  eig <- eigen(A)
  idx <- which.max(Re(eig$values))
  EC <- abs(Re(eig$vectors[, idx]))
  EC_min <- min(EC)
  EC_max <- max(EC)
  if (EC_max > EC_min) {
    EC <- (EC - EC_min)/(EC_max - EC_min)
  } else {
    EC[] <- 0
  }
  names(EC) <- rownames(A)
  EC
}
```

---

## Step 6: Calculate Participation Coefficient

The directed coefficient matrix is converted to an undirected weighted network for Louvain module detection. The participation coefficient measures how broadly each gene distributes its outgoing connections across the detected network modules. Genes connected to several modules receive higher participation coefficients.

```r
##############################
## Participation Coefficient 
##############################
calculate_PC_with_louvain <- function(A) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- abs(A)
  diag(A) <- 0
  ##############################
  ## Step 1: Louvain module detection
  ##############################
  A_undirected <- (A + t(A)) / 2
  if (sum(A_undirected) == 0) {
    module_membership <- seq_len(nrow(A))
    names(module_membership) <- rownames(A)
  } else {

    g <- graph_from_adjacency_matrix(
      A_undirected,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )

    lv <- cluster_louvain(
      g,
      weights = E(g)$weight
    )
    module_membership <- membership(lv)
    module_membership <- module_membership[rownames(A)]
  }
  ##############################
  ## Step 2: Participation Coefficient
  ##############################
  modules <- sort(unique(module_membership))
  PC_raw <- numeric(nrow(A))
  names(PC_raw) <- rownames(A)
  for (j in seq_len(nrow(A))) {
    k_j <- sum(A[j, ])
    if (k_j == 0) {
      PC_raw[j] <- 0
      next
    }
    PC_raw[j] <- 1 - sum(
      sapply(
        modules,
        function(m) {
          genes_in_module <- which(module_membership == m)
          k_jm <- sum(A[j, genes_in_module])
          (k_jm / k_j)^2
        }
      )
    )
  }
  ##############################
  ## Step 3: Min-max normalization
  ##############################
  if (max(PC_raw) > min(PC_raw)) {
    PC_normalized <-
      (PC_raw - min(PC_raw)) /
      (max(PC_raw) - min(PC_raw))
  } else {
    PC_normalized <- rep(0, length(PC_raw))
    names(PC_normalized) <- names(PC_raw)
  }
  list(
    module = module_membership,
    PC_raw = PC_raw,
    PC_normalized = PC_normalized
  )
}
```

---

## Step 7: Calculate Communicability and Topology-Weighted Regulatory Activity

Outgoing communicability is calculated from the matrix exponential of the absolute coefficient matrix. Eigenvector centrality, normalized participation coefficient, and normalized communicability are combined using a logistic transformation to produce the topology weight (`W_TOPO`). The raw regulatory effect is multiplied by this topology weight to obtain the topology-weighted regulatory activity stored in `THETA`.

```r
##############################
## Communicability 
##############################
calculate_outgoing_communicability <- function(A) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- abs(A)
  diag(A) <- 0
  if (is.null(rownames(A))) {
    rownames(A) <- paste0("G", seq_len(nrow(A)))
  }
  if (is.null(colnames(A))) {
    colnames(A) <- rownames(A)
  }
  G <- expm::expm(A)
  diag(G) <- 0
  CM_raw <- rowSums(G)
  names(CM_raw) <- rownames(A)
  ##############################
  ## Min-max normalization
  ##############################
  CM_min <- min(CM_raw, na.rm = TRUE)
  CM_max <- max(CM_raw, na.rm = TRUE)
  if (is.finite(CM_min) &&
      is.finite(CM_max) &&
      CM_max > CM_min) {
    CM_normalized <- (CM_raw - CM_min) /
                     (CM_max - CM_min)
  } else {
    CM_normalized <- rep(0, length(CM_raw))
    names(CM_normalized) <- names(CM_raw)
  }
  list(
    communicability_matrix = G,
    CM_raw = CM_raw,
    CM_normalized = CM_normalized
  )
}

##################################################
### Infile Data
##################################################
EXP<-read.table(paste("ToyDATA_CoMIEA\\EXP.csv",sep=""),sep=",")
Modulator<-read.table(paste("ToyDATA_CoMIEA\\M.csv",sep=""),sep=",")
PWgenes<-read.table(paste("ToyDATA_CoMIEA\\PWgenes.csv",sep=""),sep=",")[,1]
library(igraph)
library(expm)

## Parameters
p <- 250
n <- 50
NoPM<-201
eps <- 1e-8

SCORE<-matrix(seq(0,0,length=(NoPM+1)),ncol=(NoPM+1))
colnames(SCORE)<-c(paste("pm",c(1:NoPM),sep=""),"Pvalue")
for (i in 1:n) { # Open loof i
BETA<-as.matrix(read.table(paste("ToyDATA_CoMIEA\\BETA_Sample",i,".csv",sep=""),sep=","))
gene_names <-colnames(BETA)

THETA<-matrix(numeric(p*n),ncol=n)
rownames(THETA)<-gene_names
colnames(THETA)<-paste0("Sample",1:n)

##############################
## Regulatory Effect
##############################
library(data.table)
estNP<-c(0,0,0)
for (c in 1:p){
if (sum(BETA[,c]!=0)>0){
estNP<-rbind(estNP,cbind(rownames(BETA)[BETA[,c]!=0],colnames(BETA)[c],BETA[,c][BETA[,c]!=0]))
}
}
estNP<-estNP[-1,]
estNP<-cbind(estNP,abs(as.numeric(estNP[,3])))
colnames(estNP)<-c("RG","TG","COEF","absCOEF")
rownames(estNP)<-paste(estNP[,1],"_",estNP[,2],sep="")
estNP<-cbind(estNP,as.numeric(estNP[,"absCOEF"])*EXP[rownames(Modulator)[i],estNP[,1]])
colnames(estNP)[5]<-"RE"
library(data.table)
dt <- as.data.table(estNP)
result <- dt[, .(RE_sum = sum(as.numeric(RE), na.rm = TRUE)), by = RG]
mat_result <- as.matrix(result)
RE<-matrix(numeric(p*1),ncol=1)
rownames(RE)<-colnames(BETA)
RE[mat_result[,1],1]<-as.matrix(as.numeric(mat_result[,2]),ncol=1)

## Min-max normalization of RE
RE <- cbind(RE_raw = RE[,1],
RE_norm = (RE[,1] - min(RE[,1])) /(max(RE[,1]) - min(RE[,1])))

##############################
## Eigenvector centrality
##############################
A <- abs(BETA)
diag(A) <- 0
EC <- directed_outgoing_EC(A)
head(sort(EC, decreasing = TRUE))

##############################
## Participation Coefficient
##############################
PC_result <- calculate_PC_with_louvain(A)
module <- PC_result$module
PC <- PC_result$PC_raw
PC_normalized <- PC_result$PC_normalized

##############################
## Communicability 
##############################
CM_result <- calculate_outgoing_communicability(A)
CM <- CM_result$CM_raw
CM_normalized <- CM_result$CM_normalized

##############################
## Topology-weighted regulatory activity
##############################
W_TOPO <-
  1 / (1 + exp(-(
    EC[rownames(A)] +
    PC_normalized[rownames(A)] +
    CM_normalized[rownames(A)]
  )))
THETA[rownames(A), i] <-
  RE[rownames(A), "RE_raw"] * W_TOPO
} # Close loof i
```

---

## Step 8: Calculate the Phenotype-Coherence Score

This step compares topology-weighted regulatory activity between every pair of samples. Sample pairs with similar phenotype values receive the within-phenotype weight `wW`, whereas phenotypically distant pairs receive the between-phenotype weight `wB`. For each gene, the phenotype-coherence score is calculated as the ratio of between-phenotype weighted differences to within-phenotype weighted differences.

```r
##############################
## Phenotype-coherence score
##############################
m <- as.numeric(Modulator[,1])
THETA <- THETA[, rownames(Modulator), drop = FALSE]
h <- 2 * var(m)
sample_pairs <- combn(seq_len(ncol(THETA)), 2)
alpha <- sample_pairs[1, ]
beta  <- sample_pairs[2, ]
wW <- exp(-((m[alpha] - m[beta])^2) / h)
wB <- 1 - wW
THETA_diff2 <- (THETA[, alpha, drop = FALSE] -THETA[, beta, drop = FALSE])^2
numerator <- rowSums(sweep(THETA_diff2, 2, wB, `*`),na.rm = TRUE)
denominator <- rowSums(sweep(THETA_diff2, 2, wW, `*`),na.rm = TRUE)
S <- numerator / (denominator + eps)
S[is.na(S)] <- 0
names(S) <- rownames(THETA)
```

---

## Step 9: Perform the Gene Permutation Test and Compute the p-value

The first score is the mean phenotype-coherence score of the observed pathway genes. The remaining scores are generated from randomly sampled gene sets of the same size. The observed score is normalized relative to the permutation distribution, and an empirical permutation p-value is calculated.

```r
##############################
## Permutation pvalue
##############################
for (pm in 1:NoPM){
if (pm==1){
SCORE[,pm]<-mean(S[PWgenes])
} else if (pm>1){
SCORE[,pm]<-mean(S[sample(1:length(S),length(PWgenes))])
}
} 
nSCORE<-SCORE
nSCORE[]<-0
nSCORE[SCORE[,1]>0,1]<-SCORE[SCORE[,1]>0,1]/mean(abs(SCORE[,2:NoPM][SCORE[,2:NoPM]>0]))
nSCORE[SCORE[,1]<0,1]<-SCORE[SCORE[,1]<0,1]/mean(abs(SCORE[,2:NoPM][SCORE[,2:NoPM]<0]))
if (nSCORE[1,1]>0){Pvalue<-sum(nSCORE[1,1]<=nSCORE[-1,1][nSCORE[-1,1]>0])/(NoPM-1)};
if (nSCORE[1,1]<0){Pvalue<-sum(abs(nSCORE[1,1])<=abs(nSCORE[-1,1][nSCORE[-1,1]<0]))/(NoPM-1)};
print(Pvalue)
```

---

## Main Output Objects

| Object | Description |
|---|---|
| `BETA` | Sample-specific directed regulatory coefficient matrix |
| `estNP` | Nonzero regulator-target coefficients and edge-level regulatory effects |
| `RE` | Raw and normalized gene-level regulatory effects |
| `EC` | Normalized outgoing eigenvector centrality |
| `module` | Louvain module membership of each gene |
| `PC` | Raw participation coefficient |
| `PC_normalized` | Min-max-normalized participation coefficient |
| `CM` | Raw outgoing communicability |
| `CM_normalized` | Min-max-normalized outgoing communicability |
| `W_TOPO` | Combined topology weight obtained from EC, PC, and CM |
| `THETA` | Topology-weighted regulatory activity matrix |
| `S` | Gene-level phenotype-coherence scores |
| `SCORE` | Observed and permutation-based pathway scores |
| `nSCORE` | Normalized observed and permutation scores |
| `Pvalue` | Empirical gene-permutation p-value |

---

