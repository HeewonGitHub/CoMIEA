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




