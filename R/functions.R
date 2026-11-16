#' Loading Data
#'
#' Load data from files.
#'
#'
#' @param Candidate_gene_file_path,Candidate_peak_file_path paths to pre-selected candidate genes and peaks for the cell type
#' @param scRNA_readcount_file_path,scRNA_gene_file_path,scRNA_cellmeta_file_path paths to filtered scRNA data of the cell type
#' @param scATAC_readcount_file_path,scATAC_peak_file_path,scATAC_cellmeta_file_path paths to filtered scATAC data of the cell type
#' @param Motif_mapping_file_path,Motif_name_file_path paths to TF motif prior on all ATAC peaks
#' @param Ref_seq_file_path path to the Refseq file for transcription starting site extraction
#'
#' @return A list named "loaded_data" containing the following elements:
#' * Common_samples: a vector of common samples between scRNA and scATAC data
#' * Candidate_Genes: a data frame of candidate genes
#' * Candidate_Peaks: a data frame of candidate peaks
#' * scRNA_Genes: a data frame of scRNA genes
#' * scRNA_cells: a data frame of scRNA cells
#' * scRNA_read_count_matrix: a sparse matrix of scRNA read counts
#' * scATAC_Peaks: a data frame of scATAC peaks
#' * scATAC_cells: a data frame of scATAC cells
#' * scATAC_read_count_matrix: a sparse matrix of scATAC read counts
#' * Motifs: a data frame of TF motifs
#' * TF_peak_binding_matrix: a sparse matrix of TF binding on ATAC peaks
#' * Refseq: a data frame of Refseq
#'
#' @import Matrix
#' @importFrom Matrix colSums rowSums sparseMatrix as.matrix
#'
#' @export

Data_loading <- function(Candidate_Gene_file_path, Candidate_Peak_file_path,
                         scRNA_readcount_file_path, scRNA_Gene_file_path, scRNA_cellmeta_file_path,
                         scATAC_readcount_file_path, scATAC_Peak_file_path, scATAC_cellmeta_file_path,
                         Motif_mapping_file_path, Motif_name_file_path, Ref_seq_file_path) {
  # load candidate Genes
  Candidate_Genes <- read.table(Candidate_Gene_file_path, sep = "\t")
  colnames(Candidate_Genes) <- c("Gene_symbols")
  print("loading candidate genes ...")
  cat("\n")
  print(paste(length(Candidate_Genes$Gene_symbols), "candidate genes are provided.", sep = " "))
  cat("\n")


  # load scRNAseq data
  print("loading scRNAseq data ...")
  cat("\n")

  scRNA_Genes <- read.table(scRNA_Gene_file_path, sep = "\t")
  colnames(scRNA_Genes) <- c("Gene_index", "Gene_symbols")

  scRNA_cells <- read.table(scRNA_cellmeta_file_path, sep = "\t")
    colnames(scRNA_cells) <- c("cell_index", "cell_barcode", "cell_type", "subject_ID", "condition")

  scRNA_read_count_table <- read.table(scRNA_readcount_file_path, sep = "\t")
  scRNA_read_count_matrix <- sparseMatrix(
    i = scRNA_read_count_table[, 1], j = scRNA_read_count_table[, 2], x = scRNA_read_count_table[, 3],
    dimnames = list(levels(scRNA_read_count_table[, 1]), levels(scRNA_read_count_table[, 2])), repr = "T"
  )
  rm(scRNA_read_count_table)

  Group_conditions <- unique(scRNA_cells$condition)
  print(paste("We detected", length(Group_conditions), "conditions from the meta file", sep = " "))
  cat("\n")

  scRNA_samples <- unique(scRNA_cells$subject_ID)
  print(paste("The input scRNAseq data includes", length(scRNA_Genes$Gene_symbols), "Genes and", length(scRNA_cells$cell_barcode), "cells from", length(scRNA_samples), "samples", sep = " "))
  cat("\n")

  # print('RNA cells of each sample')
  # table(scRNA_cells$subject_ID)


  # load candidate Peaks
  print("loading candidate peaks ...")
  cat("\n")
  Candidate_Peaks <- read.table(Candidate_Peak_file_path, sep = "\t")
  colnames(Candidate_Peaks) <- c("chr", "point1", "point2")
  print(paste(nrow(Candidate_Peaks), "candidate peaks are provided.", sep = " "))
  cat("\n")


  # load scATACseq data
  print("loading scATACseq data ...")
  cat("\n")
  scATAC_Peaks <- read.table(scATAC_Peak_file_path, sep = "\t")
  scATAC_Peaks <- scATAC_Peaks[, c(1, 2, 3, 4)]
  colnames(scATAC_Peaks) <- c("Peak_index", "chr", "point1", "point2")

  scATAC_cells <- read.table(scATAC_cellmeta_file_path, sep = "\t")
  colnames(scATAC_cells) <- c("cell_index", "cell_barcode", "cell_type", "subject_ID")

  scATAC_read_count_table <- read.table(scATAC_readcount_file_path, sep = "\t")
  scATAC_read_count_matrix <- sparseMatrix(
    i = scATAC_read_count_table[, 1], j = scATAC_read_count_table[, 2], x = scATAC_read_count_table[, 3],
    dimnames = list(levels(scATAC_read_count_table[, 1]), levels(scATAC_read_count_table[, 2])), repr = "T"
  )
  rm(scATAC_read_count_table)

  scATAC_samples <- unique(scATAC_cells$subject_ID)
  print(paste("The input scATACseq data includes", length(scATAC_Peaks$Peak_index), "Peaks and", length(scATAC_cells$cell_barcode), "cells from", length(scATAC_samples), "samples", sep = " "))
  cat("\n")

  # print('ATAC cells of each sample')
  # table(scATAC_cells$subject_ID)


  Common_samples <- intersect(scRNA_samples, scATAC_samples)
  print(paste("There are sample-paired single cell multiomics data for", length(Common_samples), "samples.", sep = " "))
  print("Please check sample IDs if this number is lower than expected.")
  cat("\n")

  # load TF motif prior
  print("loading TF motif match data ...")
  cat("\n")
  Motifs <- read.table(Motif_name_file_path, sep = "\t")
  colnames(Motifs) <- c("motif_index", "name")
  Motif_mapping_table <- read.table(Motif_mapping_file_path)
  TF_Peak_binding_matrix <- sparseMatrix(
    i = Motif_mapping_table[, 1], j = Motif_mapping_table[, 2], x = Motif_mapping_table[, 3],
    dimnames = list(levels(Motif_mapping_table[, 1]), levels(Motif_mapping_table[, 2])), repr = "T"
  )
  rm(Motif_mapping_table)
  print(paste(nrow(Motifs), "motifs are screened.", sep = " "))
  cat("\n")


  # load RefSeq
  Refseq <- read.table(Ref_seq_file_path, header = TRUE, sep = "\t")
  colnames(Refseq) <- c("chr", "strand", "start", "end", "Gene_symbols")


  loaded_data <- list(
    "Common_samples" = Common_samples,
    "Candidate_Genes" = Candidate_Genes,
    "Candidate_Peaks" = Candidate_Peaks,
    "scRNA_Genes" = scRNA_Genes,
    "scRNA_cells" = scRNA_cells,
    "scRNA_read_count_matrix" = scRNA_read_count_matrix,
    "scATAC_Peaks" = scATAC_Peaks,
    "scATAC_cells" = scATAC_cells,
    "scATAC_read_count_matrix" = scATAC_read_count_matrix,
    "Motifs" = Motifs,
    "TF_Peak_binding_matrix" = TF_Peak_binding_matrix,
    "Refseq" = Refseq
  )

  return(loaded_data)
}


#' Build candidate circuits
#'
#' This function constructs candidate circuits based on the input data WITHOUT TAD prior.
#'
#' @seealso [Candidate_circuits_construction_with_TAD()] for the version WITH TAD.
#'
#' @param loaded_data The output from `Data_loading`
#' @param distance_control Bp threshold for paring peaks and genes
#' @param TF_pct_in_candidate_Peaks_cutoff Parameter to select TFs. TFs that bind to more than this percentage of all peaks will be selected. Default is 0.05.
#' @param TF_num_cutoff Parameter to select TFs. TFs that bind to more than this number of peaks will be selected. Default is 30.
#' @param TF_enrichment_FC_cutoff Parameter to select TFs. TFs that are enriched in candidate peaks with fold change (not log fold change!) higher than this cutoff will be selected. Default is 0.8.
#'
#' @return A list named "Candidate_circuits" containing the following elements:
#' * TFs: candidate TFs
#' * TF_log2Count: log2 transformed TF expression
#' * Peaks: candidate peaks
#' * Peak_log2Count: log2 transformed peak accessibility
#' * Genes: candidate genes
#' * Gene_log2Count: log2 transformed gene expression
#' * TF_peak_binding: TF-peak binding matrix
#' * Peak_Gene_looping: peak-gene looping matrix
#'
#' @import Matrix
#' @importFrom Matrix colSums rowSums sparseMatrix as.matrix
#'
#' @export
Candidate_circuits_construction_without_TAD <- function(loaded_data, distance_control, TF_pct_in_candidate_Peaks_cutoff = 0.05, TF_num_cutoff = 30, TF_enrichment_FC_cutoff = 0.8) {
  Common_samples <- loaded_data$Common_samples
  Candidate_Genes <- loaded_data$Candidate_Genes
  Candidate_Peaks <- loaded_data$Candidate_Peaks
  scRNA_Genes <- loaded_data$scRNA_Genes
  scRNA_cells <- loaded_data$scRNA_cells
  scRNA_read_count_matrix <- loaded_data$scRNA_read_count_matrix
  scATAC_Peaks <- loaded_data$scATAC_Peaks
  scATAC_cells <- loaded_data$scATAC_cells
  scATAC_read_count_matrix <- loaded_data$scATAC_read_count_matrix
  Motifs <- loaded_data$Motifs
  TF_Peak_binding_matrix <- loaded_data$TF_Peak_binding_matrix
  Refseq <- loaded_data$Refseq

  #********** load TAD prior**********
  print(paste("Naive pairing of peaks and genes if they are within", distance_control, "bps", sep = " "))
  cat("\n")
  print("MAGICAL model initialization ...")
  cat("\n")

  # TF-Peak binding
  Candidate_Peaks <- inner_join(scATAC_Peaks, Candidate_Peaks)
  index <- match(Candidate_Peaks$Peak_index, scATAC_Peaks$Peak_index)
  Candidate_TF_Peak_Binding <- TF_Peak_binding_matrix[index, ]
  TF_num <- colSums(Candidate_TF_Peak_Binding)
  TF_pct_in_candidate_Peaks <- TF_num / nrow(Candidate_Peaks)
  TF_pct_in_all_Peaks <- colSums(TF_Peak_binding_matrix) / nrow(scATAC_Peaks)
  TF_enrichment_FC <- TF_pct_in_candidate_Peaks / TF_pct_in_all_Peaks
  TF_index <- which(TF_pct_in_candidate_Peaks > TF_pct_in_candidate_Peaks_cutoff & TF_num > TF_num_cutoff & TF_enrichment_FC > TF_enrichment_FC_cutoff)
  if (length(TF_index) == 0) {
    stop("No TFs meet the criteria. Please adjust TF parameters.")
  } else {
    print(paste0(length(TF_index), " candidate TFs selected.\n
                 Please adjust TF parameters if too few/many TFs are selected."))
    Candidate_TFs <- Motifs[TF_index, ]
    Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[, TF_index]
  }
  Peak_index <- which(rowSums(Candidate_TF_Peak_Binding) > 0)
  Candidate_Peaks <- Candidate_Peaks[Peak_index, ]
  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[Peak_index, ]


  # extract Gene TSS
  Gene_symbols <- intersect(intersect(Candidate_Genes$Gene_symbols, Refseq$Gene_symbols), scRNA_Genes$Gene_symbols)
  Gene_chr <- character(length(Gene_symbols))
  Gene_TSS <- integer(length(Gene_symbols))
  Candidate_Genes <- data.frame(Gene_symbols, Gene_chr, Gene_TSS)
  colnames(Candidate_Genes) <- c("Gene_symbols", "chr", "TSS")
  for (g in 1:length(Gene_symbols)) {
    index <- which(Refseq$Gene_symbols == Candidate_Genes$Gene_symbols[g])
    Candidate_Genes$chr[g] <- Refseq$chr[index[1]]
    if (Refseq$strand[index[1]] == "+") {
      Candidate_Genes$TSS[g] <- min(Refseq$start[index])
    } else {
      Candidate_Genes$TSS[g] <- max(Refseq$end[index])
    }
  }


  # Peak-Gene looping
  Candidate_Peak_Gene_looping_distance_control <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(nrow(Candidate_Peaks), nrow(Candidate_Genes)), repr = "T")
  for (g in 1:length(Candidate_Genes$Gene_symbols)) {
    index <- which(Candidate_Peaks$chr == Candidate_Genes$chr[g] & abs((Candidate_Peaks$point1 + Candidate_Peaks$point2) / 2 - Candidate_Genes$TSS[g]) < distance_control)
    if (length(index) > 0) {
      Candidate_Peak_Gene_looping_distance_control[index, g] <- 1
    }
  }
  Candidate_Peak_Gene_looping <- Candidate_Peak_Gene_looping_distance_control
  Peak_index <- which(rowSums(Candidate_Peak_Gene_looping) > 0)
  Gene_index <- which(colSums(Candidate_Peak_Gene_looping) > 0)
  Candidate_Peak_Gene_looping <- Candidate_Peak_Gene_looping[Peak_index, Gene_index]
  Candidate_Peaks <- Candidate_Peaks[Peak_index, ]
  Candidate_Genes <- Candidate_Genes[Gene_index, ]
  TF_index <- which(colSums(Candidate_TF_Peak_Binding[Peak_index, ]) > 0)
  Candidate_TFs <- Candidate_TFs[TF_index, ]
  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[Peak_index, TF_index]


  # pseudobulk ATAC data is calculated for model initialization
  ATAC_count <- matrix(0, nrow = nrow(scATAC_read_count_matrix), ncol = length(Common_samples))
  for (s in 1:length(Common_samples)) {
    totake <- which(scATAC_cells$subject_ID == Common_samples[s])
    if(length(totake) == 1){
      ATAC_count[, s] <- scATAC_read_count_matrix[, totake]
    } else {
      ATAC_count[, s] <- rowSums(scATAC_read_count_matrix[, totake])
    }
  }
  colnames(ATAC_count) <- Common_samples
  rownames(ATAC_count) <- scATAC_Peaks$Peak_index

  total_ATAC_reads <- 5e6
  ATAC_scale_facter <- total_ATAC_reads / (colSums(ATAC_count) + 1)
  ATAC_count <- sweep(ATAC_count, 2, ATAC_scale_facter, "*")
  ATAC_log2 <- log2(ATAC_count + 1)

  index <- match(Candidate_Peaks$Peak_index, scATAC_Peaks$Peak_index)
  Candidate_Peak_log2Count <- ATAC_log2[index, ]
  Candidate_Peak_log2Count <- sweep(Candidate_Peak_log2Count, 1, rowMeans(Candidate_Peak_log2Count), "-")


  # pseudobulk RNA data is calculated for model initialization
  RNA_count <- matrix(0, nrow = nrow(scRNA_read_count_matrix), ncol = length(Common_samples))
  for (s in 1:length(Common_samples)) {
    totake <- which(scRNA_cells$subject_ID == Common_samples[s])
    if(length(totake) == 1){
      RNA_count[, s] <- scRNA_read_count_matrix[, totake]
    } else {
      RNA_count[, s] <- rowSums(scRNA_read_count_matrix[, totake])
    }
  }
  colnames(RNA_count) <- Common_samples
  rownames(RNA_count) <- scRNA_Genes$Gene_symbols

  total_RNA_reads <- 5e6
  RNA_scale_facter <- total_RNA_reads / (colSums(RNA_count) + 1)
  RNA_count <- sweep(RNA_count, 2, RNA_scale_facter, "*")
  RNA_log2 <- log2(RNA_count + 1)

  index <- match(Candidate_Genes$Gene_symbols, scRNA_Genes$Gene_symbols)
  Candidate_Gene_log2Count <- RNA_log2[index, ]
  Candidate_Gene_log2Count <- sweep(Candidate_Gene_log2Count, 1, rowMeans(Candidate_Gene_log2Count), "-")


  # TFs that are not expressed in the current data are not considered
  #selected_TFs <- intersect(Candidate_TFs$name, scRNA_Genes$Gene_symbols)
  #index <- match(selected_TFs, Candidate_TFs$name)
  #Candidate_TFs <- Candidate_TFs[index, ]
  #Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[, index]

#  index <- match(selected_TFs, scRNA_Genes$Gene_symbols)
#  Candidate_TF_log2Count <- RNA_log2[index, ]

#  TF_index <- which(rowSums(Candidate_TF_log2Count) > 0)
#  Candidate_TFs <- Candidate_TFs[TF_index, ]
#  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[, TF_index]
#  Candidate_TF_log2Count <- Candidate_TF_log2Count[TF_index, ]
#  Candidate_TF_log2Count <- sweep(Candidate_TF_log2Count, 1, rowMeans(Candidate_TF_log2Count), "-")


  # final filtering for candidate circuits
  Peak_index <- which(rowSums(Candidate_Peak_Gene_looping) > 0 & rowSums(Candidate_TF_Peak_Binding) > 0)
  Gene_index <- which(colSums(Candidate_Peak_Gene_looping[Peak_index, ]) > 0)
  TF_index <- which(colSums(Candidate_TF_Peak_Binding[Peak_index, ]) > 0)

  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[Peak_index, TF_index]
  Candidate_Peak_Gene_looping <- Candidate_Peak_Gene_looping[Peak_index, Gene_index]

  Candidate_TFs <- Candidate_TFs[TF_index, ]

#  Candidate_TF_log2Count <- Candidate_TF_log2Count[TF_index, ]

  Candidate_Peaks <- Candidate_Peaks[Peak_index, ]
  Candidate_Peak_log2Count <- Candidate_Peak_log2Count[Peak_index, ]

  Candidate_Genes <- Candidate_Genes[Gene_index, ]
  Candidate_Gene_log2Count <- Candidate_Gene_log2Count[Gene_index, ]

  Candidate_circuits <- list(
    "TFs" = Candidate_TFs,
 #   "TF_log2Count" = Candidate_TF_log2Count,
    "Peaks" = Candidate_Peaks,
    "Peak_log2Count" = Candidate_Peak_log2Count,
    "Genes" = Candidate_Genes,
    "Gene_log2Count" = Candidate_Gene_log2Count,
    "TF_Peak_Binding" = Candidate_TF_Peak_Binding,
    "Peak_Gene_looping" = Candidate_Peak_Gene_looping
  )

  return(Candidate_circuits)
}


#' MCMC model initialization
#'
#' This function initializes the MCMC model for MAGICAL.
#'
#' For MCMC method, the initial values of parameters are not deterministic for the results.
#' Use pseudo bulk data to initialize the model as only at this dimension all data are matched.
#'
#' @param loaded_data The output from `Data_loading`
#' @param Candidate_circuits The output from "Candidate_circuits" (with or without TAD)
#'
#' @return A list of initial values for the MCMC model
#'
#' @import Matrix
#' @importFrom Matrix colSums rowSums sparseMatrix as.matrix
#'
#' @export
MAGICAL_initialization <- function(
    loaded_data,
    Candidate_circuits,
    T_prior_var = 1.0,
    B_prior_var = 0.25,
    B_inclusion_prob = 0.10
) {
  Candidate_TFs <- Candidate_circuits$TFs
#  Candidate_TF_log2Count <- Candidate_circuits$TF_log2Count

  Candidate_Peaks <- Candidate_circuits$Peaks
  Candidate_Peak_log2Count <- Candidate_circuits$Peak_log2Count

  Candidate_Genes <- Candidate_circuits$Genes
  Candidate_Gene_log2Count <- Candidate_circuits$Gene_log2Count

  Candidate_TF_Peak_Binding <- Candidate_circuits$TF_Peak_Binding
  Candidate_Peak_Gene_looping <- Candidate_circuits$Peak_Gene_looping

  print(paste("T_prior_var ", T_prior_var, "\n", "B_prior_var ", B_prior_var, "\n", "B_inclusion_prob ", B_inclusion_prob, "\n"))


  S <- length(loaded_data$Common_samples)
  M <- nrow(Candidate_TFs) # total number of candidate TFs
  P <- nrow(Candidate_Peaks) # total number of candidate Peaks
  G <- nrow(Candidate_Genes) # total number of candidate Genes

  # TF activity prior: TF RNA expression
#  T_prior <- Candidate_TF_log2Count
#  T_mean <- Candidate_TF_log2Count
  T_mean <- matrix(
    0,
    nrow = M,
    ncol = S
  )
  
  T_var <- matrix(
    T_prior_var,
    nrow = M,
    ncol = S
  )

  T_prior <- matrix(
    rnorm(
      M * S,
      mean = 0,
      sd = sqrt(T_prior_var)
    ),
    nrow = M,
    ncol = S
  )
  
  
  tf_names <- Candidate_TFs$name
  sample_names <- loaded_data$Common_samples
  
  dimnames(T_mean) <- list(tf_names, sample_names)
  dimnames(T_var) <- list(tf_names, sample_names)
  dimnames(T_prior) <- list(tf_names, sample_names)


 # ------------------------------------------------------------
 # Expression-independent TF-to-peak prior
 # ------------------------------------------------------------

 # Candidate-link mask
 B_mask <- as.matrix(
     Candidate_TF_Peak_Binding > 0
 )
  storage.mode(B_mask) <- "integer"
# storage.mode(B_mask) <- "numeric"

 # Prior settings
 sigma_B <- sqrt(B_prior_var)
 pi_B <- B_inclusion_prob

 # Fixed prior mean used in every MCMC iteration
 B_mean <- matrix(
     0,
     nrow = P,
     ncol = M,
     dimnames = dimnames(B_mask)
 )

 # Fixed prior variance used in every MCMC iteration
 B_var <- matrix(
   B_prior_var,
   nrow = 1,
   ncol = M
 )

 # Fixed binary-state prior score
 B_prob <- matrix(
     0,
     nrow = P,
     ncol = M,
     dimnames = dimnames(B_mask)
 )

 B_prob[B_mask > 0] <- pi_B

 # Random starting values, drawn only once
 B_initial <- matrix(
     0,
     nrow = P,
     ncol = M,
     dimnames = dimnames(B_mask)
 )

 B_initial[B_mask > 0] <- rnorm(
     sum(B_mask > 0),
     mean = 0,
     sd = sigma_B
 )


 # B_prior <- as.matrix(Candidate_TF_Peak_Binding)
 # B_prob <- as.matrix(Candidate_TF_Peak_Binding)
 # TF_Peak_index <- which(as.matrix(Candidate_TF_Peak_Binding) > 0, arr.ind = TRUE)

 # for (b in 1:length(Candidate_TF_Peak_Binding@x)) {
 #   mdl <- summary(lm(Candidate_Peak_log2Count[TF_Peak_index[b, 1], ] ~ Candidate_TF_log2Count[TF_Peak_index[b, 2], ]))
 #   B_prior[TF_Peak_index[b, 1], TF_Peak_index[b, 2]] <- mdl$coefficients[2, 1]
 #   B_prob[TF_Peak_index[b, 1], TF_Peak_index[b, 2]] <- 1 - mdl$coefficients[2, 4]
 # }
 # B_mean <- B_prior

 # B_var <- matrix(0.5, nrow = 1, ncol = M)

 # for (m in 1:M) {
 #   if (sum(Candidate_TF_Peak_Binding[, m]) > 1) {
 #     B_var[m] <- var(B_prior[which(Candidate_TF_Peak_Binding[, m] > 0), m]) # variance of binding events weights of TF t
 #   } else {
 #     B_var[m] <- 0.5
 #   }
 # }


  # Peak-Gene looping prior
  L_prior <- as.matrix(Candidate_Peak_Gene_looping)
  L_prob <- as.matrix(Candidate_Peak_Gene_looping)
  Peak_Gene_index <- which(as.matrix(Candidate_Peak_Gene_looping) > 0, arr.ind = TRUE)

  for (l in 1:nrow(Peak_Gene_index)) {
    mdl <- summary(lm(Candidate_Gene_log2Count[Peak_Gene_index[l, 2], ] ~ Candidate_Peak_log2Count[Peak_Gene_index[l, 1], ]))
    L_prior[Peak_Gene_index[l, 1], Peak_Gene_index[l, 2]] <- mdl$coefficients[2, 1]
    L_prob[Peak_Gene_index[l, 1], Peak_Gene_index[l, 2]] <- 1 - mdl$coefficients[2, 4]
  }

  L_mean <- L_prior
  L_var <- var(L_prior[which(L_prob > 0)])

  Initial_model <- list(
    "T_prior" = T_prior,
    "T_mean" = T_mean,
    "T_var" = T_var,
   # "B_prior" = B_prior,
    "B_initial" = B_initial,
    "B_mean" = B_mean,
    "B_var" = B_var,
    "B_prob" = B_prob,
    "B_mask" = B_mask,
    "L_prior" = L_prior,
    "L_mean" = L_mean,
    "L_var" = L_var,
    "L_prob" = L_prob
  )

  return(Initial_model)
}

#************************************************************
# Joint spike-and-slab update for TF-peak states and weights
#************************************************************
TF_peak_binding_B_and_state_sampling <- function(
    A_sample,
    TFA,
    B,
    B_state,
    B_mask,
    B_prior_mean,
    B_prior_var,
    B_prior_prob,
    sigma_A_noise,
    P,
    M
) {
  # Validate dimensions
  if (!identical(dim(B), c(P, M))) {
    stop("B must have dimensions P x M.")
  }
  
  if (!identical(dim(B_state), c(P, M))) {
    stop("B_state must have dimensions P x M.")
  }
  
  if (!identical(dim(B_mask), c(P, M))) {
    stop("B_mask must have dimensions P x M.")
  }
  
  if (!identical(dim(B_prior_mean), c(P, M))) {
    stop("B_prior_mean must have dimensions P x M.")
  }
  
  if (!identical(dim(B_prior_prob), c(P, M))) {
    stop("B_prior_prob must have dimensions P x M.")
  }
  
  if (length(B_prior_var) != M) {
    stop("B_prior_var must contain one variance per TF.")
  }
  
  if (!identical(dim(TFA$T_sample), c(M, ncol(A_sample)))) {
    stop("TFA$T_sample must have dimensions M x ncol(A_sample).")
  }
  
  if (!is.finite(sigma_A_noise) || sigma_A_noise <= 0) {
    stop("sigma_A_noise must be finite and positive.")
  }
  
  # This Bayes-factor expression assumes a zero-centered slab.
  candidate_mask <- B_mask > 0
  
  if (any(abs(B_prior_mean[candidate_mask]) > 1e-12, na.rm = TRUE)) {
    stop("This implementation assumes B_prior_mean = 0 for candidate TF-peak edges.")
  }
  
  # Force binary states and zero forbidden or inactive weights.
  B_mask <- 1L * candidate_mask
  B_state <- 1L * (B_state > 0)
  
  B_state[B_mask == 0L] <- 0L
  B[B_mask == 0L] <- 0
  B[B_state == 0L] <- 0
  
  # Current ATAC prediction:
  # B is P x M, T_sample is M x S, and pred_A is P x S.
  pred_A <- B %*% TFA$T_sample
  
  # Visit only motif-permitted TF-peak edges.
  candidate_edges <- which(B_mask > 0L, arr.ind = TRUE)
  
  if (nrow(candidate_edges) == 0L) {
    stop("No candidate TF-peak edges are present.")
  }
  
  # Random-scan update order.
  edge_order <- sample.int(nrow(candidate_edges))
  sigma2 <- sigma_A_noise
  
  for (edge_number in edge_order) {
    f <- candidate_edges[edge_number, 1L]
    m <- candidate_edges[edge_number, 2L]
    
    # TF m activity across samples.
    x <- as.numeric(TFA$T_sample[m, ])
    old_B <- B[f, m]
    
    # pred_A currently includes old_B * x. Adding it back produces
    # the residual under a model that excludes the current edge.
    residual_without_edge <- as.numeric(A_sample[f, ]) - as.numeric(pred_A[f, ]) + old_B * x
    
    tau2 <- B_prior_var[m]
    pi_edge <- B_prior_prob[f, m]
    
    if (!is.finite(tau2) || tau2 <= 0) {
      stop(sprintf("B_prior_var[%d] must be finite and positive.", m))
    }
    
    if (!is.finite(pi_edge) || pi_edge <= 0 || pi_edge >= 1) {
      stop(sprintf("B_prior_prob[%d, %d] must be strictly between 0 and 1.", f, m))
    }
    
    # Sufficient statistics.
    x2 <- sum(x^2)
    xy <- sum(x * residual_without_edge)
    
    if (!is.finite(x2) || !is.finite(xy)) {
      stop(sprintf("Non-finite sufficient statistics for TF-peak edge [%d, %d].", f, m))
    }
    
    # Slab posterior conditional on an active edge.
    # There is deliberately no division by S.
    posterior_var <- 1 / (1 / tau2 + x2 / sigma2)
    posterior_mean <- posterior_var * xy / sigma2
    
    if (!is.finite(posterior_var) || posterior_var <= 0) {
      stop(sprintf("Invalid posterior variance for TF-peak edge [%d, %d].", f, m))
    }
    
    if (!is.finite(posterior_mean)) {
      stop(sprintf("Invalid posterior mean for TF-peak edge [%d, %d].", f, m))
    }
    
    # Integrated active-versus-inactive evidence.
    # Inactive model: B[f, m] = 0.
    # Active model: B[f, m] ~ Normal(0, tau2).
    log_bayes_factor <- -0.5 * log1p(tau2 * x2 / sigma2) +
      0.5 * tau2 * xy^2 / (sigma2 * (sigma2 + tau2 * x2))
    
    # Combine prior inclusion odds with evidence from the ATAC data.
    log_prior_odds <- qlogis(pi_edge)
    log_posterior_odds <- log_prior_odds + log_bayes_factor
    inclusion_probability <- plogis(log_posterior_odds)
    
    if (!is.finite(inclusion_probability)) {
      stop(sprintf("Invalid inclusion probability for TF-peak edge [%d, %d].", f, m))
    }
    
    # Sample the binary state and the corresponding weight.
    new_state <- rbinom(1L, size = 1L, prob = inclusion_probability)
    
    if (new_state == 1L) {
      new_B <- rnorm(1L, mean = posterior_mean, sd = sqrt(posterior_var))
    } else {
      new_B <- 0
    }
    
    B_state[f, m] <- new_state
    B[f, m] <- new_B
    
    # Replace the old edge contribution with the new one.
    pred_A[f, ] <- pred_A[f, ] + (new_B - old_B) * x
  }
  
  return(list(
    B = B,
    B_state = B_state,
    pred_A = pred_A
  ))
}

# # Functions within `MAGICAL_estimation`, no need to export
# #*************************** TFA estimation ***************************
TF_activity_T_sampling <- function(
    A,
    A_sample,
    ATAC_Cell_Sample_vector,
    R,
    R_sample,
    RNA_Cell_Sample_vector,
    TFA,
    T_prior_mean,
    T_prior_var,
    B,
    B_state,
    sigma_A_noise,
    P,
    G,
    M,
    S,
    clip_sample_draw = TRUE,
    clip_limit = 3
) {
  if (!is.finite(sigma_A_noise) || sigma_A_noise <= 0) {
    stop("sigma_A_noise must be finite and positive.")
  }
  
  if (!identical(dim(A_sample), c(P, S))) {
    stop("A_sample must have dimensions P x S.")
  }
  
  if (!identical(dim(TFA$T_sample), c(M, S))) {
    stop("TFA$T_sample must have dimensions M x S.")
  }
  
  if (!identical(dim(T_prior_mean), c(M, S))) {
    stop("T_prior_mean must have dimensions M x S.")
  }
  
  if (!identical(dim(T_prior_var), c(M, S))) {
    stop("T_prior_var must have dimensions M x S.")
  }
  
  if (!identical(dim(B), c(P, M))) {
    stop("B must have dimensions P x M.")
  }
  
  if (!identical(dim(B_state), c(P, M))) {
    stop("B_state must have dimensions P x M.")
  }
  
  if (any(!is.finite(T_prior_var)) || any(T_prior_var <= 0)) {
    stop("All T_prior_var values must be finite and positive.")
  }
  
  if (!is.finite(clip_limit) || clip_limit <= 0) {
    stop("clip_limit must be finite and positive.")
  }
  
  TF_index <- sample.int(M)
  fitted_A <- B %*% TFA$T_sample
  
  for (m in TF_index) {
    b <- B[, m]
    old_T <- TFA$T_sample[m, ]
    
    if (any(B_state[, m] > 0L)) {
      residual_without_tf <- A_sample - fitted_A + b %*% t(old_T)
      
      variance_T <- 1 / (1 / T_prior_var[m, ] + sum(b^2) / sigma_A_noise)
      
      mean_T <- variance_T * (
        T_prior_mean[m, ] / T_prior_var[m, ] +
          as.numeric(crossprod(b, residual_without_tf)) / sigma_A_noise
      )
    } else {
      # With no active links, the ATAC data provide no information
      # about this TF, so draw its activity from the prior.
      mean_T <- T_prior_mean[m, ]
      variance_T <- T_prior_var[m, ]
    }
    
    if (any(!is.finite(variance_T)) || any(variance_T <= 0)) {
      stop(sprintf("Invalid posterior TFA variance for TF %d.", m))
    }
    
    if (any(!is.finite(mean_T))) {
      stop(sprintf("Invalid posterior TFA mean for TF %d.", m))
    }
    
    for (s in seq_len(S)) {
      # This is the only draw clipped in the original function.
      aa <- rnorm(1L)
      
      if (isTRUE(clip_sample_draw)) {
        aa <- max(-clip_limit, min(clip_limit, aa))
      }
      
      TFA$T_sample[m, s] <- mean_T[s] + aa * sqrt(variance_T[s])
      
      ATAC_cell_index <- which(ATAC_Cell_Sample_vector == s)
      
      if (length(ATAC_cell_index) > 0L) {
        TFA$T_A[m, ATAC_cell_index] <- rnorm(
          length(ATAC_cell_index),
          mean = TFA$T_sample[m, s],
          sd = sqrt(variance_T[s])
        )
      }
      
      RNA_cell_index <- which(RNA_Cell_Sample_vector == s)
      
      if (length(RNA_cell_index) > 0L) {
        TFA$T_R[m, RNA_cell_index] <- rnorm(
          length(RNA_cell_index),
          mean = TFA$T_sample[m, s],
          sd = sqrt(variance_T[s])
        )
      }
    }
    
    # The next TF update must use the newly sampled activity for TF m.
    fitted_A <- fitted_A + b %*% t(TFA$T_sample[m, ] - old_T)
  }
  
  return(TFA)
}
# TF_activity_T_sampling <- function(A, A_sample, ATAC_Cell_Sample_vector, R, R_sample, RNA_Cell_Sample_vector, TFA, T_prior_mean, T_prior_var, B, B_state, sigma_A_noise, P, G, M, S) {
#   TF_index <- sample(M)
#   for (i in 1:M) {
#     # during the looping, TFA is iteratively updated for all M TFs and the update order is random in each round of sampling
#     m <- TF_index[i]
#     if (sum(B_state[, m]) > 0) {
#     #remove division by S (sum of states)
#      #   temp_var <- sum(B[, m]^2) * T_prior_var[m, ] / sum(B_state[, m]) + sigma_A_noise
#      # mean_T <- (t(B[, m]) %*% (A_sample - B %*% TFA$T_sample + B[, m] %*% t(TFA$T_sample[m, ])) / sum(B_state[, m]) + T_prior_mean[m, ] * sigma_A_noise) / temp_var
#       temp_var <- sum(B[, m]^2) * T_prior_var[m, ] + sigma_A_noise
#       mean_T <- (t(B[, m]) %*% (A_sample - B %*% TFA$T_sample + B[, m] %*% t(TFA$T_sample[m, ])) * T_prior_var[m, ] + T_prior_mean[m, ] * sigma_A_noise) / temp_var
#       
#       variance_T <- T_prior_var[m, ] * sigma_A_noise / temp_var
# 
#       for (s in 1:S) {
#         aa <- rnorm(1)
#         if (aa - 3 > 0) {
#           aa <- 3
#         }
#         if (aa + 3 < 0) {
#           aa <- -3
#         }
# 
#         TFA$T_sample[m, s] <- aa * sqrt(abs(variance_T[s])) + mean_T[s]
# 
#         ATAC_cell_index <- which(ATAC_Cell_Sample_vector == s)
#         TFA$T_A[m, ATAC_cell_index] <- rnorm(length(ATAC_cell_index)) * sqrt(abs(variance_T[s])) + TFA$T_sample[m, s]
# 
#         RNA_cell_index <- which(RNA_Cell_Sample_vector == s)
#         TFA$T_R[m, RNA_cell_index] <- rnorm(length(RNA_cell_index)) * sqrt(abs(variance_T[s])) + TFA$T_sample[m, s]
#       }
#     }
#   }
#   return(TFA)
# }



# #*************************** TF-peak binding confidence estimation ***************************
# TF_peak_binding_B_sampling <- function(A, A_sample, ATAC_Cell_Sample_vector, TFA, B, B_state, B_prior_mean, B_prior_var, sigma_A_noise, P, G, M, S) {
#   T_sample <- matrix(0, nrow = M, ncol = S)
#   for (s in 1:S) {
#   totake <- which(ATAC_Cell_Sample_vector == s)
#   if (length(totake) == 1) {
#     T_sample[, s] <- TFA$T_A[, totake]
#   }else{
#     T_sample[, s] <- rowMeans(TFA$T_A[, totake])
#   }
#   }
# 
#   TF_index <- sample(M)
# 
#   for (i in 1:M) {
#     m <- TF_index[i]
#     temp_var <- sum(TFA$T_sample[m, ]^2) * B_prior_var[m] / S + sigma_A_noise
#     mean_B <- ((A_sample - B %*% TFA$T_sample + B[, m] %*% t(TFA$T_sample[m, ])) %*% TFA$T_sample[m, ] * B_prior_var[m] / S + B_prior_mean[, m] * sigma_A_noise) / temp_var
#     vairance_B <- B_prior_var[m] * sigma_A_noise / temp_var
# 
#     bb <- rnorm(P)
#     bb[which(bb - 3 > 0)] <- 3
#     bb[which(bb + 3 < 0)] <- -3
# 
#     B[, m] <- as.vector((bb * sqrt(abs(vairance_B)) + mean_B) * B_state[, m])
#   }
# 
#   return(B)
# }

#*************************** Peak-gene looping confidence estimation ***************************

Peak_gene_looping_L_samping <- function(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, sigma_R_noise, P, G, M, S) {
  T_sample <- matrix(0, nrow = M, ncol = S)
  for (s in 1:S) {
    totake <- which(RNA_Cell_Sample_vector == s)
    if (length(totake) == 1) {
      T_sample[, s] <- TFA$T_R[, totake]}
    else{
      T_sample[, s] <- rowMeans(TFA$T_R[, totake])
    }
  }

  A_estimate <- B %*% TFA$T_sample

  Peak_index <- sample(P)
  for (i in 1:P) {
    f <- Peak_index[i]
    ## remove division by S
   # temp_var <- sum(A_estimate[f, ]^2) * L_prior_var / S + sigma_R_noise
   # mean_L <- (A_estimate[f, ] %*% t(R_sample - t(L) %*% A_estimate + L[f, ] %*% t(A_estimate[f, ])) * L_prior_var / S + L_prior_mean[f, ] * sigma_R_noise) / temp_var
    temp_var <- sum(A_estimate[f, ]^2) * L_prior_var + sigma_R_noise
    mean_L <- (A_estimate[f, ] %*% t(R_sample - t(L) %*% A_estimate + L[f, ] %*% t(A_estimate[f, ])) * L_prior_var + L_prior_mean[f, ] * sigma_R_noise) / temp_var
    
    vairance_L <- L_prior_var * sigma_R_noise / temp_var

    ll <- rnorm(G)
    ll[which(ll - 3 > 0)] <- 3
    ll[which(ll + 3 < 0)] <- -3

    L[f, ] <- as.vector((ll * sqrt(vairance_L) + mean_L) * L_state[f, ])
  }
  return(L)
}

# #*************************** TF-peak binding state update ***************************
# TF_peak_binary_binding_B_state_sampling <- function(A, A_sample, ATAC_Cell_Sample_vector, TFA,
#                                                     B, B_state, B_prior_mean, B_prior_var, B_prior_prob, sigma_A_noise, P, G, M, S) {
#   T_sample <- matrix(0, nrow = M, ncol = S)
#   for (s in 1:S) {
#     totake <- which(ATAC_Cell_Sample_vector == s)
#     if (length(totake) == 1) {
#       T_sample[, s] <- TFA$T_A[, totake]
#     }else{
#       T_sample[, s] <- rowMeans(TFA$T_A[, which(ATAC_Cell_Sample_vector == s)])
#     }
#   }
# 
#   TF_index <- sample(M)
#   Peak_index <- sample(P)
# 
#   for (i in 1:P) {
#     f <- Peak_index[i]
#     temp <- t(B[f, ]) %*% TFA$T_sample
# 
#     for (j in 1:M) {
#       m <- TF_index[j]
#       old_B <- B[f, m]
#       temp_var <- sum(TFA$T_sample[m, ]^2) * B_prior_var[m] / S + sigma_A_noise
# 
#       if (B_state[f, m] > 0) {
#         mean_B <- (sum((A_sample[f, ] - temp + B[f, m] * TFA$T_sample[m, ]) * TFA$T_sample[m, ]) * B_prior_var[m] / S + B_prior_mean[f, m] * sigma_A_noise) / temp_var
# 
#         vairance_B <- B_prior_var[m] * sigma_A_noise / temp_var
# 
#         post_b1 <- exp(-(B[f, m] * 1 - mean_B)^2 / (2 * vairance_B)) * (B_prior_prob[f, m] + 0.25) + 1e-6
# 
#         post_b0 <- exp(-(B[f, m] * 0 - mean_B)^2 / (2 * vairance_B)) * (1 - B_prior_prob[f, m] + 0.25) + 1e-6
# 
#         P1 <- post_b1 / (post_b1 + post_b0)
# 
#         threshold_c <- runif(1)
# 
#         if (is.na(P1) || is.nan(P1)) {
#           P1 <- 0.5
#         }
# 
#         if (P1 < threshold_c) {
#           B[f, m] <- 0
#           B_state[f, m] <- 0
#         } else {
#           B[f, m] <- B[f, m]
#           B_state[f, m] <- 1
#         }
#       } #### THIS IS WHERE WE CHANGE FROM ORIGINAL CODE. PER COPILOT: 
#         ####      "If an active link is switched off by the first block, the second block immediately runs and can switch the link back on during the same iteration."
#        else  if (B_state[f, m] == 0 && B_prior_prob[f, m] > 0) {
#         mean_B <- (sum((A_sample[f, ] - temp) * TFA$T_sample[m, ]) * B_prior_var[m] / S + B_prior_mean[f, m] * sigma_A_noise) / temp_var
#         vairance_B <- B_prior_var[m] * sigma_A_noise / temp_var
# 
#         bb <- rnorm(1)
# 
#         if (bb - 3 > 0) {
#           bb <- 3
#         }
# 
#         if (bb + 3 < 0) {
#           bb <- -3
#         }
# 
#         B_temp <- bb * sqrt(abs(vairance_B)) + mean_B
# 
#         post_b1 <- exp(-(bb^2) / 2) * (B_prior_prob[f, m] + 0.25) + 1e-6
# 
#         post_b0 <- exp(-(mean_B)^2 / (2 * vairance_B)) * (1 - B_prior_prob[f, m] + 0.25) + 1e-6
# 
#         P1 <- post_b1 / (post_b1 + post_b0)
# 
#         threshold_c <- runif(1)
# 
#         if (is.na(P1) || is.nan(P1)) {
#           P1 <- 0.5
#         }
# 
#         if (P1 < threshold_c) {
#           B[f, m] <- 0
#           B_state[f, m] <- 0
#         } else {
#           B[f, m] <- B_temp
#           B_state[f, m] <- 1
#         }
#        }
#       temp <- temp + (B[f, m] - old_B) * TFA$T_sample[m, ]
#     }
#   }
#   Updated_B_and_State <- list("B" = B, "B_state" = B_state)
#   return(Updated_B_and_State)
# }

#*************************** Peak-gene looping state update ***************************
Peak_gene_binary_looping_L_state_samping <- function(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, L_prior_prob, sigma_R_noise, P, G, M, S) {
  Peak_index <- sample(P)
  Gene_index <- sample(G)

  T_sample <- matrix(0, nrow = M, ncol = S)
  for (s in 1:S) {
    totake <- which(RNA_Cell_Sample_vector == s)
    if (length(totake) == 1) {
      T_sample[, s] <- TFA$T_R[, totake]
    } else {
      T_sample[, s] <- rowMeans(TFA$T_R[, totake])
    }
  }


  A_estimate <- B %*% TFA$T_sample


  for (i in 1:G) {
    g <- Gene_index[i]
    temp <- L[, g] %*% A_estimate

    for (j in 1:P) {
      f <- Peak_index[j]
      old_L <- L[f, g]
      #remove division by S
  #      temp_var <- sum(A_estimate[f, ]^2) * L_prior_var / S + sigma_R_noise
      temp_var <- sum(A_estimate[f, ]^2) * L_prior_var  + sigma_R_noise
      
      if (L_state[f, g] > 0) {
        #remove division by S
       # mean_L <- (sum((R_sample[g, ] - temp + L[f, g] * A_estimate[f, ]) * A_estimate[f, ]) * L_prior_var / S + L_prior_mean[f, g] * sigma_R_noise) / temp_var
        mean_L <- (sum((R_sample[g, ] - temp + L[f, g] * A_estimate[f, ]) * A_estimate[f, ]) * L_prior_var  + L_prior_mean[f, g] * sigma_R_noise) / temp_var
        
        vairance_L <- L_prior_var * sigma_R_noise / temp_var
   
        post_l1 <- exp(-(L[f, g] * 1 - mean_L)^2 / (2 * vairance_L)) * (L_prior_prob[f, g] + 0.25) + 1e-6

        post_l0 <- exp(-(L[f, g] * 0 - mean_L)^2 / (2 * vairance_L)) * (1 - L_prior_prob[f, g] + 0.25) + 1e-6

        P1 <- post_l1 / (post_l1 + post_l0)

        threshold_c <- runif(1)

        if (is.na(P1) || is.nan(P1)) {
          P1 <- 0.5
        }

        if (P1 < threshold_c) {
          L[f, g] <- 0
          L_state[f, g] <- 0
        } else {
          L[f, g] <- L[f, g]
          L_state[f, g] <- 1
        }
      }
      #### THIS IS WHERE WE CHANGE FROM ORIGINAL CODE. PER COPILOT: 
      ####      "If an active link is switched off by the first block, the second block immediately runs and can switch the link back on during the same iteration."
      else if (L_state[f, g] == 0 && L_prior_prob[f, g] > 0) {
      #remove division by S
        #  mean_L <- (sum((R_sample[g, ] - temp) * A_estimate[f, ]) * L_prior_var / S + L_prior_mean[f, g] * sigma_R_noise) / temp_var
        mean_L <- (sum((R_sample[g, ] - temp) * A_estimate[f, ]) * L_prior_var  + L_prior_mean[f, g] * sigma_R_noise) / temp_var
        
        vairance_L <- L_prior_var * sigma_R_noise / temp_var

        ll <- rnorm(1)

        if (ll - 3 > 0) {
          ll <- 3
        }
        if (ll + 3 < 0) {
          ll <- -3
        }

        L_temp <- ll * sqrt(abs(vairance_L)) + mean_L

        post_l1 <- exp(-ll^2 / 2) * (L_prior_prob[f, g] + 0.1) + 1e-6

        post_l0 <- exp(-(mean_L)^2 / (2 * vairance_L)) * (1 - L_prior_prob[f, g] + 0.1) + 1e-6

        P1 <- post_l1 / (post_l1 + post_l0)

        threshold_c <- runif(1)

        if (is.na(P1) || is.nan(P1)) {
          P1 <- 0.5
        }

        if (P1 < threshold_c) {
          L[f, g] <- 0
          L_state[f, g] <- 0
        } else {
          L[f, g] <- L_temp
          L_state[f, g] <- 1
        }
      }
      temp <- temp + (L[f, g] - old_L) * A_estimate[f, ]
    }
  }
  Updated_L_and_State <- list("L" = L, "L_state" = L_state)
  return(Updated_L_and_State)
}




choose_diagnostic_intervals <- function(
    iteration_num,
    burnin,
    target_checkpoints = 30L,
    target_windows = 6L,
    minimum_check_every = 10L,
    minimum_window_size = 50L
) {
  n_postburnin <- iteration_num - burnin
  
  if (n_postburnin <= 0L) {
    stop("`burnin` must be smaller than `iteration_num`.")
  }
  
  check_every <- max(
    minimum_check_every,
    floor(n_postburnin / target_checkpoints)
  )
  
  window_size <- max(
    minimum_window_size,
    floor(n_postburnin / target_windows)
  )
  
  list(
    n_postburnin = n_postburnin,
    check_every = as.integer(check_every),
    window_size = as.integer(window_size),
    expected_checkpoints =
      floor(n_postburnin / check_every),
    expected_complete_windows =
      floor(n_postburnin / window_size)
  )
}

sample_residual_variance <- function(
    residuals,
    alpha = 1,
    beta = 1
) {
  valid <- is.finite(residuals)
  N <- sum(valid)
  
  if (N == 0L) {
    stop("No finite residuals available.")
  }
  
  SSE <- sum(residuals[valid]^2)
  
  1 / rgamma(
    1,
    shape = alpha + N / 2,
    rate = beta + SSE / 2
  )
}

#' MCMC sampling
#'
#' This function estimates the parameters.
#'
#' @param loaded_data The output from `Data_loading`
#' @param Candidate_circuits The output from "Candidate_circuits" (with or without TAD)
#' @param Initial_model The output from `MAGICAL_initialization`
#' @param iteration_num The number of iterations
#'
#' @return A list named "Circuits_linkage_posterior" containing the following elements:
#' * TF_Peak_Binding_prob: The posterior probability of TF binding to peaks
#' * Peak_Gene_Looping_prob: The posterior probability of looping between peaks and genes
#' * Noise_parameters: The posterior noise parameters
#'
#' @import Matrix
#' @importFrom Matrix colSums rowSums sparseMatrix as.matrix
#'
#' @export

q99 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(quantile(x, probs = 0.99, names = FALSE))
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

choose_diagnostic_intervals <- function(iteration_num, burnin, target_checkpoints = 30L, target_windows = 6L, minimum_check_every = 10L, minimum_window_size = 50L) {
  n_postburnin <- iteration_num - burnin
  if (n_postburnin <= 0L) stop("`burnin` must be smaller than `iteration_num`.")
  check_every <- max(minimum_check_every, floor(n_postburnin / target_checkpoints))
  window_size <- max(minimum_window_size, floor(n_postburnin / target_windows))
  list(n_postburnin = n_postburnin, check_every = as.integer(check_every), window_size = as.integer(window_size), expected_checkpoints = floor(n_postburnin / check_every), expected_complete_windows = floor(n_postburnin / window_size))
}

initialize_mcmc_diagnostics <- function(B_state, L_state, B_candidate_index, L_candidate_index, window_size) {
  d <- new.env(parent = emptyenv())
  
  d$n_kept <- 0L
  d$stable_check_count <- 0L
  
  d$B_candidate_index <- B_candidate_index
  d$L_candidate_index <- L_candidate_index
  d$window_size <- window_size
  
  d$previous_B_state <- B_state
  d$previous_L_state <- L_state
  d$previous_B_prob <- NULL
  d$previous_L_prob <- NULL
  d$previous_B_window_prob <- NULL
  d$previous_L_window_prob <- NULL
  
  d$B_sum <- matrix(0, nrow = nrow(B_state), ncol = ncol(B_state), dimnames = dimnames(B_state))
  d$L_sum <- matrix(0, nrow = nrow(L_state), ncol = ncol(L_state), dimnames = dimnames(L_state))
  d$B_switches <- matrix(0L, nrow = nrow(B_state), ncol = ncol(B_state), dimnames = dimnames(B_state))
  d$L_switches <- matrix(0L, nrow = nrow(L_state), ncol = ncol(L_state), dimnames = dimnames(L_state))
  
  d$window_n <- 0L
  d$window_start_iteration <- NA_integer_
  d$B_window_sum <- matrix(0, nrow = nrow(B_state), ncol = ncol(B_state), dimnames = dimnames(B_state))
  d$L_window_sum <- matrix(0, nrow = nrow(L_state), ncol = ncol(L_state), dimnames = dimnames(L_state))
  d$B_window_switches <- matrix(0L, nrow = nrow(B_state), ncol = ncol(B_state), dimnames = dimnames(B_state))
  d$L_window_switches <- matrix(0L, nrow = nrow(L_state), ncol = ncol(L_state), dimnames = dimnames(L_state))
  
  d$convergence_log <- data.frame(iteration = integer(), n_kept = integer(), max_B_probability_change = numeric(), q99_B_probability_change = numeric(), max_L_probability_change = numeric(), q99_L_probability_change = numeric(), mean_B_switch_rate = numeric(), mean_L_switch_rate = numeric(), sigma_A_recent_mean = numeric(), sigma_R_recent_mean = numeric(), stable = logical())
  
  d$window_log <- data.frame(window_start = integer(), window_end = integer(), window_n = integer(), q99_B_change_from_previous_window = numeric(), max_B_change_from_previous_window = numeric(), q99_L_change_from_previous_window = numeric(), max_L_change_from_previous_window = numeric(), mean_B_probability = numeric(), mean_L_probability = numeric(), mean_B_switch_rate = numeric(), mean_L_switch_rate = numeric(), median_ATAC_noise = numeric(), q25_ATAC_noise = numeric(), q75_ATAC_noise = numeric(), q95_ATAC_noise = numeric(), q99_ATAC_noise = numeric(), max_ATAC_noise = numeric(), median_RNA_noise = numeric(), q25_RNA_noise = numeric(), q75_RNA_noise = numeric(), q95_RNA_noise = numeric(), q99_RNA_noise = numeric(), max_RNA_noise = numeric())
  
  d
}

save_diagnostic_window <- function(d, iteration, Noise_parameters) {
  B_window_prob <- d$B_window_sum / d$window_n
  L_window_prob <- d$L_window_sum / d$window_n
  B_window_switch_rate <- d$B_window_switches / d$window_n
  L_window_switch_rate <- d$L_window_switches / d$window_n
  
  if (is.null(d$previous_B_window_prob)) {
    q99_B_change <- NA_real_
    max_B_change <- NA_real_
    q99_L_change <- NA_real_
    max_L_change <- NA_real_
  } else {
    B_change <- abs(B_window_prob[d$B_candidate_index] - d$previous_B_window_prob[d$B_candidate_index])
    L_change <- abs(L_window_prob[d$L_candidate_index] - d$previous_L_window_prob[d$L_candidate_index])
    q99_B_change <- q99(B_change)
    max_B_change <- safe_max(B_change)
    q99_L_change <- q99(L_change)
    max_L_change <- safe_max(L_change)
  }
  
  window_noise <- Noise_parameters[d$window_start_iteration:iteration, , drop = FALSE]
  atac_noise <- window_noise[, "sigma_A_noise"]
  rna_noise <- window_noise[, "sigma_R_noise"]
  
  d$window_log <- rbind(d$window_log, data.frame(
    window_start = d$window_start_iteration,
    window_end = iteration,
    window_n = d$window_n,
    q99_B_change_from_previous_window = q99_B_change,
    max_B_change_from_previous_window = max_B_change,
    q99_L_change_from_previous_window = q99_L_change,
    max_L_change_from_previous_window = max_L_change,
    mean_B_probability = mean(B_window_prob[d$B_candidate_index], na.rm = TRUE),
    mean_L_probability = mean(L_window_prob[d$L_candidate_index], na.rm = TRUE),
    mean_B_switch_rate = mean(B_window_switch_rate[d$B_candidate_index], na.rm = TRUE),
    mean_L_switch_rate = mean(L_window_switch_rate[d$L_candidate_index], na.rm = TRUE),
    median_ATAC_noise = median(atac_noise, na.rm = TRUE),
    q25_ATAC_noise = as.numeric(quantile(atac_noise, 0.25, na.rm = TRUE, names = FALSE)),
    q75_ATAC_noise = as.numeric(quantile(atac_noise, 0.75, na.rm = TRUE, names = FALSE)),
    q95_ATAC_noise = as.numeric(quantile(atac_noise, 0.95, na.rm = TRUE, names = FALSE)),
    q99_ATAC_noise = as.numeric(quantile(atac_noise, 0.99, na.rm = TRUE, names = FALSE)),
    max_ATAC_noise = max(atac_noise, na.rm = TRUE),
    median_RNA_noise = median(rna_noise, na.rm = TRUE),
    q25_RNA_noise = as.numeric(quantile(rna_noise, 0.25, na.rm = TRUE, names = FALSE)),
    q75_RNA_noise = as.numeric(quantile(rna_noise, 0.75, na.rm = TRUE, names = FALSE)),
    q95_RNA_noise = as.numeric(quantile(rna_noise, 0.95, na.rm = TRUE, names = FALSE)),
    q99_RNA_noise = as.numeric(quantile(rna_noise, 0.99, na.rm = TRUE, names = FALSE)),
    max_RNA_noise = max(rna_noise, na.rm = TRUE)
  ))
  
  d$previous_B_window_prob <- B_window_prob
  d$previous_L_window_prob <- L_window_prob
  
  d$B_window_sum[,] <- 0
  d$L_window_sum[,] <- 0
  d$B_window_switches[,] <- 0L
  d$L_window_switches[,] <- 0L
  d$window_n <- 0L
  d$window_start_iteration <- NA_integer_
}

update_mcmc_diagnostics <- function(d, iteration, burnin, B_state, L_state, Noise_parameters, check_every, stability_tolerance, required_stable_checks) {
  B_switch <- abs(B_state - d$previous_B_state)
  L_switch <- abs(L_state - d$previous_L_state)
  
  d$previous_B_state <- B_state
  d$previous_L_state <- L_state
  
  if (iteration <= burnin) return(invisible(FALSE))
  
  d$n_kept <- d$n_kept + 1L
  
  if (d$window_n == 0L) d$window_start_iteration <- iteration
  
  d$window_n <- d$window_n + 1L
  d$B_sum <- d$B_sum + B_state
  d$L_sum <- d$L_sum + L_state
  d$B_switches <- d$B_switches + B_switch
  d$L_switches <- d$L_switches + L_switch
  d$B_window_sum <- d$B_window_sum + B_state
  d$L_window_sum <- d$L_window_sum + L_state
  d$B_window_switches <- d$B_window_switches + B_switch
  d$L_window_switches <- d$L_window_switches + L_switch
  
  if (d$window_n == d$window_size) save_diagnostic_window(d, iteration, Noise_parameters)
  
  if (d$n_kept < check_every || d$n_kept %% check_every != 0L) return(invisible(FALSE))
  
  current_B_prob <- d$B_sum / d$n_kept
  current_L_prob <- d$L_sum / d$n_kept
  
  if (is.null(d$previous_B_prob)) {
    message("Iteration ", iteration, ": saving first post-burn-in checkpoint as baseline.")
    d$previous_B_prob <- current_B_prob
    d$previous_L_prob <- current_L_prob
    return(invisible(FALSE))
  }
  
  B_change <- abs(current_B_prob[d$B_candidate_index] - d$previous_B_prob[d$B_candidate_index])
  L_change <- abs(current_L_prob[d$L_candidate_index] - d$previous_L_prob[d$L_candidate_index])
  
  q99_B_change <- q99(B_change)
  q99_L_change <- q99(L_change)
  max_B_change <- safe_max(B_change)
  max_L_change <- safe_max(L_change)
  
  B_switch_rate <- d$B_switches / d$n_kept
  L_switch_rate <- d$L_switches / d$n_kept
  mean_B_switch_rate <- mean(B_switch_rate[d$B_candidate_index], na.rm = TRUE)
  mean_L_switch_rate <- mean(L_switch_rate[d$L_candidate_index], na.rm = TRUE)
  
  recent_start <- max(burnin + 1L, iteration - check_every + 1L)
  sigma_A_recent_mean <- mean(Noise_parameters[recent_start:iteration, "sigma_A_noise"], na.rm = TRUE)
  sigma_R_recent_mean <- mean(Noise_parameters[recent_start:iteration, "sigma_R_noise"], na.rm = TRUE)
  
  stable_now <- is.finite(q99_B_change) && is.finite(q99_L_change) && q99_B_change <= stability_tolerance && q99_L_change <= stability_tolerance
  d$stable_check_count <- if (stable_now) d$stable_check_count + 1L else 0L
  
  d$convergence_log <- rbind(d$convergence_log, data.frame(iteration = iteration, n_kept = d$n_kept, max_B_probability_change = max_B_change, q99_B_probability_change = q99_B_change, max_L_probability_change = max_L_change, q99_L_probability_change = q99_L_change, mean_B_switch_rate = mean_B_switch_rate, mean_L_switch_rate = mean_L_switch_rate, sigma_A_recent_mean = sigma_A_recent_mean, sigma_R_recent_mean = sigma_R_recent_mean, stable = stable_now))
  
  message(sprintf("Iteration %d: B q99 change = %.5f; L q99 change = %.5f; stable checks = %d/%d", iteration, q99_B_change, q99_L_change, d$stable_check_count, required_stable_checks))
  
  d$previous_B_prob <- current_B_prob
  d$previous_L_prob <- current_L_prob
  
  invisible(d$stable_check_count >= required_stable_checks)
}

finalize_mcmc_diagnostics <- function(d, final_iteration, Noise_parameters) {
  if (d$n_kept == 0L) stop("No post-burn-in iterations were retained.")
  if (d$window_n > 0L) save_diagnostic_window(d, final_iteration, Noise_parameters)
  
  list(
    B_prob = d$B_sum / d$n_kept,
    L_prob = d$L_sum / d$n_kept,
    B_switch_rate = d$B_switches / d$n_kept,
    L_switch_rate = d$L_switches / d$n_kept,
    convergence_log = d$convergence_log,
    window_log = d$window_log,
    n_kept = d$n_kept,
    stable_check_count = d$stable_check_count
  )
}

MAGICAL_estimation <- function(
    loaded_data,
    Candidate_circuits,
    Initial_model,
    iteration_num,
    burnin = round(iteration_num * 0.25),
    check_every = NULL,
    window_size = NULL,
    target_checkpoints = 30L,
    target_windows = 6L,
    stability_tolerance = 0.005,
    required_stable_checks = 4L,
    stop_when_stable = FALSE
) {
  iteration_num <- as.integer(iteration_num)
  burnin <- as.integer(burnin)
  required_stable_checks <- as.integer(required_stable_checks)

  if (iteration_num < 1L) stop("`iteration_num` must be at least 1.")
  if (burnin < 0L || burnin >= iteration_num) stop("`burnin` must be non-negative and smaller than `iteration_num`.")
  if (required_stable_checks < 1L) stop("`required_stable_checks` must be at least 1.")
  
  diagnostic_intervals <- choose_diagnostic_intervals(iteration_num, burnin, target_checkpoints, target_windows)
  n_postburnin_planned <- diagnostic_intervals$n_postburnin
  
  if (is.null(check_every)) check_every <- diagnostic_intervals$check_every
  if (is.null(window_size)) window_size <- diagnostic_intervals$window_size
  
  check_every <- as.integer(check_every)
  window_size <- as.integer(window_size)
  
  if (length(check_every) != 1L || !is.finite(check_every) || check_every < 1L) stop("`check_every` must be one positive integer.")
  if (length(window_size) != 1L || !is.finite(window_size) || window_size < 2L) stop("`window_size` must be one integer of at least 2.")
  if (check_every > n_postburnin_planned) stop("`check_every` is larger than the number of planned post-burn-in iterations.")
  if (window_size > n_postburnin_planned) stop("`window_size` is larger than the number of planned post-burn-in iterations.")
  
  planned_checkpoints <- floor(n_postburnin_planned / check_every)
  planned_complete_windows <- floor(n_postburnin_planned / window_size)
  
  message(sprintf(
    "MCMC schedule: iterations=%d; burn-in=%d; post-burn-in=%d; check every=%d; window size=%d",
    iteration_num, burnin, n_postburnin_planned, check_every, window_size
  ))
  
  # Input dimensions and sample assignments
  Common_samples <- loaded_data$Common_samples
  S <- length(Common_samples)
  
  Peak_index <- match(Candidate_circuits$Peaks$Peak_index, loaded_data$scATAC_Peaks$Peak_index)
  if (anyNA(Peak_index)) stop("Some candidate peaks were not found in `loaded_data$scATAC_Peaks`.")
  
  A <- loaded_data$scATAC_read_count_matrix[Peak_index, ]
  A_sample <- Candidate_circuits$Peak_log2Count
  P <- length(Peak_index)
  
  ATAC_Cell_Sample_vector <- matrix(0L, nrow = 1L, ncol = nrow(loaded_data$scATAC_cells))
  for (s in seq_len(S)) ATAC_Cell_Sample_vector[which(loaded_data$scATAC_cells$subject_ID == Common_samples[s])] <- s
  
  Gene_index <- match(Candidate_circuits$Genes$Gene_symbols, loaded_data$scRNA_Genes$Gene_symbols)
  if (anyNA(Gene_index)) stop("Some candidate genes were not found in `loaded_data$scRNA_Genes`.")
  
  R <- loaded_data$scRNA_read_count_matrix[Gene_index, ]
  R_sample <- Candidate_circuits$Gene_log2Count
  G <- length(Gene_index)
  
  RNA_Cell_Sample_vector <- matrix(0L, nrow = 1L, ncol = nrow(loaded_data$scRNA_cells))
  for (s in seq_len(S)) RNA_Cell_Sample_vector[1L, which(loaded_data$scRNA_cells$subject_ID == Common_samples[s])] <- s
  
  # TF-peak parameters

  B_prior_mean <- Initial_model$B_mean
  B_prior_var <- Initial_model$B_var
  B_prior_prob <- Initial_model$B_prob
  B <- Initial_model$B_initial
  
  B_mask <- Initial_model$B_mask
  if (is.null(B_mask)) {
    stop("Initial_model$B_mask is missing.")
  }
  B_mask <- 1L * (B_mask > 0)
  #B_state <- as.matrix(Candidate_circuits$TF_Peak_Binding)

  
  B_state <- matrix(
    0L,
    nrow = nrow(B_mask),
    ncol = ncol(B_mask),
    dimnames = dimnames(B_mask)
  )
  
  allowed_B_edges <- B_mask > 0
  
  B_state[allowed_B_edges] <- rbinom(
    n = sum(allowed_B_edges),
    size = 1L,
    prob = B_prior_prob[allowed_B_edges]
  )
  
  B <- Initial_model$B_initial * B_state
  
  # Peak-gene parameters
  L_prior_mean <- Initial_model$L_mean
  L_prior_var <- Initial_model$L_var
  L_prior_prob <- Initial_model$L_prob
  L <- Initial_model$L_prior
  L_state <- as.matrix(Candidate_circuits$Peak_Gene_looping)
  
  # TF activity parameters
  T_prior_mean <- Initial_model$T_mean
  T_prior_var <- Initial_model$T_var
  T_sample <- Initial_model$T_prior

  M <- nrow(T_sample)
  
  # Basic dimension checks
  if (!identical(dim(A_sample), c(P, S))) stop(sprintf("`A_sample` must have dimensions %d x %d.", P, S))
  if (!identical(dim(R_sample), c(G, S))) stop(sprintf("`R_sample` must have dimensions %d x %d.", G, S))
  if (!identical(dim(B), c(P, M))) stop(sprintf("Initial `B` must have dimensions %d x %d.", P, M))
  if (!identical(dim(B_state), c(P, M))) stop(sprintf("Initial `B_state` must have dimensions %d x %d.", P, M))
  if (!identical(dim(L), c(P, G))) stop(sprintf("Initial `L` must have dimensions %d x %d.", P, G))
  if (!identical(dim(L_state), c(P, G))) stop(sprintf("Initial `L_state` must have dimensions %d x %d.", P, G))
  if (!identical(dim(T_sample), c(M, S))) stop(sprintf("Initial `T_sample` must have dimensions %d x %d.", M, S))
  
  # Initialize cell-level TF activity
  T_A <- matrix(0, nrow = M, ncol = nrow(loaded_data$scATAC_cells))
  T_R <- matrix(0, nrow = M, ncol = nrow(loaded_data$scRNA_cells))
  
  for (s in seq_len(S)) {
    ATAC_index <- which(ATAC_Cell_Sample_vector == s)
    RNA_index <- which(RNA_Cell_Sample_vector == s)
    
    for (m in seq_len(M)) {
      T_A[m, ATAC_index] <- rnorm(length(ATAC_index), mean = T_prior_mean[m, s], sd = sqrt(abs(T_prior_var[m, s])))
      T_R[m, RNA_index] <- rnorm(length(RNA_index), mean = T_prior_mean[m, s], sd = sqrt(abs(T_prior_var[m, s])))
    }
  }
  
  TFA <- list(T_A = T_A, T_R = T_R, T_sample = T_sample)
  
  # Residual-variance priors
  alpha_A <- 1
  beta_A <- 1
  alpha_R <- 1
  beta_R <- 1
  
  # Initial predictions, residuals, and residual-variance draws
  pred_A <- B %*% TFA$T_sample
  pred_R <- t(L) %*% pred_A
  
  resid_A <- A_sample - pred_A
  resid_R <- R_sample - pred_R
  
  N_A <- sum(is.finite(resid_A))
  N_R <- sum(is.finite(resid_R))
  
  if (N_A == 0L) stop("No finite ATAC residuals are available.")
  if (N_R == 0L) stop("No finite RNA residuals are available.")
  
  SSE_A <- sum(resid_A^2, na.rm = TRUE)
  SSE_R <- sum(resid_R^2, na.rm = TRUE)
  
  sigma_A_noise <- 1 / rgamma(1, shape = alpha_A + N_A / 2, rate = beta_A + SSE_A / 2)
  sigma_R_noise <- 1 / rgamma(1, shape = alpha_R + N_R / 2, rate = beta_R + SSE_R / 2)
  
  # Fixed null-model quantities
  null_MSE_A <- mean(A_sample^2, na.rm = TRUE)
  null_MSE_R <- mean(R_sample^2, na.rm = TRUE)
  observed_sd_A <- sd(as.vector(A_sample), na.rm = TRUE)
  observed_sd_R <- sd(as.vector(R_sample), na.rm = TRUE)
  
  if (!is.finite(null_MSE_A) || null_MSE_A <= 0) stop("ATAC null MSE is not finite and positive.")
  if (!is.finite(null_MSE_R) || null_MSE_R <= 0) stop("RNA null MSE is not finite and positive.")
  if (!is.finite(observed_sd_A) || observed_sd_A <= 0) stop("ATAC observations have no finite positive standard deviation.")
  if (!is.finite(observed_sd_R) || observed_sd_R <= 0) stop("RNA observations have no finite positive standard deviation.")
  
  cat(sprintf("Null MSE: ATAC = %.5g; RNA = %.5g\n", null_MSE_A, null_MSE_R))
  cat(sprintf("Observed SD: ATAC = %.5g; RNA = %.5g\n", observed_sd_A, observed_sd_R))
  
  # Store variance, fit, scale, and correlation diagnostics
  Noise_parameters <- matrix(
    NA_real_,
    nrow = iteration_num,
    ncol = 10L,
    dimnames = list(
      NULL,
      c(
        "sigma_A_noise",
        "sigma_R_noise",
        "MSE_A",
        "MSE_R",
        "model_null_ratio_A",
        "model_null_ratio_R",
        "predicted_sd_A",
        "predicted_sd_R",
        "prediction_cor_A",
        "prediction_cor_R"
      )
    )
  )
  
  B_candidate_index <- which(as.matrix(Candidate_circuits$TF_Peak_Binding) > 0)
  L_candidate_index <- which(as.matrix(Candidate_circuits$Peak_Gene_looping) > 0)
  
  if (length(B_candidate_index) == 0L) stop("No candidate TF-peak links are present.")
  if (length(L_candidate_index) == 0L) stop("No candidate peak-gene links are present.")
  
  diagnostics <- initialize_mcmc_diagnostics(B_state, L_state, B_candidate_index, L_candidate_index, window_size)
  iteration_seg <- max(1L, floor(iteration_num / 10L))

  cat("Common_samples:\n")
  print(Common_samples)
  
  cat("A_sample columns:\n")
  print(colnames(A_sample))
  
  cat("R_sample columns:\n")
  print(colnames(R_sample))
  
  cat("T_sample columns:\n")
  print(colnames(TFA$T_sample))

  if (!is.null(colnames(A_sample))) A_sample <- A_sample[, Common_samples, drop = FALSE]
  if (!is.null(colnames(R_sample))) R_sample <- R_sample[, Common_samples, drop = FALSE]
  if (!is.null(colnames(TFA$T_sample))) TFA$T_sample <- TFA$T_sample[, Common_samples, drop = FALSE]

  if (!is.null(colnames(A_sample)) && !identical(colnames(A_sample), Common_samples)) stop("A_sample columns are not aligned to Common_samples.")
  if (!is.null(colnames(R_sample)) && !identical(colnames(R_sample), Common_samples)) stop("R_sample columns are not aligned to Common_samples.")
  if (!is.null(colnames(TFA$T_sample)) && !identical(colnames(TFA$T_sample), Common_samples)) stop("TFA$T_sample columns are not aligned to Common_samples.")
  
  
  #****** print diagnostics *****#
  #*B_candidate_index <- which(B_mask > 0L)
  L_candidate_index <- which(as.matrix(Candidate_circuits$Peak_Gene_looping) > 0)
  
  if (length(B_candidate_index) == 0L) {
    stop("No candidate TF-peak links are present.")
  }
  
  if (length(L_candidate_index) == 0L) {
    stop("No candidate peak-gene links are present.")
  }
  
  cat("Candidate B edges:", length(B_candidate_index), "\n")
  cat("Candidate L edges:", length(L_candidate_index), "\n")
  
  diagnostics <- initialize_mcmc_diagnostics(
    B_state,
    L_state,
    B_candidate_index,
    L_candidate_index,
    window_size
  )
  
  
  
  
  print("MAGICAL integration starts ...")
  
  for (i in seq_len(iteration_num)) {
    # Step 1: sample TF activity
    TFA <- TF_activity_T_sampling(A, A_sample, ATAC_Cell_Sample_vector, R, R_sample, RNA_Cell_Sample_vector, TFA, T_prior_mean, T_prior_var, B, B_state, sigma_A_noise, P, G, M, S)
    
    ## Step 2: sample TF-peak weights
    #B <- TF_peak_binding_B_sampling(A, A_sample, ATAC_Cell_Sample_vector, TFA, B, B_state, B_prior_mean, B_prior_var, sigma_A_noise, P, G, M, S)
    
  #  # Step 3: sample TF-peak states
  #  Updated_B_and_State <- TF_peak_binary_binding_B_state_sampling(A, A_sample, ATAC_Cell_Sample_vector, TFA, B, B_state, B_prior_mean, B_prior_var, B_prior_prob, sigma_A_noise, P, G, M, S)
  #  B <- Updated_B_and_State$B
  #  B_state <- Updated_B_and_State$B_state
 
    #New Step 2: jointly sample TF-peak states and weights
    
    Updated_B_and_State <-
      TF_peak_binding_B_and_state_sampling(
        A_sample = A_sample,
        TFA = TFA,
        B = B,
        B_state = B_state,
        B_mask = B_mask,
        B_prior_mean = B_prior_mean,
        B_prior_var = B_prior_var,
        B_prior_prob = B_prior_prob,
        sigma_A_noise = sigma_A_noise,
        P = P,
        M = M
      )
    
    B <- Updated_B_and_State$B
    B_state <- Updated_B_and_State$B_state
    
    
    # Step 4: sample peak-gene weights
    L <- Peak_gene_looping_L_samping(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, sigma_R_noise, P, G, M, S)
    
    # Step 5: sample peak-gene states
    Updated_L_and_State <- Peak_gene_binary_looping_L_state_samping(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, L_prior_prob, sigma_R_noise, P, G, M, S)
    L <- Updated_L_and_State$L
    L_state <- Updated_L_and_State$L_state
    
    # Step 6: calculate current predictions and residuals
    pred_A <- B %*% TFA$T_sample
    pred_R <- t(L) %*% pred_A
    
    resid_A <- A_sample - pred_A
    resid_R <- R_sample - pred_R
    
    N_A <- sum(is.finite(resid_A))
    N_R <- sum(is.finite(resid_R))
    
    if (N_A == 0L) stop(sprintf("No finite ATAC residuals at iteration %d.", i))
    if (N_R == 0L) stop(sprintf("No finite RNA residuals at iteration %d.", i))
    
    SSE_A <- sum(resid_A^2, na.rm = TRUE)
    SSE_R <- sum(resid_R^2, na.rm = TRUE)
    
    MSE_A <- SSE_A / N_A
    MSE_R <- SSE_R / N_R
    
    # Step 7: sample residual variances
    sigma_A_noise <- 1 / rgamma(1, shape = alpha_A + N_A / 2, rate = beta_A + SSE_A / 2)
    sigma_R_noise <- 1 / rgamma(1, shape = alpha_R + N_R / 2, rate = beta_R + SSE_R / 2)
    
    # Step 8: calculate current fit diagnostics
    model_null_ratio_A <- MSE_A / null_MSE_A
    model_null_ratio_R <- MSE_R / null_MSE_R
    predicted_sd_A <- sd(as.vector(pred_A), na.rm = TRUE)
    predicted_sd_R <- sd(as.vector(pred_R), na.rm = TRUE)
    
    prediction_cor_A <- if (is.finite(predicted_sd_A) && predicted_sd_A > 0) cor(as.vector(A_sample), as.vector(pred_A), use = "complete.obs") else NA_real_
    prediction_cor_R <- if (is.finite(predicted_sd_R) && predicted_sd_R > 0) cor(as.vector(R_sample), as.vector(pred_R), use = "complete.obs") else NA_real_
    
    peak_cor <- vapply(
      seq_len(nrow(A_sample)),
      function(p) {
        if (sd(A_sample[p, ], na.rm = TRUE) == 0 ||
            sd(pred_A[p, ], na.rm = TRUE) == 0) {
          return(NA_real_)
        }
        
        cor(
          A_sample[p, ],
          pred_A[p, ],
          use = "complete.obs"
        )
      },
      numeric(1)
    )
    
    median_peak_cor <- median(peak_cor, na.rm = TRUE)
    q25_peak_cor <- as.numeric(quantile(peak_cor, 0.25, na.rm = TRUE, names = FALSE))
    q75_peak_cor <- as.numeric(quantile(peak_cor, 0.75, na.rm = TRUE, names = FALSE))
    
    Noise_parameters[i, ] <- c(
      sigma_A_noise,
      sigma_R_noise,
      MSE_A,
      MSE_R,
      model_null_ratio_A,
      model_null_ratio_R,
      predicted_sd_A,
      predicted_sd_R,
      prediction_cor_A,
      prediction_cor_R
    )
    
    # Step 9: update post-burn-in convergence diagnostics
    reached_stability <- update_mcmc_diagnostics(diagnostics, i, burnin, B_state, L_state, Noise_parameters, check_every, stability_tolerance, required_stable_checks)
    
    # Print progress ten times per run
    if (i %% iteration_seg == 0L || i == iteration_num) {
      cat(sprintf("MAGICAL finished %.0f percent\n", 100 * i / iteration_num))
      cat(sprintf("Iteration %d\n", i))
      cat(sprintf("  ATAC: variance = %.5g, MSE = %.5g, variance/MSE = %.3f, model/null = %.3f\n", sigma_A_noise, MSE_A, sigma_A_noise / MSE_A, model_null_ratio_A))
      cat(sprintf("  RNA:  variance = %.5g, MSE = %.5g, variance/MSE = %.3f, model/null = %.3f\n", sigma_R_noise, MSE_R, sigma_R_noise / MSE_R, model_null_ratio_R))
      cat(sprintf("  ATAC SD: observed = %.4f, predicted = %.4f\n", observed_sd_A, predicted_sd_A))
      cat(sprintf("  RNA SD:  observed = %.4f, predicted = %.4f\n", observed_sd_R, predicted_sd_R))
      cat(sprintf("  ATAC observed/predicted correlation = %.3f\n", prediction_cor_A))
      cat(sprintf("  RNA observed/predicted correlation  = %.3f\n", prediction_cor_R))
      
      cat(sprintf(
        "  Peak-wise ATAC correlation: median = %.3f, q25 = %.3f, q75 = %.3f\n",
        median_peak_cor,
        q25_peak_cor,
        q75_peak_cor
      ))
    }
    
    if (stop_when_stable && reached_stability) {
      message(sprintf("Stopping at iteration %d after %d consecutive stable checkpoints.", i, diagnostics$stable_check_count))
      break
    }
  }
  
  diagnostic_results <- finalize_mcmc_diagnostics(diagnostics, i, Noise_parameters)
  
  Circuits_linkage_posterior <- list(
    TF_Peak_Binding_prob = diagnostic_results$B_prob,
    Peak_Gene_Looping_prob = diagnostic_results$L_prob,
    TF_Peak_Binding_switch_rate = diagnostic_results$B_switch_rate,
    Peak_Gene_Looping_switch_rate = diagnostic_results$L_switch_rate,
    Noise_parameters = Noise_parameters[seq_len(i), , drop = FALSE],
    convergence_log = diagnostic_results$convergence_log,
    window_log = diagnostic_results$window_log,
    window_size = window_size,
    burnin = burnin,
    n_iterations_completed = i,
    n_post_burnin_iterations = diagnostic_results$n_kept,
    reached_stability_criterion = diagnostic_results$stable_check_count >= required_stable_checks,
    null_MSE_A = null_MSE_A,
    null_MSE_R = null_MSE_R,
    observed_sd_A = observed_sd_A,
    observed_sd_R = observed_sd_R,
    final_B = B,
    final_B_state = B_state,
    final_L = L,
    final_L_state = L_state,
    final_TFA = TFA,
    diagnostic_schedule = list(
      n_postburnin_planned = n_postburnin_planned,
      planned_checkpoints = planned_checkpoints,
      planned_complete_windows = planned_complete_windows
    )
  )
  
  return(Circuits_linkage_posterior)
}

#' MAGICAL output
#'
#' @param Output_file_path The output file path
#' @param Candidate_circuits The output from "Candidate_circuits" (with or without TAD)
#' @param Circuits_linkage_posterior The output from `MAGICAL_estimation`
#' @param prob_threshold_TF_peak_binding The threshold of TF-peak binding probability. Default is 0.7
#' @param prob_threshold_peak_gene_looping The threshold of peak-gene looping probability. Default is 0.85
#'
#' @return A data.frame with following columns: Gene_symbol, Gene_chr, Gene_TSS, Peak_chr, Peak_start, Peak_end,  Looping_prob, TFs(binding prob)
#'
#' @import Matrix
#' @importFrom Matrix colSums rowSums sparseMatrix as.matrix
#'
#' @export
MAGICAL_circuits_output <- function(Output_file_path, Candidate_circuits, Circuits_linkage_posterior,
                                    prob_threshold_TF_peak_binding = 0.7, prob_threshold_peak_gene_looping = 0.85) {
  # Initialize results data.frame
  results <- data.frame(
    Gene_symbol = character(),
    Gene_chr = character(),
    Gene_TSS = numeric(),
    Peak_chr = character(),
    Peak_start = numeric(),
    Peak_end = numeric(),
    Looping_prob = numeric(),
    TFs_binding_prob = character(),
    stringsAsFactors = FALSE
  )
  
  Peak_Gene_index <- which(Circuits_linkage_posterior$Peak_Gene_Looping_prob > prob_threshold_peak_gene_looping, arr.ind = TRUE)
  circuit_flag <- matrix(0, nrow = nrow(Peak_Gene_index), ncol = 1)
  
  TF_vector <- matrix(0, nrow(Candidate_circuits$TFs), ncol = 1)
  for (i in 1:nrow(Peak_Gene_index)) {
    TF_prob <- sort.int(Circuits_linkage_posterior$TF_Peak_Binding_prob[Peak_Gene_index[i, 1], ], decreasing = TRUE, index.return = TRUE)
    index <- which(TF_prob$x > prob_threshold_TF_peak_binding)
    
    if (length(index) > 0) {
      circuit_flag[i] <- 1
      
      TFs_binding <- c()
      for (j in 1:length(index)) {
        TF_binding <- paste(Candidate_circuits$TFs[TF_prob$ix[j], 2], "(", TF_prob$x[j], ")", sep = "")
        TFs_binding <- c(TFs_binding, TF_binding)
        TF_vector[TF_prob$ix[j]] <- TF_vector[TF_prob$ix[j]] + 1
      }
      
      # Add the row to the results data.frame
      results <- rbind(
        results,
        data.frame(
          Gene_symbol = Candidate_circuits$Genes$Gene_symbols[Peak_Gene_index[i, 2]],
          Gene_chr = Candidate_circuits$Genes$chr[Peak_Gene_index[i, 2]],
          Gene_TSS = Candidate_circuits$Genes$TSS[Peak_Gene_index[i, 2]],
          Peak_chr = Candidate_circuits$Peaks$chr[Peak_Gene_index[i, 1]],
          Peak_start = Candidate_circuits$Peaks$point1[Peak_Gene_index[i, 1]],
          Peak_end = Candidate_circuits$Peaks$point2[Peak_Gene_index[i, 1]],
          Looping_prob = Circuits_linkage_posterior$Peak_Gene_Looping_prob[Peak_Gene_index[i, 1], Peak_Gene_index[i, 2]],
          TFs_binding_prob = paste(TFs_binding, collapse = ", "),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  # Write the data.frame to the output file
  write.table(
    results,
    file = Output_file_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  print(paste("MAGICAL selected regulatory circuits with", length(which(TF_vector > 1)), "TFs,", length(unique(Peak_Gene_index[circuit_flag > 0, 1])), "peaks, and", length(unique(Peak_Gene_index[circuit_flag > 0, 2])), "genes", sep = " "))
  
  # Return the results data.frame
  return(results)
}
plot_noise <- function(Circuits_linkage_posterior){
  noise <- Circuits_linkage_posterior$Noise_parameters
  post <- seq.int(
    Circuits_linkage_posterior$burnin + 1L,
    nrow(noise)
  )
  
  pdf(
    file = paste0(DIR2, "/", t, "_noise_vs_MSE.pdf"),
    width = 10,
    height = 8
  )
  
  par(mfrow = c(2, 2))
  
  plot(
    post,
    noise[post, "sigma_A_noise"],
    type = "l",
    xlab = "Iteration",
    ylab = "Residual variance",
    main = "ATAC: sampled variance versus MSE"
  )
  lines(post, noise[post, "MSE_A"], col = "red")
  legend(
    "topright",
    legend = c("Sampled variance", "MSE"),
    col = c("black", "red"),
    lty = 1,
    bty = "n"
  )
  
  plot(
    post,
    noise[post, "sigma_R_noise"],
    type = "l",
    xlab = "Iteration",
    ylab = "Residual variance",
    main = "RNA: sampled variance versus MSE"
  )
  lines(post, noise[post, "MSE_R"], col = "red")
  legend(
    "topright",
    legend = c("Sampled variance", "MSE"),
    col = c("black", "red"),
    lty = 1,
    bty = "n"
  )
  
  plot(
    post,
    noise[post, "sigma_A_noise"] / noise[post, "MSE_A"],
    type = "l",
    xlab = "Iteration",
    ylab = "Variance / MSE",
    main = "ATAC ratio"
  )
  abline(h = 1, col = "red", lty = 2)
  
  plot(
    post,
    noise[post, "sigma_R_noise"] / noise[post, "MSE_R"],
    type = "l",
    xlab = "Iteration",
    ylab = "Variance / MSE",
    main = "RNA ratio"
  )
  abline(h = 1, col = "red", lty = 2)
  
  dev.off()  

}
plot_convergence <- function(
    fit,
    stability_tolerance = 0.005,
    output_file =
      "MAGICAL_convergence_plots/MAGICAL_convergence_diagnostics.pdf",
    full_windows_only = TRUE
) {
  
  # --------------------------------------------------------
  # Validate input
  # --------------------------------------------------------
  
  if (!is.list(fit)) {
    stop(
      "`fit` must be the list returned by ",
      "MAGICAL_estimation()."
    )
  }
  
  required_elements <- c(
    "Noise_parameters",
    "convergence_log",
    "window_log",
    "burnin"
  )
  
  missing_elements <- setdiff(
    required_elements,
    names(fit)
  )
  
  if (length(missing_elements) > 0L) {
    stop(
      "The following elements are missing from `fit`: ",
      paste(missing_elements, collapse = ", ")
    )
  }
  
  convergence_log <- fit$convergence_log
  window_log <- fit$window_log
  noise <- fit$Noise_parameters
  
  if (!is.matrix(noise) && !is.data.frame(noise)) {
    stop("`fit$Noise_parameters` must be a matrix or data frame.")
  }
  
  if (ncol(noise) < 2L) {
    stop("`fit$Noise_parameters` must have at least two columns.")
  }
  
  # --------------------------------------------------------
  # Prepare output directory
  # --------------------------------------------------------
  
  output_dir <- dirname(output_file)
  
  if (!dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  # --------------------------------------------------------
  # Select full windows if requested
  # --------------------------------------------------------
  
  if (!is.null(window_log) && nrow(window_log) > 0L) {
    
    if (full_windows_only &&
        "window_n" %in% names(window_log) &&
        !is.null(fit$window_size)) {
      
      window_plot_data <- window_log[
        window_log$window_n == fit$window_size,
        ,
        drop = FALSE
      ]
      
    } else {
      
      window_plot_data <- window_log
    }
    
  } else {
    
    window_plot_data <- data.frame()
  }
  
  # --------------------------------------------------------
  # Print diagnostic information
  # --------------------------------------------------------
  
  cat("\n========================================\n")
  cat("MAGICAL convergence diagnostics\n")
  cat("========================================\n\n")
  
  cat("Output file:\n")
  cat(output_file, "\n\n")
  
  cat("Iterations completed: ",
      fit$n_iterations_completed,
      "\n",
      sep = ""
  )
  
  cat("Burn-in: ",
      fit$burnin,
      "\n",
      sep = ""
  )
  
  cat("Post-burn-in iterations: ",
      fit$n_post_burnin_iterations,
      "\n",
      sep = ""
  )
  
  if (!is.null(fit$window_size)) {
    cat("Window size: ",
        fit$window_size,
        "\n",
        sep = ""
    )
  }
  
  cat("Stability tolerance: ",
      stability_tolerance,
      "\n",
      sep = ""
  )
  
  if (!is.null(fit$reached_stability_criterion)) {
    cat(
      "Reached cumulative stability criterion: ",
      fit$reached_stability_criterion,
      "\n",
      sep = ""
    )
  }
  
  cat("\nCumulative convergence log:\n")
  
  if (!is.null(convergence_log) &&
      nrow(convergence_log) > 0L) {
    
    print(convergence_log)
    
    cat("\nLast cumulative checkpoint:\n")
    
    print(
      tail(convergence_log, 1L)
    )
    
  } else {
    
    cat(
      "No cumulative checkpoint comparisons available.\n"
    )
  }
  
  cat("\nWindow diagnostics:\n")
  
  if (!is.null(window_log) &&
      nrow(window_log) > 0L) {
    
    print(window_log)
    
    cat("\nLast window:\n")
    
    print(
      tail(window_log, 1L)
    )
    
    if (full_windows_only &&
        nrow(window_plot_data) < nrow(window_log)) {
      
      cat(
        "\nNote: incomplete final windows are printed, ",
        "but excluded from the plots.\n",
        sep = ""
      )
    }
    
  } else {
    
    cat("No window summaries available.\n")
  }
  
  # Print post-burn-in noise summaries
  post_start <- fit$burnin + 1L
  post_end <- nrow(noise)
  
  if (post_start <= post_end) {
    
    post_noise <- noise[
      post_start:post_end,
      ,
      drop = FALSE
    ]
    
    cat("\nPost-burn-in ATAC noise summary:\n")
    
    print(
      quantile(
        post_noise[, 1],
        probs = c(
          0,
          0.25,
          0.50,
          0.75,
          0.95,
          0.99,
          1
        ),
        na.rm = TRUE
      )
    )
    
    cat("\nPost-burn-in RNA noise summary:\n")
    
    print(
      quantile(
        post_noise[, 2],
        probs = c(
          0,
          0.25,
          0.50,
          0.75,
          0.95,
          0.99,
          1
        ),
        na.rm = TRUE
      )
    )
  }
  
  # --------------------------------------------------------
  # Determine which panels contain usable data
  # --------------------------------------------------------
  
  has_B_cumulative <- !is.null(convergence_log) &&
    nrow(convergence_log) > 0L &&
    all(
      c(
        "iteration",
        "q99_B_probability_change"
      ) %in% names(convergence_log)
    ) &&
    any(
      is.finite(
        convergence_log$
          q99_B_probability_change
      )
    )
  
  has_L_cumulative <- !is.null(convergence_log) &&
    nrow(convergence_log) > 0L &&
    all(
      c(
        "iteration",
        "q99_L_probability_change"
      ) %in% names(convergence_log)
    ) &&
    any(
      is.finite(
        convergence_log$
          q99_L_probability_change
      )
    )
  
  has_B_window <- nrow(window_plot_data) > 0L &&
    all(
      c(
        "window_end",
        "q99_B_change_from_previous_window"
      ) %in% names(window_plot_data)
    ) &&
    any(
      is.finite(
        window_plot_data$
          q99_B_change_from_previous_window
      )
    )
  
  has_L_window <- nrow(window_plot_data) > 0L &&
    all(
      c(
        "window_end",
        "q99_L_change_from_previous_window"
      ) %in% names(window_plot_data)
    ) &&
    any(
      is.finite(
        window_plot_data$
          q99_L_change_from_previous_window
      )
    )
  
  has_ATAC_noise <- nrow(noise) > 0L &&
    any(is.finite(noise[, 1]))
  
  has_RNA_noise <- nrow(noise) > 0L &&
    any(is.finite(noise[, 2]))
  
  has_ATAC_window_noise <- nrow(window_plot_data) > 0L &&
    all(
      c(
        "window_end",
        "median_ATAC_noise",
        "q95_ATAC_noise",
        "q99_ATAC_noise"
      ) %in% names(window_plot_data)
    ) &&
    any(
      is.finite(
        window_plot_data$median_ATAC_noise
      )
    )
  
  has_RNA_window_noise <- nrow(window_plot_data) > 0L &&
    all(
      c(
        "window_end",
        "median_RNA_noise",
        "q95_RNA_noise",
        "q99_RNA_noise"
      ) %in% names(window_plot_data)
    ) &&
    any(
      is.finite(
        window_plot_data$median_RNA_noise
      )
    )
  
  # --------------------------------------------------------
  # Helper for an empty panel
  # --------------------------------------------------------
  
  empty_panel <- function(title_text, message_text) {
    plot.new()
    title(title_text)
    text(
      0.5,
      0.5,
      message_text,
      cex = 0.9
    )
  }
  
  # --------------------------------------------------------
  # Open PDF
  # --------------------------------------------------------
  
  pdf(
    file = output_file,
    width = 12,
    height = 16,
    onefile = TRUE
  )
  
  on.exit(dev.off(), add = TRUE)
  
  par(
    mfrow = c(4, 2),
    mar = c(4.5, 4.5, 3.5, 1.5) + 0.1,
    oma = c(0, 0, 2, 0)
  )
  
  # ========================================================
  # Row 1: cumulative link-probability stabilization
  # ========================================================
  
  if (has_B_cumulative) {
    
    keep <- is.finite(
      convergence_log$iteration
    ) &
      is.finite(
        convergence_log$
          q99_B_probability_change
      )
    
    plot(
      convergence_log$iteration[keep],
      convergence_log$
        q99_B_probability_change[keep],
      type = "l",
      xlab = "Iteration",
      ylab = "99th percentile change",
      main = "TF-peak cumulative stabilization"
    )
    
    abline(
      h = stability_tolerance,
      col = "red",
      lty = 2
    )
    
  } else {
    
    empty_panel(
      "TF-peak cumulative stabilization",
      paste(
        "No cumulative checkpoint comparisons.",
        "At least two checkpoints are required.",
        sep = "\n"
      )
    )
  }
  
  if (has_L_cumulative) {
    
    keep <- is.finite(
      convergence_log$iteration
    ) &
      is.finite(
        convergence_log$
          q99_L_probability_change
      )
    
    plot(
      convergence_log$iteration[keep],
      convergence_log$
        q99_L_probability_change[keep],
      type = "l",
      xlab = "Iteration",
      ylab = "99th percentile change",
      main = "Peak-gene cumulative stabilization"
    )
    
    abline(
      h = stability_tolerance,
      col = "red",
      lty = 2
    )
    
  } else {
    
    empty_panel(
      "Peak-gene cumulative stabilization",
      paste(
        "No cumulative checkpoint comparisons.",
        "At least two checkpoints are required.",
        sep = "\n"
      )
    )
  }
  
  # ========================================================
  # Row 2: adjacent-window link-probability differences
  # ========================================================
  
  if (has_B_window) {
    
    keep <- is.finite(
      window_plot_data$window_end
    ) &
      is.finite(
        window_plot_data$
          q99_B_change_from_previous_window
      )
    
    plot(
      window_plot_data$window_end[keep],
      window_plot_data$
        q99_B_change_from_previous_window[keep],
      type = "b",
      pch = 16,
      xlab = "Window end iteration",
      ylab = "99th percentile difference",
      main = "TF-peak adjacent-window differences"
    )
    
    abline(
      h = stability_tolerance,
      col = "red",
      lty = 2
    )
    
  } else {
    
    empty_panel(
      "TF-peak adjacent-window differences",
      paste(
        "No adjacent-window comparison available.",
        "At least two completed windows are required.",
        sep = "\n"
      )
    )
  }
  
  if (has_L_window) {
    
    keep <- is.finite(
      window_plot_data$window_end
    ) &
      is.finite(
        window_plot_data$
          q99_L_change_from_previous_window
      )
    
    plot(
      window_plot_data$window_end[keep],
      window_plot_data$
        q99_L_change_from_previous_window[keep],
      type = "b",
      pch = 16,
      xlab = "Window end iteration",
      ylab = "99th percentile difference",
      main = "Peak-gene adjacent-window differences"
    )
    
    abline(
      h = stability_tolerance,
      col = "red",
      lty = 2
    )
    
  } else {
    
    empty_panel(
      "Peak-gene adjacent-window differences",
      paste(
        "No adjacent-window comparison available.",
        "At least two completed windows are required.",
        sep = "\n"
      )
    )
  }
  
  # ========================================================
  # Row 3: complete noise traces
  # ========================================================
  
  if (has_ATAC_noise) {
    
    keep <- is.finite(noise[, 1])
    
    plot(
      which(keep),
      noise[keep, 1],
      type = "l",
      log = "y",
      xlab = "Iteration",
      ylab = "ATAC residual variance",
      main = "ATAC noise trace, log scale"
    )
    
    if (length(fit$burnin) == 1L &&
        is.finite(fit$burnin)) {
      
      abline(
        v = fit$burnin,
        col = "red",
        lty = 2
      )
    }
    
  } else {
    
    empty_panel(
      "ATAC noise trace",
      "No finite ATAC noise values"
    )
  }
  
  if (has_RNA_noise) {
    
    keep <- is.finite(noise[, 2])
    
    plot(
      which(keep),
      noise[keep, 2],
      type = "l",
      log = "y",
      xlab = "Iteration",
      ylab = "RNA residual variance",
      main = "RNA noise trace, log scale"
    )
    
    if (length(fit$burnin) == 1L &&
        is.finite(fit$burnin)) {
      
      abline(
        v = fit$burnin,
        col = "red",
        lty = 2
      )
    }
    
  } else {
    
    empty_panel(
      "RNA noise trace",
      "No finite RNA noise values"
    )
  }
  
  # ========================================================
  # Row 4: window-based noise summaries
  # ========================================================
  
  if (has_ATAC_window_noise) {
    
    keep <- is.finite(
      window_plot_data$window_end
    ) &
      is.finite(
        window_plot_data$median_ATAC_noise
      ) &
      is.finite(
        window_plot_data$q95_ATAC_noise
      ) &
      is.finite(
        window_plot_data$q99_ATAC_noise
      )
    
    x <- window_plot_data$window_end[keep]
    
    y_values <- c(
      window_plot_data$median_ATAC_noise[keep],
      window_plot_data$q95_ATAC_noise[keep],
      window_plot_data$q99_ATAC_noise[keep]
    )
    
    y_values <- y_values[
      is.finite(y_values) &
        y_values > 0
    ]
    
    plot(
      x,
      window_plot_data$median_ATAC_noise[keep],
      type = "b",
      pch = 16,
      log = "y",
      ylim = range(y_values),
      xlab = "Window end iteration",
      ylab = "ATAC residual variance",
      main = "ATAC noise by window"
    )
    
    lines(
      x,
      window_plot_data$q95_ATAC_noise[keep],
      type = "b",
      pch = 17,
      col = "blue"
    )
    
    lines(
      x,
      window_plot_data$q99_ATAC_noise[keep],
      type = "b",
      pch = 15,
      col = "red"
    )
    
    legend(
      "topright",
      legend = c(
        "Median",
        "95th percentile",
        "99th percentile"
      ),
      col = c(
        "black",
        "blue",
        "red"
      ),
      lty = 1,
      pch = c(16, 17, 15),
      bty = "n",
      cex = 0.8
    )
    
  } else {
    
    empty_panel(
      "ATAC noise by window",
      "No completed window-noise summaries"
    )
  }
  
  if (has_RNA_window_noise) {
    
    keep <- is.finite(
      window_plot_data$window_end
    ) &
      is.finite(
        window_plot_data$median_RNA_noise
      ) &
      is.finite(
        window_plot_data$q95_RNA_noise
      ) &
      is.finite(
        window_plot_data$q99_RNA_noise
      )
    
    x <- window_plot_data$window_end[keep]
    
    y_values <- c(
      window_plot_data$median_RNA_noise[keep],
      window_plot_data$q95_RNA_noise[keep],
      window_plot_data$q99_RNA_noise[keep]
    )
    
    y_values <- y_values[
      is.finite(y_values) &
        y_values > 0
    ]
    
    plot(
      x,
      window_plot_data$median_RNA_noise[keep],
      type = "b",
      pch = 16,
      log = "y",
      ylim = range(y_values),
      xlab = "Window end iteration",
      ylab = "RNA residual variance",
      main = "RNA noise by window"
    )
    
    lines(
      x,
      window_plot_data$q95_RNA_noise[keep],
      type = "b",
      pch = 17,
      col = "blue"
    )
    
    lines(
      x,
      window_plot_data$q99_RNA_noise[keep],
      type = "b",
      pch = 15,
      col = "red"
    )
    
    legend(
      "topright",
      legend = c(
        "Median",
        "95th percentile",
        "99th percentile"
      ),
      col = c(
        "black",
        "blue",
        "red"
      ),
      lty = 1,
      pch = c(16, 17, 15),
      bty = "n",
      cex = 0.8
    )
    
  } else {
    
    empty_panel(
      "RNA noise by window",
      "No completed window-noise summaries"
    )
  }
  
  mtext(
    "MAGICAL MCMC convergence and window diagnostics",
    outer = TRUE,
    line = 0.5,
    cex = 1.3,
    font = 2
  )
  
  cat("\nDiagnostic PDF written to:\n")
  cat(normalizePath(output_file, mustWork = FALSE), "\n")
  
  invisible(output_file)
}
