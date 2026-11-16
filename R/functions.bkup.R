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
#' This function constructs candidate circuits based on the input data WITH TAD prior.
#'
#' @seealso [Candidate_circuits_construction_without_TAD()] for the version WITHOUT TAD.
#'
#' @param loaded_data The output from `Data_loading`
#' @param TAD_file_path The path to the TAD prior file
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

Candidate_circuits_construction_with_TAD <- function(loaded_data, TAD_file_path, TF_pct_in_candidate_Peaks_cutoff = 0.05, TF_num_cutoff = 30, TF_enrichment_FC_cutoff = 0.8) {
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
  print("loading TAD data ...")
  cat("\n")
  TAD <- read.table(TAD_file_path, sep = "\t")
  colnames(TAD) <- c("chr", "left_boundary", "right_boundary")
  print(paste(nrow(TAD), "TAD segments provided", sep = " "))
  cat("\n")
  print("MAGICAL model initialization ...")
  cat("\n")

  # TF-Peak binding
  library(dplyr)
  Candidate_Peaks <- inner_join(scATAC_Peaks, Candidate_Peaks, c("chr", "point1", "point2"))
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
  Candidate_Peak_Gene_looping_TAD <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(nrow(Candidate_Peaks), nrow(Candidate_Genes)), repr = "T")
  for (t in 1:nrow(TAD)) {
    Peak_index <- which(Candidate_Peaks$chr == TAD$chr[t] & Candidate_Peaks$point1 > TAD$left_boundary[t] & Candidate_Peaks$point2 < TAD$right_boundary[t])
    Gene_index <- which(Candidate_Genes$chr == TAD$chr[t] & Candidate_Genes$TSS > TAD$left_boundary[t] & Candidate_Genes$TSS < TAD$right_boundary[t])
    if (length(Peak_index) > 0 && length(Gene_index) > 0) {
      Candidate_Peak_Gene_looping_TAD[Peak_index, Gene_index] <- 1
    }
  }


  # some TAD could be extra wide, based on HiC interaction scale, we limit the Peak-Gene distance upto 1M bps
  Candidate_Peak_Gene_looping_distance_control <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(nrow(Candidate_Peaks), nrow(Candidate_Genes)), repr = "T")
  for (g in 1:length(Candidate_Genes$Gene_symbols)) {
    index <- which(Candidate_Peaks$chr == Candidate_Genes$chr[g] & abs((Candidate_Peaks$point1 + Candidate_Peaks$point2) / 2 - Candidate_Genes$TSS[g]) < 1e6)
    if (length(index) > 0) {
      Candidate_Peak_Gene_looping_distance_control[index, g] <- 1
    }
  }

  Candidate_Peak_Gene_looping <- as(Candidate_Peak_Gene_looping_TAD * Candidate_Peak_Gene_looping_distance_control, "TsparseMatrix")


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
    if(length(totake) == 1) {
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
   if(length(totake) == 1) {
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
  selected_TFs <- intersect(Candidate_TFs$name, scRNA_Genes$Gene_symbols)

  index <- match(selected_TFs, Candidate_TFs$name)
  Candidate_TFs <- Candidate_TFs[index, ]
  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[, index]

  index <- match(selected_TFs, scRNA_Genes$Gene_symbols)
  Candidate_TF_log2Count <- RNA_log2[index, ]

  TF_index <- which(rowSums(Candidate_TF_log2Count) > 0)
  Candidate_TFs <- Candidate_TFs[TF_index, ]
  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[, TF_index]
  Candidate_TF_log2Count <- Candidate_TF_log2Count[TF_index, ]
  Candidate_TF_log2Count <- sweep(Candidate_TF_log2Count, 1, rowMeans(Candidate_TF_log2Count), "-")


  # final filtering for candidate circuits
  Peak_index <- which(rowSums(Candidate_Peak_Gene_looping) > 0 & rowSums(Candidate_TF_Peak_Binding) > 0)
  Gene_index <- which(colSums(Candidate_Peak_Gene_looping[Peak_index, ]) > 0)
  TF_index <- which(colSums(Candidate_TF_Peak_Binding[Peak_index, ]) > 0)

  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[Peak_index, TF_index]
  Candidate_Peak_Gene_looping <- Candidate_Peak_Gene_looping[Peak_index, Gene_index]

  Candidate_TFs <- Candidate_TFs[TF_index, ]
  Candidate_TF_log2Count <- Candidate_TF_log2Count[TF_index, ]

  Candidate_Peaks <- Candidate_Peaks[Peak_index, ]
  Candidate_Peak_log2Count <- Candidate_Peak_log2Count[Peak_index, ]

  Candidate_Genes <- Candidate_Genes[Gene_index, ]
  Candidate_Gene_log2Count <- Candidate_Gene_log2Count[Gene_index, ]

  Candidate_circuits <- list(
    "TFs" = Candidate_TFs,
    "TF_log2Count" = Candidate_TF_log2Count,
    "Peaks" = Candidate_Peaks,
    "Peak_log2Count" = Candidate_Peak_log2Count,
    "Genes" = Candidate_Genes,
    "Gene_log2Count" = Candidate_Gene_log2Count,
    "TF_Peak_Binding" = Candidate_TF_Peak_Binding,
    "Peak_Gene_looping" = Candidate_Peak_Gene_looping
  )

  return(Candidate_circuits)
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
  selected_TFs <- intersect(Candidate_TFs$name, scRNA_Genes$Gene_symbols)
  index <- match(selected_TFs, Candidate_TFs$name)
  Candidate_TFs <- Candidate_TFs[index, ]
  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[, index]

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
MAGICAL_initialization <- function(loaded_data, Candidate_circuits, T_VAR = 10) {
  Candidate_TFs <- Candidate_circuits$TFs
#  Candidate_TF_log2Count <- Candidate_circuits$TF_log2Count

  Candidate_Peaks <- Candidate_circuits$Peaks
  Candidate_Peak_log2Count <- Candidate_circuits$Peak_log2Count

  Candidate_Genes <- Candidate_circuits$Genes
  Candidate_Gene_log2Count <- Candidate_circuits$Gene_log2Count

  Candidate_TF_Peak_Binding <- Candidate_circuits$TF_Peak_Binding
  Candidate_Peak_Gene_looping <- Candidate_circuits$Peak_Gene_looping


  S <- length(loaded_data$Common_samples)
  M <- nrow(Candidate_TFs) # total number of candidate TFs
  P <- nrow(Candidate_Peaks) # total number of candidate Peaks
  G <- nrow(Candidate_Genes) # total number of candidate Genes

  # TF activity prior: TF RNA expression
#  T_prior <- Candidate_TF_log2Count
#  T_mean <- Candidate_TF_log2Count
  T_MEAN <- 0
  T_mean <- matrix(T_MEAN, nrow = M, ncol = S)
  T_var <- matrix(T_VAR, nrow = M, ncol = S)

  T_prior <- matrix(
    rnorm(
        M * S,
        mean = T_MEAN,
        sd = sqrt(T_VAR)
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

 storage.mode(B_mask) <- "numeric"

 # Prior settings
 sigma_B <- sqrt(10)
 pi_B <- 0.5

 # Fixed prior mean used in every MCMC iteration
 B_mean <- matrix(
     0,
     nrow = P,
     ncol = M,
     dimnames = dimnames(B_mask)
 )

 # Fixed prior variance used in every MCMC iteration
 B_var <- matrix(
     sigma_B^2,
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
    "L_prior" = L_prior,
    "L_mean" = L_mean,
    "L_var" = L_var,
    "L_prob" = L_prob
  )

  return(Initial_model)
}


# Functions within `MAGICAL_estimation`, no need to export
#*************************** TFA estimation ***************************
TF_activity_T_sampling <- function(A, A_sample, ATAC_Cell_Sample_vector, R, R_sample, RNA_Cell_Sample_vector, TFA, T_prior_mean, T_prior_var, B, B_state, sigma_A_noise, P, G, M, S) {
  TF_index <- sample(M)
  for (i in 1:M) {
    # during the looping, TFA is iteratively updated for all M TFs and the update order is random in each round of sampling
    m <- TF_index[i]
    if (sum(B_state[, m]) > 0) {
      temp_var <- sum(B[, m]^2) * T_prior_var[m, ] / sum(B_state[, m]) + sigma_A_noise
      mean_T <- (t(B[, m]) %*% (A_sample - B %*% TFA$T_sample + B[, m] %*% t(TFA$T_sample[m, ])) / sum(B_state[, m]) + T_prior_mean[m, ] * sigma_A_noise) / temp_var
      variance_T <- T_prior_var[m, ] * sigma_A_noise / temp_var

      for (s in 1:S) {
        aa <- rnorm(1)
        if (aa - 3 > 0) {
          aa <- 3
        }
        if (aa + 3 < 0) {
          aa <- -3
        }

        TFA$T_sample[m, s] <- aa * sqrt(abs(variance_T[s])) + mean_T[s]

        ATAC_cell_index <- which(ATAC_Cell_Sample_vector == s)
        TFA$T_A[m, ATAC_cell_index] <- rnorm(length(ATAC_cell_index)) * sqrt(abs(variance_T[s])) + TFA$T_sample[m, s]

        RNA_cell_index <- which(RNA_Cell_Sample_vector == s)
        TFA$T_R[m, RNA_cell_index] <- rnorm(length(RNA_cell_index)) * sqrt(abs(variance_T[s])) + TFA$T_sample[m, s]
      }
    }
  }
  return(TFA)
}

#*************************** TF-peak binding confidence estimation ***************************
TF_peak_binding_B_sampling <- function(A, A_sample, ATAC_Cell_Sample_vector, TFA, B, B_state, B_prior_mean, B_prior_var, sigma_A_noise, P, G, M, S) {
  T_sample <- matrix(0, nrow = M, ncol = S)
  for (s in 1:S) {
  totake <- which(ATAC_Cell_Sample_vector == s)
  if (length(totake) == 1) {
    T_sample[, s] <- TFA$T_A[, totake]
  }else{
    T_sample[, s] <- rowMeans(TFA$T_A[, totake])
  }
  }

  TF_index <- sample(M)

  for (i in 1:M) {
    m <- TF_index[i]
    temp_var <- sum(TFA$T_sample[m, ]^2) * B_prior_var[m] / S + sigma_A_noise
    mean_B <- ((A_sample - B %*% TFA$T_sample + B[, m] %*% t(TFA$T_sample[m, ])) %*% TFA$T_sample[m, ] * B_prior_var[m] / S + B_prior_mean[, m] * sigma_A_noise) / temp_var
    vairance_B <- B_prior_var[m] * sigma_A_noise / temp_var

    bb <- rnorm(P)
    bb[which(bb - 3 > 0)] <- 3
    bb[which(bb + 3 < 0)] <- -3

    B[, m] <- as.vector((bb * sqrt(abs(vairance_B)) + mean_B) * B_state[, m])
  }

  return(B)
}

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
    temp_var <- sum(A_estimate[f, ]^2) * L_prior_var / S + sigma_R_noise
    mean_L <- (A_estimate[f, ] %*% t(R_sample - t(L) %*% A_estimate + L[f, ] %*% t(A_estimate[f, ])) * L_prior_var / S + L_prior_mean[f, ] * sigma_R_noise) / temp_var
    vairance_L <- L_prior_var * sigma_R_noise / temp_var

    ll <- rnorm(G)
    ll[which(ll - 3 > 0)] <- 3
    ll[which(ll + 3 < 0)] <- -3

    L[f, ] <- as.vector((ll * sqrt(vairance_L) + mean_L) * L_state[f, ])
  }
  return(L)
}

#*************************** TF-peak binding state update ***************************
TF_peak_binary_binding_B_state_sampling <- function(A, A_sample, ATAC_Cell_Sample_vector, TFA,
                                                    B, B_state, B_prior_mean, B_prior_var, B_prior_prob, sigma_A_noise, P, G, M, S) {
  T_sample <- matrix(0, nrow = M, ncol = S)
  for (s in 1:S) {
    totake <- which(ATAC_Cell_Sample_vector == s)
    if (length(totake) == 1) {
      T_sample[, s] <- TFA$T_A[, totake]
    }else{
      T_sample[, s] <- rowMeans(TFA$T_A[, which(ATAC_Cell_Sample_vector == s)])
    }
  }

  TF_index <- sample(M)
  Peak_index <- sample(P)

  for (i in 1:P) {
    f <- Peak_index[i]
    temp <- t(B[f, ]) %*% TFA$T_sample

    for (j in 1:M) {
      m <- TF_index[j]

      temp_var <- sum(TFA$T_sample[m, ]^2) * B_prior_var[m] / S + sigma_A_noise

      if (B_state[f, m] > 0) {
        mean_B <- (sum((A_sample[f, ] - temp + B[f, m] * TFA$T_sample[m, ]) * TFA$T_sample[m, ]) * B_prior_var[m] / S + B_prior_mean[f, m] * sigma_A_noise) / temp_var

        vairance_B <- B_prior_var[m] * sigma_A_noise / temp_var

        post_b1 <- exp(-(B[f, m] * 1 - mean_B)^2 / (2 * vairance_B)) * (B_prior_prob[f, m] + 0.25) + 1e-6

        post_b0 <- exp(-(B[f, m] * 0 - mean_B)^2 / (2 * vairance_B)) * (1 - B_prior_prob[f, m] + 0.25) + 1e-6

        P1 <- post_b1 / (post_b1 + post_b0)

        threshold_c <- runif(1)

        if (is.na(P1) || is.nan(P1)) {
          P1 <- 0.5
        }

        if (P1 < threshold_c) {
          B[f, m] <- 0
          B_state[f, m] <- 0
        } else {
          B[f, m] <- B[f, m]
          B_state[f, m] <- 1
        }
      } #### THIS IS WHERE WE CHANGE FROM ORIGINAL CODE. PER COPILOT: 
        ####      "If an active link is switched off by the first block, the second block immediately runs and can switch the link back on during the same iteration."
       else  if (B_state[f, m] == 0 && B_prior_prob[f, m] > 0) {
        mean_B <- (sum((A_sample[f, ] - temp) * TFA$T_sample[m, ]) * B_prior_var[m] / S + B_prior_mean[f, m] * sigma_A_noise) / temp_var
        vairance_B <- B_prior_var[m] * sigma_A_noise / temp_var

        bb <- rnorm(1)

        if (bb - 3 > 0) {
          bb <- 3
        }

        if (bb + 3 < 0) {
          bb <- -3
        }

        B_temp <- bb * sqrt(abs(vairance_B)) + mean_B

        post_b1 <- exp(-(bb^2) / 2) * (B_prior_prob[f, m] + 0.25) + 1e-6

        post_b0 <- exp(-(mean_B)^2 / (2 * vairance_B)) * (1 - B_prior_prob[f, m] + 0.25) + 1e-6

        P1 <- post_b1 / (post_b1 + post_b0)

        threshold_c <- runif(1)

        if (is.na(P1) || is.nan(P1)) {
          P1 <- 0.5
        }

        if (P1 < threshold_c) {
          B[f, m] <- 0
          B_state[f, m] <- 0
        } else {
          B[f, m] <- B_temp
          B_state[f, m] <- 1
        }
      }
    }
  }
  Updated_B_and_State <- list("B" = B, "B_state" = B_state)
  return(Updated_B_and_State)
}
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
      temp_var <- sum(A_estimate[f, ]^2) * L_prior_var / S + sigma_R_noise

      if (L_state[f, g] > 0) {
        mean_L <- (sum((R_sample[g, ] - temp + L[f, g] * A_estimate[f, ]) * A_estimate[f, ]) * L_prior_var / S + L_prior_mean[f, g] * sigma_R_noise) / temp_var
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
        mean_L <- (sum((R_sample[g, ] - temp) * A_estimate[f, ]) * L_prior_var / S + L_prior_mean[f, g] * sigma_R_noise) / temp_var
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

  # ------------------------------------------------------------
  # Validate MCMC and convergence-diagnostic settings
  # ------------------------------------------------------------
  
  iteration_num <- as.integer(iteration_num)
  burnin <- as.integer(burnin)
  required_stable_checks <- as.integer(
    required_stable_checks
  )
  
  if (iteration_num < 1L) {
    stop("`iteration_num` must be at least 1.")
  }
  
  if (burnin < 0L || burnin >= iteration_num) {
    stop(
      "`burnin` must be non-negative and smaller ",
      "than `iteration_num`."
    )
  }
  
  if (required_stable_checks < 1L) {
    stop("`required_stable_checks` must be at least 1.")
  }
  
  diagnostic_intervals <- choose_diagnostic_intervals(
    iteration_num = iteration_num,
    burnin = burnin,
    target_checkpoints = target_checkpoints,
    target_windows = target_windows
  )
  
  n_postburnin_planned <-
    diagnostic_intervals$n_postburnin
  
  # Use automatically selected values only when the user did not
  # explicitly supply a value.
  if (is.null(check_every)) {
    check_every <-
      diagnostic_intervals$check_every
  }
  
  if (is.null(window_size)) {
    window_size <-
      diagnostic_intervals$window_size
  }
  
  check_every <- as.integer(check_every)
  window_size <- as.integer(window_size)
  
  if (length(check_every) != 1L ||
      !is.finite(check_every) ||
      check_every < 1L) {
    
    stop("`check_every` must be one positive integer.")
  }
  
  if (length(window_size) != 1L ||
      !is.finite(window_size) ||
      window_size < 2L) {
    
    stop("`window_size` must be one integer of at least 2.")
  }
  
  if (check_every > n_postburnin_planned) {
    stop(
      "`check_every` is larger than the number of ",
      "planned post-burn-in iterations."
    )
  }
  
  if (window_size > n_postburnin_planned) {
    stop(
      "`window_size` is larger than the number of ",
      "planned post-burn-in iterations."
    )
  }
  
  planned_checkpoints <- floor(
    n_postburnin_planned / check_every
  )
  
  planned_complete_windows <- floor(
    n_postburnin_planned / window_size
  )
  
  # The first cumulative checkpoint only creates a baseline.
  # At least one additional checkpoint is needed to calculate change.
  minimum_checkpoints <-
    required_stable_checks + 1L
  
  if (planned_checkpoints < minimum_checkpoints) {
    warning(
      "The current settings provide only ",
      planned_checkpoints,
      " cumulative checkpoints. At least ",
      minimum_checkpoints,
      " are required to satisfy ",
      "`required_stable_checks = ",
      required_stable_checks,
      "`."
    )
  }
  
  # Adjacent-window comparison requires at least two complete windows.
  if (planned_complete_windows < 2L) {
    warning(
      "The current settings provide only ",
      planned_complete_windows,
      " complete post-burn-in window(s). ",
      "At least two are needed for an adjacent-window comparison."
    )
  }
  
  message(
    paste0(
      "MCMC diagnostic schedule:\n",
      "  total iterations: ", iteration_num, "\n",
      "  burn-in: ", burnin, "\n",
      "  planned post-burn-in iterations: ",
      n_postburnin_planned, "\n",
      "  cumulative check interval: ",
      check_every, "\n",
      "  expected cumulative checkpoints: ",
      planned_checkpoints, "\n",
      "  window size: ", window_size, "\n",
      "  expected complete windows: ",
      planned_complete_windows
    )
  )
  
   
  Common_samples <- loaded_data$Common_samples
  S <- length(Common_samples) # number of samples

  # ATAC data
  Peak_index <- match(Candidate_circuits$Peaks$Peak_index, loaded_data$scATAC_Peaks$Peak_index)
  A <- loaded_data$scATAC_read_count_matrix[Peak_index, ]
  ATAC_Cell_Sample_vector <- matrix(0, nrow = 1, ncol = nrow(loaded_data$scATAC_cells))
  for (s in 1:S) {
    ATAC_Cell_Sample_vector[which(loaded_data$scATAC_cells$subject_ID == Common_samples[s])] <- s
  } # ATAC cells of each sample
  A_sample <- Candidate_circuits$Peak_log2Count
  P <- length(Peak_index) # number of peaks



  # RNA data
  Gene_index <- match(Candidate_circuits$Genes$Gene_symbols, loaded_data$scRNA_Genes$Gene_symbols)
  R <- loaded_data$scRNA_read_count_matrix[Gene_index, ]
  RNA_Cell_Sample_vector <- matrix(0, nrow = 1, ncol = nrow(loaded_data$scRNA_cells))
  for (s in 1:S) {
    RNA_Cell_Sample_vector[1, which(loaded_data$scRNA_cells$subject_ID == Common_samples[s])] <- s
  } # RNA cells of each sample
  R_sample <- Candidate_circuits$Gene_log2Count
  G <- length(Gene_index) # number of genes


  # TF-peak binding initial values
  B_prior_mean <- Initial_model$B_mean
  B_prior_var <- Initial_model$B_var
  B_prior_prob <- Initial_model$B_prob
 # B <- Initial_model$B_prior
  B <- Initial_model$B_initial
  B_state <- as.matrix(Candidate_circuits$TF_Peak_Binding)



  # Peak-gene looping initial values
  L_prior_mean <- Initial_model$L_mean
  L_prior_var <- Initial_model$L_var
  L_prior_prob <- Initial_model$L_prob
  L <- Initial_model$L_prior
  L_state <- as.matrix(Candidate_circuits$Peak_Gene_looping)



  # TFA initial values
  T_prior_mean <- Initial_model$T_mean
  T_prior_var <- Initial_model$T_var
  T_sample <- Initial_model$T_prior
  M <- nrow(T_sample) # number of TFs

  T_A <- matrix(0, nrow = nrow(Candidate_circuits$TFs), ncol = nrow(loaded_data$scATAC_cells))
  for (s in 1:S) {
    for (m in 1:M) {
      index <- which(ATAC_Cell_Sample_vector == s)
      T_A[m, index] <- rnorm(length(index), mean = T_prior_mean[m, s], sd = sqrt(abs(T_prior_var[m, s])))
    }
  } # initial values of the hidden TF activity of ATAC cells for each sample

  T_R <- matrix(0, nrow = nrow(Candidate_circuits$TFs), ncol = nrow(loaded_data$scRNA_cells))
  for (s in 1:S) {
    for (m in 1:M) {
      index <- which(RNA_Cell_Sample_vector == s)
      T_R[m, index] <- rnorm(length(index), mean = T_prior_mean[m, s], sd = sqrt(abs(T_prior_var[m, s])))
    }
  } # initial values of the hidden TF activity of RNA cells for each sample

  TFA <- list("T_A" = T_A, "T_R" = T_R, "T_sample" = T_sample)



  # ATAC fitting residual variance initial values
  alpha_A <- 1
  beta_A <- 1
  sigma_A_noise <- 1 / rgamma(1, shape = alpha_A + 1 / 2, rate = (beta_A + sum((A_sample - B %*% TFA$T_sample)^2)) / (2 * P * S))



  # RNA fitting residual variance initial values
  alpha_R <- 1
  beta_R <- 1
  sigma_R_noise <- 1 / rgamma(1, shape = alpha_R + 1 / 2, scale = 1 / ((beta_R + sum((R_sample - t(L) %*% (B %*% TFA$T_sample))^2)) / (2 * G * S)))



  # MCMC sampling initial state
  iteration_seg <- max(
    1L,
    floor(iteration_num / 10L)
  )
  B_state_frq <- B_state
  L_state_frq <- L_state
    

    n_kept <- 0L

    
  print("MAGICAL integration starts ...")
  cat("\n")

  #0 - set up tracking of state switching/convergence
  n_kept <- 0L
  stable_check_count <- 0L
  
  previous_B_prob <- NULL
  previous_L_prob <- NULL
  
  previous_B_state <- B_state
  previous_L_state <- L_state
  
  B_state_frq_diagnostic <- matrix(
    0,
    nrow = nrow(B_state),
    ncol = ncol(B_state),
    dimnames = dimnames(B_state)
  )
  
  L_state_frq_diagnostic <- matrix(
    0,
    nrow = nrow(L_state),
    ncol = ncol(L_state),
    dimnames = dimnames(L_state)
  )
  
  B_state_switches <- matrix(
    0L,
    nrow = nrow(B_state),
    ncol = ncol(B_state),
    dimnames = dimnames(B_state)
  )
  
  L_state_switches <- matrix(
    0L,
    nrow = nrow(L_state),
    ncol = ncol(L_state),
    dimnames = dimnames(L_state)
  )
  
  convergence_log <- data.frame(
    iteration = integer(),
    n_kept = integer(),
    max_B_probability_change = numeric(),
    q99_B_probability_change = numeric(),
    max_L_probability_change = numeric(),
    q99_L_probability_change = numeric(),
    mean_B_switch_rate = numeric(),
    mean_L_switch_rate = numeric(),
    sigma_A_recent_mean = numeric(),
    sigma_R_recent_mean = numeric(),
    stable = logical()
  )
  
B_candidate_index <- which(
    as.matrix(Candidate_circuits$TF_Peak_Binding) > 0
)

L_candidate_index <- which(
    as.matrix(Candidate_circuits$Peak_Gene_looping) > 0
)
    
# ------------------------------------------------------------
# Fixed-window diagnostic tracking
# ------------------------------------------------------------

window_n <- 0L
window_start_iteration <- NA_integer_

B_window_sum <- matrix(
  0,
  nrow = nrow(B_state),
  ncol = ncol(B_state),
  dimnames = dimnames(B_state)
)

L_window_sum <- matrix(
  0,
  nrow = nrow(L_state),
  ncol = ncol(L_state),
  dimnames = dimnames(L_state)
)

B_window_switches <- matrix(
  0L,
  nrow = nrow(B_state),
  ncol = ncol(B_state),
  dimnames = dimnames(B_state)
)

L_window_switches <- matrix(
  0L,
  nrow = nrow(L_state),
  ncol = ncol(L_state),
  dimnames = dimnames(L_state)
)

previous_B_window_prob <- NULL
previous_L_window_prob <- NULL

window_log <- data.frame(
  window_start = integer(),
  window_end = integer(),
  window_n = integer(),
  
  q99_B_change_from_previous_window = numeric(),
  max_B_change_from_previous_window = numeric(),
  q99_L_change_from_previous_window = numeric(),
  max_L_change_from_previous_window = numeric(),
  
  mean_B_probability = numeric(),
  mean_L_probability = numeric(),
  mean_B_switch_rate = numeric(),
  mean_L_switch_rate = numeric(),
  
  median_ATAC_noise = numeric(),
  q25_ATAC_noise = numeric(),
  q75_ATAC_noise = numeric(),
  q95_ATAC_noise = numeric(),
  q99_ATAC_noise = numeric(),
  max_ATAC_noise = numeric(),
  
  median_RNA_noise = numeric(),
  q25_RNA_noise = numeric(),
  q75_RNA_noise = numeric(),
  q95_RNA_noise = numeric(),
  q99_RNA_noise = numeric(),
  max_RNA_noise = numeric()
)


  Noise_parameters <- matrix(0, nrow = iteration_num, ncol = 2)
  
  for (i in 1:iteration_num) {
    
    # Step 1: TF activity sampling
    TFA <- TF_activity_T_sampling(
      A, A_sample, ATAC_Cell_Sample_vector, R, R_sample, RNA_Cell_Sample_vector,
      TFA, T_prior_mean, T_prior_var, B, B_state, sigma_A_noise, P, G, M, S
    )


    # Step 2: TF-peak binding weight sampling
    B <- TF_peak_binding_B_sampling(
      A, A_sample, ATAC_Cell_Sample_vector, TFA,
      B, B_state, B_prior_mean, B_prior_var, sigma_A_noise, P, G, M, S
    )


    # Step 3: TF-peak binding state update
    Updated_B_and_State <- TF_peak_binary_binding_B_state_sampling(
      A, A_sample, ATAC_Cell_Sample_vector,
      TFA, B, B_state, B_prior_mean, B_prior_var, B_prior_prob,
      sigma_A_noise, P, G, M, S
    )
    B <- Updated_B_and_State$B
    B_state <- Updated_B_and_State$B_state
    B_state_frq <- B_state_frq + B_state


    # Step 4: Peak-Gene looping weight sampling
    L <- Peak_gene_looping_L_samping(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, sigma_R_noise, P, G, M, S)


    # Step 5: Peak-Gene looping state update
    Updated_L_and_State <- Peak_gene_binary_looping_L_state_samping(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, L_prior_prob, sigma_R_noise, P, G, M, S)
    L <- Updated_L_and_State$L
    L_state <- Updated_L_and_State$L_state
    L_state_frq <- L_state_frq + L_state


    # Step 6: fitting residue variance control
    sigma_A_noise <- 1 / rgamma(1, shape = alpha_A + 1 / 2, rate = (beta_A + sum((A_sample - B %*% TFA$T_sample)^2)) / (2 * P * S))

    sigma_R_noise <- 1 / rgamma(1, shape = alpha_R + 1 / 2, rate = (beta_R + sum((R_sample - t(L) %*% (B %*% TFA$T_sample))^2)) / (2 * G * S))

    Noise_parameters[i, 1] <- sigma_A_noise
    Noise_parameters[i, 2] <- sigma_R_noise
    
    
   # ------------------------------------------------------------
    # Post-burn-in diagnostic accumulation
    # ------------------------------------------------------------
    
    if (i > burnin) {
      
      n_kept <- n_kept + 1L
      
      if (window_n == 0L) {
        window_start_iteration <- i
      }
      
      window_n <- window_n + 1L
      
      B_switch <- abs(B_state - previous_B_state)
      L_switch <- abs(L_state - previous_L_state)
      
      # Cumulative post-burn-in diagnostics
      B_state_frq_diagnostic <-
        B_state_frq_diagnostic + B_state
      
      L_state_frq_diagnostic <-
        L_state_frq_diagnostic + L_state
      
      B_state_switches <-
        B_state_switches + B_switch
      
      L_state_switches <-
        L_state_switches + L_switch
      
      # Current-window diagnostics
      B_window_sum <-
        B_window_sum + B_state
      
      L_window_sum <-
        L_window_sum + L_state
      
      B_window_switches <-
        B_window_switches + B_switch
      
      L_window_switches <-
        L_window_switches + L_switch
    }
    
    previous_B_state <- B_state
    previous_L_state <- L_state
    
    # ------------------------------------------------------------
    # Save a completed fixed-window summary
    # ------------------------------------------------------------
    
    if (i > burnin && window_n == window_size) {
      
      B_window_prob <- B_window_sum / window_n
      L_window_prob <- L_window_sum / window_n
      
      B_window_switch_rate <-
        B_window_switches / window_n
      
      L_window_switch_rate <-
        L_window_switches / window_n
      
      B_candidate_window_prob <-
        B_window_prob[B_candidate_index]
      
      L_candidate_window_prob <-
        L_window_prob[L_candidate_index]
      
      mean_B_window_probability <- mean(
        B_candidate_window_prob,
        na.rm = TRUE
      )
      
      mean_L_window_probability <- mean(
        L_candidate_window_prob,
        na.rm = TRUE
      )
      
      mean_B_window_switch_rate <- mean(
        B_window_switch_rate[B_candidate_index],
        na.rm = TRUE
      )
      
      mean_L_window_switch_rate <- mean(
        L_window_switch_rate[L_candidate_index],
        na.rm = TRUE
      )
      
      if (is.null(previous_B_window_prob)) {
        
        q99_B_window_change <- NA_real_
        max_B_window_change <- NA_real_
        q99_L_window_change <- NA_real_
        max_L_window_change <- NA_real_
        
      } else {
        
        B_window_change <- abs(
          B_window_prob[B_candidate_index] -
            previous_B_window_prob[B_candidate_index]
        )
        
        L_window_change <- abs(
          L_window_prob[L_candidate_index] -
            previous_L_window_prob[L_candidate_index]
        )
        
        q99_B_window_change <- as.numeric(
          quantile(
            B_window_change,
            probs = 0.99,
            na.rm = TRUE,
            names = FALSE
          )
        )
        
        max_B_window_change <- max(
          B_window_change,
          na.rm = TRUE
        )
        
        q99_L_window_change <- as.numeric(
          quantile(
            L_window_change,
            probs = 0.99,
            na.rm = TRUE,
            names = FALSE
          )
        )
        
        max_L_window_change <- max(
          L_window_change,
          na.rm = TRUE
        )
      }
      
      window_noise <- Noise_parameters[
        window_start_iteration:i,
        ,
        drop = FALSE
      ]
      
      window_log <- rbind(
        window_log,
        data.frame(
          window_start = window_start_iteration,
          window_end = i,
          window_n = window_n,
          
          q99_B_change_from_previous_window =
            q99_B_window_change,
          
          max_B_change_from_previous_window =
            max_B_window_change,
          
          q99_L_change_from_previous_window =
            q99_L_window_change,
          
          max_L_change_from_previous_window =
            max_L_window_change,
          
          mean_B_probability =
            mean_B_window_probability,
          
          mean_L_probability =
            mean_L_window_probability,
          
          mean_B_switch_rate =
            mean_B_window_switch_rate,
          
          mean_L_switch_rate =
            mean_L_window_switch_rate,
          
          median_ATAC_noise =
            median(window_noise[, 1], na.rm = TRUE),
          
          q25_ATAC_noise = as.numeric(
            quantile(
              window_noise[, 1],
              0.25,
              na.rm = TRUE,
              names = FALSE
            )
          ),
          
          q75_ATAC_noise = as.numeric(
            quantile(
              window_noise[, 1],
              0.75,
              na.rm = TRUE,
              names = FALSE
            )
          ),
          
          q95_ATAC_noise = as.numeric(
            quantile(
              window_noise[, 1],
              0.95,
              na.rm = TRUE,
              names = FALSE
            )
          ),
          
          q99_ATAC_noise = as.numeric(
            quantile(
              window_noise[, 1],
              0.99,
              na.rm = TRUE,
              names = FALSE
            )
          ),
          
          max_ATAC_noise =
            max(window_noise[, 1], na.rm = TRUE),
          
          median_RNA_noise =
            median(window_noise[, 2], na.rm = TRUE),
          
          q25_RNA_noise = as.numeric(
            quantile(
              window_noise[, 2],
              0.25,
              na.rm = TRUE,
              names = FALSE
            )
          ),
          
          q75_RNA_noise = as.numeric(
            quantile(
              window_noise[, 2],
              0.75,
              na.rm = TRUE,
              names = FALSE
            )
          ),
          
          q95_RNA_noise = as.numeric(
            quantile(
              window_noise[, 2],
              0.95,
              na.rm = TRUE,
              names = FALSE
            )
          ),
          
          q99_RNA_noise = as.numeric(
            quantile(
              window_noise[, 2],
              0.99,
              na.rm = TRUE,
              names = FALSE
            )
          ),
          
          max_RNA_noise =
            max(window_noise[, 2], na.rm = TRUE)
        )
      )
      
      # Save this window for comparison with the next window
      previous_B_window_prob <- B_window_prob
      previous_L_window_prob <- L_window_prob
      
      # Reset the window accumulators
      B_window_sum[,] <- 0
      L_window_sum[,] <- 0
      B_window_switches[,] <- 0L
      L_window_switches[,] <- 0L
      
      window_n <- 0L
      window_start_iteration <- NA_integer_
    }
    
    
    
    # ------------------------------------------------------------
    # Convergence checkpoint
    # ------------------------------------------------------------
    
    if (i > burnin &&
        n_kept >= check_every &&
        n_kept %% check_every == 0L) {
      
      current_B_prob <-
        B_state_frq_diagnostic / n_kept
      
      current_L_prob <-
        L_state_frq_diagnostic / n_kept
      
      if (is.null(previous_B_prob) ||
          is.null(previous_L_prob)) {
        
        message(
          "Iteration ", i,
          ": saving first post-burn-in checkpoint as baseline."
        )
        
      } else {
        
        B_change <- abs(
          current_B_prob - previous_B_prob
        )
        
        L_change <- abs(
          current_L_prob - previous_L_prob
        )
        
        # Restrict diagnostics to candidate links
        B_change_candidates <-
          B_change[B_candidate_index]
        
        L_change_candidates <-
          L_change[L_candidate_index]
        
        max_B_change <- max(
          B_change_candidates,
          na.rm = TRUE
        )
        
        max_L_change <- max(
          L_change_candidates,
          na.rm = TRUE
        )
        
        q99_B_change <- as.numeric(
          quantile(
            B_change_candidates,
            probs = 0.99,
            na.rm = TRUE,
            names = FALSE
          )
        )
        
        q99_L_change <- as.numeric(
          quantile(
            L_change_candidates,
            probs = 0.99,
            na.rm = TRUE,
            names = FALSE
          )
        )
        
        B_switch_rate <-
          B_state_switches / n_kept
        
        L_switch_rate <-
          L_state_switches / n_kept
        
        mean_B_switch_rate <- mean(
          B_switch_rate[B_candidate_index],
          na.rm = TRUE
        )
        
        mean_L_switch_rate <- mean(
          L_switch_rate[L_candidate_index],
          na.rm = TRUE
        )
        
        recent_start <- max(
          burnin + 1L,
          i - check_every + 1L
        )
        
        sigma_A_recent_mean <- mean(
          Noise_parameters[recent_start:i, 1],
          na.rm = TRUE
        )
        
        sigma_R_recent_mean <- mean(
          Noise_parameters[recent_start:i, 2],
          na.rm = TRUE
        )
        
        stable_now <-
          is.finite(q99_B_change) &&
          is.finite(q99_L_change) &&
          q99_B_change <= stability_tolerance &&
          q99_L_change <= stability_tolerance
        
        if (stable_now) {
          stable_check_count <-
            stable_check_count + 1L
        } else {
          stable_check_count <- 0L
        }
        
        convergence_log <- rbind(
          convergence_log,
          data.frame(
            iteration = i,
            n_kept = n_kept,
            max_B_probability_change =
              max_B_change,
            q99_B_probability_change =
              q99_B_change,
            max_L_probability_change =
              max_L_change,
            q99_L_probability_change =
              q99_L_change,
            mean_B_switch_rate =
              mean_B_switch_rate,
            mean_L_switch_rate =
              mean_L_switch_rate,
            sigma_A_recent_mean =
              sigma_A_recent_mean,
            sigma_R_recent_mean =
              sigma_R_recent_mean,
            stable = stable_now
          )
        )
        
        message(sprintf(
          paste0(
            "Iteration %d: ",
            "B q99 change = %.5f; ",
            "L q99 change = %.5f; ",
            "stable checks = %d/%d"
          ),
          i,
          q99_B_change,
          q99_L_change,
          stable_check_count,
          required_stable_checks
        ))
        
        if (stop_when_stable &&
            stable_check_count >=
            required_stable_checks) {
          
          message(sprintf(
            paste0(
              "Stopping at iteration %d after %d ",
              "consecutive stable checkpoints."
            ),
            i,
            stable_check_count
          ))
          
          break
        }
      }
      
      previous_B_prob <- current_B_prob
      previous_L_prob <- current_L_prob
    }
    
    if (i %% iteration_seg == 0) {
      print(paste("MAGICAL finished", 10 * i / iteration_seg, "percent", sep = " "))
      cat("\n")
    }
  }

  # Preserve original MAGICAL probabilities
Candidate_TF_Peak_Binding_prob <-
    B_state_frq / (i + 1L)

Candidate_Peak_Gene_Looping_prob <-
    L_state_frq / (i + 1L)

if (n_kept == 0L) {
    stop("No post-burn-in iterations were retained.")
}

B_state_switch_rate <-
    B_state_switches / n_kept

L_state_switch_rate <-
    L_state_switches / n_kept



# ------------------------------------------------------------
# Save the final incomplete window, if present
# ------------------------------------------------------------

if (window_n > 0L) {
  
  B_window_prob <- B_window_sum / window_n
  L_window_prob <- L_window_sum / window_n
  
  B_window_switch_rate <-
    B_window_switches / window_n
  
  L_window_switch_rate <-
    L_window_switches / window_n
  
  if (is.null(previous_B_window_prob)) {
    
    q99_B_window_change <- NA_real_
    max_B_window_change <- NA_real_
    q99_L_window_change <- NA_real_
    max_L_window_change <- NA_real_
    
  } else {
    
    B_window_change <- abs(
      B_window_prob[B_candidate_index] -
        previous_B_window_prob[B_candidate_index]
    )
    
    L_window_change <- abs(
      L_window_prob[L_candidate_index] -
        previous_L_window_prob[L_candidate_index]
    )
    
    q99_B_window_change <- as.numeric(
      quantile(
        B_window_change,
        0.99,
        na.rm = TRUE,
        names = FALSE
      )
    )
    
    max_B_window_change <- max(
      B_window_change,
      na.rm = TRUE
    )
    
    q99_L_window_change <- as.numeric(
      quantile(
        L_window_change,
        0.99,
        na.rm = TRUE,
        names = FALSE
      )
    )
    
    max_L_window_change <- max(
      L_window_change,
      na.rm = TRUE
    )
  }
  
  window_noise <- Noise_parameters[
    window_start_iteration:i,
    ,
    drop = FALSE
  ]
  
  window_log <- rbind(
    window_log,
    data.frame(
      window_start = window_start_iteration,
      window_end = i,
      window_n = window_n,
      
      q99_B_change_from_previous_window =
        q99_B_window_change,
      
      max_B_change_from_previous_window =
        max_B_window_change,
      
      q99_L_change_from_previous_window =
        q99_L_window_change,
      
      max_L_change_from_previous_window =
        max_L_window_change,
      
      mean_B_probability = mean(
        B_window_prob[B_candidate_index],
        na.rm = TRUE
      ),
      
      mean_L_probability = mean(
        L_window_prob[L_candidate_index],
        na.rm = TRUE
      ),
      
      mean_B_switch_rate = mean(
        B_window_switch_rate[B_candidate_index],
        na.rm = TRUE
      ),
      
      mean_L_switch_rate = mean(
        L_window_switch_rate[L_candidate_index],
        na.rm = TRUE
      ),
      
      median_ATAC_noise =
        median(window_noise[, 1], na.rm = TRUE),
      
      q25_ATAC_noise = as.numeric(
        quantile(
          window_noise[, 1],
          0.25,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      
      q75_ATAC_noise = as.numeric(
        quantile(
          window_noise[, 1],
          0.75,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      
      q95_ATAC_noise = as.numeric(
        quantile(
          window_noise[, 1],
          0.95,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      
      q99_ATAC_noise = as.numeric(
        quantile(
          window_noise[, 1],
          0.99,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      
      max_ATAC_noise =
        max(window_noise[, 1], na.rm = TRUE),
      
      median_RNA_noise =
        median(window_noise[, 2], na.rm = TRUE),
      
      q25_RNA_noise = as.numeric(
        quantile(
          window_noise[, 2],
          0.25,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      
      q75_RNA_noise = as.numeric(
        quantile(
          window_noise[, 2],
          0.75,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      
      q95_RNA_noise = as.numeric(
        quantile(
          window_noise[, 2],
          0.95,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      
      q99_RNA_noise = as.numeric(
        quantile(
          window_noise[, 2],
          0.99,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      
      max_RNA_noise =
        max(window_noise[, 2], na.rm = TRUE)
    )
  )
}


Circuits_linkage_posterior <- list(
    TF_Peak_Binding_prob =
        Candidate_TF_Peak_Binding_prob,

    Peak_Gene_Looping_prob =
        Candidate_Peak_Gene_Looping_prob,

    TF_Peak_Binding_switch_rate =
        B_state_switch_rate,

    Peak_Gene_Looping_switch_rate =
        L_state_switch_rate,

    Noise_parameters =
        Noise_parameters[
            seq_len(i),
            ,
            drop = FALSE
        ],

    convergence_log = convergence_log,

    burnin = burnin,

    n_iterations_completed =  i,

    n_post_burnin_iterations = n_kept,

    reached_stability_criterion = 
      stable_check_count >= required_stable_checks,
    window_log = window_log,
    window_size = window_size,
    diagnostic_schedule = list(
      n_postburnin_planned =
        n_postburnin_planned,
      
      planned_checkpoints =
        planned_checkpoints,
      
      planned_complete_windows =
        planned_complete_windows
    ),
    
    burnin =
      burnin,
    
    n_iterations_completed =
      i,
    
    n_post_burnin_iterations =
      n_kept,
    
    reached_stability_criterion =
      stable_check_count >=
      required_stable_checks
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
