# scATACseq pipeline
library(Signac)
library(Seurat)
library(GenomicRanges)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(AnnotationHub)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)

### ----- Processing started ----- ###
import_atac <- function(peak_file,frag_file,sc_file){
  counts <- Read10X_h5(filename = peak_file)
  metadata <- read.csv(
    file = sc_file,
    header = TRUE,
    row.names = 1
  )
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = frag_file,
    min.cells = 10,
    min.features = 200
  )
  seurat_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  

  # Annotation
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- 'UCSC'

  Annotation(seurat_obj) <- annotations
  

  # QC
  seurat_obj <- NucleosomeSignal(object = seurat_obj) #fragment ratio 147-294: <147
  seurat_obj <- TSSEnrichment(object = seurat_obj, fast = FALSE)
  seurat_obj$blacklist_ratio <- seurat_obj$blacklist_region_fragments / seurat_obj$peak_region_fragments
  
  # Add fraction of reads in peaks
  seurat_obj$pct_reads_in_peaks <- seurat_obj$peak_region_fragments / seurat_obj$passed_filters * 100 
  
  # Calculate Boundery
  low_prf <- quantile(seurat_obj[["peak_region_fragments"]]$peak_region_fragments, probs = 0.02)
  hig_prf <- quantile(seurat_obj[["peak_region_fragments"]]$peak_region_fragments, probs = 0.98)
  low_prp <- quantile(seurat_obj[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
  
  high_blr <- quantile(seurat_obj[["blacklist_ratio"]]$blacklist_ratio, probs = 0.98)
  
  hig_ns <- quantile(seurat_obj[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
  
  low_ts <- quantile(seurat_obj[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)
  
  # Subset the data
  seurat_obj <- subset(
    x = seurat_obj,
    subset = peak_region_fragments > low_prf &
      peak_region_fragments < hig_prf &
      pct_reads_in_peaks > low_prp &
      blacklist_ratio < high_blr &
      nucleosome_signal < hig_ns &
      TSS.enrichment > low_ts
  )
  
  return(seurat_obj)
}


# Create Seurat objects for Experiment and Control group data
young_mice<-import_atac("./mice_data/GSM5723631_Young_HSC_filtered_peak_bc_matrix.h5","./mice_data/GSM5723631_Young_HSC_fragments.tsv.gz","./mice_data/GSM5723631_Young_HSC_singlecell.csv.gz")
aged_mice<-import_atac("./mice_data/GSM5723632_Aged_HSC_filtered_peak_bc_matrix.h5","./mice_data/GSM5723632_Aged_HSC_fragments.tsv.gz","./mice_data/GSM5723632_Aged_HSC_singlecell.csv.gz")

# Label sample
young_mice@meta.data$dataset <- "young"
aged_mice@meta.data$dataset <- "aged"
data <- merge(young_mice,aged_mice)


# Normalization, dimension reduction
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q0')
data <- RunSVD(data)

data <- RunUMAP(object = data, reduction = 'lsi', dims = 2:30)
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3, resolution = .4)


# Analysis
gene.activities <- GeneActivity(data)
data[['RNA']] <- CreateAssayObject(counts = gene.activities)

data <- NormalizeData(
  object = data,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)

data[['RNA']]

DefaultAssay(data) <- 'RNA'

FeaturePlot(
  object = data,
  features = c('Kit', 'Pecam1', 'Itgam'),
  max.cutoff = 'q95'
)

DefaultAssay(data) <- 'peaks'

# saveRDS(data, file = "processed_data.rds")


##### ------ Making Scaled Data for testing  ------ #######

downsample_seurat <- function(seurat_obj, fraction = 1/3) {
  set.seed(123)
  cells <- Cells(seurat_obj)
  total_cells <- length(cells)

  if (total_cells == 0) {
    stop("Error: No cells found in the Seurat object.")
  }

  sampled_cells <- sample(cells, size = round(total_cells * fraction), replace = FALSE)

  seurat_obj <- subset(seurat_obj, cells = sampled_cells)

  return(seurat_obj)
}

filter_low_variance_features <- function(seurat_obj, nfeatures = 2000) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Error: Input is not a Seurat object.")
  }

  DefaultAssay(seurat_obj) <- "peaks"

  total_features <- nrow(seurat_obj[["peaks"]]@counts)

  if (total_features < nfeatures) {
    warning("Fewer than ", nfeatures, " features available; keeping all features.")
    var_features <- rownames(seurat_obj[["peaks"]])
  } else {
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
    var_features <- VariableFeatures(seurat_obj)
  }

  seurat_obj <- subset(seurat_obj, features = var_features)

  return(seurat_obj)
}

save_reduced_seurat <- function(seurat_obj, output_path = "small_seurat.rds") {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Error: The object to save is not a Seurat object.")
  }

  saveRDS(seurat_obj, file = output_path)
  print(paste("Saved reduced Seurat object to:", output_path))
}

seurat_obj <- readRDS("processed_data.rds")

seurat_obj <- downsample_seurat(seurat_obj, fraction = 1/3)

seurat_obj <- filter_low_variance_features(seurat_obj, nfeatures = 2000)

save_reduced_seurat(seurat_obj, "small_seurat.rds")

##### ------ Making Scaled Data for testing  ------ #######


### ----- Processing finished ----- ###