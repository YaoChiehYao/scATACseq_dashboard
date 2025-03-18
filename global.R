library(gdata)
library(tools)
library(Seurat)
library(dplyr)
library(shinyjs)
library(shiny)
library(DT)
library(shinydashboard)
library(shinydashboardPlus)
library(ggplot2)
library(shinybusy)
library(glue)
library(markdown)
library(ggthemes)
library(bslib)


load_seurat_obj<- function(path) {
  errors <- c()
  # Check if the ext is rds
  if (tolower(tools::file_ext(path)) != "rds") {
    errors <- c(errors, "Invalid file extension")
    return(errors)
  }
  # try to read file
  obj <- tryCatch({
    readRDS(path)
  }, error = function(e) {
    errors <- c(errors, "Invalid RDS file")
    return(errors)
  })
  # Check if it is seurat object
  if (!inherits(obj, "Seurat")) {
    errors <- c(errors, "Input is not a Seurat object")
    return(errors)
  }
  return(obj)
}

create_da_table <- function(data,region) {
  da_peaks <- FindMarkers(
    object = data,
    ident.1 = rownames(data[[]][data$dataset == "aged",]),
    ident.2 = rownames(data[[]][data$dataset == "young",]),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'peak_region_fragments'
  )
  closest_features <- ClosestFeature(data, regions = rownames(da_peaks))
  da_peaks$closet_gene <- ClosestFeature(data, regions = rownames(da_peaks))$gene_name
  da_peaks$distance <- ClosestFeature(data, regions = rownames(da_peaks))$distance
  return(da_peaks)
}


create_accessibility_plot <- function(obj, region) {
  req(region)  # Ensure the region is selected
  
  # Extract Data from Seurat object
  plot_data <- obj@assays$peaks@data[region, , drop = FALSE]
  
  # Check if region exist
  if (nrow(plot_data) == 0) {
    print(paste("Warning: Selected region", region, "not found in dataset."))
    return(NULL)
  }
  
  # Check for NA and 0 variant
  if (all(is.na(plot_data)) || length(unique(as.numeric(plot_data))) == 1) {
    print(paste("Warning: No variation in selected region", region))
    return(NULL)
  }
  
  print(paste("Plotting Violin Plot for:", region))
  
  plot1 <- VlnPlot(
    object = obj,
    features = region,
    group.by = "dataset"
  )
  
  plot2 <- FeaturePlot(
    object = obj,
    features = region,
    max.cutoff = "q95"
  )
  
  return(plot1 | plot2)
}

create_coverage_plot <- function(obj, region){
  coverage_plot<-CoveragePlot(
    object = obj,
    region = region,
    extend.upstream = 10000,
    extend.downstream = 5000,
    group.by = "dataset"
  )
  return(coverage_plot)
}
