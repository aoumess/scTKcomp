## Functions contributing to filling the gaps in the singleCellTK modules

## Launch singleCellTK
scTK_launch <- function() {
  library(Matrix)
  ## Force using web browser rather than Rstudio
  options('Matrix.warnDeprecatedCoerce' = 0, "shiny.launch.browser" = .rs.invokeShinyWindowExternal)
  message(options("Matrix.warnDeprecatedCoerce"))
  singleCellTK::singleCellTK()
  options('Matrix.warnDeprecatedCoerce' = 0)
  message(options("Matrix.warnDeprecatedCoerce"))
}

## Loading data into a Seurat object
# 1) Loading data
# 2) Filtering duplicated cell barcodes
# 3) Rename ensembl genes id by genes symbols
# 4) Remove empty droplets
# 5) Plotting saturation and Kneeplot
# 6) Creation of the Seurat object

## Loading SC data to SCE object
scTK_load <- function(data_path = NULL, sample_name = NULL, assay = 'RNA', out_dir = getwd(), return_data = FALSE) {
  
  ## Checks
  if (is.null(sample_name)) stop('A sample name is required !')
  sample_name <- as.character(sample_name)
  if (!dir.exists(data_path)) stop('Input data path does not exist !')
  if (!dir.exists(out_dir)) stop('Output directory does not exist !')
  if (!is.character(assay)) stop('Assay name should be a character (string)')
  if (!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  
  message("Loading data ...")
  
  ## Loading data
  source.format <- ''
  if(file.exists(paste0(data_path, "/matrix.mtx")) | file.exists(paste0(data_path, "/matrix.mtx.gz"))) {
    ### Cell Ranger
    source.format <- "CellRanger"
    scmat <- Seurat::Read10X(data_path)
    if ('Gene Expression' %in% names (scmat)) {
      message("Keeping gene expression only")
      scmat <- scmat[['Gene Expression']]
    }
  } else if(file.exists(paste0(data_path, "/", sample_name, ".mtx"))) {
    ### BUStools
    source.format <- "BUStools"
    library(Matrix)
    options("Matrix.warnDeprecatedCoerce" = 1)
    scmat <- BUSpaRse::read_count_output(dir = data_path, name = sample_name, tcc = FALSE)
  } else if (file.exists(paste0(data_path, "/quants_mat.gz"))) {
    ### Alevin
    source.format <- "Alevin"
    scmat <- Seurat::ReadAlevin(data_path)
  } else if (file.exists(paste0(data_path, "/", sample_name, "_counts.tsv.gz"))) {
    ### UMI-tools
    source.format <- "UMIt-ools"
    scmat <- read.table(file = paste0(data_path, "/", sample_name, "_counts.tsv.gz"), header = TRUE, sep = "\t", quote = "", check.names = FALSE, row.names = 1)
  } else {
    stop(paste0("No data found in [", data_path, "] (wrong path ?), or unsupported format ..."))
  }
  message(paste0("Found ", source.format, " data"))
  
  scmat <- scmat[,order(colnames(scmat))]
  
  message('Full droplets matrix dimensions :')
  droplets.nb <- ncol(scmat)
  print(dim(scmat))
  
  message('Total UMIs :')
  umi.total.nb <- sum(scmat)
  print(umi.total.nb)
  
  ## Computing some metrics
  numi_drop <- Matrix::colSums(x = scmat)
  scmat.bin <- scmat
  attr(scmat.bin, "x")[attr(scmat.bin, "x") >= 1] <- 1
  scmat.bin <- Matrix::drop0(scmat.bin)
  ngen_drop <- Matrix::colSums(x = scmat.bin)
  rm(scmat.bin)
  
  ## Creation of the SCE object
  sobj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=scmat), mainExpName = assay)
  rm(scmat)
  ## Adding metrics
  sobj$sample <- as.factor(sample_name)
  sobj[[paste0('nCount_', assay)]] <- numi_drop
  sobj[[paste0('nFeature_', assay)]] <- ngen_drop
  sobj[[paste0('log_nCount_', assay)]] <- log(sobj[[paste0('nCount_', assay)]]+1)
  sobj[[paste0('log_nFeature_', assay)]] <- log(sobj[[paste0('nFeature_', assay)]]+1)
  
  ## Fill misc metadata : metrics
  mymeta <- list(misc=list())
  mymeta$misc$params$Droplet_Quality$captured_droplet <- droplets.nb
  mymeta$misc$params$Droplet_Quality$total_number_UMI <- umi.total.nb
  mymeta$misc$samplename <- sample_name
  sobj@metadata <- mymeta
  
  ## File output
  saveRDS(object = sobj, file = paste0(out_dir, '/', sample_name, '_raw_SCE.RDS'), compress = 'bzip2')
  
  ## Return data ?
  if (return_data) return(sobj)
}

## Example
# scTK_load(data_path = "/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix", sample_name = "PBMC1K10X", out_dir = "/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/")


## Filtering of Empty droplets
scTK_edf <- function(sobj = NULL, assay = 'RNA', droplets_min = 1E+04, emptydrops_fdr = 1E-03, emptydrops_retain = NULL, my_seed = 1337L, out_dir = getwd(), draw_plots = TRUE, return_data = FALSE) {
  
  ## Checks
  if (is.null(sobj)) stop('the sobj parameter should not be NULL !')
  if (!dir.exists(out_dir)) stop('Output directory does not exist !')
  if (!is.character(assay)) stop('Assay name should be a character (string)')
  if (!is.logical(draw_plots)) stop('The return_data parameter should be a logical (boolean)')
  if (!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  
  ## Checking sobj format
  in.format <- is(sobj)[1]
  if (tolower(in.format) != 'singlecellexperiment') stop ("Input object is in unknown format ! Should be a 'SingleCellExperiment' object.")
  
  message('Full droplets matrix dimensions :')
  droplets.nb <- ncol(sobj@assays@data@listData$counts)
  print(dim(sobj))
  
  message('Total UMIs :')
  umi.total.nb <- sum(sobj@assays@data@listData$counts)
  print(umi.total.nb)
  
  if (!is.null(droplets_min) && ncol(sobj@assays@data@listData$counts) > droplets_min && !is.null(emptydrops_fdr)) {
    ## Removing empty droplets
    message("Removing empty droplets with emptyDrops ...")
    bc_rank <- DropletUtils::barcodeRanks(sobj@assays@data@listData$counts)
    set.seed(my_seed)
    bc_rank2 <- DropletUtils::emptyDrops(sobj@assays@data@listData$counts, retain = emptydrops_retain)
    keep.bc <- bc_rank2$FDR < emptydrops_fdr
    rm(bc_rank2)
    keep.bc[is.na(keep.bc)] <- FALSE
    ## Assessing if some barcodes are kept
    if (!any(keep.bc, na.rm = TRUE)) {
      ### Nope ...
      message("WARNING : emptyDrops estimated all droplets as empty! emptyDrops won't get used ...")
      emptydrops_fdr <- NULL
      umi.kept.nb <- umi.total.nb
    } else {
      ### Yerp ! Filtering empty droplets
      # scmat <- sobj@assays@data@listData$counts[, which(keep.bc)]
      message('Filtered droplets matrix dimensions :')
      print(dim(sobj@assays@data@listData$counts[,keep.bc]))
      message('Total UMIs (filtered) :')
      umi.kept.nb <- sum(sobj@assays@data@listData$counts[,keep.bc])
      print(umi.kept.nb)
      message('Fraction of UMIs in cells :')
      print(umi.kept.nb / umi.total.nb)
      
      ## Plotting when requested
      if(draw_plots){
        ### Kneeplot
        png(filename = paste0(out_dir, '/', sobj@metadata$misc$samplename, "_kneeplot.png"), width = 1000, height = 700)
        plot(bc_rank$rank, bc_rank$total+1, log = "xy", xlab = "Rank", ylab = "Total", col = ifelse(keep.bc, "red", "black"), pch = 20, cex = ifelse(keep.bc, 1, 2), main = paste0(sobj@metadata$misc$samplename, ' kneeplot (', length(which(keep.bc)), ' barcodes kept)'))
        o <- order(bc_rank$rank)
        lines(bc_rank$rank[o], bc_rank$fitted[o], col = "red")
        abline(h = S4Vectors::metadata(bc_rank)$knee, col = "dodgerblue", lty = 2)
        abline(h = S4Vectors::metadata(bc_rank)$inflection, col = "forestgreen", lty = 2)
        legend("bottomleft", lty = c(2, 2, NA, NA), pch = c(NA, NA, 20, 20), col = c("dodgerblue", "forestgreen", "red", "black"), legend = c("knee", "inflection", "cell", "empty"))
        dev.off()
        rm(bc_rank)
        ### Saturation plot
        png(filename = paste0(out_dir, '/', sobj@metadata$misc$samplename, "_satplot.png"), width = 1000, height = 700)
        suppressWarnings(plot(sobj[[paste0('nCount_', assay)]], sobj[[paste0('nFeature_', assay)]], pch = 20, log = "xy", col = ifelse(keep.bc, "red", "black"), xlab = 'Nb of UMIs in droplet (log)', ylab = 'Nb of genes with at least 1 UMI count in droplets (log)', main = paste0(sample_name, ' saturation plot (', length(which(keep.bc)), ' barcodes kept)')))
        legend("topleft", pch = c(20, 20), col = c("red", "black"), legend = c("cell", "empty"))
        dev.off()
        ### Saturation distrib plot
        # png(filename = paste0(out_dir, sample_name, "_satdistplot.png"), width = 1000, height = 700)
        # umipgen <- numi_drop / ngen_drop
        # upg.range <- range(umipgen)
        # plot(density(umipgen[keep.bc], col = "red", xlim = c(0, upg.range[2], )))
      }
      
      ## Applying filter
      sobj <- sobj[,keep.bc]
    }
  }
  
  ## Fill misc metadata : metrics
  sobj@metadata$misc$params$Droplet_Quality$captured_droplet <- droplets.nb
  sobj@metadata$misc$params$Droplet_Quality$total_number_UMI <- umi.total.nb
  sobj@metadata$misc$params$Droplet_Quality$estimated_cells <- ncol(sobj)
  sobj@metadata$misc$params$Droplet_Quality$estimated_UMI <- umi.kept.nb
  sobj@metadata$misc$params$Droplet_Quality$fraction_read_in_cells <- umi.kept.nb / umi.total.nb
  
  ## Fill misc metadata : tool parameters
  sobj@metadata$misc$params$sobj_creation$emptydrops_fdr <- emptydrops_fdr
  sobj@metadata$misc$params$sobj_creation$droplets_min <- droplets_min
  sobj@metadata$misc$params$sobj_creation$emptydrops_fdr <- emptydrops_fdr
  sobj@metadata$misc$params$sobj_creation$emptydrops_retain <- emptydrops_retain
  sobj@metadata$misc$params$seed <- my_seed
  
  ## File output
  saveRDS(object = sobj, file = paste0(out_dir, '/', sobj@metadata$misc$samplename, '_EDf.RDS'), compress = 'bzip2')
  
  ## Return data ?
  if (return_data) return(sobj)
}

## Example
# scTK_edf(sobj = readRDS("/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/PBMC1K10X_raw.RDS"), out_dir = "/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/")


## Estimate cell cycle : TRICYCLE
scTK_cc_tricycle <- function(sobj = NULL, assay = 'RNA', species = 'human', feature_type = 'SYMBOL', out_dir = getwd(), return_data = FALSE) {
  
  ## Checks
  if (is.null(sobj)) stop('the sobj parameter should not be NULL !')
  if (!dir.exists(out_dir)) stop('Output directory does not exist !')
  if (!is.character(assay)) stop('Assay name should be a character (string)')
  if (!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  if (!tolower(species) %in% c('human', 'mouse')) stop('Tricycle is only compatible with human and mouse.')
  if(!tolower(feature_type) %in% c('ensembl', 'symbol')) stop("Feature type should be  'ENSEMBL' or 'SYMBOL'")
  
  ## Checking sobj format
  in.format <- is(sobj)[1]
  if (tolower(in.format) != 'singlecellexperiment') stop ("Input object is in unknown format ! Should be a 'SingleCellExperiment' object.")
  
  ## Converting to Seurat obj
  sobj_seu <- Seurat::CreateSeuratObject(counts = sobj@assays@data@listData$counts, project = sobj@metadata$misc$samplename, assay = assay, meta.data = as.data.frame(sobj@colData@listData))
  ## Quick log-normalization
  sobj_seu <- Seurat::NormalizeData(object = sobj_seu, normalization_method = 'LogNormalize')
  ## Running Tricycle
  sobj_seu <- SeuratWrappers::Runtricycle(object = sobj_seu, assay = assay, gname.type = 'SYMBOL', species = 'human', reduction.name = 'TRIC', reduction.key = 'TRIC_')
  sobj$cc_tricycle <- sobj_seu$tricyclePosition
  rm(sobj_seu)
  
  ## File output
  # saveRDS(object = sobj, file = paste0(out_dir, '/', sobj@metadata$misc$samplename, '_CC.tricycle.RDS'), compress = 'bzip2')
  saveRDS(object = sobj, file = paste0(out_dir, '/', paste(c(sobj@metadata$misc$samplename, SingleCellExperiment::mainExpName(sobj), assay, 'CC.tricycle.RDS'), collapse = '_')), compress = 'bzip2')
  
  ## Return data ?
  if (return_data) return(sobj)
}

## Example
# scTK_cc_tricycle(sobj = readRDS("/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/PBMC1K10X_edf.RDS"), out_dir = "/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/", assay = 'RNA', species = 'human', feature_type = 'SYMBOL')


## Estimate cell cycle : SCHWABE
scTK_cc_schwabe <- function(sobj = NULL, assay = 'RNA', species = 'human', feature_type = 'SYMBOL', out_dir = getwd(), return_data = FALSE) {
  
  ## Checks
  if (is.null(sobj)) stop('the sobj parameter should not be NULL !')
  if (!dir.exists(out_dir)) stop('Output directory does not exist !')
  if (!is.character(assay)) stop('Assay name should be a character (string)')
  if (!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  if (!tolower(species) %in% c('human', 'mouse')) stop('Tricycle is only compatible with human and mouse.')
  if(!tolower(feature_type) %in% c('ensembl', 'symbol')) stop("Feature type should be  'ENSEMBL' or 'SYMBOL'")
  
  ## Checking sobj format
  in.format <- is(sobj)[1]
  if (tolower(in.format) != 'singlecellexperiment') stop ("Input object is in unknown format ! Should be a 'SingleCellExperiment' object.")
  
  ## Converting to Seurat obj
  sobj_seu <- Seurat::CreateSeuratObject(counts = sobj@assays@data@listData$counts, project = sobj@metadata$misc$samplename, assay = assay, meta.data = as.data.frame(sobj@colData@listData))
  ## Quick log-normalization
  sobj_seu <- Seurat::NormalizeData(object = sobj_seu, normalization_method = 'LogNormalize')
  ## Running Schwabe method
  sobj$cc_schwabe <- suppressMessages(tricycle::estimate_Schwabe_stage(sobj_seu@assays$RNA@data, gname.type = toupper(feature_type), species = tolower(species)))
  message("Cell cycle phases according to Schwabe method : ")
  print(table(sobj$cc_schwabe, useNA = "always"))
  
  ## File output
  # saveRDS(object = sobj, file = paste0(out_dir, '/', sobj@metadata$misc$samplename, '_CC.schwabe.RDS'), compress = 'bzip2')
  saveRDS(object = sobj, file = paste0(out_dir, '/', paste(c(sobj@metadata$misc$samplename, SingleCellExperiment::mainExpName(sobj), assay, 'CC.schwabe.RDS'), collapse = '_')), compress = 'bzip2')
  ## Return data ?
  if (return_data) return(sobj)
}

## Example
# scTK_cc_schwabe(sobj = readRDS("/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/PBMC1K10X_cc.tricycle.RDS"), out_dir = "/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/", assay = 'RNA', species = 'human', feature_type = 'SYMBOL')


## Estimate cell cycle : Seurat
scTK_cc_seurat <- function(sobj = NULL, exp.name = NULL, assay = 'counts', cc_seurat_file = NULL, nbin = 24, out_dir = getwd(), return_data = FALSE) {
  
  ## Checks
  message('Checks ...')
  if (is.null(sobj)) stop('The sobj parameter should not be NULL !')
  if (is.null(cc_seurat_file)) stop('The RDS containing Seurat gene lists is required !')
  if (!dir.exists(out_dir)) stop('Output directory does not exist !')
  if (!is.character(assay)) stop('Assay name should be a character (string)')
  if (!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  
  ## Checking sobj format
  in.format <- is(sobj)[1]
  if (tolower(in.format) != 'singlecellexperiment') stop ("Input object is in unknown format ! Should be a 'SingleCellExperiment' object.")
  
  ## Loading cc file
  message('Loading Seurat cell cyle genes ...')
  cc_data <- readRDS(cc_seurat_file)
  
  ## Converting to Seurat obj
  message('Converting to Seurat object ...')
  if (exp.name == SingleCellExperiment::mainExpName(sobj)) {
    seu.counts <- SummarizedExperiment::assay(x = sobj, i = assay)
    # seu.meta <- SummarizedExperiment::colData(x = sobj)
  } else if (exp.name %in% SingleCellExperiment::altExpNames(sobj)) {
    alt.exp <- SingleCellExperiment::altExp(x = sobj, e = exp.name)
    seu.counts <- SummarizedExperiment::assay(x = alt.exp, i = assay)
    # seu.meta <- SummarizedExperiment::colData(x = alt.exp)
  } else stop('Unknown Experiment name and/or Assay name !')
  
  sobj_seu <- Seurat::CreateSeuratObject(counts = seu.counts, project = sobj@metadata$misc$samplename, assay = exp.name, meta.data = as.data.frame(sobj@colData@listData))
  
  ## Quick log-normalization
  sobj_seu <- Seurat::NormalizeData(object = sobj_seu, normalization_method = 'LogNormalize')
  
  ### Seurat
  message('Cell cycle scoring with Seurat ...')
  sobj_seu <- Seurat::CellCycleScoring(object = sobj_seu, s.features = cc_data$s.genes, g2m.features = cc_data$g2m.genes, assay = exp.name, nbin = nbin, seed = sobj@metadata$misc$params$seed)
  sobj$cc_seurat.S.Score <- sobj_seu$S.Score
  sobj$cc_seurat.G2M.Score <- sobj_seu$G2M.Score
  sobj$cc_seurat.SmG2M.Score <- sobj_seu$S.Score - sobj_seu$G2M.Score
  sobj$cc_seurat.Phase <- sobj_seu$Phase
  rm(sobj_seu)
  message("Cell cycle phases according to Seurat: ")
  print(table(sobj$cc_seurat.Phase, useNA = 'always'))
  
  ## File output
  # saveRDS(object = sobj, file = paste0(out_dir, '/', sobj@metadata$misc$samplename, '_CC.seurat.RDS'), compress = 'bzip2')
  saveRDS(object = sobj, file = paste0(out_dir, '/', paste(c(sobj@metadata$misc$samplename, SingleCellExperiment::mainExpName(sobj), assay, 'CC.seurat.RDS'), collapse = '_')), compress = 'bzip2')
  ## Return data ?
  if (return_data) return(sobj)
}

## Example
# scTK_cc_seurat(sobj = readRDS("/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/PBMC1K10X_cc.schwabe.RDS"), out_dir = "/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/", assay = 'counts', cc_seurat_file = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/RESOURCES/GENELISTS/homo_sapiens_Seurat_cc.genes_20191031.rds')


## Estimate cell cycle : Cyclone (scran)
scTK_cc_cyclone <- function(sobj = NULL, assay = 'counts', cc_cyclone_file = NULL, out_dir = getwd(), return_data = FALSE) {
  
  ## Checks
  if (is.null(sobj)) stop('The sobj parameter should not be NULL !')
  if (is.null(cc_cyclone_file)) stop('The RDS containing Cyclone gene pairs is required !')
  if (!dir.exists(out_dir)) stop('Output directory does not exist !')
  if (!is.character(assay)) stop('Assay name should be a character (string)')
  if (!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  if(!assay %in% names(SummarizedExperiment::assays(sobj))) stop('Provided assay not found in the provided SCE object.')
  
  ## Checking sobj format
  in.format <- is(sobj)[1]
  if (tolower(in.format) != 'singlecellexperiment') stop ("Input object is in unknown format ! Should be a 'SingleCellExperiment' object.")
  
  ## Loading cc file
  cc_data <- readRDS(cc_cyclone_file)
  
  ## Cyclone
  set.seed(sobj@metadata$misc$params$seed)
  # doParallel::registerDoParallel(nthreads)
  # BPPARAM <- BiocParallel::SerialParam()
  cycres <- scran::cyclone(sobj, pairs = cc_data, assay.type = assay, verbose = FALSE)
  sobj$cc_cyclone.Phase <- as.factor(cycres$phases)
  sobj$cc_cyclone.G1.Score <- cycres$normalized.scores$G1
  sobj$cc_cyclone.S.Score <- cycres$normalized.scores$S
  sobj$cc_cyclone.G2M.Score <- cycres$normalized.scores$G2M
  sobj$cc_cyclone.SmG2M.Score <- cycres$normalized.scores$S - cycres$normalized.scores$G2M
  message("Cell cycle phases according to Cyclone: ")
  print(table(sobj$cc_cyclone.Phase, useNA = 'always'))
  
  ## File output
  saveRDS(object = sobj, file = paste0(out_dir, '/', paste(c(sobj@metadata$misc$samplename, SingleCellExperiment::mainExpName(sobj), assay, 'CC.cyclone.RDS'), collapse = '_')), compress = 'bzip2')
  
  ## Return data ?
  if (return_data) return(sobj)
}

## Example
# scTK_cc_cyclone(sobj = readRDS("/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/PBMC1K10X_cc.seurat.RDS"), out_dir = "/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/", assay = 'RNA', cc_cyclone_file = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/RESOURCES/GENELISTS/homo_sapiens_cyclone_pairs_symbols_20191001.rds')


## Function to assess the weight of annotation covariates in a (sample x annotations) data.frame, on a (feature x sample) data matrix, through correlation (for a continuous covariate) or Kruskal-Wallis statistic (for factors)
## . sobj               [SCE object]            An SCE object with at least one assay
## . assay              [char]                  Name of an assay present in sobj
## . factor_names       [vec(char)]             Column name(s) from the colData() of sobj, corresponding to factor covariate(s)
## . conti.colnames     [vec(char)]             CColumn name(s) from the colData() of sobj, corresponding to contiunous covariates
## . ctrl_features      [vec(char)]             Name of feature(s) of sobj, to consider as control gene(s)
## . marker_features    [vec(char)]             Name of feature(s) of sobj, to consider as marker gene(s)
## . red_method         [char]                  Dimension reduction method ['PCA', 'MDS.euc', 'MDS.spear']
## . ndim_max           [int>0]                 Number of dimensions to compute and plot
## . center             [bool]                  Center the matrix ?
## . scale              [bool]                  Scale the matrix ?
## . coef_cut           [0<=num<1]              Do not display (as colors) coefficients inferior to this value on the heatmap
## . color_palette      [vec(col)]              Color vector (length 2 : start,end) for the heatmap
## . out_dir            [char]                  Path to the output image
scTK_assess_covar <- function(sobj = NULL, assay = 'counts', factor_names = NULL, conti_names = NULL, ctrl_features, marker_features, red_method = 'pca', ndim_max = 10, center = TRUE, scale = TRUE, coef_cut = 0, color_palette = c("white", "orangered3"), out_dir = getwd()) {
  ## Checks
  ### Mandatory
  if (!dir.exists(out_dir)) stop('Output directory does not exist !')
  if (!is.character(assay)) stop('Assay name should be a character (string)')
  if (is.null(sobj)) stop('The sobj parameter should not be NULL !') else if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if (all(is.null(c(factor_names, conti_names, ctrl_features, marker_features)))) stop('At least one of [factor_names], [conti.colnames], [ctrl_features] or [marker_features] should not be NULL.')
  if (ndim_max <= 0) stop('[ndim_max] should be a non-null positive integer (and <= s samples).')
  if (!tolower(red_method) %in% c('pca', 'mds.euc', 'mds.spear')) stop('Unknown or unsupported dimension reduction method')
  
  ## Loading data 
  mat <- SummarizedExperiment::assay(x = sobj, i = assay)
  annot.df <- as.data.frame(sobj@colData)
  
  ### Compatibility
  #### Number of requested dimensions
  if (ndim_max > ncol(mat)) {
    message('\tWARNING : requested [ndim_max] is higher than columns in [mat]. Reducing it to [mat] samples.')
    ndim_max <- ncol(mat)-1
  }
  #### Existing covariates
  if (nrow(annot.df) != ncol(mat)) stop('There should be the same number of samples in [mat] (columns) and [annot.df] (rows)')
  if(!is.null(factor_names)) {
    if(!all(factor_names %in% colnames(annot.df))) stop('All [factor_names] should be in colnames of [annot.df].')
  }
  if(!is.null(conti_names)) {
    if(!all(conti_names %in% colnames(annot.df))) stop('All [conti_names] should be in colnames of [annot.df].')
  }
  #### Filtering features
  cf.check <- ctrl_features %in% rownames(sobj)
  mf.check <- marker_features %in% rownames(sobj)
  if(!all(cf.check)) {
    message('\tWARNING : [', paste(ctrl_features[!cf.check], collapse = ', '), '] CONTROL feature(s) not found in sobj. It/They will be discarded.')
    ctrl_features <- ctrl_features[cf.check]
    if(length(ctrl_features) == 0) message('\t\tWARNING : No CONTROL feature left !')
  }
  if(!all(mf.check)) {
    message('\tWARNING : [', paste(marker_features[!mf.check], collapse = ', '), '] MARKER feature(s) not found in sobj. It/They will be discarded.')
    marker_features <- marker_features[mf.check]
    if(length(marker_features) == 0) message('\t\tWARNING : No MARKER feature left !')
  }
  if (length(c(ctrl_features, marker_features)) == 0) message('WARNING : There is no CONTROL nor MARKER feature to assess !')
  
  ## RUN
  
  ## Get experiment name
  exp.name <- SingleCellExperiment::mainExpName(sobj)
  
  ## Center / scale ?
  if (any(c(center, scale))) {
    message('Centering and/or scaling ...')
    mat <- base::scale(x = mat, center = center, scale = scale)
  }
  ## Dimension reduction
  message('Performing dimension reduction ...')
  # if (tolower(red_method) == 'pca') norm.red <- base::svd(x = mat, nv = ndim_max)$v
  if (tolower(red_method) == 'pca') norm.red <- irlba::prcomp_irlba(x = mat, n = ndim_max, center = FALSE, scale. = FALSE)$rotation
  if (tolower(red_method) == 'mds.euc') norm.red <- stats::cmdscale(d = dist(x = t(mat), method = 'euclidean'), k = ndim_max)
  if (tolower(red_method) == 'mds.spear') norm.red <- stats::cmdscale(d = as.dist(1-cor(mat, method = 'spearman')), k = ndim_max)
  ## Prepping covariates
  covar.names <- c(factor_names, conti_names)
  covar.types <- c(rep('C.factor', length(factor_names)), rep('C.continuous', length(conti_names)))
  ## Prepping features
  feat.names <- c(ctrl_features, marker_features)
  feat.types <- c(rep('F.control', length(ctrl_features)), rep('F.marker', length(marker_features)))
  ## All
  cf.names <- c(covar.names, feat.names)
  cf.types <- c(covar.types, feat.types)
  ## Setting output matrix
  bc.mat <- matrix(NA, ncol = length(cf.names), nrow = ndim_max, dimnames = list(paste0(toupper(red_method), seq_len(ndim_max)), cf.names))
  ## Filling matrix
  ### Covariates ...
  if(length(covar.names) > 0) {
    message('Assessing covariate(s) ...')
    for (cn in seq_along(covar.names)) {
      message("\t ...", covar.names[cn])
      if (covar.names[cn] %in% conti_names) {
        cv2cor <- annot.df[[covar.names[cn]]]
        nona <- !is.na(cv2cor)
        bc.mat[, cn] <-  abs(cor(x = cv2cor[nona], y = norm.red[nona,], method = 'spearman'))
      } else if (covar.names[cn] %in% factor_names & length(unique(annot.df[[covar.names[cn]]])) > 1) {
        b2kw <- annot.df[[covar.names[cn]]]
        nona <- !is.na(b2kw)
        for (si in seq_len(ndim_max)) {
          bc.mat[si,cn] <- kruskal.test(x = norm.red[nona,si], g = as.factor(b2kw[nona]))$statistic / nrow(norm.red)
        }
      }
    }
  }
  ### Covariates ...
  if(length(feat.names) > 0) {
    message('Assessing features ...')
    for (fn in seq_along(feat.names)) {
      message("\t ...", feat.names[fn])
      cv2cor <- mat[feat.names[fn],]
      nona <- !is.na(cv2cor)
      bc.mat[, fn + length(covar.names)] <-  abs(cor(x = cv2cor[nona], y = norm.red[nona,], method = 'spearman'))
    }
  }
  ## Cutting values if requested
  bc.mat[bc.mat < coef_cut] <- 0
  ## Heatmap
  myRamp.col <- circlize::colorRamp2(c(0, 1), color_palette)
  BC.hm <- ComplexHeatmap::Heatmap(matrix = bc.mat,
                                   name = 'Weight',
                                   col = myRamp.col,
                                   na_col = 'grey75',
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   rect_gp = grid::gpar(col = "darkgrey", lwd=0.5),
                                   column_title = 'Batch factors and covariates weight on dataset',
                                   row_title = 'Dimensions',
                                   column_split = cf.types,
                                   top_annotation = ComplexHeatmap::HeatmapAnnotation(Type = cf.types, col = list(Type = setNames(object = c('blue', 'lightblue', 'pink', 'red'), nm = c('C.factor', 'C.continuous', 'F.control', 'F.marker')))))
  png(filename = paste0(out_dir, '/', paste(c(sobj@metadata$misc$samplename, exp.name, assay, 'covar.png'), collapse = '_')), width = 400+(50*length(cf.names)), height = 1000)
  ComplexHeatmap::draw(BC.hm)
  dev.off()
}

## EXAMPLE
# scTK_assess_covar(sobj = sobj, assay = 'SoupX', factor_names = c('hidden', 'cc_seurat.Phase', 'cc_cyclone.Phase'), conti_names = c('log_nCount_RNA', 'log_nFeature_RNA', 'soupX_nUMIs', 'subsets_mm.ribo_percent', 'subsets_mm.stress_percent', 'subsets_Mito_percent', 'cc_seurat.SmG2M.Score', 'cc_cyclone.SmG2M.Score'), ctrl_features = c('Gapdh'), marker_features = c('Il2ra', 'Cd8b1', 'Cd8a', 'Cd4', 'Ccr7', 'Itm2a', 'Aif1', 'Hba-a1'), ndim_max = 25, out_dir = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/TEST_DATASET/ANALYSIS/')


## List assay names in an SCE object
scTK_list_assays <- function(sobj = NULL) {
  ## Checks
  if (is.null(sobj)) stop('The sobj parameter should not be NULL !')
  message('Experiment [', SingleCellExperiment::mainExpName(sobj), '] (main)')
  message('\tAssay(s) : ', paste0('[', paste(names(sobj@assays), collapse = '], ['), ']'))
  alt.names <- SingleCellExperiment::altExpNames(sobj)
  for (en in seq_along(alt.names)) {
    message('Experiment [', alt.names[en], '] (alt)')
    message('\tAssay(s) : ', paste0('[', paste(names(SingleCellExperiment::altExp(x = sobj, e = alt.names[en])@assays), collapse = '], ['), ']'))
  }
}

## EXAMPLE
# scTK_list_assays(sobj = sobj)


## Normalization with covariate regression, using Seurat methods (SCTransform / LogNormalize / CLR))
## Expected value for 'normalization_method' : 'SCT', 'LN', 'CLR'
scTK_norm_covar <- function(sobj = NULL, assay = 'counts', normalization_method = 'SCT', features = 3000, covariates = NULL, seed = 12345, out_dir = getwd(), return_data = FALSE) {
  
  valid.norm <- c('SCT', 'LN', 'CLR', 'SCBFA')
  
  ## Checks
  message('Checks ...')
  if (!dir.exists(out_dir)) stop('Output directory does not exist !')
  if (!is.character(assay)) stop('Assay name should be a character (string)')
  if (is.null(sobj)) stop('The sobj parameter should not be NULL !') else if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if (!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  if (!tolower(normalization_method) %in% tolower(valid.norm)) stop('Unknown normalization method !')

  ## Cleaning covariates
  message('Cleaning covariates ...')
  if(!is.null(covariates)) covariates <- sort(unique(covariates))
  
  ## Save command
  sobj@metadata$misc$pipeline_commands=c(sobj@metadata$misc$pipeline_commands, paste0("scTK_norm_covar(sobj = sobj, assay = ", assay, ", normalization_method = ", normalization_method, ", features = ", features, ", covariates = ", if(is.null(covariates)) "NULL" else paste0("c(", paste(covariates, collapse = ","),")"), ")"))
  
  ## Convert to Seurat object
  message('Converting to Seurat object ...')
  sobj_seu <- Seurat::CreateSeuratObject(counts = SummarizedExperiment::assay(x = sobj, i = assay), project = sobj@metadata$misc$samplename, assay = assay, meta.data = as.data.frame(sobj@colData))
  
  ## Backing metadata
  mymeta <- sobj@metadata
  
  ## Normalize
  message('Performing normalization ...')
  
  ## tesT scbfa
  if (tolower(normalization_method) == tolower('SCBFA')) {
    ## Dealing with covariates
    minimeta <- sobj@colData[,colnames(sobj@colData) %in% covariates, drop = FALSE]
    X <- matrix(ncol = 0, nrow = nrow(minimeta))
    for(v in covariates) {
      if(is.character(minimeta[[v]])) minimeta[[v]] <- as.factor(minimeta[[v]])
      if(any(is.na(minimeta[[v]]))) stop(paste0("Covariate '", v, "' contains NA value(s) !"))
      if(is.factor(minimeta[[v]])) {
        message(paste0("Converting '", v, "' factor into model matrix and adding to the regression..."))
        mm <- model.matrix(~minimeta[[v]])[,-1, drop = FALSE]
        # mm <- model.matrix(~minimeta[[v]])
        X <- cbind(X, mm)
      } else {
        message(paste0("Adding '", v, "' covariate", if(vtr.scale) ' (scaled)' else NULL,  ' to the regression ...'))
        X <- cbind(X, if(vtr.scale) scale(minimeta[[v]]) else minimeta[[v]])
      }
    }
    
    system.time(sobj_scbfa <- scBFA::scBFA(scData = sobj, numFactors = max.dims, X = X))
  }
  
  if (tolower(normalization_method) == tolower('SCT')) {
    ## SCT
    sobj_seu <- suppressWarnings(Seurat::SCTransform(object = sobj_seu, assay = assay, new.assay.name = paste0('SCT.', assay), vars.to.regress = covariates, seed.use = seed, variable.features.n = features, return.only.var.genes = FALSE))
  } else if (tolower(normalization_method) == tolower('LN')) {
    ## LogNorm
    sobj <- Seurat::NormalizeData(sobj_seu, normalization_method = 'LogNormalize', assay = assay, )
    sobj <- Seurat::FindVariableFeatures(sobj, assay = assay, nfeatures = features)
    sobj <- Seurat::ScaleData(object = sobj, vars.to.regress = covariates, do.scale = TRUE, scale.max = 10, block.size = 1000)
  } else if (tolower(normalization_method) == tolower('CLR')) {
    ## CLR
    sobj <- Seurat::NormalizeData(sobj, normalization_method = 'CLR', assay = assay)
    sobj <- Seurat::FindVariableFeatures(sobj, assay = assay, nfeatures = features)
    sobj <- Seurat::ScaleData(object = sobj, vars.to.regress = covariates, do.scale = TRUE, scale.max = 10, block.size = 1000)
  } else stop('Unknown or unsupported normalization method !')
  
  ## Converting back to SCE
  sobj <- Seurat::as.SingleCellExperiment(sobj_seu)
  sobj@metadata <- mymeta
  rm(sobj_seu)
  
  ## File output
  saveRDS(object = sobj, file = paste0(out_dir, '/', paste(c(sobj@metadata$misc$samplename, assay, toupper(normalization_method)), collapse = '_'), '.rds'), compress = 'bzip2')
  
  ## Return data ?
  if (return_data) return(sobj)
}

## EXAMPLE
# scTK_norm_covar(sobj = readRDS('/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/TEST_DATASET/ANALYSIS/BACON_CC.seurat.RDS'), covariates = c('subsets_mm.ribo_percent', 'cc_seurat.SmG2M.Score'), assay = 'SoupX', normalization_method = 'SCT', features = 2000, out_dir = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/TEST_DATASET/ANALYSIS/')

