##  FUNCTIONS CONTRIBUTING TO FILLING THE GAPS IN THE singleCellTK MODULES

## REQUIRED PACKAGES
# 'singleCellTK'            BioConductor    BiocManager::install('singleCellTK')    >= 2.4.0
# 'BUSpaRse'                BioConductor    BiocManager::install('BUSpaRse')    1.10.0
# 'scSensitiveGeneDefine'   github          remotes::install_github('Zechuan-Chen/scSensitiveGeneDefine')   0.1
# 'entropy'                 CRAN            install.packages('entropy')
# 'TInGa'                   github          remotes::install_github('Helena-todd/TInGa')    0.0.0.9000
# cerebroApp                github          remotes::install_github('romanhaa/cerebroApp')   1.3.1

## To get python dependencies : singleCellTK::sctkPythonInstallVirtualEnv(envname = 'sctkpy)

## LAUNCH singleCellTK (this function won't actually work in every place for singleCelTK v2.4, as some parts of the workflow will overload the Matrix hack)
scTK_launch <- function(external_browser = TRUE) {
  # library(reticulate)
  ## Install venv
  # sctkPythonInstallConda(envname = 'sctk', packages = c("pip", "scipy", "numpy", "astroid", "six"), pipIgnoreInstalled = FALSE, pythonVersion = '3.9')
  # sctkPythonInstallVirtualEnv(envname = 'sctkvenv', packages = c("pip", "scipy", "numpy", "astroid", "six"))
  ## Select python environment
  # selectSCTKVirtualEnvironment('sctkvenv')
  # selectSCTKConda('sctk')
  
  # options("shiny.launch.browser" = .rs.invokeShinyWindowExternal)
  # singleCellTK::singleCellTK()

  library(Matrix)
  ## Override Matrix deprecation error status
  options('Matrix.warnDeprecatedCoerce' = 0)
  ## Force using web browser rather than Rstudio
  if(external_browser) options("shiny.launch.browser" = .rs.invokeShinyWindowExternal)
  ## Launch scTK
  singleCellTK::singleCellTK()
}

scTK_recompress <- function(in_sobj = NULL, compress = 'bzip2') {
  sobj <- readRDS(in_sobj)
  message('Recompressing [', in_sobj, '] ...')
  saveRDS(object = sobj, file = in_sobj, compress = 'bzip2')
}
## LOAD 10X / BUStools / Alevin / UMI-tools DATA TO SingleCellExperiment OBJECT, SAVE IT ON DISK AND/OR RETURN IT
### data_path       [char]    Path to SC data (10X, BUStools, Alevin, or UMI-tools). If input data come from BUStools or UMI-tools, this path should also contain the rootname of the generated files. Default [NULL]
### sample_name     [char]    Name to give for the loaded data. Default ['SAMPLE']
### exp_name        [char]    The main Experiment name to set for the output SCE object (see ?SingleCellExperiment::mainExpName). Default ['RNA'].
### out_rds         ['auto'|char|NULL]    Path+name.rds of the disk output of the generated SCE object, save as a bzip2-compressed RDS archive. If 'auto', the RDS filename will be generated automatically, and written in the same folder as data_path. IF NULL, nothing will be written on disk. Default [auto].
### return_data     [logic]   Should the SCE object be returned by the function ? Default [FALSE].
scTK_load <- function(data_path = NULL, sample_name = 'SAMPLE', exp_name = 'RNA', out_rds = 'auto', return_data = FALSE) {
  
  ## Checks
  if(is.null(sample_name)) stop('A sample name is required !')
  sample_name <- as.character(sample_name)
  if(!dir.exists(dirname(data_path))) stop('Input data path does not exist !')
  if(is.null(out_rds)) message('No data will be written on disk')
  if(!out_rds == 'auto' & !dir.exists(dirname(out_rds))) stop('Output directory does not exist !')
  # if(file.exists(out_rds)) stop('The requested output rds file name already exists !')
  if(!is.character(exp_name)) stop('Experiment name should be a character (string)')
  if(!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  if(is.null(out_rds) & !return_data) stop('No "out_rds" ste to NULL and "return_data" set to FALSE : Nothing to do !')
  
  if(out_rds != 'auto') out_dir <- dirname(out_rds)
  
  ## Loading data
  message("Loading data ...")
  if(file.exists(paste0(data_path, "/matrix.mtx")) | file.exists(paste0(data_path, "/matrix.mtx.gz"))) {
    ### Cell Ranger
    message('Found CellRanger data.')
    scmat <- Seurat::Read10X(data_path)
    if ('Gene Expression' %in% names (scmat)) {
      message("Multiple matrices found. Keeping gene expression only")
      scmat <- scmat[['Gene Expression']]
    }
    if(out_rds == 'auto') out_dir <- data_path
  # } else if(file.exists(paste0(data_path, "/", sample_name, ".mtx"))) {
  } else if(file.exists(paste0(data_path, ".mtx"))) {
    ### BUStools
    message('Found BUStools data.')
    library(Matrix)
    options("Matrix.warnDeprecatedCoerce" = 1)
    scmat <- BUSpaRse::read_count_output(dir = dirname(data_path), name = basename(data_path), tcc = FALSE)
    if(out_rds == 'auto') out_dir <- dirname(data_path)
  } else if (file.exists(paste0(data_path, "/quants_mat.gz"))) {
    ### Alevin
    message('Found Alevin data.')
    scmat <- Seurat::ReadAlevin(data_path)
    if(out_rds == 'auto') out_dir <- data_path
  } else if (file.exists(paste0(data_path, '_counts.tsv.gz'))) {
    ### UMI-tools
    message('Found UMI-tools data.')
    scmat <- read.table(file = paste0(data_path, '_counts.tsv.gz'), header = TRUE, sep = "\t", quote = "", check.names = FALSE, row.names = 1)
    if(out_rds == 'auto') out_dir <- dirname(data_path)
  } else {
    stop(paste0("No data found in [", data_path, "]. Wrong path, or unsupported format ..."))
  }
  
  ## Reorder scmat
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
  sobj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=scmat), mainExpName = exp_name)
  rm(scmat)
  ## Adding metrics
  sobj$sample <- as.factor(sample_name)
  sobj[[paste0('nCount_', exp_name)]] <- numi_drop
  sobj[[paste0('nFeature_', exp_name)]] <- ngen_drop
  sobj[[paste0('log_nCount_', exp_name)]] <- log(sobj[[paste0('nCount_', exp_name)]]+1)
  sobj[[paste0('log_nFeature_', exp_name)]] <- log(sobj[[paste0('nFeature_', exp_name)]]+1)
  
  ## Fill misc metadata : metrics
  mymeta <- list(misc=list())
  mymeta$misc$params$Droplet_Quality$captured_droplet <- droplets.nb
  mymeta$misc$params$Droplet_Quality$total_number_UMI <- umi.total.nb
  mymeta$misc$samplename <- sample_name
  mymeta$misc$id <- 0
  sobj@metadata <- mymeta
  
  ## File output
  if(is.null(out_rds)) {
    message('Saving RDS ...')
    if (out_rds == 'auto') out_rds <- paste0(out_dir, '/', paste(c(sample_name, paste0(sprintf('%02d', sobj@metadata$misc$id), 'a'), 'raw_SCE.rds'), collapse = '_'))
  saveRDS(object = sobj, file = out_rds, compress = 'bzip2')
  }
  ## Return data ?
  if (return_data) return(sobj)
}

## scTK_load() EXAMPLE
# scTK_load(data_path = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix', sample_name = 'PBMC1K10X')


## SCE OBJECT DESCRIPTION
## in_rds         [char]      Path to a SCE object saved as a RDS
## max_levels     [int>0]     Maximal number of unique values to consider a barcode annotation as a factor rather than a continuous numeric vector
## describe       ['all', 'assays', 'dimred', 'coldata']    Type of entries to describe
scTK_descriptor <- function(in_rds = NULL, describe = 'all', sparse_level = TRUE, max_levels = 100, assay_plot = FALSE, out_dir = 'auto', return_data = FALSE) {
  ## Checks
  if (is.null(in_rds)) stop('A RDS containing a SingleCellExperiment object is required !')
  if(!file.exists(in_rds)) stop('Provided RDS not found !')
  
  describe <- tolower(describe)
  desc.valids <- c('all', 'assay', 'dimred', 'coldata')
  if(!all(describe %in% desc.valids)) stop('At least one requested item type to describe is not valid. Expecting any combination of ["', paste(c(desc.valids), collapse = '", "'), '] ("all" supersedes any other).')
  if('all' %in% describe) describe <- desc.valids[-1]
  
  ## Loading SCE
  message('Loading SCE object ...')
  sobj <- readRDS(in_rds)
  
  ## Additional checks on sobj
  if(!is(sobj, 'SingleCellExperiment')) stop('Provided RDS is not a proper SingleCellExperiment object !')
  
  if (assay_plot) {
    if(out_dir == 'auto') out_dir <- dirname(in_rds)
    rootname <- sub(pattern = '.rds', replacement = '', x = basename(in_rds), ignore.case = TRUE)
  } 
  retlist <- list()
  ## Handling merged case
  if(length(sobj@metadata$misc) > 1 & all(vapply(sobj@metadata$misc, is.list, TRUE)))
 sobj@metadata$misc <- list(samplename = paste(vapply(sobj@metadata$misc, function(x) x$samplename, "a"), collapse = '.'), id = max(vapply(sobj@metadata$misc, function(x) x$id, 1)))
  
  ## Getting samplename
  samplename <- sobj@metadata$misc$samplename
  
  ## MAIN EXPERIMENT
  if(!return_data) {
    ### EXPERIMENT NAME
    exp_name <- SingleCellExperiment::mainExpName(sobj)
    message('EXPERIMENT : [', exp_name, '] (main)')
    ### ASSAYS
    if('assay' %in% describe) {
      expassays <- SummarizedExperiment::assays(x = sobj)
      for (ea in seq_along(expassays)) {
        assay_name <- names(expassays)[ea]
        message('\tASSAY ', ea, ' : [', assay_name, ']  Dims:[', nrow(expassays[[ea]]), ' x ', ncol(expassays[[ea]]), ']  Range:[', paste(sprintf('%.2f', range(expassays[[ea]], na.rm = TRUE)), collapse = '-'), ']  Type:[', sobj@metadata$assayType$assayTag[sobj@metadata$assayType$assayName == assay_name], ']')
        if (sparse_level & is(expassays[[ea]], 'dgCMatrix')) {
          splev <- sum(sparseMatrixStats::colCounts(x = expassays[[ea]], value = 0)) / prod(dim(expassays[[ea]]))
          message('\t\tSparsity level : ', sprintf('%.2f', splev * 100), '%')
          message('\t\tCounts : ', round(sum(expassays[[ea]])))
        }
        if(assay_plot) {
          png(paste0(out_dir, '/', paste(c(rootname, exp_name, assay_name))))
          plot3D::persp3D(z=log(as.matrix(SummarizedExperiment::assay(x = sobj, i = assay_name))+1), xlab = 'Features', ylab = 'Cells', zlab = "Value", main = assay_name)
          dev.off()
        }
      }
    }
    ### DIMRED
    if ('dimred' %in% describe) {
      dimreds <- SingleCellExperiment::reducedDims(x = sobj)
      for (dr in seq_along(dimreds)) message('\tDIMRED ', dr, ' : [', names(dimreds)[dr], ']  Dims:[', nrow(dimreds[[dr]]), ' x ', ncol(dimreds[[dr]]), ']')
    } else retlist[['main']] <- list(name = SingleCellExperiment::mainExpName(sobj), assay = names(sobj@assays))
  }
  
  ## Alt exp
  alt.names <- SingleCellExperiment::altExpNames(sobj)
  for (en in seq_along(alt.names)) {
    if(!return_data) {
      message('EXPERIMENT : [', alt.names[en], '] (alt)')
      expassays <- names(SingleCellExperiment::altExp(x = sobj, e = alt.names[en])@assays)
      if('assay' %in% describe) {
        for (ea in seq_along(expassays)) {
          message('\tASSAY ', ea, ' : [', expassays[ea], ']  Dims:[', nrow(SummarizedExperiment::assay(x = SingleCellExperiment::altExp(x = sobj, e = alt.names[en]), i = expassays[ea])), ' x ', ncol(SummarizedExperiment::assay(x = SingleCellExperiment::altExp(x = sobj, e = alt.names[en]), i = expassays[ea])), ']  Range:[', paste(sprintf('%.2f', range(SummarizedExperiment::assay(x = SingleCellExperiment::altExp(x = sobj, e = alt.names[en]), i = expassays[ea]), na.rm = TRUE)), collapse = '-'), ']  Type:[', sobj@metadata$assayType$assayTag[sobj@metadata$assayType$assayName == expassays[ea]], ']')
          if(assay_plot) {
            png(paste0(out_dir, '/', paste(c(rootname, exp_name, assay_name))))
            plot3D::persp3D(z=log(as.matrix(SingleCellExperiment::altExp(x = sobj, e = SummarizedExperiment::assay(x = sobj, i = assay_name)))+1), xlab = 'Features', ylab = 'Cells', zlab = "Value", main = assay_name)
            dev.off()
          }
        }
      }
    } else retlist[['alt']][[alt.names[en]]] <- names(SingleCellExperiment::altExp(x = sobj, e = alt.names[en])@assays)
  }
  
  ### BARCODE META
  if(!return_data & 'coldata' %in% describe) {
    bc.df <- as.data.frame(sobj@colData)
    message('BARCODES METADATA :')
    for (b in seq_len(ncol(bc.df))) {
      b.fac <- as.factor(bc.df[,b])
      if (nlevels(b.fac) < max_levels) {
        b.tbl <- as.data.frame(table(b.fac, useNA = 'always'))
        colnames(b.tbl)[1] <- colnames(bc.df)[b]
        print(knitr::kable(b.tbl, escape = FALSE))
      } else {
        # message(colnames(bc.df)[b], ' : CONTI')
        message('\n',colnames(bc.df)[b])
        print(summary(bc.df[,b]))
        # message('\t\t', colnames(bc.df)[b], ' : ', print(summary(b.fac)))
      }
    }
  }
  if(return_data) return(retlist)
}


## FILTER EMPTY DROPLETS USING DropletUtils::emptyDrops
### in_rds              [char]        Path to the SCE object save as RDS to read. Default [NULL]
### droplets_min        [int>>0]      Minimum number of barcodes found in the input SCE object. This is a lose flag to avoid filtering empty droplets on a dataset already filtered for those. The value should be fat superior to the expected amount of cells. Default [1E+04]
### emptydrops_fdr      [0<num<<1]    Maximum p-value from the emptyDrops test to consider a droplet as non-empty. Default [1E-03].
### emptydrops_retain   [char]        Minimum UMI count above which all barcodes are assumed to contain cells (see ?DropletUtils::emptyDrops). Default [NULL] (ie, no minimum value)
### my_seed             [int>0]       Seed used to start the RNG. Default [1337]
### draw_plots          [logic]       Should metric plots be drawn ? Default [TRUE]
### out_rds             ['auto'|char|NULL]    Path+name.rds of the disk output of the generated SCE object, save as a bzip2-compressed RDS archive. If 'auto', the RDS filename will be generated automatically, and written in the same folder as data_path. IF NULL, nothing will be written on disk. Default [auto].
### return_data         [logical]     Should the SCE object be returned by the function ? Default [FALSE].
scTK_edf <- function(in_rds = NULL, assay = 'counts', droplets_min = 1E+04, emptydrops_fdr = 1E-03, emptydrops_retain = NULL, my_seed = 1337, draw_plots = TRUE, out_rds = 'auto', return_data = FALSE) {
  
  ## Checks
  if(is.null(in_rds)) stop('A RDS containing a SingleCellExperiment object is required !')
  if(!file.exists(in_rds)) stop('Provided RDS not found !')
  if(is.null(out_rds)) message('No data will be written on disk')
  if(out_rds != 'auto' & !dir.exists(dirname(out_rds))) stop('Output directory does not exist !')
  if(!is.character(assay)) stop('Assay name should be a character (string)')
  if(!is.logical(draw_plots)) stop('The draw_plots parameter should be a logical (boolean)')
  if(!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  if(is.null(out_rds) & !return_data) stop('No "out_rds" ste to NULL and "return_data" set to FALSE : Nothing to do !')
  
  ## Loading sobj
  message('Loading SCE object ...')
  sobj <- readRDS(in_rds)
  
  ## Additional checks on sobj
  if(!is(sobj, 'SingleCellExperiment')) stop('Provided RDS is not a proper SingleCellExperiment object !')
  if(!assay %in% names(SummarizedExperiment::assays(sobj))) stop('Provided assay not found in the provided SCE object.')
  
  ## Setting out_dir
  out_dir <- if(out_rds == 'auto') dirname(in_rds) else dirname(out_rds)
  
  ## Load sc matrix
  scmat <- SummarizedExperiment::assay(x = sobj, i = assay)
  
  message('Full droplets matrix dimensions :')
  droplets.nb <- ncol(scmat)
  print(dim(sobj))
  
  message('Total UMIs :')
  umi.total.nb <- sum(scmat)
  print(umi.total.nb)
  
  ## Handling merged case
  if(length(sobj@metadata$misc) > 1 & all(vapply(sobj@metadata$misc, is.list, TRUE)))
    sobj@metadata$misc <- list(samplename = paste(vapply(sobj@metadata$misc, function(x) x$samplename, "a"), collapse = '.'), id = max(vapply(sobj@metadata$misc, function(x) x$id, 1)))
  
  ## Incrementing id
  sobj@metadata$misc$id <- sobj@metadata$misc$id + 1
  
  if (!is.null(droplets_min) && droplets.nb > droplets_min && !is.null(emptydrops_fdr)) {
    ## Removing empty droplets
    message("Removing empty droplets with emptyDrops ...")
    bc_rank <- DropletUtils::barcodeRanks(scmat)
    set.seed(my_seed)
    bc_rank2 <- DropletUtils::emptyDrops(scmat, retain = emptydrops_retain)
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
      print(dim(scmat[,keep.bc]))
      message('Total UMIs (filtered) :')
      umi.kept.nb <- sum(scmat[,keep.bc])
      print(umi.kept.nb)
      message('Fraction of UMIs in cells :')
      print(umi.kept.nb / umi.total.nb)
      
      ## Plotting when requested
      if(draw_plots){
        
        ## Computing some metrics
        numi_drop <- Matrix::colSums(x = scmat[,keep.bc])
        scmat.bin <- scmat[,keep.bc]
        attr(scmat.bin, "x")[attr(scmat.bin, "x") >= 1] <- 1
        scmat.bin <- Matrix::drop0(scmat.bin)
        ngen_drop <- Matrix::colSums(x = scmat.bin)
        rm(scmat.bin)
        
        ### Kneeplot
        png(filename = paste0(out_dir, '/', paste(c(sobj@metadata$misc$samplename, sprintf('%02d', sobj@metadata$misc$id), "kneeplot.png"), collapse = '_')), width = 1000, height = 700)
        plot(bc_rank$rank, bc_rank$total+1, log = "xy", xlab = "Rank", ylab = "Total", col = ifelse(keep.bc, "red", "black"), pch = 20, cex = ifelse(keep.bc, 1, 2), main = paste0(sobj@metadata$misc$samplename, ' kneeplot (', length(which(keep.bc)), ' barcodes kept)'))
        o <- order(bc_rank$rank)
        lines(bc_rank$rank[o], bc_rank$fitted[o], col = "red")
        abline(h = S4Vectors::metadata(bc_rank)$knee, col = "dodgerblue", lty = 2)
        abline(h = S4Vectors::metadata(bc_rank)$inflection, col = "forestgreen", lty = 2)
        legend("bottomleft", lty = c(2, 2, NA, NA), pch = c(NA, NA, 20, 20), col = c("dodgerblue", "forestgreen", "red", "black"), legend = c("knee", "inflection", "cell", "empty"))
        dev.off()
        rm(bc_rank)
        ### Saturation plot
        png(filename = paste0(out_dir, '/', paste(c(sobj@metadata$misc$samplename, sprintf('%02d', sobj@metadata$misc$id), "satplot.png"), collapse = '_')), width = 1000, height = 700)
        suppressWarnings(plot(numi_drop, ngen_drop, pch = 20, log = "xy", col = ifelse(keep.bc, "red", "black"), xlab = 'Nb of UMIs in droplet (log)', ylab = 'Nb of genes with at least 1 UMI count in droplets (log)', main = paste0(sobj@metadata$misc$samplename, ' saturation plot (', length(which(keep.bc)), ' barcodes kept)')))
        legend("topleft", pch = c(20, 20), col = c("red", "black"), legend = c("cell", "empty"))
        dev.off()
        ### Saturation distrib plot
        # png(filename = paste0(out_dir, '/', paste(c(sample_name, sprintf('%02d', sobj@metadata$misc$id), "satdistplot.png"), collapse = '_')), width = 1000, height = 700)
        # umipgen <- numi_drop / ngen_drop
        # upg.range <- range(umipgen)
        # plot(density(umipgen[keep.bc], col = "red", xlim = c(0, upg.range[2], )))
      }
      
      ## Applying filter
      sobj <- sobj[,keep.bc]
    }
    rm(scmat)
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
  if(!is.null(out_rds)) {
    samplename <- if(is.list(sobj@metadata$misc[[1]])) paste(vapply(sobj@metadata$misc, function(x) x$samplename, "a"), collapse = '.') else sobj@metadata$misc$samplename
    message('Saving RDS ...')
    if(out_rds == 'auto') out_rds <- paste0(out_dir, '/', paste(c(samplename, paste0(sprintf('%02d', sobj@metadata$misc$id), 'a'), 'EDf.rds'), collapse = '_'))
    saveRDS(object = sobj, file = out_rds, compress = 'bzip2')
  }
  
  ## Return data ?
  if (return_data) return(sobj)
}

## scTK_edf() EXXAMPLE
# scTK_edf(in_rds = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/PBMC1K10X_raw.rds')


## CREATE A QUICK AND DIRTY 2D UMAP
scTK_QnDuMAP <- function(in_rds = NULL, exp_name = NULL, assay = 'SLN', pca_comp = 10, my_seed = 1337, return_data = FALSE) {
  
  message('Checks ...')
  ### Mandatory
  if (is.null(in_rds)) stop('A RDS containing a SingleCellExperiment object is required !')
  if(!file.exists(in_rds)) stop('Provided RDS not found !')
  if(!is.character(assay)) stop('Assay name should be a character (string)')
  if(pca_comp <= 0) stop('[pca_comp] should be a non-null positive integer (and <= N cells).')
  
  ## Loading sobj
  message('Loading SCE object ...')
  sobj <- readRDS(in_rds)
  
  ## Additional checks on sobj
  if(!is(sobj, 'SingleCellExperiment')) stop('Provided RDS is not a proper SingleCellExperiment object !')
  
  ## Checking if requested exp and assay (and feature subset) exist
  expassay.check <- suppressMessages(d2(sobj))
  ## Setting type of experience
  exp_type = NULL
  if(is.null(exp_name)) {
    if (is.null(expassay.check$main$name)) exp_type <- 'main'
  } else {
    if (exp_name %in% names(expassay.check$alt)) exp_type <- 'alt'
  }
  if(is.null(exp_type)) stop('Could not find the requested experiment in neither main nor alternate experiments !')
  if(exp_type == 'main') {
    if (!assay %in% expassay.check[['main']][['assay']]) stop('Requested assay does not exist for the main experiment !')
  }
  if(exp_type == 'alt') {
    if(!assay %in% expassay.check[['alt']][[exp_name]]) stop('Requested assay does not exist for the requested alternate experiment !')
  }
  
  message('Exp type : ', exp_type)
  ## Setting out_dir
  out_dir <- dirname(in_rds)
  
  ## Loading data 
  scmat <- if(exp_type == 'main') SummarizedExperiment::assay(x = sobj, i = assay) else SummarizedExperiment::assay(x = SingleCellExperiment::altExp(x = sobj, e = exp_name), i = assay)
  
  ## Compute quick UMPA with scater
  message("Generating a quick'n dirty uMAP ...")
  set.seed(my_seed)
  umap <- scater::calculateUMAP(x = scmat, pca = pca_comp, ntop = nrow(scmat))
  rm(scmat)
  
  ## Adding to sobj
  SingleCellExperiment::reducedDim(x = sobj, type = paste0('QnD_', assay)) <- umap
  rm(umap)
  
  ## Saving obj
  message('Saving to RDS ...')
  saveRDS(object = sobj, file = in_rds, compress = 'bzip2')
  
  if(return_data) return(sobj)
}

## ESTIMATE CELL CYCLE USING THE SEURAT METHOD
### in_rds            [char]      Path to the SCE object save as RDS to read. Default [NULL]
### exp_name          [char]      Experiment name to use. Default [NULL] will use the mainExpName
### assay             [char]      Name of the SCE assay to use (ie, the matrix level). Default ['counts'].
### cc_seurat_file    [char]      Path+Name to a RDS file that contains gene lists used to score the cell cycle phases. Default [NULL]
### nbin              [int>0]     Number of bins of aggregate expression levels for all analyzed features (see ?Seurat::AddModuleScore). Default [24]
### out_rds           ['auto'|char|NULL]    Path+name.rds of the disk output of the generated SCE object, save as a bzip2-compressed RDS archive. If 'auto', the RDS filename will be generated automatically, and written in the same folder as data_path. IF NULL, nothing will be written on disk. Default [auto].
### return_data       [logical]   Should the SCE object be returned by the function ? Default [FALSE].
scTK_cc_seurat <- function(in_rds = NULL, assay = 'counts', cc_seurat_file = NULL, nbin = 24, my_seed = 1337, out_rds = 'auto', return_data = FALSE) {
  
  ## Checks
  message('Checks ...')
  if(is.null(in_rds)) stop('A RDS containing a SingleCellExperiment object is required !')
  if(!file.exists(in_rds)) stop('Provided RDS not found !')
  if(is.null(out_rds)) message('No data will be written on disk')
  if(!out_rds == 'auto' & !dir.exists(dirname(out_rds))) stop('Output directory does not exist !')
  if(is.null(cc_seurat_file)) stop('The RDS containing Seurat gene lists is required !')
  if(!is.character(assay)) stop('Assay name should be a character (string)')
  if(!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  if(is.null(out_rds) & !return_data) stop('No "out_rds" ste to NULL and "return_data" set to FALSE : Nothing to do !')
  
  ## Loading sobj
  message('Loading SCE object ...')
  sobj <- readRDS(in_rds)
  
  ## Additional checks on sobj
  if(!is(sobj, 'SingleCellExperiment')) stop('Provided RDS is not a proper SingleCellExperiment object !')
  if(!assay %in% names(SummarizedExperiment::assays(sobj))) stop('Provided assay not found in the provided SCE object.')
  
  ## Setting out_dir
  out_dir <- if(out_rds == 'auto') dirname(in_rds) else dirname(out_rds)
  
  ## Loading cc file
  message('Loading Seurat cell cyle genes ...')
  cc_data <- readRDS(cc_seurat_file)
  
  ## Converting to Seurat obj
  message('Converting to Seurat object ...')
  seu.counts <- SummarizedExperiment::assay(x = sobj, i = assay)
  
  seu.meta <- as.data.frame(sobj@colData)
  rownames(seu.meta) <- colnames(seu.counts)
  sobj_seu <- Seurat::CreateSeuratObject(counts = seu.counts, project = sobj@metadata$misc$samplename, assay = 'temp', meta.data = seu.meta)
  rm(seu.meta, seu.counts)
  
  ## Quick log-normalization
  sobj_seu <- Seurat::NormalizeData(object = sobj_seu, normalization.method = 'LogNormalize')
  
  ### Seurat
  message('Cell cycle scoring with Seurat ...')
  sobj_seu <- Seurat::CellCycleScoring(object = sobj_seu, s.features = cc_data$s.genes, g2m.features = cc_data$g2m.genes, assay = 'temp', nbin = nbin, seed = my_seed)
  sobj$cc_seurat.S.Score <- sobj_seu$S.Score
  sobj$cc_seurat.G2M.Score <- sobj_seu$G2M.Score
  sobj$cc_seurat.SmG2M.Score <- sobj_seu$S.Score - sobj_seu$G2M.Score
  sobj$cc_seurat.Phase <- sobj_seu$Phase
  rm(sobj_seu)
  message("Cell cycle phases according to Seurat: ")
  print(table(sobj$cc_seurat.Phase, useNA = 'always'))
  
  ## Handling merged case
  if(length(sobj@metadata$misc) > 1 & all(vapply(sobj@metadata$misc, is.list, TRUE)))
    sobj@metadata$misc <- list(samplename = paste(vapply(sobj@metadata$misc, function(x) x$samplename, "a"), collapse = '.'), id = max(vapply(sobj@metadata$misc, function(x) x$id, 1)))
  
  ## Adding id
  sobj@metadata$misc$id <- sobj@metadata$misc$id+1
  
  ## File output
  message('Saving results ...')
  if (!is.null(out_rds)) {
    message('Saving RDS ...')
    samplename <- sobj@metadata$misc$samplename
    out_rds <- if(out_rds == 'auto') paste0(out_dir, '/', paste(c(samplename, paste0(sprintf('%02d', sobj@metadata$misc$id), 'a'), assay, 'CC.seurat.rds'), collapse = '_'))
    saveRDS(object = sobj, file = out_rds, compress = 'bzip2')
  }
  ## Return data ?
  if (return_data) return(sobj)
}

## scTK_cc_seurat() EXAMPLE
# scTK_cc_seurat(in_rds = readRDS('/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/PBMC1K10X.rds'), cc_seurat_file = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/RESOURCES/GENELISTS/homo_sapiens_Seurat_cc.genes_20191031.rds')


## scTK_cc_cyclone() EXAMPLE
# scTK_cc_cyclone(in_rds = readRDS('/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/10X_DATASET_1kPBMC_CR3v3/COUNT_MATRIX/pbmc_1k_v3_raw_feature_bc_matrix/PBMC1K10X_cc.seurat.rds'), cc_cyclone_file = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/RESOURCES/GENELISTS/homo_sapiens_cyclone_pairs_symbols_20191001.rds')


## GET SENSITIVE GENES INSTEAD OF SEURAT HVG
scTK_scSG <- function(in_rds = NULL, exp_name = NULL, raw_assay = 'counts', norm_assay = 'SLN', n_features = 2000, n_PC = 10, resolution = 0.6, out_rds = 'auto', return_data = FALSE) {
  
  ## Installing
  # remotes::install_github('Zechuan-Chen/scSensitiveGeneDefine')
  # install.packages('entropy')
  
  message('Checks ...')
  ### Mandatory
  if (is.null(in_rds)) stop('A RDS containing a SingleCellExperiment object is required !')
  if(!file.exists(in_rds)) stop('Provided RDS not found !')
  if(!is.character(raw_assay)) stop('Raw assay name should be a character (string)')
  if(!is.character(norm_assay)) stop('Normalized assay name should be a character (string)')
  
  ## Loading sobj
  message('Loading SCE object ...')
  sobj <- readRDS(in_rds)
  
  ## Additional checks on sobj
  if(!is(sobj, 'SingleCellExperiment')) stop('Provided RDS is not a proper SingleCellExperiment object !')
  
  ## Checking if requested exp and assay (and feature subset) exist
  expassay.check <- scTK_descriptor(in_rds = in_rds, return_data = TRUE)
  ## Setting type of experience
  exp_type = NULL
  if(is.null(exp_name)) {
    if (is.null(expassay.check$main$name)) exp_type <- 'main'
  } else {
    if (exp_name %in% names(expassay.check$alt)) exp_type <- 'alt'
  }
  if(is.null(exp_type)) stop('Could not find the requested experiment in neither main nor alternate experiments !')
  if(exp_type == 'main') {
    if (!assay %in% expassay.check[['main']][['assay']]) stop('Requested raw assay does not exist for the main experiment !')
    if (!norm_assay %in% expassay.check[['main']][['assay']]) stop('Requested normalized assay does not exist for the main experiment !')
  }
  if(exp_type == 'alt') {
    if(!assay %in% expassay.check[['alt']][[exp_name]]) stop('Requested raw assay does not exist for the requested alternate experiment !')
    if(!assay %in% expassay.check[['alt']][[exp_name]]) stop('Requested normalized assay does not exist for the requested alternate experiment !')
  }
  
  ## Setting out_dir
  # out_dir <- if(out_rds == 'auto') dirname(in_rds) else dirname(out_rds)
  
  ## Loading data 
  raw_scmat <- if(exp_type == 'main') SummarizedExperiment::assay(x = sobj, i = raw_assay) else SummarizedExperiment::assay(x = SingleCellExperiment::altExp(x = sobj, e = exp_name), i = raw_assay)
  norm_scmat <- if(exp_type == 'main') SummarizedExperiment::assay(x = sobj, i = norm_assay) else SummarizedExperiment::assay(x = SingleCellExperiment::altExp(x = sobj, e = exp_name), i = norm_assay)
  
  ## Handling merged case
  if(length(sobj@metadata$misc) > 1 & all(vapply(sobj@metadata$misc, is.list, TRUE)))
    sobj@metadata$misc <- list(samplename = paste(vapply(sobj@metadata$misc, function(x) x$samplename, "a"), collapse = '.'), id = max(vapply(sobj@metadata$misc, function(x) x$id, 1)))
  
  ## Getting samplename
  samplename <- sobj@metadata$misc$samplename
  
  ## Converting to Seurat Object
  sobj_seu <- Seurat::CreateSeuratObject(counts = raw_scmat, data = norm_data, project = samplename, assay = 'RNA', meta.data = as.data.frame(sobj@colData))
  sobj_seu <- Seurat::FindVariableFeatures(object = sobj_seu, nfeatures = n_features)
  sobj_seu <- Seurat::ScaleData(object = sobj_seu)
  sobj_seu <- Seurat::RunPCA(object = sobj_seu, npcs = n_PC+1, verbose = FALSE)
  sobj_seu <- Seurat::FindNeighbors(object = sobj_seu)
  sobj_seu <- Seurat::FindClusters(object = sobj_seu, resolution = .6)
  
  ## Incrementing id
  sobj@metadata$misc$id <- sobj@metadata$misc$id +1
  
  if(!is.null(out_rds)) {
    out_root <- if(out_rds == 'auto') paste0(dirname(in_rds), '/', paste(c(samplename, paste0(sprintf('%02d', sobj@metadata$misc$id), 'a'), paste0('SG', n_features)), collapse = '_')) else out_root <- sub(pattern = '.rds$', replacement = '', x = out_rds, ignore.case = TRUE)
  }
  
  ## Variance plot for Seurat HVGs
  if (!is.null(out_rds)) {
    png(paste0(out_root, '_VarPlot_HVGs.png'), width = 800, height = 600)
    print(Seurat::VariableFeaturePlot(sobj_seu))
    dev.off()
  }
  
  ## Getting sensitive features
  HVG_Anno <- HVG_Statistic2(object = sobj_seu , First_time_unsupervised_clustering_label = paste0('RNA_snn_res.', resolution), nfeatures = n_features)
  SensitiveGene <- GetSensitivegene2(object = sobj_seu, HVG_Anno = HVG_Anno)
  
  ## getting filtered features
  sobj_seu <- ReSelectVariableFeatures2(object = sobj_seu, SensitiveGene = SensitiveGene, nfeatures = n_features)
  ok_sg <- sobj_seu@assays$RNA@var.features
  ## Variance plot for SGs
  if (!is.null(out_rds)) {
    png(paste0(out_root, '_VarPlot_SGs.png'), width = 800, height = 600)
    print(Seurat::VariableFeaturePlot(sobj_seu))
    dev.off()
  }
  
  rm(sobj_seu)
  
  ## Inserting into the SCE
  newname <- paste0('SG', n_features)
  assaylist <- list(assay = norm_scmat[ok_sg,])
  names(assaylist) <- newname
  SingleCellExperiment::altExp(x = sobj, e = newname) <- SingleCellExperiment::SingleCellExperiment(assays = assaylist, mainExpName = newname)
  rm(raw_scmat, norm_scmat)
  
  ## Updating metadata
  sobj@metadata$assayType <- rbind(sobj@metadata$assayType, c('hvg', newname))
  
  ## File output
  if(!is.null(out_rds)) {
    message('Saving RDS ...')
    if(out_rds == 'auto') out_rds <- paste0(out_root, '.rds')
    saveRDS(object = sobj, file = out_rds, compress = "bzip2")
  }
  
  ## Return data ?
  if (return_data) return(sobj)
}


## ASSES THE WEIGHT OF PUTATIVE COVARIATES FROM A [SAMPLES x VARIABLES] ANNOTATION DATAFRAME, ON A [FEATURES x SAMPLES] DATA MATRIX, THROUGH CORRELATION TEST (FOR CONTINUOUS VARIABLES), OR KRUSKAL-WALLIS TEST (FOR FACTORS)
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
## . out_png            ['auto'|char]           Path+name to the output image. If 'auto', the image will be named thanks to the 'in_rds' value. Default ['auto']
scTK_assess_covar <- function(in_rds = NULL, exp_name = NULL, assay = 'counts', factor_names = NULL, conti_names = NULL, ctrl_features = NULL, marker_features = NULL, red_method = 'pca', ndim_max = 10, center = TRUE, scale = TRUE, coef_cut = 0, out_png = 'auto', color_palette = c("white", "orangered3")) {
  
  message('Checks ...')
  ### Mandatory
  if (is.null(in_rds)) stop('A RDS containing a SingleCellExperiment object is required !')
  if(!file.exists(in_rds)) stop('Provided RDS not found !')
  if(!is.character(assay)) stop('Assay name should be a character (string)')
  if(all(is.null(c(factor_names, conti_names, ctrl_features, marker_features)))) stop('At least one of [factor_names], [conti.colnames], [ctrl_features] or [marker_features] should not be NULL.')
  if(ndim_max <= 0) stop('[ndim_max] should be a non-null positive integer (and <= N cells).')
  if(!tolower(red_method) %in% c('pca', 'mds.euc', 'mds.spear')) stop('Unknown or unsupported dimension reduction method')
  if(tolower(out_png) != 'auto' & !dir.exists(dirname(out_png))) stop('Provided path for "out_png" does not exist !')
  
  ## Loading sobj
  message('Loading SCE object ...')
  sobj <- readRDS(in_rds)
  
  ## Additional checks on sobj
  if(!is(sobj, 'SingleCellExperiment')) stop('Provided RDS is not a proper SingleCellExperiment object !')
  
  ## Checking if requested exp and assay (and feature subset) exist
  expassay.check <- suppressMessages(scTK_descriptor(in_rds = in_rds, return_data = TRUE))
  ## Setting type of experience
  exp_type = NULL
  if(is.null(exp_name)) {
    if (is.null(expassay.check$main$name)) exp_type <- 'main'
  } else {
    if (exp_name %in% names(expassay.check$alt)) exp_type <- 'alt'
  }
  if(is.null(exp_type)) stop('Could not find the requested experiment in neither main nor alternate experiments !')
  if(exp_type == 'main') {
    if (!assay %in% expassay.check[['main']][['assay']]) stop('Requested assay does not exist for the main experiment !')
  }
  if(exp_type == 'alt') {
    if(!assay %in% expassay.check[['alt']][[exp_name]]) stop('Requested assay does not exist for the requested alternate experiment !')
  }
  
  message('Exp type : ', exp_type)
  
  ## Setting out_dir
  out_dir <- dirname(in_rds)
  
  ## Loading data 
  mat <- if(exp_type == 'main') SummarizedExperiment::assay(x = sobj, i = assay) else SummarizedExperiment::assay(x = SingleCellExperiment::altExp(x = sobj, e = exp_name), i = assay)
  annot.df <- as.data.frame(sobj@colData)
  
  ## Filtering absent feature(s)
  feat.names <- c(ctrl_features, marker_features)
  feat.types <- c(rep('F.control', length(ctrl_features)), rep('F.marker', length(marker_features)))
  if(!is.null(feat.names)) {
    out.feat.idx <- !feat.names %in% rownames(mat)
    if(all(out.feat.idx)) stop('All requested features are absent of the assay !')
    if(any(out.feat.idx)) {
      message('Feature(s) [', paste(feat.names[out.feat.idx], collapse = ', '), '] is/are absent from the assay, thus is/are discarded.')
      feat.names <- feat.names[!out.feat.idx]
      feat.types <- feat.types[!out.feat.idx]
    }
  }
  
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
  
  ## Center / scale ?
  if (any(c(center, scale))) {
    message('Centering and/or scaling data ...')
    mat <- base::scale(x = mat, center = center, scale = scale)
  }
  ## Dimension reduction
  message('Performing dimension reduction ...')
  
  norm.red <- if (tolower(red_method) == 'pca') irlba::prcomp_irlba(x = t(mat), n = ndim_max, center = FALSE, scale. = FALSE)$x else if (tolower(red_method) == 'mds.euc') stats::cmdscale(d = dist(x = t(mat), method = 'euclidean'), k = ndim_max) else if (tolower(red_method) == 'mds.spear') stats::cmdscale(d = as.dist(1-cor(mat, method = 'spearman')), k = ndim_max)
  
  ## Prepping covariates
  covar.names <- c(factor_names, conti_names)
  covar.types <- c(rep('C.factor', length(factor_names)), rep('C.continuous', length(conti_names)))
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
  ### Features ...
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
  
  if(tolower(out_png) == 'auto'){
    samplename <- if(is.list(sobj@metadata$misc[[1]])) paste(vapply(sobj@metadata$misc, function(x) x$samplename, "a"), collapse = '.') else sobj@metadata$misc$samplename
    out_png <- paste0(out_dir, '/', paste(c(samplename, sprintf('%02d', sobj@metadata$misc$id), exp_name, assay, 'covar.png'), collapse = '_'))
  }
  png(filename = out_png, width = 400+(50*length(cf.names)), height = 1000)
  ComplexHeatmap::draw(BC.hm)
  dev.off()
}

## scTK_assess_covar() EXAMPLE
# scTK_assess_covar(rds_in = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/TEST_DATASET/ANALYSIS/testRDS.rds', assay = 'SoupX', factor_names = c('hsn', 'cc_seurat.Phase'), conti_names = c('log_nCount_RNA', 'log_nFeature_RNA', 'soupX_nUMIs', 'subsets_Ribo_percent', 'subsets_Stress_percent', 'subsets_Mito_percent', 'cc_seurat.SmG2M.Score'), ctrl_features = c('Gapdh'), marker_features = c('Il2ra', 'Cd8b1', 'Cd8a', 'Cd4', 'Ccr7', 'Itm2a', 'Aif1', 'Hba-a1'), ndim_max = 25)

# scTK_assess_covar(in_rds = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/DATASETS/PAIVA/20221107/PAIVA_02d_SLN2K.rds', exp_name = 'SLN2K', assay = 'SLN2K', factor_names = c('cc_seurat.Phase'), conti_names = c('subsets_Ribo_percent', 'subsets_Stress_percent', 'subsets_Mito_percent', 'cc_seurat.SmG2M.Score'), ctrl_features = c('Gapdh'), marker_features = c('Il2ra', 'Plac8', 'Ly6d', 'Ccr7', 'Itm2a', 'C1qb', 'Hmgn2'), ndim_max = 25)


# system.time(scTK_assess_covar(in_rds = '/home/job/WORKSPACE/ENCADREMENT/2022/EBAII_n1_2022/ATELIER_SC/DATASETS/PAIVA/GSE160135_RAW/G3_02c_SLNst2K.rds', , exp_name = 'SLNst2K', assay = 'SLNst2K', factor_names = c('cc_seurat.Phase'), conti_names = c('log_nCount', 'log_nFeature', 'subsets_Ribo_percent', 'subsets_StressA_percent', 'subsets_Mito_percent', 'cc_seurat.SmG2M.Score'), ctrl_features = c('Gapdh'), marker_features = c('Il2ra', 'Plac8', 'Ly6d', 'Ccr7', 'Itm2a', 'C1qb', 'Hmgn2', 'Cd28', 'Slc3a2', 'Cd69'), ndim_max = 10))


## REGRESS COVARIATES
## Function to scale and / or regress covariates on a normalized/transformed matrix (ie NOT on counts nor dimred).
### in_rds            [char]        Path to the SCE object save as RDS to read. Default [NULL]
### exp               [char|NULL]   Name of the SCE experiment to use. If NULL, the default main experiment will be used. Default ['counts'].
### assay             [char]        Name of the SCE assay to use (ie, the matrix level). Default ['counts'].
### model_use         [char]        Model to use for regression. Sould be one of 'linear', 'poisson', 'negbinom'. Default ['linear']
### scale_residuals   [bool]        Perform scaling of the regressed matrix. Default [TRUE]
### scale_limit       [int|NULL]    Perform a trimming of the regressed matrix (erasing outlier). Default [NULL]
### center_residuals  [bool]        Perform centering of the regressed matrix. Default [TRUE]
### out_rds           ['auto'|char|NULL]    Path+name.rds of the disk output of the generated SCE object, save as a bzip2-compressed RDS archive. If 'auto', the RDS filename will be generated automatically, and written in the same folder as data_path. IF NULL, nothing will be written on disk. Default [auto].
### return_data       [logical]   Should the SCE object be returned by the function ? Default [FALSE].
scTK_regress_covar <- function(in_rds = NULL, exp_name = NULL, assay = 'counts', vars_to_regress = NULL, model_use, scale_residuals = TRUE, scale_limit = 10, center_residuals = TRUE, out_rds = 'auto', return_data = FALSE) {
  
  ## Parameters checks
  message('Checks ...')
  if (is.null(in_rds)) stop('A RDS containing a SingleCellExperiment object is required !')
  if(!file.exists(in_rds)) stop('Provided RDS not found !')
  if (!is.character(assay)) stop('Assay name should be a character (string)')
  if (!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  
  ## Loading sobj
  message('Loading SCE object ...')
  sobj <- readRDS(in_rds)
  if(!is(sobj, 'SingleCellExperiment')) stop('Provided RDS is not a proper SingleCellExperiment object !')
  
  ## Checking if requested exp and assay (and feature subset) exist
  expassay.check <- suppressMessages(d2(sobj))
  exp_type = NULL
  if(is.null(exp_name)) {
    if (is.null(expassay.check$main$name)) exp_type <- 'main'
  } else {
    if (exp_name %in% names(expassay.check$alt)) exp_type <- 'alt'
  }
  if(is.null(exp_type)) stop('Could not find the requested experiment in neither main nor alternate experiments !')
  if(exp_type == 'main') {
    if (!assay %in% expassay.check[['main']][['assay']]) stop('Requested assay does not exist for the main experiment !')
  }
  if(exp_type == 'alt') {
    if(!assay %in% expassay.check[['alt']][[exp_name]]) stop('Requested assay does not exist for the requested alternate experiment !')
  }
  
  ## Setting out_dir
  out_dir <- if(out_rds == 'auto') dirname(in_rds) else dirname(out_rds)
  
  ## Loading data
  scmat <- if(exp_type == 'main') SummarizedExperiment::assay(x = sobj, i = assay) else SummarizedExperiment::assay(x = SingleCellExperiment::altExp(x = sobj, e = exp_name), i = assay)
  
  ## Handling multi misc
  
  ## Handling merged case
  if(length(sobj@metadata$misc) > 1 & all(vapply(sobj@metadata$misc, is.list, TRUE)))
    sobj@metadata$misc <- list(samplename = paste(vapply(sobj@metadata$misc, function(x) x$samplename, "a"), collapse = '.'), id = max(vapply(sobj@metadata$misc, function(x) x$id, 1)))
  
  ## Getting samplename
  samplename <- sobj@metadata$misc$samplename
  
  ## Converting to Seurat
  sobj_seu <- Seurat::CreateSeuratObject(counts = SummarizedExperiment::assay(x = sobj, i ='counts'), project = samplename, assay = 'temp', meta.data = as.data.frame(sobj@colData))
  sobj_seu@assays$temp@data <- as.matrix(scmat)
  ## Perform scaling + regression
  sobj_seu <- Seurat::ScaleData(object = sobj_seu, vars.to.regress = vars_to_regress, do.scale = scale_residuals, do.center = center_residuals, scale.max = scale_limit, model.use = model_use)
  regmat <- sobj_seu@assays$temp@scale.data
  rm(sobj_seu)
  ## Inserting into the SCE
  newname <- paste(c(paste(c(assay, if(!is.null(feature_subset)) 'sub' else NULL, 'R'), collapse = '_'), if(scale_residuals) 'S' else NULL), collapse = '')
  assaylist = list(regmat)
  names(assaylist) <- newname
  SingleCellExperiment::altExp(x = sobj, e = newname) <- SingleCellExperiment::SingleCellExperiment(assays = assaylist, mainExpName = newname)
  ## Adding metadata
  sobj@metadata$assayType <- rbind(sobj@metadata$assayType, c('transformed', newname))
  sobj@metadata$misc$id <- sobj@metadata$misc$id + 1
  
  ## File output
  if(!is.null(out_rds)) {
    message('Saving RDS ...')
    out_name <- paste(c(samplename, paste0(sprintf('%02d', sobj@metadata$misc$id), 'a'), model_use, exp_name, assay, 'REG'), collapse = '_')
    if(out_rds == 'auto') out_rds <- paste0(out_dir, '/', out_name, '.rds')
    saveRDS(object = sobj, file = out_rds, compress = 'bzip2')
    reg_table <- data.frame(Regressed.covariate = vars_to_regress, Type = unname(vapply(vars_to_regress, function(x) is(sobj[[x]])[1], 'a')))
    write.table(reg_table, file = sub(pattern = '.rds$', replacement = '.tsv', x = out_rds), sep = '\t', quote = FALSE, row.names = FALSE)
  }
  
  ## Return data ?
  if (return_data) return(sobj)
}


## PERFORM SEURAT UMAP
## Function to perform a Seurat UMAP dimension reduction. This is actually an external hack, as scTK performs weirdly when doing a UMAP (be it Seurat or Scater) on a FMNN reduction : the resulting object is abnormaly huge, can be used in the current session he was generated in, but can't be reloaded when saved to a RDS !
### in_rds            [char]      Path to the SCE object save as RDS to read. Default [NULL]
### in_dimred         [char]      Name of a dimension reduction to use as input. Default [NULL].
### ndim_max          [int>0]     Number of dimension from in_dimred to use. Default [10]
### out_dimred        [char]      Name of the generated UMAP in the SCE. Default ['UMAP'].
### out_rds           ['auto'|char|NULL]    Path+name.rds of the disk output of the generated SCE object, save as a bzip2-compressed RDS archive. If 'auto', the RDS filename will be generated automatically, and written in the same folder as data_path. IF NULL, nothing will be written on disk. Default [auto].
### return_data       [logical]   Should the SCE object be returned by the function ? Default [FALSE].
### ...               [...]       Any parameter to give to Seurat::runUMAP()
scTK_SeuratUMAP <- function(in_rds = NULL, in_dimred = NULL, ndim_max = 10, out_dimred = 'UMAP', out_rds = 'auto', return_data = FALSE, ...) {
  
  ## Parameters checks
  message('Checks ...')
  if (is.null(in_rds)) stop('A RDS containing a SingleCellExperiment object is required !')
  if(!file.exists(in_rds)) stop('Provided RDS not found !')
  if (!is.character(in_dimred)) stop('Assay name should be a character (string)')
  if (!is.character(out_dimred)) stop('Assay name should be a character (string)')
  if (!is.logical(return_data)) stop('The return_data parameter should be a logical (boolean)')
  
  ## Loading sobj
  message('Loading SCE object ...')
  sobj <- readRDS(in_rds)
  if(!is(sobj, 'SingleCellExperiment')) stop('Provided RDS is not a proper SingleCellExperiment object !')
  
  ## Checking if requested exp and assay (and feature subset) exist
  expassay.check <- suppressMessages(d2(sobj))
  exp_type = NULL
  if(is.null(exp_name)) {
    if (is.null(expassay.check$main$name)) exp_type <- 'main'
  } else {
    if (exp_name %in% names(expassay.check$alt)) exp_type <- 'alt'
  }
  if(is.null(exp_type)) stop('Could not find the requested experiment in neither main nor alternate experiments !')
  if(exp_type == 'main') {
    if (!assay %in% expassay.check[['main']][['assay']]) stop('Requested assay does not exist for the main experiment !')
  }
  if(exp_type == 'alt') {
    if(!assay %in% expassay.check[['alt']][[exp_name]]) stop('Requested assay does not exist for the requested alternate experiment !')
  }
  
  ## Setting out_dir
  out_dir <- if(out_rds == 'auto') dirname(in_rds) else dirname(out_rds)
  
  ## Loading data
  scmat <- if(exp_type == 'main') SummarizedExperiment::assay(x = sobj, i = assay) else SummarizedExperiment::assay(x = SingleCellExperiment::altExp(x = sobj, e = exp_name), i = assay)
  
  ## Handling merged case
  if(length(sobj@metadata$misc) > 1 & all(vapply(sobj@metadata$misc, is.list, TRUE)))
    sobj@metadata$misc <- list(samplename = paste(vapply(sobj@metadata$misc, function(x) x$samplename, "a"), collapse = '.'), id = max(vapply(sobj@metadata$misc, function(x) x$id, 1)))
  
  ## Getting samplename
  samplename <- sobj@metadata$misc$samplename
  
  ## Converting to Seurat
  sobj_seu <- Seurat::CreateSeuratObject(counts = SummarizedExperiment::assay(x = sobj, i ='counts'), project = samplename, assay = 'temp', meta.data = as.data.frame(sobj@colData))
  sobj_seu@assays$temp@data <- as.matrix(scmat)
  ## Perform scaling + regression
  sobj_seu <- Seurat::ScaleData(object = sobj_seu, vars.to.regress = vars_to_regress, do.scale = scale_residuals, do.center = center_residuals, scale.max = scale_limit, model.use = model_use, ...)
  regmat <- sobj_seu@assays$temp@scale.data
  rm(sobj_seu)
  ## Inserting into the SCE
  newname <- paste(c(paste(c(assay, if(!is.null(feature_subset)) 'sub' else NULL, 'R'), collapse = '_'), if(scale_residuals) 'S' else NULL), collapse = '')
  
  assaylist = list(regmat)
  names(assaylist) <- newname
  SingleCellExperiment::altExp(x = sobj, e = newname) <- SingleCellExperiment::SingleCellExperiment(assays = assaylist, mainExpName = newname)
  ## Adding metadata
  sobj@metadata$assayType <- rbind(sobj@metadata$assayType, c('transformed', newname))
  sobj@metadata$misc$id <- sobj@metadata$misc$id + 1
  
  ## File output
  if(!is.null(out_rds)) {
    message('Saving RDS ...')
    out_name <- paste(c(samplename, paste0(sprintf('%02d', sobj@metadata$misc$id), 'a'), model_use, exp_name, assay, 'REG'), collapse = '_')
    if(out_rds == 'auto') out_rds <- paste0(out_dir, '/', out_name, '.rds')
    saveRDS(object = sobj, file = out_rds, compress = 'bzip2')
    reg_table <- data.frame(Regressed.covariate = vars_to_regress, Type = unname(vapply(vars_to_regress, function(x) is(sobj[[x]])[1], 'a')))
    write.table(reg_table, file = sub(pattern = '.rds$', replacement = '.tsv', x = out_rds), sep = '\t', quote = FALSE, row.names = FALSE)
  }
  
  ## Return data ?
  if (return_data) return(sobj)
}

## (Unexported function) Hack for scSensitiveGenes to support any number of features (was hardcoded to 2000 to compare to Seurat default)
HVG_Statistic2 <- function (object, First_time_unsupervised_clustering_label = "First_time_unsupervised_clustering", nfeatures = 2000) {
  require("entropy")
  require("Seurat")
  Idents(object) <- First_time_unsupervised_clustering_label
  label1 <- levels(object)
  subtype_cells <- list()
  Variable_list <- list()
  Common_HVG <- c()
  meta.features <- list()
  for (i in 1:length(label1)) {
    subtype_cells[[label1[i]]] <- subset(object, idents = label1[i])
    subtype_cells[[label1[i]]] <- FindVariableFeatures(subtype_cells[[label1[i]]], selection.method = "vst", nfeatures = nfeatures)
    Variable_list[[label1[i]]] <- VariableFeatures(subtype_cells[[label1[i]]])
    meta.features[[label1[i]]] <- subtype_cells[[label1[i]]]@assays$RNA@meta.features
    if (i != 1) Common_HVG <- intersect(Common_HVG, Variable_list[[label1[i]]]) else Common_HVG = Variable_list[[label1[i]]]
  }
  Common_HVG.statistic <- data.frame(row.names = label1, nHVG = rep(NA, length(label1)))
  emtropu_value <- c()
  for (j in 1:length(label1)) emtropu_value <- c(emtropu_value, length(Variable_list[[label1[j]]]))
  Common_HVG.statistic$nHVG <- emtropu_value
  Total_HVG <- c()
  for (k in 1:length(label1)) Total_HVG <- c(Total_HVG, Variable_list[[label1[k]]])
  HVG_Type <- unique(Total_HVG)
  HVG_Type.statistic <- data.frame(row.names = HVG_Type, nCluster = rep(NA, length(HVG_Type)), nCell = rep(NA, length(HVG_Type)))
  x <- as.data.frame(table(Total_HVG))
  rownames(x) <- x$Total_HVG
  HVG_Type.statistic$nCluster <- x[HVG_Type, ]$Freq
  for (l in 1:length(HVG_Type)) {
    nCell = 0
    rm(i)
    for (i in 1:length(label1)) {
      if (HVG_Type[l] %in% Variable_list[[label1[i]]]) nCell = nCell + dim(subtype_cells[[label1[i]]]@meta.data)[1]
    }
    HVG_Type.statistic[l, 2] <- nCell
  }
  HVG_Anno <- list(Variable_list = Variable_list, Common_HVG.statistic = Common_HVG.statistic, HVG_Type.statistic = HVG_Type.statistic, meta.features = meta.features)
  return(HVG_Anno)
}

## (Unexported function) Hack for scSensitiveGenes to support any number of features (was hardcoded to 2000 to compare to Seurat default)
GetSensitivegene2 <- function (object = NULL, min_nClusters = "Default", min_nCell = 0, HVG_Anno = NULL) {
  require(entropy)
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(glue)
  if (min_nClusters == "Default") {
    min_nClusters <- floor(length(levels(object))/2)
  }
  pre_Senstive_gene.statistic <- HVG_Anno$HVG_Type.statistic[HVG_Anno$HVG_Type.statistic$nCluster > min_nClusters & HVG_Anno$HVG_Type.statistic$nCell > min_nCell,]
  pre_Senstive_gene <- rownames(pre_Senstive_gene.statistic)
  avg_expression <- AverageExpression(object, assays = "RNA", features = pre_Senstive_gene)
  data <- apply(avg_expression$RNA, 1, function(x) entropy(x))
  data <- as.data.frame(data)
  data$gene <- rownames(data)
  data <- data[order(data$data), ]
  emtropu_value <- data
  emtropu_value$nCluster <- pre_Senstive_gene.statistic[rownames(data),]$nCluster
  emtropu_value$nCell <- pre_Senstive_gene.statistic[rownames(data), ]$nCell
  colnames(emtropu_value) <- c("entropy_value", "gene", "nCluster", "nCell")
  SensitiveGene <- data[data$data > median(data$data, na.rm = TRUE), ]
  colnames(SensitiveGene) <- c("entropy_value", "gene")
  return(SensitiveGene)
}

## (Unexported function) Hack for scSensitiveGenes to support any number of features (was hardcoded to 2000 to compare to Seurat default)
ReSelectVariableFeatures2 <- function (object = NULL, SensitiveGene = NULL, nfeatures = 2000) {
  matrix <- object@assays$RNA@meta.features
  matrix <- matrix[order(matrix$vst.variance.standardized, decreasing = T), ]
  matrix <- matrix[setdiff(rownames(matrix), SensitiveGene$gene),]
  new_Features <- rownames(matrix)[1:nfeatures]
  object@assays$RNA@meta.features$vst.variable <- FALSE
  object@assays$RNA@meta.features[new_Features,]$vst.variable <- TRUE
  VariableFeatures(object) <- new_Features
  return(object)
}

## (Unexported function) Returns a list of experiment names and their assay names
d2 <- function(sobj = NULL) {
  ## Checks
  if(!is(sobj, 'SingleCellExperiment')) stop('Provided RDS is not a proper SingleCellExperiment object !')
  retlist <- list()
  ## MAIN EXP
  expassays <- SummarizedExperiment::assayNames(x = sobj)
  retlist[['main']] <- list(name = SingleCellExperiment::mainExpName(sobj), assay = SummarizedExperiment::assayNames(x = sobj))
  ## ALT EXP
  alt.names <- SingleCellExperiment::altExpNames(sobj)
  for (en in seq_along(alt.names)) retlist[['alt']][[alt.names[en]]] <- names(SingleCellExperiment::altExp(x = sobj, e = alt.names[en])@assays)
  ## RETURN
  return(retlist)
}
