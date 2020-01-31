    #' Carries out WGCNA with default settings or custom settings
    #' @import WGCNA
    #' @import ape
    #' @importFrom grDevices dev.off recordPlot
    #' @importFrom graphics abline plot text
    #' @importFrom stats as.dist hclust setNames
    #' @param colname_correct a logical value. If TRUE (default), "." in gene 
    #' names will be replaced
    #' with "-". This corrects a name change that is induced by R when creating 
    #' a data.frame. If FALSE,
    #' no changes will be made.
    #' @param sft_RsquaredCut desired minimum scale free topology fitting 
    #' index R^2.
    #' Default is 0.80.
    #' @inheritParams WGCNA::blockwiseModules
    #' @inheritParams WGCNA::pickSoftThreshold
    #' @seealso \code{\link[WGCNA]{blockwiseModules}}
    #' @inheritParams WGCNA::adjacency
    #' @seealso \code{\link[WGCNA]{adjacency}}
    #' @note
    #' This is a wrapper for WGCNA.
    #' @export
    #' @return Returns a lists containing network input parameters used 
    #' for WGCNA,
    #' WGCNA module information, and quality control plots.
    #' @examples 
    #' sample_dat_dir<-system.file("extdata", "sample_dat.Rdata", 
    #' package = "GmicR", mustWork = TRUE)
    #' load(sample_dat_dir)
    #' GMIC_Builder<-Auto_WGCNA(sample_dat, mergeCutHeight = 0.35, 
    #' minModuleSize = 10)
    
    Auto_WGCNA<-function(datExpr, colname_correct = TRUE, minModuleSize = 10, 
    deepSplit = 4, networkType = "signed hybrid", TOMType = "unsigned", 
    corFnc = "bicor", mergeCutHeight = 0.25, sft_RsquaredCut = 0.85,
    removeFirst = FALSE,
    reassignThreshold=1e-6, maxBlockSize = 25000,nThreads=NULL){
    
    # replacing . with - in columnn names
    if(isTRUE(colname_correct)){
    colnames(datExpr)<-gsub(".", "-", colnames(datExpr), fixed = TRUE)}
    
    # checking columns names
    message("verify that colnames contain official gene symbols")
    print(colnames(datExpr[seq(1,10)])) #seq_len
    
    # getting ref_genes for analysis
    ref_genes = colnames(datExpr)
    
    # soft power --------------------------------------------------------------
    powers = c(seq(1,10), seq(from = 12, to=20, by=2))
    allowWGCNAThreads(nThreads = nThreads)
    sft = pickSoftThreshold(datExpr, networkType = networkType,
    corFnc = corFnc, RsquaredCut = sft_RsquaredCut,
    removeFirst=removeFirst,
    powerVector = powers, verbose = 5)
    disableWGCNAThreads()
    
    # soft thresholds
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab='Soft Threshold (power)',
    ylab='Scale Free Topology Model Fit,signed R^2',
    type='n', main = paste('Soft Threshold seletion'));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=1,col='red'); abline(h=sft_RsquaredCut,col='red');
    
    soft_threshold_plot <- recordPlot()
    dev.off() ## clean up device
    
    # softpower estimation ----------------------------------------------------
    softPower<-  sft$powerEstimate
    
    {allowWGCNAThreads(nThreads=nThreads)
    
    net = blockwiseModules(datExpr, power = softPower,
    TOMType = TOMType, networkType = networkType,
    corType = corFnc, minModuleSize = minModuleSize, deepSplit = deepSplit,
    mergeCutHeight = mergeCutHeight, reassignThreshold = reassignThreshold,
    maxBlockSize = maxBlockSize, numericLabels = TRUE,
    saveTOMs = FALSE, verbose = 3)}
    disableWGCNAThreads()
    modules <- net$colors
    module.colours = labels2colors(net$colors)
    
    # modules -----------------------------------------------------------------
    
    plotDendroAndColors(net$dendrograms[[1]], 
    module.colours[net$blockGenes[[1]]],
    "Module colors", main = "Gene dendrogram and module colors",
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    
    net_dendrogram <- recordPlot()
    dev.off() ## clean up device
    
    # modules and plot.phylo --------------------------------------------------
    
    # calculate eigengenes
    MEs = net$MEs
    rownames(MEs)<-rownames(datExpr)
    
    #calculate dissimilarity of module eigengenes
    allowWGCNAThreads(nThreads=nThreads)
    MEDiss = 1-cor(MEs);
    disableWGCNAThreads()
    
    #cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = 'average');
    
    # ME plots ----------------------------------------------------------------
    plot(METree, main = "Clustering of module eigengenes", xlab = "", 
    sub = "", cex = 0.8)
    
    module_clustering <- recordPlot()
    dev.off() ## clean up device
    
    # network input parameters report
    Input_Parameters<-list(networkType=networkType,
    TOMType=TOMType, corFnc=corFnc,
    sft_RsquaredCut=sft_RsquaredCut,
    softPower=softPower,
    minModuleSize=minModuleSize,
    deepSplit=deepSplit,
    mergeCutHeight=mergeCutHeight,
    reassignThreshold=reassignThreshold,
    maxBlockSize=maxBlockSize)
    
    # Network Output
    Network_Output<-list(softPower=softPower,
    ref_genes=ref_genes,
    modules=modules,
    MEs=MEs,
    datExpr=datExpr)
    
    # Output plots
    
    Output_plots<-list(soft_threshold_plot=soft_threshold_plot,
    module_clustering=module_clustering,
    net_dendrogram=net_dendrogram)
    
    Final_Output<-list(Input_Parameters=Input_Parameters,
    Network_Output=Network_Output,
    Output_plots=Output_plots)
    
    message("Table for modules and gene counts")
    print(table(Final_Output$Network_Output$modules))
    return(Final_Output)
    
    }
    
    #' Scales and centers data by sample/row in preparation for discretization
    #' @param File the name of the text file generated by xCell that 
    #' contains the
    #' cell signature scores.
    #' @return xCell signatures scaled and centered by sample. For GMIC,
    #' ImmuneScore, StromaScore, and MicroenvironmentScore are removed.
    #' @importFrom utils read.delim
    #' @examples file_dir<-system.file("extdata", "IRIS_xCell_sig.txt", 
    #' package = "GmicR", mustWork = TRUE)
    #' Xcell_sig<-xCell_loader(file_dir)
    #' @export
    
    xCell_loader<-function(File=NULL){
    
    # Function to check if file directory is defined --------------------------
    is.not.null <- function(x) !is.null(x)
    # loading file ------------------------------------------------------------
    
    if(is.not.null(File)){
    xCell_Sigs <- data.frame(read.delim(File),check.names = FALSE,
    stringsAsFactors = FALSE)
    
    # setting rownames
    rownames(xCell_Sigs)<-xCell_Sigs$X
    xCell_Sigs$X<-NULL
    
    # removing ImmuneScore and stromal scores
    clean_xCell_Sigs<-xCell_Sigs[seq(1,64),] 
    message("ImmuneScore, StromaScore, and MicroenvironmentScore removed")
    
    # Scaling
    Zclean_xCell_Sigs<-data.frame(scale(clean_xCell_Sigs), check.names = FALSE)
    xCell_df<-data.frame(t(Zclean_xCell_Sigs), check.names = FALSE)
    }
    
    # fixing names for R
    colnames(xCell_df)<-gsub("+", "", colnames(xCell_df), fixed = TRUE)
    colnames(xCell_df)<-gsub(" ", "_", colnames(xCell_df), fixed = TRUE)
    colnames(xCell_df)<-gsub("-", "", colnames(xCell_df), fixed = TRUE)
    
    return(xCell_df)
    }
    
    
    #' Discretizes biological assay data in preparation for bayensian network 
    #' learning
    #' @param Auto_WGCNA_OUTPUT R object generated from Auto_WGCNA function. 
    #' @param Remove_ME0 a logical value. If FALSE (default), ME0 is not 
    #' removed.
    #' If TRUE the eigengene for module 0 is removed prior to analysis.
    #' @param xCell_Signatures the name of the text file generated by 
    #' xCell that 
    #' contains the cell signature scores. If NULL (default) the only module 
    #' eigenegnes will be processed. If not NULL and if 
    #' Auto_WGCNA_OUTPUT is NULL,
    #' cell signature scores will be discretized.
    #' @param ibreaks an integer that indicates the number of ibreaks 
    #' used for 
    #' discretization.
    #' The default value is 60.
    #' @param Numeric_Pheno_scores a data.frame with rows indicating sample ID 
    #' and columns representing additional phenotype data to be included in 
    #' BN learning. If NULL (default) no data will be included. If provided, the
    #' data.frame will be merged with MEs and discretized into three levels.
    #' @return a list containing a data.frame with 
    #' module eigenegnes merged with Xcell signature scores and 
    #' discretized into
    #' three levels: L, M, H. If Auto_WGCNA_OUTPUT is NULL, both scaled and 
    #' discretized cell signatures will be return.
    #' @note Please verify that the sample name formatting is 
    #' consistent between 
    #' both datasets. Rownames in the module eigengenes data.frame 
    #' and the column 
    #' names of xCell signatures scores text file are matched for 
    #' merging. Only 
    #' samples that are present in both will be processed!
    #' @examples file_dir<-system.file("extdata", "IRIS_xCell_sig.txt", 
    #' package = "GmicR", mustWork = TRUE)
    #' Disc_Xcell_sig<-Data_Prep(xCell_Signatures=file_dir, ibreaks = 10)
    #' Disc_Xcell_sig$disc_data
    #' @export
    
    Data_Prep<-function(Auto_WGCNA_OUTPUT=NULL, Remove_ME0=FALSE, 
    Numeric_Pheno_scores = NULL,
    xCell_Signatures=NULL,
    ibreaks = 60){
    
    
    if(is.null(Numeric_Pheno_scores)){
    MEs<-Auto_WGCNA_OUTPUT$Network_Output$MEs
    }else if(!is.null(Numeric_Pheno_scores)){
    MEs<-Auto_WGCNA_OUTPUT$Network_Output$MEs
    merged_dat<-merge(Numeric_Pheno_scores, MEs, by = "row.names")
    rownames(merged_dat)<-merged_dat$Row.names
    merged_dat$Row.names<-NULL
    MEs<-merged_dat
    }
    
    
    
    if(!is.null(Auto_WGCNA_OUTPUT)){
    MEs<-MEs
    # removing ME0 ------------------------------------------------------------
    if(isTRUE(Remove_ME0)){
    MEs$ME0<-NULL
    }
    
    # checking for cell sigs --------------------------------------------------
    if(is.null(xCell_Signatures)){
    disc_data<-discretize(MEs, method = "hartemink", breaks = 3,
    ibreaks = 60, idisc = "quantile")
    
    # setting levels
    for (i in names(disc_data))
    levels(disc_data[, i]) = c("L","M", "H")
    } else {
    
    xCell_Signatures<-xCell_loader(xCell_Signatures)
    merged_data<-merge(MEs, xCell_Signatures, by = "row.names",
    all = FALSE)
    
    rownames(merged_data)<-merged_data$Row.names
    merged_data$Row.names<-NULL
    
    # discretizing data  ----------------------------------------------
    disc_data<-discretize(merged_data, method = "hartemink", breaks = 3,
    ibreaks = ibreaks, idisc = "quantile")
    for (i in names(disc_data))
    levels(disc_data[, i]) = c("L","M", "H")
    }
    
    rownames(disc_data)<-rownames(MEs)
    Auto_WGCNA_OUTPUT$disc_data<-disc_data
    return(Auto_WGCNA_OUTPUT)
    } else if(is.null(Auto_WGCNA_OUTPUT)&!is.null(xCell_Signatures)){
    message("discretizing cell signature data")
    xCell_Signatures<-xCell_loader(xCell_Signatures)
    
    # discretizing data  ----------------------------------------------
    disc_data<-discretize(xCell_Signatures, method = "hartemink", breaks = 3,
    ibreaks = ibreaks, idisc = "quantile")
    for (i in names(disc_data))
    levels(disc_data[, i]) = c("L","M", "H")
    
    xCell_Signatures_output<-list(xCell_Signatures=xCell_Signatures,
    disc_data=disc_data)
    return(xCell_Signatures_output)
    }
    
    }
    
