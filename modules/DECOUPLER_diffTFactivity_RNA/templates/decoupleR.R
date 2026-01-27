#!/opt/conda/bin/R

# usage: Rscript decoupleR.R --dds <dds> --outdir <outdir_integration> --group_order <group_order> --currentCovariate <currentCovariate>

args = commandArgs(trailingOnly = TRUE)
# Simple argument parser
parse_args = function(args) {
  res = list()
  for (i in seq(1, length(args), by = 2)) {
    key = gsub("^--", "", args[i])
    res[[key]] = args[i + 1]
  }
  res
}
opt = parse_args(args)
dds = opt$dds
logNormCount  = opt$logNormCount
design = opt$design
outdir = opt$outdir
group_order = opt$group_order # comma separated string, e.g. "Control,HFrEF,HFpEF"
currentCovariate = opt$currentCovariate # e.g. "age"
dir.create(outdir, recursive=T, showWarning=F)

# singularity shell /proj/nobackup/sens2025644/wharf/ekolyal/ekolyal-sens2025644/nf-atacco/conts/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-deseq2_pruned_a1ceaff5bb99ce9d.sif
# outdir = "/proj/sens2025644/TF/results/Dorothea_VIPER_TF_RNA/"
# dds = "/proj/sens2025644/nbis8271/results/RNAseq/deseq2.dds.RData"
# group_order = "Control,HFrEF,HFpEF"


# Set up functions
VisualizeTopTFactivities_VIPER = function(
    tfs=NULL, n_tfs = 25, sample_acts=NULL, design, outdir, feature_meta=NULL,target_expression=NULL) {
    print("Plot heatmap of top variable TF activities across samples...")
    # Transform to wide matrix
    if (!is.null(sample_acts)) {
    sample_acts_mat = sample_acts %>%
        tidyr::pivot_wider(id_cols = 'condition',
        names_from = 'source',
        values_from = 'score') %>%
        tibble::column_to_rownames('condition') %>%
        as.matrix()
    } else {sample_acts_mat = t(target_expression)}
    design2 = design %>% filter(sample %in% rownames(sample_acts_mat))
    sample_acts_mat = sample_acts_mat[design2$sample,]
    if (!is.null(feature_meta)) {
        feature_meta = filter(feature_meta, target %in% colnames(sample_acts_mat))
        sample_acts_mat = sample_acts_mat[,feature_meta$target]
    }
    if (is.null(tfs)) {
        # Get top tfs with more variable means across samples
        tfs = sample_acts %>%
            dplyr::group_by(source) %>%
            dplyr::summarise(std = stats::sd(score)) %>%
            dplyr::arrange(-abs(std)) %>%
            head(n_tfs) %>%
            dplyr::pull(source)
        file = paste0(outdir,"/TF_activity_heatmap_topVariable.pdf")
    } else if (is.null(feature_meta)) {
        file = paste0(outdir,"/TF_activity_heatmap_topDifferential.pdf")
        } else {file = paste0(outdir,"/TF_target_heatmap_topDifferentialTargets.pdf")}
    sample_acts_mat = sample_acts_mat[,intersect(tfs, colnames(sample_acts_mat))]
    # Scale per gene
    sample_acts_mat = scale(sample_acts_mat)
    # Choose color palette
    colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
    colors.use = grDevices::colorRampPalette(colors = colors)(100)
    my_breaks = c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
                    seq(0.05, 2, length.out = floor(100 / 2)))
    # Plot
    if (!is.null(feature_meta)) {
        cluster_cols = FALSE
        ann_col = data.frame(
            TFs = feature_meta$TF,
            Targets = feature_meta$target
            ) %>% tibble::column_to_rownames(var="Targets") 
        ann_col$TFs = factor(ann_col$TFs, 
                levels=unique(ann_col$TFs))
        h = (nrow(sample_acts_mat) * 20)/100 + 5.5
        w = (length(tfs) * 20)/100 + 6
        main = "Top TF target expression across samples"
    } else {
        cluster_cols = TRUE
        ann_col = NA
        h = (nrow(sample_acts_mat) * 20)/100 + 3.5
        w = (length(tfs) * 20)/100 + 4
        main = "Top TF activities across samples"
    }
    ann_row = data.frame(
        groups = design2$condition,
        sample = design2$sample
    ) %>% tibble::column_to_rownames(var="sample") 
    ann_row$groups = factor(ann_row$groups, 
            levels=unique(ann_row$groups))
    pdf(file, h=h, w=w)
    print(pheatmap::pheatmap(mat = sample_acts_mat,
        annotation_row = ann_row,
        annotation_col = ann_col,
        cluster_rows = FALSE,
        cluster_cols = cluster_cols,
        color = colors.use,
        border_color = "white",
        breaks = my_breaks,
        cellwidth = 15,
        cellheight = 15,
        treeheight_row = 20,
        treeheight_col = 20,
        main = main))
    dev.off()
}

#' PCA plot function
#' @param data_matrix matrix with features as rows and samples as columns
#' @param HVG whether to use highly variable genes for PCA
#' @param topHVGnr number of top highly variable genes to use for PCA
#' @param meta metadata dataframe for samples
#' @param metaVal vector of metadata column names to color/shape the PCA plots
#' @param joinVariable variable name to join data_matrix and meta, usually sample names 
#' 
makePCAplot = 
  function(data_matrix, HVG=TRUE, topHVGnr=500, meta, metaVal, joinVariable="sample"){
      log_sd = apply(data_matrix, 1, sd)
      HVGs = names(sort(log_sd, T))[1:topHVGnr] # top 500 HVGs
      if (HVG) {data_matrix=data_matrix[HVGs,]}
      c_pca = prcomp(t(na.omit(data_matrix)),
                     scale=TRUE, center = TRUE)
      frac_var = function(x) x^2/sum(x^2)
      var_explained = c_pca$sdev %>% 
        tibble::as_tibble() %>% 
        frac_var() %>% 
        round(digits = 2)
      c_pca_all = c_pca$x %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = joinVariable) %>%
        left_join(meta, by=joinVariable)
      pcs = list(c("PC1","PC2"),c("PC1","PC3"),c("PC3","PC2"),c("PC2","PC4"))
      vars = list()
      for (i in metaVal) {
        vars=c(vars, list(c("Group", i)))}
      lapply(vars, function(V) {
        plots = paste0("p2_", 1:4)
        for (i in 1:4) {
          pc = pcs[[i]]
          pc_nr = gsub("PC", "", pc) %>% as.numeric()
          assign(plots[i], ggplot(aes_string(pc[1], pc[2]), data = c_pca_all) +
          geom_point(size = 3, aes_string(shape = V[1], col=V[2])) +
          labs(x = paste(pc[1], " (", var_explained[pc_nr[1],"value"]*100, "% explained var.)", sep = ""),
               y = paste(pc[2], " (", var_explained[pc_nr[2],"value"]*100, "% explained var.)", sep = "")))}
        print(cowplot::plot_grid(p2_1, p2_2, p2_3, p2_4, ncol=2))  
      })
      plot(var_explained$value, type = "l")
}


## Limma differential analysis
#' @param dat the dataset needs to have features as rownames
#' @param meta meta data needs to have at least two columns named "Samples", "Group", for the moment, we only support Group with two categories
#' 
limmaTest = function(dat, paired=TRUE, currentCovariate=NULL, 
                     checkCovariate=FALSE, meta=meta,
                     paired_variable="Patient", groups=c("Control", "HFpEF")) {
  meta = meta %>% filter(Group %in% groups) %>%
        filter(Samples%in%colnames(dat))
  # Make sure the samples are in the same order in both meta and dat:
  dat = dat[, meta$Samples]
  meta = meta %>% mutate(across(where(is.factor), droplevels))
  # print(dim(dat))
  #if (!is.null(currentCovariate)) {currentCovariate=NULL}
  # Only keep features that have at least 2 values in each group
    keep = sapply(1:nrow(dat), function(i) {
        row = dat[i, ]
        grp1 = row[meta$Group==groups[1]] %>% na.omit()
        grp2 = row[meta$Group==groups[2]] %>% na.omit()
        return((length(grp1)>=2) & (length(grp2)>=2))
    })
    dat = dat[keep, ]
    print(paste0("After filtering, ", nrow(dat), " features are kept for differential analysis."))
  baseFormula = ~ Group
  if (is.null(currentCovariate)) {
    currentFormula = baseFormula
  } else {
    currentFormula = as.formula(paste0("~ ", currentCovariate, " + Group"))}
  
  design = model.matrix(currentFormula, data=meta) # each column is for coefficient
  group.categories = meta$Group %>% unique() %>% length()
  colnames(design) = gsub("Group", "", colnames(design))
  
  if (paired) {
    # Fit with correlated arrays
    dupcor = duplicateCorrelation(dat, design, 
                            block=meta[,paired_variable,drop=T])
    fit = lmFit(dat, design, block=meta[,paired_variable,drop=T],
                correlation=dupcor$consensus)
  } else {
    fit = lmFit(dat, design)
  }
  # fit3 = eBayes(fit, robust=TRUE)
  # print("Check whether we should include a factor in the model by checking for how many metabolites this factor is significant in the model: ")
  if (checkCovariate) {
    fit.con_co = eBayes(fit)
    rlst_co = lapply(2:(ncol(design)-group.categories+1), function(i) {
        topTable(fit.con_co, n=Inf, coef=colnames(design)[i])}) %>%
        setNames(colnames(design)[2:(ncol(design)-group.categories+1)])
    # n refers to the number of top-ranked genes to be returned.
  } else {rlst_co = NA}
  
  fit.con = eBayes(fit, robust = TRUE, trend=TRUE)
  rlst_interest = topTable(fit.con, n=Inf, coef=levels(meta$Group)[2]) %>%
      arrange(-t) 
  
  return(list(rslt_of_interest=rlst_interest, rslt_co=rlst_co, data=dat, fitted.model=fit))
}


## Volcano plot
#' @param sig.cutoff4label significance (i.e. sig.statistic) value quantile cutoff for selecting labels to show, between 0-1.
#' 
makeVolcano = function(rslt_list, sig.statistic="adj.P.Val", 
                        # sig.cutoff=0.05, logFC.cutoff=0.3,
                        sig.cutoff4label=0.1, OUTDIR_pairDiff) {
    doVolcano = function(resi, NAME, sig.statistic="adj.P.Val", 
                        # sig.cutoff=0.05, logFC.cutoff=0.3,
                        sig.cutoff4label=0.1) {
        selectedCols = c("Feature", "logFC", "P.Value", "adj.P.Val", "t")
        toExport = resi[, selectedCols]
        toExport = toExport[order(toExport$P.Value),]

        tmp = as.data.frame(toExport) %>% filter(!is.na(adj.P.Val)) 
        library(ggrepel)
        # Categorize Results based on P-value & FDR for plotting
        tmp$sig = "NS"
        tmp$sig[tmp$P.Value < 0.05] = "P<0.05"
        tmp$sig[tmp$adj.P.Val < 0.05] = "FDR<0.05"
        tmp$sig = factor(tmp$sig,
            levels = c("NS", "P<0.05", "FDR<0.05"))

        tmpSignif = tmp %>%
            filter(tmp[[sig.statistic]]<quantile(abs(tmp[[sig.statistic]]),sig.cutoff4label)) %>%
            filter(sig!="NS")
        if (nrow(tmpSignif)==0) {
            tmpSignif = tmp %>% filter(tmp[[sig.statistic]]<quantile(abs(tmp[[sig.statistic]]),0.2)) %>%
                filter(sig!="NS")
        }
        p = ggplot(data = tmp, aes(logFC, -log10(P.Value), col = sig)) + 
            geom_point() + 
            xlab("Log2 Fold Change") +
            ggrepel::geom_text_repel(data = tmpSignif, aes(logFC, -log10(P.Value), 
                                            label = Feature), color="cornsilk4") +
            geom_hline(yintercept = -log10(1), lty = 2) +
            geom_vline(xintercept = 0, lty = 2) + 
            # scale_colour_manual(values = c("down"="blue", "nonSig"="black", "up"="red")) + 
            scale_color_manual(values = c(
                `FDR<0.05` = "purple",
                `P<0.05` = "orange2",
                `NS` = "gray"),
                guide = guide_legend(override.aes = list(size = 4))) +
            #theme(legend.position="none") +
            ggtitle(NAME)
        return(p)
    }

    invisible(capture.output(lapply(names(rslt_list), function(Pair) {
        OUT = paste0(OUTDIR_pairDiff, "/", Pair, "/")
        dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
        pdf(paste0(OUT, "VolcanoPlot_", Pair, "_TFactivity.pdf"), h=7, w=8.5)
        print(doVolcano(rslt_list[[Pair]] %>% as.data.frame() %>%
            tibble::rownames_to_column(var="Feature"), 
            NAME=Pair, sig.statistic=sig.statistic, 
            # sig.cutoff=sig.cutoff, logFC.cutoff=logFC.cutoff,
            sig.cutoff4label=sig.cutoff4label))
        dev.off()
    })))
}

inferTFactivity = function(outdir, logNormCount, design) {
    print("# Infer TF activities using Dorothea regulons...")    
    print("## Extract dorothea regulons for human (hsapiens) with confidence A, B, and C")
    conf_map = c(A = 1, B = 0.8, C = 0.6, D = 0.4, E = 0.2)
    net = dorothea::dorothea_hs %>%
        dplyr::filter(confidence %in% c("A", "B", "C")) %>%
        dplyr::mutate(
        likelihood = conf_map[confidence]) %>%
        dplyr::rename(
            source = tf,
            target = target
        ) %>%
        dplyr::select(source, target, mor)
    saveRDS(net, file=paste0(outdir,"/dorothea_regulons_ABC.RDS"))
    write.table(
        net,
        file=paste0(outdir,"/dorothea_regulons_ABC.tsv"), 
        row.names=FALSE, sep="\t")

    print("## Run VIPER to infer TF activities...")
    sample_acts = decoupleR::run_viper(
        mat = logNormCount,
        network = net
    )

    print("In the result table, the p_value indicates how statistically significant that inferred activity score is, 
    often based on permutations or null models used in VIPER to assess 
    if the activity is different from random expectation.")
    saveRDS(sample_acts, file=paste0(outdir,"/TFactivities_VIPER.RDS"))
    write.table(
        sample_acts,
        file=paste0(outdir,"/TFactivities_VIPER.tsv"), 
        row.names=FALSE, sep="\t")

    print("## Filtering: only use confidently inferred TF activities (p_value < 0.05) for downstream analysis...")
    sample_acts_filtered = sample_acts %>% filter(p_value < 0.05)
    saveRDS(sample_acts_filtered, file=paste0(outdir,"/TFactivities_filtered_VIPER.RDS"))
    write.table(
        sample_acts_filtered,
        file=paste0(outdir,"/TFactivities_filtered_VIPER.tsv"), 
        row.names=FALSE, sep="\t")
    write.table(
        net %>% filter(source %in% unique(sample_acts_filtered$source)) %>%
            filter(target %in% unique(rownames(logNormCount))) %>%
            as.data.frame() %>%
            mutate(TF=source) %>%
            dplyr::select(TF, target, mor),
        file=paste0(outdir,"/dorothea_regulons_ABC_usedbythisstudy.tsv"), 
        row.names=FALSE, sep="\t")

    print("## Visualize the most variable TF activities across samples")
    VisualizeTopTFactivities_VIPER(n_tfs = 25, 
        sample_acts=sample_acts_filtered, design=design, outdir=outdir)

    print("## PCA plot of samples based on TF activities")
    pdf(paste0(outdir,"/PCA_TFactivities.pdf"), h=6, w=9)
        dat4PCA = sample_acts_filtered %>%
            tidyr::pivot_wider(id_cols = 'source',
            names_from = 'condition',
            values_from = 'score') %>%
            tibble::column_to_rownames('source') %>%
            as.matrix()
        makePCAplot(data_matrix=dat4PCA, 
            HVG=FALSE,
            meta=mutate(design, Group=condition), 
            metaVal=c("sex", "age"), 
            joinVariable="sample")
        rm(dat4PCA)
    dev.off()

    return(list(sample_acts_filtered=sample_acts_filtered, net=net))
}

runLimmaDiffTFactivity = function(
    sample_acts_filtered, currentCovariate,
    design, logNormCount, outdir, group_order) {
    print("# Limma differential analysis of TF activities between conditions...")
    print("## Prepare data matrix for limma...")
    dat4limma = sample_acts_filtered %>%
        tidyr::pivot_wider(id_cols = 'source',
        names_from = 'condition',
        values_from = 'score') %>%
        tibble::column_to_rownames('source')
    write.table(
        dat4limma %>% as.data.frame() %>% 
            tibble::rownames_to_column(var="TF"), 
        file=paste0(outdir,"/TFactivities_matrix_for_limma.tsv"), 
        row.names=FALSE, sep="\t")

    print("## Perform pairwise differential analysis of TF activities between groups...")
    Group_order = stringr::str_split(group_order, ",")[[1]]
    # c("Control", "HFrEF", "HFpEF") # Reference first
    pairs = combn(Group_order, 2, simplify = FALSE)
    design$condition = factor(design$condition, levels=Group_order)
    tf_activity_diff_list = lapply(pairs, function(gp) {
        Pair = paste0(gp[2], "vs", gp[1])
        dir.create(paste0(outdir,"/", Pair), showWarnings=F)
        print(paste0("Comparing ", Pair, "..."))
        tf_activity_diff = limmaTest(dat=dat4limma, 
            currentCovariate=currentCovariate,
            paired=FALSE, 
            meta=mutate(design, Group=condition, Samples=sample),
            groups=gp)
        saveRDS(tf_activity_diff, 
            file=paste0(outdir,"/", Pair, "/TFactivity_diff_", Pair, ".RDS"))
        write.table(
            tf_activity_diff$rslt_of_interest %>%
                tibble::rownames_to_column(var="TF") %>%
                left_join(tf_activity_diff$data %>% 
                        tibble::rownames_to_column(var="TF"), by="TF"),
            file=paste0(outdir,"/", Pair, "/TFactivity_diff_", Pair, ".tsv"), 
            row.names=FALSE, sep="\t")   
        print("Make volcano plot of TF activities...")
        makeVolcano(
            rslt_list=
                list(tf_activity_diff$rslt_of_interest) %>% setNames(Pair), 
            sig.statistic="P.Value", 
            sig.cutoff4label=0.1, 
            OUTDIR_pairDiff=outdir)
        return(tf_activity_diff)
    }) %>% setNames(sapply(pairs, function(gp) {paste0(gp[2], "vs", gp[1])}))
    return(tf_activity_diff_list)
}

visualizeDiffTFactivity = function(
    tf_activity_diff_list, sample_acts_filtered, design, logNormCount, net, outdir) {
    print("## Visualize the differential TF activities (P.Value < 0.05) across samples")
    Tfs = sapply(tf_activity_diff_list, function(diff) {
        temp = diff$rslt_of_interest %>% filter(P.Value < 0.05) 
        rownames(temp)
    }) %>% unlist() %>% unique()
    VisualizeTopTFactivities_VIPER(tfs = Tfs, 
        sample_acts=sample_acts_filtered, 
        design=design, outdir=outdir)

    print("## Visualize the target expression of differential TFs (P.Value< 0.005) across samples")
    Tfs2 = sapply(tf_activity_diff_list, function(diff) {
        temp = diff$rslt_of_interest %>% filter(P.Value < 0.005) 
        rownames(temp)
    }) %>% unlist() %>% unique()
    Tfs2_targets_meta = net %>% 
        filter(source %in% unique(sample_acts_filtered$source)) %>%
        filter(target %in% unique(rownames(logNormCount))) %>%
        as.data.frame() %>%
        mutate(TF=source) %>%
        dplyr::select(TF, target, mor) %>% 
        filter(TF %in% Tfs2) %>%
        dplyr::group_by(TF)
    Tfs2_targets_meta = Tfs2_targets_meta[!duplicated(Tfs2_targets_meta$target),]
    logNormCount_Tfs2_targets = logNormCount[unique(Tfs2_targets_meta$target), ]
    VisualizeTopTFactivities_VIPER(
        tfs=pull(Tfs2_targets_meta,target) %>% unique(), 
        design=design, 
        outdir=outdir, 
        feature_meta=Tfs2_targets_meta,
        target_expression=as.data.frame(logNormCount_Tfs2_targets))
}


suppressPackageStartupMessages({
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(DESeq2)
library(org.Hs.eg.db)
library(dorothea)
library(limma)
# library(OmnipathR)
})


print("Infer TF activities using decoupleR (VIPER) based on RNA-seq data and Dorothea regulons...")
print("# Prepare input data for decoupleR...")
if (!is.null(dds)) {
    print(paste0("Processing dds file: ", dds))
    load(dds) # loads dds object 
    rownames(dds) = mapIds(org.Hs.eg.db,
                        keys = rownames(dds),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
    ### Remove rows with NA gene symbols
    dds = dds[!is.na(rownames(dds)), ]

    ### Log transformed normalized counts
    logNormCount = log2(counts(dds, normalized = TRUE)+1)
    ### Remove NAs and set row names
    logNormCount = logNormCount %>% as.data.frame() %>%
            dplyr::mutate_if(~ any(is.na(.x)), 
                            ~ dplyr::if_else(is.na(.x), 0, .x)) %>% 
            as.matrix()
    write.table(logNormCount %>% as.data.frame() %>% 
        tibble::rownames_to_column(var="TF"), 
        file=paste0(outdir,"/Normalized_counts_log2_matrix.tsv"), 
        row.names=FALSE, sep="\t")  

    ## Prepare design matrix
    design = colData(dds) %>% 
        as.data.frame() %>% tibble::rownames_to_column(var="sample")
    write.table(design, 
        file=paste0(outdir,"/Design_matrix.tsv"), 
        row.names=FALSE, sep="\t")
} else if (!is.null(logNormCount)&!is.null(design)) {
    logNormCount = read.table(logNormCount) %>%
        filter(!is.na(TF)) %>%
        tibble::column_to_rownames(var="TF") %>%
        as.matrix()
    copy.file(from=logNormCount, 
        to=paste0(outdir, "/Normalized_counts_log2_matrix.tsv"))
    design = read.table(design) 
    copy.file(from=design, 
        to=paste0(outdir, "/Design_matrix.tsv"))
} else {
    stop("Please provide either dds or logNormCount and design matrix!")
}

# Infer TF activity using VIPER with the Dorothea TF database
TFactivity_rslt = inferTFactivity(
    outdir, logNormCount, design)

# Perform differential TF activity analysis between groups ...
tf_activity_diff_list = runLimmaDiffTFactivity(
    sample_acts_filtered=TFactivity_rslt$sample_acts_filtered, 
    currentCovariate=currentCovariate,
    design, 
    logNormCount, 
    outdir, 
    group_order)

# Visualize differential TF activities...
visualizeDiffTFactivity(
    tf_activity_diff_list, 
    TFactivity_rslt$sample_acts_filtered, 
    design, logNormCount, 
    TFactivity_rslt$net, 
    outdir)

print("## Save session info...")
sessionInfo = capture.output(sessionInfo())
writeLines(sessionInfo, paste0(outdir,"/sessionInfo_decoupleR_TFactivitys.txt"))
print("All done!")
