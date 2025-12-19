#!/opt/conda/bin/R

# Usage: Rscript Integration.R --outdir <outdir_integration> --diffTFexpr_files <diffTFexpr_files> \
    # --diffTFbinding_file <diffTFbinding_file> --diffTFactivity_files <diffTFactivity_files> \
    # --monalisa_files <monalisa_files> --group_order <group_order> \
    # --human_mouse_symbol_ortholog <human_mouse_symbol_ortholog_file> \
    # --ggvenn_script <ggvenn_script>

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
library(tidyr)
# library(OmnipathR)
})


outdir_integration = opt$outdir
TFexprFiles = opt$diffTFexpr_files %>% stringr::str_split(.,";") %>% .[[1]]
TOBIAS_file = opt$diffTFbinding_file
diffTFactivityFiles = opt$diffTFactivity_files %>% stringr::str_split(.,";") %>% .[[1]]
monalisa_files = opt$monalisa_files %>% stringr::str_split(., ";") %>% .[[1]]
Group_order = opt$group_order %>% stringr::str_split(.,",") %>% .[[1]]
human_mouse_symbol_ortholog_file = opt$human_mouse_symbol_ortholog
ggvenn_script = opt$ggvenn_script


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

GetVennParts = function(vennlist) {
    if (length(vennlist)==3) {
        A = vennlist[[1]]
        B = vennlist[[2]]
        C = vennlist[[3]]
        A_only = setdiff(A, union(B, C))
        B_only = setdiff(B, union(A, C))
        C_only = setdiff(C, union(A, B))
        A_B_only = setdiff(intersect(A, B), C)
        A_C_only = setdiff(intersect(A, C), B)
        B_C_only = setdiff(intersect(B, C), A)
        A_B_C = Reduce(intersect, list(A, B, C))

        venn_parts = list(
            AAAAAAA_only   = A_only,
            BBBBBBB_only   = B_only,
            CCCCCCC_only   = C_only,
            AAAAAAA_BBBBBBB      = A_B_only,
            AAAAAAA_CCCCCCC      = A_C_only,
            BBBBBBB_CCCCCCC      = B_C_only,
            AAAAAAA_BBBBBBB_CCCCCCC    = A_B_C)
        names(venn_parts) =
            gsub("AAAAAAA", names(vennlist)[1], names(venn_parts)) %>%
            gsub("BBBBBBB", names(vennlist)[2], .) %>%
            gsub("CCCCCCC", names(vennlist)[3], .)
    }
    if (length(vennlist)==4) {
        A = vennlist[[1]]
        B = vennlist[[2]]
        C = vennlist[[3]]
        D = vennlist[[4]]
        A_only = setdiff(A, Reduce(union, list(D, B, C)))
        B_only = setdiff(B, Reduce(union, list(A, D, C)))
        C_only = setdiff(C, Reduce(union, list(A, B, D)))
        D_only = setdiff(D, Reduce(union, list(A, B, C)))
        A_B_only = setdiff(intersect(A, B), union(C,D))
        A_C_only = setdiff(intersect(A, C), union(B,D))
        A_D_only = setdiff(intersect(A, D), union(B,C))
        B_C_only = setdiff(intersect(B, C), union(A,D))
        B_D_only = setdiff(intersect(B, D), union(A,C))
        C_D_only = setdiff(intersect(C, D), union(A,B))
        A_B_C_only = setdiff(Reduce(intersect,list(A,B,C)), D)
        A_B_D_only = setdiff(Reduce(intersect,list(A,B,D)), C)
        A_C_D_only = setdiff(Reduce(intersect,list(A,C,D)), B)
        B_C_D_only = setdiff(Reduce(intersect,list(B,C,D)), A)
        A_B_C_D = Reduce(intersect,list(A,B,C,D))

        venn_parts = list(
            AAAAAAA_only   = A_only,
            BBBBBBB_only   = B_only,
            CCCCCCC_only   = C_only,
            DDDDDDD_only   = D_only,
            AAAAAAA_BBBBBBB_only      = A_B_only,
            AAAAAAA_CCCCCCC_only      = A_C_only,
            AAAAAAA_DDDDDDD_only      = A_D_only,
            BBBBBBB_CCCCCCC_only      = B_C_only,
            BBBBBBB_DDDDDDD_only      = B_D_only,
            CCCCCCC_DDDDDDD_only      = C_D_only,
            AAAAAAA_BBBBBBB_CCCCCCC_only    = A_B_C_only,
            AAAAAAA_BBBBBBB_DDDDDDD_only    = A_B_D_only,
            AAAAAAA_CCCCCCC_DDDDDDD_only    = A_C_D_only,
            BBBBBBB_CCCCCCC_DDDDDDD_only    = B_C_D_only,
            AAAAAAA_BBBBBBB_CCCCCCC_DDDDDDD    = A_B_C_D)
        names(venn_parts) =
            gsub("AAAAAAA", names(vennlist)[1], names(venn_parts)) %>%
            gsub("BBBBBBB", names(vennlist)[2], .) %>%
            gsub("CCCCCCC", names(vennlist)[3], .) %>%
            gsub("DDDDDDD", names(vennlist)[4], .)
    }
    return(venn_parts)
}

getVennDiagram = function(
    OUT_DIR=paste0(outdir_integration, "/diffTFexpr_RNA/"), 
    up_list=degs_up_list, 
    dn_list=degs_dn_list) {
    dir.create(OUT_DIR, showWarnings=F, recursive=TRUE)
    source(ggvenn_script)
    pdf(paste0(OUT_DIR, "/VenDiagram.pdf"), h=6,w=10)
        venn1_grob = ggvenn(up_list, , show_percentage = FALSE)
        venn2_grob = ggvenn(dn_list, , show_percentage = FALSE)
        print(cowplot::plot_grid(venn1_grob, venn2_grob, 
            ncol = 2, labels=c("          Upregulated TFs", "       Downregulated TFs")))
    dev.off()

    VennParts_degs_up_list = GetVennParts(up_list)
    VennParts_degs_dn_list = GetVennParts(dn_list)
    out = c(VennParts_degs_up_list, VennParts_degs_dn_list) %>%
        setNames(c(
            sapply(names(VennParts_degs_up_list), function(x) {paste0(x, "_Up")}),
            sapply(names(VennParts_degs_dn_list), function(x) {paste0(x, "_Down")})
        )) %>%
        lapply(function(x) {paste0(x, collapse = ";")}) %>%
        as.data.frame() %>% t() %>% as.data.frame() %>% 
        tibble::rownames_to_column(var="Comparison_Direction") %>%
        setNames(c("Comparison_Direction", "TFs"))
    write.table(out,
        file=paste0(OUT_DIR, "/VenDiagramInput.tsv"), 
        row.names=FALSE, sep="\t")
} 

convertMouse2HumanGeneSymbol = function(gene_list, human_mouse_symbol_map) {
    library(stringr)
    library(data.table)
    # Example input vector
    names_vec = gene_list
    # 1. Split all at once (list of vectors)
    split_names = str_split(names_vec, pattern = fixed("::"))
    # 2. Flatten all symbols into one vector with indexes to map back
    all_symbols = unlist(split_names)
    group_ids = rep(seq_along(split_names), lengths(split_names))
    # 3. Vectorized check: which symbols are not in human SYMBOL keys
    human_symbols = keys(org.Hs.eg.db, keytype = "SYMBOL")
    not_in_human = !(all_symbols %in% human_symbols)
    # 4. Vectorized map mouse symbols to human using your mapping table/vector
    # Assume human_mouse_symbol_map is a named vector: names are mouse symbols, values human symbols
    all_symbols2 = all_symbols
    all_symbols2[not_in_human] = human_mouse_symbol_map[all_symbols2[not_in_human]]
    all_symbols2[is.na(all_symbols2)] = all_symbols[is.na(all_symbols2)] 
    # 5. Recombine symbols back into the original groups
    dt = data.table(group = group_ids, symbol = all_symbols)
    result = dt[, paste(symbol, collapse = "::"), by = group][order(group), V1]
    # 6. Convert to character vector
    return(as.character(result))
}

addHeterodimer = function(df, heterodimerTFs_flatten, heterodimerTFs, TF_col="Feature", effect_col="log2FoldChange") {
    temp = df[df[[TF_col]] %in% heterodimerTFs_flatten, TF_col] %>% 
        as.character()
    if (length(temp)>0) {
        temp2 = sapply(heterodimerTFs, function(h){
            sapply(temp, function(t) grepl(t,h)) %>% sum()
        })
        temp2 = names(temp2[temp2==2])
        if (length(temp2)>0) {
            temp3 = lapply(temp2, function(t) {
                tl = strsplit(t, "::")[[1]]
                # 1. Extract log2FoldChange for ETV7 and NFE2L3
                l1 = df[[effect_col]][df[[TF_col]] == tl[1]]
                l2 = df[[effect_col]][df[[TF_col]] == tl[2]]
                # 2. Calculate mean log2FoldChange
                mean_log2FC = mean(c(l1, l2))
                # 3. Create a new row with Feature "ETV7::NFE2L3" and log2FoldChange = mean_log2FC
                new_row = df[1,,drop=F]
                new_row[[TF_col]] = t
                new_row[[effect_col]] = mean_log2FC
                new_row
            }) %>% Reduce(rbind, .)
            df = rbind(df, temp3)
        }
    }
    return(df)
}



# singularity shell /proj/nobackup/sens2025644/wharf/ekolyal/ekolyal-sens2025644/nf-atacco/conts/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-deseq2_pruned_a1ceaff5bb99ce9d.sif
# outdir_integration = "/proj/sens2025644/TF/results/Integration/"
# dir.create(outdir_integration, showWarnings = FALSE, recursive = TRUE)
# files = list.files("/proj/sens2025644/nbis8271/results/RNAseq/DEGs/", full.names = TRUE)
# # Now filter out anything matching Males or Females
# TFexprFiles = files[!grepl("_M_|_F_|Males|Females", files)]
# outdir_TOBIAS = "/proj/sens2025644/TF/results/TOBIAS_TFfootprinting/BINDetect/DiffTFBinding/"
# TOBIAS_file = paste0(outdir_TOBIAS, "bindetect_results.xlsx")
# dirs = list.dirs("/proj/sens2025644/TF/results/diffTFactivity_RNA/", 
#         full.names = TRUE) %>% setdiff("/proj/sens2025644/TF/results/diffTFactivity_RNA/")
# sapply(dirs, function(dir) {
#     list.files(dir, pattern="TFactivity_diff_.*.tsv", 
#     full.names = TRUE)
# }) %>% as.character() -> diffTFactivityFiles
# monalisa_files = "/proj/sens2025644/TF/results/MonaLisa/HFpEF_vs_Control/result_files/binnedEnr.RDS;/proj/sens2025644/TF/results/MonaLisa/HFpEF_vs_HFrEF/result_files/binnedEnr.RDS;/proj/sens2025644/TF/results/MonaLisa/HFrEF_vs_Control/result_files/binnedEnr.RDS" %>%
#    stringr::str_split(., ";") %>% .[[1]]
# human_mouse_symbol_ortholog_file = "/proj/nobackup/sens2025644/wharf/ekolyal/ekolyal-sens2025644/nf-atacco/data/human_mouse_symbol_ortholog.csv"
# Group_order = c("Control", "HFrEF", "HFpEF") # Reference first


print("# Integration TF inference across methods and comparisons")



print("## Venn diagram across comparisons for each method")
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------


print("### Venn diagram for differential TF binding from TOBIAS across comparisons")
#-------------------------------------------------------------------------------------
TOBIAS_rslt = readxl::read_excel(TOBIAS_file)
print("#### But first let's convert mouse gene symbols to human symbols, 
    TOBIAS and Monalisa use TF motif data from differet species, 
    which are ususally similar across species...")
human_mouse_symbol_ortholog = 
    read.table(human_mouse_symbol_ortholog_file, header=T, sep="\t") 
human_mouse_symbol_map = 
    setNames(human_mouse_symbol_ortholog$X9606, human_mouse_symbol_ortholog$X10090)
TOBIAS_rslt$name = convertMouse2HumanGeneSymbol(
    gene_list = TOBIAS_rslt$name, 
    human_mouse_symbol_map = human_mouse_symbol_map)
pairs = combn(Group_order, 2, simplify = FALSE)
TOBIAS_rslt_list = lapply(pairs, function(gp) {
    Pair = paste0(gp[2], "vs", gp[1])
    cols = colnames(TOBIAS_rslt)[!grepl(setdiff(Group_order, gp), colnames(TOBIAS_rslt))] 
    TOBIAS_rslt[, cols] %>% 
        setNames(gsub(".*_pvalue", "pvalue", cols)) %>%
        setNames(gsub(".*_change", "change", colnames(.))) %>%
        na.omit() %>% mutate(pvalue=as.numeric(pvalue)) %>%
        mutate(padj = p.adjust(pvalue, method="fdr"))
}) %>% setNames(sapply(pairs, function(gp) {paste0(gp[2], "vs", gp[1])}))

TOBIAS_diffTFbinding_up_table_list = lapply(names(TOBIAS_rslt_list), function(Pair) {
    TOBIAS_rslt_list[[Pair]] %>%
        filter(padj<0.05 & change>0) 
}) %>% setNames(names(TOBIAS_rslt_list))
TOBIAS_diffTFbinding_up_list = lapply(TOBIAS_diffTFbinding_up_table_list, function(x) {
    x %>% pull(name) %>% as.character()
})
TOBIAS_diffTFbinding_dn_table_list = lapply(names(TOBIAS_rslt_list), function(Pair) {
    TOBIAS_rslt_list[[Pair]] %>%
        filter(padj<0.05 & change<0)
}) %>% setNames(names(TOBIAS_rslt_list))
TOBIAS_diffTFbinding_dn_list = lapply(TOBIAS_diffTFbinding_dn_table_list, function(x) {
    x %>% pull(name) %>% as.character()
})

getVennDiagram(
    OUT_DIR = paste0(outdir_integration, "/diffTFbinding_ATAC/"), 
    up_list = TOBIAS_diffTFbinding_up_list, 
    dn_list = TOBIAS_diffTFbinding_dn_list)



print("### Venn diagram for TF motif enrichments from Monalisa across comparisons")
#-------------------------------------------------------------------------------------
print("We will only keep motifs with pvalue < 0.005 in at least one of top 2 negative logFC or top 2 positive bins")
monalisa_table_list = lapply(monalisa_files, function(file) {
    monalis_out = readRDS(file)
    temp = lapply(names(assays(monalis_out)), function(col) {
        temp = assay(monalis_out, col) %>% as.data.frame() %>%
            tibble::rownames_to_column(var="Motif") %>%
            gather(key="Bin", value=col, -Motif)
        colnames(temp)[3] = col
        temp
    }) %>% Reduce(merge,.) %>%
    as.data.frame() %>%
    # filter(negLog10Padj > -log10(0.05)) %>% dim()
    filter(negLog10P > -log10(0.005)) 

    lower_bounds = as.numeric(sub("^.?(.*?),.*$", "\\1", unique(temp$Bin)))
    # Order bins by the extracted lower bound
    extreme_bins_down = 
        unique(temp$Bin)[order(lower_bounds)] %>% 
        .[c(1,2)]
    extreme_bins_up = 
        unique(temp$Bin)[order(lower_bounds)] %>% 
        .[c(length(unique(temp$Bin))-1,length(unique(temp$Bin)))] 
    dn = filter(temp, Bin %in% c(extreme_bins_down)) %>% 
        filter(log2enr>0) %>% 
        group_by(Motif) %>%
        summarise(meanlog2enr = -mean(log2enr)) %>%
        mutate(TF = rowData(monalis_out)[Motif,"motif.name"]) %>%
        mutate(TF = convertMouse2HumanGeneSymbol(
                gene_list = TF, 
                human_mouse_symbol_map = human_mouse_symbol_map))
    up = filter(temp, Bin %in% c(extreme_bins_up)) %>% 
        filter(log2enr>0) %>% 
        group_by(Motif) %>%
        summarise(meanlog2enr = mean(log2enr)) %>%
        mutate(TF = rowData(monalis_out)[Motif,"motif.name"]) %>%
        mutate(TF = convertMouse2HumanGeneSymbol(
                gene_list = TF, 
                human_mouse_symbol_map = human_mouse_symbol_map))

    list(dn = dn, up = up)
}) %>% setNames(
    sapply(monalisa_files, function(file) {
        stringr::str_extract(file, "[^/]+_vs_[^/]+")}) %>%
        as.character() %>% gsub("_","",.))

monalisa_table_list_up = lapply(monalisa_table_list, function(table) table$up)
monalisa_table_list_dn = lapply(monalisa_table_list, function(table) table$dn)
monalisa_list_up = lapply(monalisa_table_list, 
    function(table) table$up %>% pull(TF) %>% unique()) 
monalisa_list_dn = lapply(monalisa_table_list, 
    function(table) table$dn %>% pull(TF) %>% unique())

paste0(outdir_integration, "/TFmotifEnr_Monalisa/") %>%
    dir.create(., recursive=T, showWarnings=FALSE)
getVennDiagram(
    OUT_DIR = paste0(outdir_integration, "/TFmotifEnr_Monalisa/"), 
    up_list = monalisa_list_up, 
    dn_list = monalisa_list_dn)



print("### Venn diagram for differential TF expression from RNA-seq across comparisons")
#-------------------------------------------------------------------------------------
TF_list = c(dorothea::dorothea_hs %>% pull(tf), TOBIAS_rslt$name) %>% unique()
heterodimerTFs = TF_list %>% grep("::",.,value=T)
heterodimerTFs_flatten = 
    sapply(heterodimerTFs, function(x) {strsplit(x,"::")}) %>% unlist() %>% as.character()
degs_up_table_list = lapply(TFexprFiles, function(file) {
    df = read.table(file, header=TRUE, sep=",") %>%
        as.data.frame() %>%
        na.omit() %>% mutate(Ensemble = X) %>%
        dplyr::select(-X) %>%
        mutate(Feature = mapIds(org.Hs.eg.db,
            keys = Ensemble,
            column = "SYMBOL",
            keytype = "ENSEMBL",
            multiVals = "first")) %>%
        filter(Feature %in% TF_list) %>%
        filter(pvalue<0.005&log2FoldChange>0) 
    addHeterodimer(df=df, heterodimerTFs_flatten, heterodimerTFs) 
}) %>% setNames(sapply(TFexprFiles, function(file) {
    gsub(".csv|_", "", gsub("DEGs", "", basename(file)))}))

degs_up_list = lapply(degs_up_table_list, function(x){
    x %>% pull(Feature) %>% as.character()
})
#degs_up_list = degs_up_list[!grepl("::", degs_up_list)]

degs_dn_table_list = lapply(TFexprFiles, function(file) {
    df = read.table(file, header=TRUE, sep=",") %>%
        as.data.frame() %>%
        na.omit() %>% mutate(Ensemble = X) %>%
        dplyr::select(-X) %>%
        mutate(Feature = mapIds(org.Hs.eg.db,
            keys = Ensemble,
            column = "SYMBOL",
            keytype = "ENSEMBL",
            multiVals = "first")) %>%
        filter(Feature %in% TF_list) %>%
        filter(pvalue<0.005&log2FoldChange<0) 
    addHeterodimer(df=df, heterodimerTFs_flatten, heterodimerTFs) 
}) %>% setNames(sapply(TFexprFiles, function(file) {
    gsub(".csv|_", "", gsub("DEGs", "", basename(file)))}))

degs_dn_list = lapply(degs_dn_table_list, function(x){
    x %>% pull(Feature) %>% as.character()
})
getVennDiagram(
    OUT_DIR=paste0(outdir_integration, "/diffTFexpr_RNA/"), 
    up_list=degs_up_list, 
    dn_list=degs_dn_list)


print("### Venn diagram for differential TF activities from decoupleR across comparisons")
#-------------------------------------------------------------------------------------
decoupleR_diffTFactivity_up_table_list = 
    lapply(diffTFactivityFiles, function(file) {
        df = read.table(file, header=TRUE, sep="\t") %>%
            as.data.frame() %>%
            na.omit() %>%
            filter(P.Value<0.05 & logFC>0)
        addHeterodimer(df=df, 
            heterodimerTFs_flatten, heterodimerTFs,
            TF_col="TF", effect_col="logFC") 
}) %>% setNames(sapply(diffTFactivityFiles, function(file) {
    gsub(".tsv|TFactivity_diff_", "", basename(file))}))  
decoupleR_diffTFactivity_up_list = 
    lapply(decoupleR_diffTFactivity_up_table_list, function(x) {
    x %>% pull(TF) %>% as.character()
})
decoupleR_diffTFactivity_dn_table_list = 
    lapply(diffTFactivityFiles, function(file) {
        df = read.table(file, header=TRUE, sep="\t") %>%
            as.data.frame() %>%
            na.omit() %>%
            filter(P.Value<0.05 & logFC<0) 
        addHeterodimer(df=df, 
            heterodimerTFs_flatten, heterodimerTFs,
            TF_col="TF", effect_col="logFC") 
}) %>% setNames(sapply(diffTFactivityFiles, function(file) {
    gsub(".tsv|TFactivity_diff_", "", basename(file))}))   
decoupleR_diffTFactivity_dn_list = 
    lapply(decoupleR_diffTFactivity_dn_table_list, function(x) {
    x %>% pull(TF) %>% as.character()
})

getVennDiagram(
    OUT_DIR = paste0(outdir_integration, "/diffTFactivity_RNA_Dorothea_VIPER/"), 
    up_list = decoupleR_diffTFactivity_up_list, 
    dn_list = decoupleR_diffTFactivity_dn_list)





print("## Venn diagram across methods per comparison")
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
print("Venn diagram of differential TF activities from decoupleR vs differential TF binding 
        from TOBIAS vs differential TF expression from RNA-seq, per comparison")
lapply(names(degs_up_list), function(pair){
    TFs_degs_up = degs_up_list[[pair]]
    TFs_degs_dn = degs_dn_list[[pair]]
    TOBIAS_diffTFbinding_up = TOBIAS_diffTFbinding_up_list[[pair]]
    TOBIAS_diffTFbinding_dn = TOBIAS_diffTFbinding_dn_list[[pair]]
    decoupleR_diffTFactivity_up = decoupleR_diffTFactivity_up_list[[pair]]
    decoupleR_diffTFactivity_dn = decoupleR_diffTFactivity_dn_list[[pair]]
    monalisa_up = monalisa_list_up[[pair]]
    monalisa_dn = monalisa_list_dn[[pair]]       
    source(ggvenn_script)
    up_list = list(
                diffTFexpr_RNA = TFs_degs_up,
                diffTFbinding_ATAC = TOBIAS_diffTFbinding_up,
                diffTFactivity_RNA = decoupleR_diffTFactivity_up,
                TFmotifEnr_ATAC = monalisa_up
            )
    names(up_list) = gsub("_.*", "", names(up_list))
    dn_list = list(
                diffTFexpr_RNA = TFs_degs_dn,
                diffTFbinding_ATAC = TOBIAS_diffTFbinding_dn,
                diffTFactivity_RNA = decoupleR_diffTFactivity_dn,
                TFmotifEnr_ATAC = monalisa_dn
            )
    names(dn_list) = gsub("_.*", "", names(dn_list))
    dir.create(paste0(outdir_integration, "/VennDiagram_crossMethods/"), 
        showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(outdir_integration, 
        "/VennDiagram_crossMethods/VennDiagram_TFs_across_methods_", pair, ".pdf"), h=6,w=12)
        venn1_grob = ggvenn(
            up_list, show_percentage = FALSE)
        venn2_grob = ggvenn(
            dn_list, show_percentage = FALSE)
        print(cowplot::plot_grid(venn1_grob, venn2_grob, 
            ncol = 2, labels=c("          Upregulated TFs", "       Downregulated TFs")))
    dev.off()

    VennParts_across_methods_up = GetVennParts(up_list)
    VennParts_across_methods_dn = GetVennParts(dn_list)
    out = c(VennParts_across_methods_up, VennParts_across_methods_dn) %>%
        setNames(c(
            sapply(names(VennParts_across_methods_up), function(x) {paste0(x, "_Up")}),
            sapply(names(VennParts_across_methods_dn), function(x) {paste0(x, "_Down")})
        )) %>%
        lapply(function(x) {paste0(x, collapse = ";")}) %>%
        as.data.frame() %>% t() %>% as.data.frame() %>% 
        tibble::rownames_to_column(var="Comparison_Direction") %>%
        setNames(c("Method_Direction", "TFs"))
    write.table(out,
        file=paste0(outdir_integration, 
            "/VennDiagram_crossMethods/TFs_cross_methods_VenDiagram", pair, ".tsv"), 
        row.names=FALSE, sep="\t")
})





print("## Scatter plots across methods per comparison")
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
print("Scatter plots of differential TF activities from decoupleR vs differential TF binding 
        from TOBIAS vs differential TF expression from RNA-seq, per comparison")
lapply(names(degs_up_table_list), function(pair){
    methods_list =list(
        diffTFexpr_RNA = 
            rbind(degs_up_table_list[[pair]], degs_dn_table_list[[pair]]) %>%
            dplyr::select(Feature, log2FoldChange) %>%
            setNames(c("TF", "log2FC_diffTFexpr_RNA")),
        diffTFbinding_ATAC = rbind(TOBIAS_diffTFbinding_up_table_list[[pair]],
            TOBIAS_diffTFbinding_dn_table_list[[pair]]) %>%
            dplyr::select(name, change) %>%
            setNames(c("TF", "change_diffTFbinding_ATAC")),
        diffTFactivity_RNA = rbind(decoupleR_diffTFactivity_up_table_list[[pair]],
            decoupleR_diffTFactivity_dn_table_list[[pair]]) %>%
            dplyr::select(TF, logFC) %>%
            setNames(c("TF", "logFC_diffTFactivity_RNA")),
        TFmotifEnr_ATAC = rbind(monalisa_table_list_up[[pair]],
            monalisa_table_list_dn[[pair]]) %>%
            dplyr::select(TF, meanlog2enr) %>%
            setNames(c("TF", "log2enr_TFmotifEnr_ATAC"))
    )
    # Merge tables and plot scatter plots
    for (i in 1:length(methods_list)) {
        if (i < length(methods_list)) {
            for (j in (i+1):length(methods_list)) {
                method1_name = names(methods_list)[i]
                method2_name = names(methods_list)[j]
                merged_table = methods_list[[i]] %>%
                    inner_join(methods_list[[j]], by="TF")
                merged_table$sync = ifelse(rowSums(merged_table[,2:3]>0)==1, "Discordant", "Concordant")
                p = ggplot(merged_table, 
                    aes_string(x=colnames(merged_table)[2],
                            y=colnames(merged_table)[3])) +
                    geom_point(aes(col=sync)) +
                    ggrepel::geom_text_repel(
                        data = merged_table, aes_string(colnames(merged_table)[2], colnames(merged_table)[3], 
                        label = "TF"), color="cornsilk4") +
                    theme_bw() +
                    ggtitle(paste0("Scatter Plot of TFs - ", pair, "\n",
                        method1_name, " vs ", method2_name)) +
                    scale_colour_manual(values = c("Discordant"="orange", "Concordant"="green4")) 
                dir.create(paste0(outdir_integration, "/ScatterPlot_crossMethods/"), 
                    showWarnings = FALSE, recursive = TRUE)
                pdf(paste0(outdir_integration, 
                    "/ScatterPlot_crossMethods/ScatterPlot_sigTFs_", 
                        method1_name, "_vs_", method2_name, "_", pair, ".pdf"), h=5,w=7)
                    print(p)
                dev.off()
            }
        }
    }
})


print("## Save session info...")
sessionInfo = capture.output(sessionInfo())
writeLines(sessionInfo, paste0(outdir_integration,"/sessionInfo_integration.txt"))
print("All done!")

## Integration of multiple TF activity inference methods using Robust Rank Aggregation (RRA), per comparison
# library(RobustRankAggreg)

# # Suppose you have 3 ranked lists of TFs:
# ranks_monalisa = c("TF1", "TF3", "TF2", "TF5", "TF4")
# ranks_tobias = c("TF3", "TF2", "TF1", "TF4", "TF6")
# ranks_decoupleR = c("TF2", "TF1", "TF4", "TF5", "TF7")

# rank_lists = list(monalisa = ranks_monalisa, tobias = ranks_tobias, decoupleR = ranks_decoupleR)

# # Run RRA
# rra_res = aggregateRanks(rank_lists)

# print(rra_res)


