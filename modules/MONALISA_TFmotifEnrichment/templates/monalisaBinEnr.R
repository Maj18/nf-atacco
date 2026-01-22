#!/opt/conda/bin/R

# Usage: Rscript monalisaBinEnr.R <difftable> <peakAnnotation> <pfm_file> <outdir> <coreNum>

args = commandArgs(trailingOnly = TRUE)
difftables = args[1] # ; separated if multiple files
peakAnnotation = args[2]
pfm_file = args[3]
out_dir = args[4]
coreNum = args[5]

# singularity shell  /proj/nobackup/sens2025644/wharf/ekolyal/ekolyal-sens2025644/nf-atacco/conts/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-genomicranges_pruned_55fb052e06404cfb.sif
# difftables = "/proj/sens2025644/nbis8271/results/nf-atac_drop_103_114_121_results/dar/HFpEF_vs_HFrEF.tsv"
# peakAnnotation = "/proj/sens2025644/TF/results/topAR/AnnotatedPeaks_all.bed"
# pfm_file = "/proj/nobackup/sens2025644/wharf/ekolyal/ekolyal-sens2025644/nf-atacco/data/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
# out_dir = "test" # "monalisa_output"

# List all objects with their class and size
ls_obj_info = function(pos = 1, order.by = "Size", decreasing = TRUE, head.n = 10) {
  napply = function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))
  
  obj_names = ls(pos = pos)
  obj_class = napply(obj_names, class)
  obj_size  = napply(obj_names, object.size)
  
  obj_df = data.frame(
    Name = obj_names,
    Class = as.character(obj_class),
    Size = obj_size
  )
  
  obj_df = obj_df[order(obj_df[[order.by]], decreasing = decreasing), ]
  obj_df$Size = format(obj_df$Size, units = "auto")
  head(obj_df, head.n)
}

prepareData = function(dat, outdir) {
    print("----Creating GRanges object...")
    gr = GRanges(
        seqnames = dat$seqnames,
        ranges = IRanges(start = dat$start, end = dat$end),
        strand = dat$strand
    )
    ## Add metadata columns (all others except these three)
    mcols(gr) = dat[, setdiff(colnames(dat), c("seqnames", "start", "end", "strand"))]

    print("----Calculating GC content for each peak...")
    # ff <- FaFile("${ref}")
    # peakSeqs = getSeq(x = Sequence, gr)
    gcContentPeaks = letterFrequency(DNAStringSet(gr$Sequence), "GC", as.prob=TRUE)[,1]
    mcols(gr)$gc = gcContentPeaks

    # print("----Keep enhancers at least 1kb away from any TSS and not in any gene, this is where most difference occur!")
    # keep = gr$annotation %in% c("Distal Intergenic", "Promoter (1-2kb)", "Promoter (2-3kb)")
    # gr = gr[keep]
    # table(gr$annotation)

    print("----Only keep enhancers on autosomes and sex chromosomes...")
    ## subset autosomal enhancers
    # seqlevels(gr) = paste0("chr", c(1:22, "X", "Y"))
    keep = seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y"))
    gr = gr[keep]
    # seqnames(gr) = droplevels(seqnames(gr))
    table(gr$annotation)

    ## Fix names
    names(gr) = gr$peakID

    print("----Plot logFC vs GC content...")
    paste0(outdir,"/data/") %>% dir.create(., recursive=T, showWarnings=F)
    pdf(paste0(outdir,"/data/logFCvsGC.pdf"))
    # par(mfrow=c(1,2))
    plot(gr$gc, gr$log2FoldChange, pch = ".")
    abline(h = 0, col = "red", lty = 5)
    # smoothScatter(gr$gc, gr$logFC)
    # abline(h = 0, col = "red", lty = 5)
    dev.off()
    saveRDS(gr, paste0(outdir,"/data/GRanges_object.RDS"))
    return(gr)
}


runBinnedMotifEnr = function(gr, pfm_file, outdir) {
    print("Now let's perform binned Motif Enrichment analysis...")
    col2rm = intersect(
        c("seqnames", "ranges", "strand", "seqlevels", 
            "seqlengths", "isCircular", "start", "end", "width", "element"),
        colnames(mcols(gr)))
    colnames(mcols(gr))[colnames(mcols(gr)) %in% col2rm] = paste0(col2rm, "_2")
    ### region size distribution
    summary(width(gr))
    print("----Resize the regions and trim out-of bounds ranges...")
    gr = trim(resize(gr, width = median(width(gr)), fix = "center"))
    summary(width(gr))
    # width(gr) returns a vector of widths of each range in the GRanges object gr.
    # median(width(gr)) calculates the median width of all these ranges.
    # resize(gr, width = median(width(gr)), fix = "center") resizes all ranges to this median width, keeping the center of each range fixed.
    # trim() ensures that any resized ranges do not extend beyond the chromosome boundaries (i.e., trims them to valid coordinates).
    print("----Plot log2FC histogram...")
    paste0(outdir,"/data/") %>% dir.create(., recursive=T, showWarnings=F)
    pdf(paste0(outdir, "/data/logFC_histgram.pdf"))
        print(ggplot(data = data.frame(logFC = gr$log2FoldChange)) +
            geom_histogram(aes(x = logFC), bins = 100, fill = "steelblue") +
            xlab("Batf cKO vs Wt logFC") +
            theme_bw())
    dev.off()

    print("----Bin the histogram...")
    bins = bin(x = gr$log2FoldChange, 
        binmode = "equalN", nElement = 1400, minAbsX = 0.3)
    table(bins)
    # binmode = "equalN": Bins are made so that each bin contains approximately the same number of elements (nElement specifies the bin size).
    # nElement = 800: Each bin contains about 800 values.
    # minAbsX = 0.3: Possibly a threshold to ignore or filter very small absolute values before binning.
    print("--------Plot binned histogram...")
    pdf(paste0(outdir, "/data/binnedHistogram.pdf"), w=8, h=6)
        print(plotBinDensity(x = gr$log2FoldChange, b = bins) +
            xlab("logFC"))
    dev.off()

    print("-----Before proceeding with the enrichment analysis, 
            let us  check if there is any sequence bias associated")
    # with the bins. oﬀers some plo%ng func!ons for this purpose.
    # extract DNA sequences of the enhancers
    seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
    pdf(paste0(outdir, "/data/GCfraction.pdf"), h=6, w=6)
        # by GC fraction
        print(plotBinDiagnostics(seqs = seqs, bins = bins, aspect = "GCfrac"))
        # by dinucleotide frequency
        print(plotBinDiagnostics(seqs = seqs, bins = bins, aspect = "dinucfreq"))
    dev.off()

    print("----Get PWMs from Jaspar...")
    ## extract PWMs of vertebrate TFs from JASPAR2024
    ## Read the PFMs from the file
    pfm = readJASPARMatrix(pfm_file, matrixClass=c("PFM"))
    pwm = toPWM(pfm)
    names(bins) = names(seqs)
    names(pwm) = ID(pwm)
    list(seqs=seqs, bins=bins, pwm=pwm) %>% saveRDS(
        paste0(outdir, "/data/monalisa_input_seqs_bins_pwm.RDS"))

    print("----Run binned enrichment")
    ## motif enrichment using 4 cores
    # library(BiocParallel)
    bp = SnowParam(workers = coreNum, type = "SOCK", 
        exportglobals = TRUE, timeout = 3600)
    # bp = BiocParallel::MulticoreParam(1)
    on.exit({
        bpstop(bp)
        closeAllConnections()
        gc()
    }, add = TRUE)
    se = calcBinnedMotifEnrR(seqs = seqs,
        bins = bins,
        pwmL = pwm,
        background = "otherBins",
        BPPARAM = bp)

    print("----Save SummarizedExperiment object...")
    paste0(outdir, "/result_files") %>%
        dir.create(., recursive=T, showWarnings=FALSE)
    saveRDS(se, paste0(outdir, "/result_files/binnedEnr.RDS"))
    rowData(se)[, c("motif.id", "motif.name", "motif.percentGC")] %>%
        as.data.frame() %>%
        write.csv(file = paste0(outdir, "/result_files/binnedEnr_motifInfo.csv"),
            row.names = FALSE)
    write.csv(colData(se),
        file = paste0(outdir, "/result_files/binnedEnr_binInfo.csv"),
        row.names = TRUE)
    for (assay_name in names(assays(se))) {
        write.csv(as.data.frame(assays(se)[[assay_name]]),
            file = paste0(outdir, "/result_files/binnedEnr_assay_", assay_name, ".csv"),
            row.names = TRUE)
    }
    # assays(se)
    # head(assays(se)$log2enr)
    return(se)
}


visualizeEnrichment = function(se, outdir) {
    print("----Visualization...")
    print("-------select strongly enriched motifs")
    TOProwMax = apply(assay(se, "negLog10Padj"), 1,
        function(x) max(abs(x), 0, na.rm = TRUE)) %>%
        as.numeric() %>%
        quantile(., probs=0.915)
    sel = apply(assay(se, "negLog10Padj"), 1,
        function(x) max(abs(x), 0, na.rm = TRUE)) > TOProwMax
    sum(sel)
    seSel = se[sel, ]

    print("-------plot heatmap of selected motifs")
    pdf(paste0(outdir,"/MotifHeatmap.pdf"), w=8, h=2+sum(sel)/8)
    print(plotMotifHeatmaps(
        x = seSel, which.plots = c("log2enr", "negLog10Padj"),
        width = 2.0, cluster = TRUE, maxEnr = 2, maxSig = 10,
        show_motif_GC = TRUE))

    print("-------plot with motif sequence logos")
    SimMatSel = motifSimilarity(rowData(seSel)$motif.pfm)
    paste0(outdir, "/result_files") %>%
        dir.create(., recursive=T, showWarnings=FALSE)
    saveRDS(SimMatSel, paste0(outdir,"/result_files/Motif_Similarity_Matrix.RDS"))
    range(SimMatSel)
    # Create hclust object, similarity defined by 1-Person correlation
    hcl = hclust(as.dist(1-SimMatSel), method = "average")
    print(plotMotifHeatmaps(
        x = seSel, which.plots = c("log2enr", "negLog10Padj"),
        width = 1.8, cluster = hcl, maxEnr = 2, maxSig = 10,
        show_dendrogram = TRUE, show_seqlogo = TRUE,
        width.seqlogo = 1.2))
dev.off()
}

getTFmotifHits = function(peakannotation, pfm_file, out_dir, difftables) {
    # peakannotation = readr::read_tsv(peakAnnotation)
    pfm = readJASPARMatrix(pfm_file, matrixClass=c("PFM"))
    pwm = toPWM(pfm)
    names(pwm) = ID(pwm)
    hits = lapply(1:nrow(peakannotation), function(i) {
        seq = getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
            peakannotation[i, "seqnames", drop=TRUE], 
            start=peakannotation[i, "start", drop=TRUE], 
            end=peakannotation[i, "end", drop=TRUE])
        seq_set = DNAStringSet(seq)
        names(seq_set) = peakannotation[i, "Geneid", drop=TRUE]
        bp = SnowParam(workers = coreNum, type = "SOCK", 
                    exportglobals = TRUE, timeout = 3600)
        on.exit({
            bpstop(bp)
            closeAllConnections()
            gc()
        }, add = TRUE)
        hit = findMotifHits(query = pwm,
            subject = seq_set,
            min.score = 10.0,
            method = "matchPWM",
            BPPARAM = bp)
        hit
    }) 

    saveRDS(hits, paste0(out_dir, "TFmotif_hits_list_allPeaks.RDS"))
    combined_hits = do.call(c, hits)
    rm(hits)
    gc()
    # we can summarize the number of predicted hits per promoter in matrix format
    hitsMatrix = table(
        factor(seqnames(combined_hits), levels = unique(peakannotation$Geneid)),
        factor(combined_hits$pwmname, levels = unique(name(pwm)))) 
    hitsMatrix_dat = hitsMatrix %>% as.data.frame.matrix() %>% 
        tibble::rownames_to_column(var="peakID") 
    rm(hitsMatrix)
    gc()
    saveRDS(hitsMatrix_dat, paste0(out_dir, "TFmotif_hitsMatrix_allPeaks.RDS"))
    write.table(hitsMatrix_dat, 
        file = paste0(out_dir, "TFmotif_hitsMatrix_allPeaks.tsv"), 
        sep = "\t", row.names = FALSE, quote = FALSE)

    lapply(difftables, function(difftable) {
        dat = hitsMatrix_dat %>%
            dplyr::right_join(read.table(difftable) %>% 
                tibble::rownames_to_column(var="peakID") %>%
                select(-c(Chr, Start, End)) %>%
                dplyr::left_join(peakannotation %>% 
                mutate(peakID=Geneid) %>%
                dplyr::select(-Geneid)) %>% na.omit() %>%
                arrange(., pvalue))
        outdir = paste0(out_dir, "/", 
                gsub(".tsv$", "", basename(difftable)))
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        write.table(dat, 
            file = paste0(outdir, "TFmotif_hitsMatrix_diff.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
    })
}




suppressPackageStartupMessages({
library(monaLisa)
library(BiocParallel)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2024)
library(TFBSTools)
library(SummarizedExperiment)
library(RSQLite)
library(dplyr)
library(Biostrings)
library(TFBSTools)
library(stringr)
library(tibble)
})

print("Loading peak annotation...")
print(peakAnnotation)
peakannotation = readr::read_tsv(peakAnnotation)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

print("Loading differential accessibility tables ...")
file.copy(from=pfm_file, 
    to=paste0(out_dir,"/", basename(pfm_file)))
print("----Please be aware that Monalia requires 0-based peak coordinate!!!")
difftables = strsplit(difftables, ";")[[1]]

# Get and save TF motif hits for all peaks
print("Get TF motif hits for all peaks...")
getTFmotifHits(peakannotation, pfm_file, out_dir, difftables)

# For each differential table, run binned motif enrichment analysis
print("Run binned motif enrichment analysis for each differential table...")
lapply(difftables, function(difftable) {
    print(paste0("Processing ", difftable, "..."))
    dat = read.table(difftable) %>% 
        tibble::rownames_to_column(var="peakID") %>%
        select(-c(Chr, Start, End)) %>%
        dplyr::left_join(peakannotation %>% 
        mutate(peakID=Geneid) %>%
        dplyr::select(-Geneid)) %>% na.omit()
    outdir = paste0(out_dir, "/", 
        gsub(".tsv$", "", basename(difftable)))
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    print("Prepare data...")
    gr = prepareData(dat, outdir)
    
    print("Run binned motif enrichment analysis...")
    se = runBinnedMotifEnr(gr, pfm_file, outdir)
    
    print("Visualize enrichment results...")
    visualizeEnrichment(se, outdir)
})

# ls_obj_info()

# print("Binned k-mer enrichment analysis...")
# print("----Calculating k-mer enrichment...")
# seKmer = calcBinnedKmerEnr(seqs = seqs,
#     bins = bins,
#     kmerLen = 6,
#     includeRevComp = TRUE,
#     BPPARAM = bp)
# print("----Select the most enriched kmers (negLog10Padj>4)...")
# selKmer = apply(assay(seKmer, "negLog10Padj"), 1
#     function(x) max(abs(x), 0, na.rm = TRUE)) > 4.0
# sum(selKmer)
# seKmerSel = seKmer[selKmer, ]
# print("----Calculate similarity between enriched kmers and enriched motifs...")
# pfmSel = rowData(seSel)$motif.pfm
# sims = motifKmerSimilarity(
#     x=pfmSel, kmers=rownames(seKmerSel), 
#     includeRevComp=TRUE)
# dim(sims)
# print("----Plot heatmap of similarity between selected motifs and enriched k-mers...")
# maxwidth = max(sapply(TFBSTools::Matrix(pfmSel), ncol))
# seqlogoGrobs = lapply(pfmSel, seqLogoGrob, xmax = maxwidth)
# hmSeqlogo = rowAnnotation(logo=annoSeqlogo(seqlogoGrobs, which="row"),
#     annotation_width = unit(1.5, "inch"),
#     show_annotation_name = FALSE)
# pdf(paste0(outdir,"/Kmer_Motif_SimilarityHeatmap.pdf"), 
#     w=3+dim(sims)[2]/8, h=2+dim(sims)[1]/8)
# Heatmap(sims, show_row_names=TRUE, row_names_gp=gpar(fontsize=8),
#     show_column_names=TRUE, column_names_gp=gpar(fontsize=8),
#     name="Similarity", 
#     column_title="Selected TFs and enriched k-mers",
#     right_annotation=hmSeqlogo)
# dev.off()

out = capture.output({sessionInfo()})
writeLines(out, paste0(out_dir,"/Monalisa_sessionInfo.txt"))



