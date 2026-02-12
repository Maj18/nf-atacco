#!/opt/conda/bin/R

# Usage: Rscript monalisaBinEnr.R <difftable> <peakAnnotation> <pfm_file> <outdir> <coreNum>

args = commandArgs(trailingOnly = TRUE)
difftables = args[1] # ; separated if multiple files
peakAnnotation = args[2]
pfm_file = args[3]
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

getTFmotifHits = function(peakannotation, pfm_file, difftables) {
    # peakannotation = readr::read_tsv(peakAnnotation)
    pfm = readJASPARMatrix(pfm_file, matrixClass=c("PFM"))
    pwm = toPWM(pfm)
    names(pwm) = ID(pwm)
    bp = SnowParam(workers = as.numeric(coreNum)*2, type = "SOCK", 
                    exportglobals = TRUE, timeout = 7200)
    print("Now this version of peak annotation is 0-based coordinate, but getSeq requires 1-based coordinate, so we need to convert the data while running getSeq")
    # seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    #            names = peakannotation$seqnames,
    #            start = peakannotation$start+1,
    #            end = peakannotation$end)
    # names(seqs) = peakannotation$Geneid
    # hits = findMotifHits(query = pwm,
    #                   subject = seqs,
    #                   min.score = 10.0,
    #                   method = "matchPWM",
    #                   BPPARAM = bp)

    # on.exit({
    #     bpstop(bp)
    #     closeAllConnections()
    #     gc()
    # }, add = TRUE)
    hits = bplapply(1:nrow(peakannotation), function(i, pwm_local) {
        seq = getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
            peakannotation[i, "seqnames", drop=TRUE], 
            start=peakannotation[i, "start", drop=TRUE] + 1,
            end=peakannotation[i, "end", drop=TRUE])
        seq_set = DNAStringSet(seq)
        names(seq_set) = peakannotation[i, "Geneid", drop=TRUE]
        hit = findMotifHits(query = pwm_local,
            subject = seq_set,
            min.score = 10.0,
            method = "matchPWM")
        hit
    }, BPPARAM = bp, pwm_local=pwm) 

    # Clean up
    bpstop(bp)
    closeAllConnections()
    gc()

    dir.create("PeakTFmotifHits/", recursive=TRUE, showWarnings=FALSE)
    saveRDS(hits, paste0("PeakTFmotifHits/TFmotif_hits_list_allPeaks.RDS"))
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
    saveRDS(hitsMatrix_dat, paste0("PeakTFmotifHits/TFmotif_hitsMatrix_allPeaks.RDS"))
    write.table(hitsMatrix_dat, 
        file = paste0("PeakTFmotifHits/TFmotif_hitsMatrix_allPeaks.tsv"), 
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
        outdir = paste0("PeakTFmotifHits/", 
                    gsub(".tsv$", "", basename(difftable)))
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        write.table(dat, 
            file = paste0(outdir, "/TFmotif_hitsMatrix_diff.tsv"),
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

print("----Please be aware that Monalia requires 0-based peak coordinate!!!")
difftables = strsplit(difftables, ";")[[1]]

# Get and save TF motif hits for all peaks
print("Get TF motif hits for all peaks...")
getTFmotifHits(peakannotation, pfm_file, difftables)

out = capture.output({sessionInfo()})
writeLines(out, paste0("/PeakTFmotifHits/Monalisa_sessionInfo.txt"))



