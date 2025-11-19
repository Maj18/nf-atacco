#!/usr/bin/env bash

# Read env vars
bedgraphs=$1
peakCallBed=$2
geneModelGTF=$3
trackIniNameSuffix=$4
# convert string to array by splitting on ','
IFS=',' read -r -a bedgraphs <<< "$bedgraphs"
    # Write header
    cat <<HEADER > tracks${trackIniNameSuffix}.ini
[x-axis]
where = top

[spacer]
height = 0.3
HEADER

    # Loop through all bedgraph files
    for f in ${bedgraphs[@]}; do
        fname=$(basename "$f" .bedgraph); name=${fname#P35002_}; name=${name%%.filtered.bedgraph}; name=${name/_REP1.mLb.clN.sorted.__/}; name=${name//__/}

        cat <<TRACK >> tracks${trackIniNameSuffix}.ini
[bedgraph]
file = $f
height = 4
title = $name
color = green
number_of_bins = 1000
min_value = 0
max_value = 0.5

[spacer]
height = 0.3
TRACK
    done

    # Append BED and genes tracks
    cat <<TAIL >> tracks${trackIniNameSuffix}.ini
[bed]
file = ${peakCallBed}
height = 1
labels = true
max_labels = 1000
title = consensus Genrich peaks merged
color = crimson

[spacer]
height = 0.5

[genes]
file = ${geneModelGTF}
file_type = gtf
title = gene models
style = UCSC
labels = true
max_labels = 1000
arrow_interval = 3
display = stacked
fontsize = 10
gene_rows = 10
height = 7
all_labels_inside = true
merge_transcripts = true
TAIL
