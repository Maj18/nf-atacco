#!/usr/bin/env bash

# Read env vars
bedgraphs="$1"
peakCallBed="$2"
geneModelGTF="$3"
groups="$4"
trackIniNameSuffix="$5"

# echo $geneModelGTF
# echo $groups

# convert string to array by splitting on ','
IFS=',' read -r -a bedgraphs <<< "$bedgraphs"
IFS=',' read -r -a groups <<< "$groups"
# echo ${bedgraphs[@]}
# echo ${groups[@]}

uniqueG=($(printf "%s\n" "${groups[@]}" | sort -u))
# echo ${uniqueG[@]}

colors=("#9467bd" "#ff7f0e" "#2ca02c" "#d62728" "#8c564b" "#e377c2" "#7f7f7f" "#bcbd22" "#17becf" "#1f77b4")
# echo ${colors[0]}
# echo ${colors[1]}
    # Write header
    cat <<HEADER > tracks${trackIniNameSuffix}.ini
[x-axis]
where = top

[spacer]
height = 0.3
HEADER

    # Loop through all bedgraph files
    # for f in ${bedgraphs[@]}; do
    for ((i=0; i<${#bedgraphs[@]}; i++)); do
        f=${bedgraphs[$i]}
        g=${groups[$i]}
        # echo "$g"
        for ((j=0; j<${#uniqueG[@]}; j++)); do
	    # echo "$j"
            if [[ "${uniqueG[$j]}" == "$g" ]]; then
                color="${colors[$j]}"
                break
            fi
        done
	# echo "$color"
        fname=$(basename "$f" .bedgraph); name=${fname#P35002_}; name=${name%%.filtered}; name=${name/_REP1.mLb.clN.sorted.__/}; name=${name//__/}

        cat <<TRACK >> tracks${trackIniNameSuffix}.ini
[bedgraph]
file = $f
height = 4
title = $name
color = $color
# number_of_bins = 1000
min_value = 0
max_value = 15
binSize = 10

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

