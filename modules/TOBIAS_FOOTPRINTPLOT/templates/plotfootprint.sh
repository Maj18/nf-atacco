#!/usr/bin/env bash

# Read env vars
TF="$1"
groups="$2"
# ftscores="$3"
corrected_dirs="$3"
ls ${corrected_dirs}

echo "Generate pairwise footprint plot for ${TF}..."
# groups=${groups.split(',').join(" ")}
# signals=${ftscores.split(',').join(" ")}
IFS=',' read -r -a groups <<< "$groups"
IFS=',' read -r -a corrected_dirs <<< "$corrected_dirs"
echo $groups
echo $corrected_dirs
for ((i=0; i<${#groups[@]}-1; i++)); do
    for ((j=i+1; j<${#groups[@]}; j++)); do
        g1=${groups[$i]}
        g2=${groups[$j]}
	signal1=${corrected_dirs[$i]}/${g1}_merged_filtered_corrected.bw
        signal2=${corrected_dirs[j]}/${g2}_merged_filtered_corrected.bw
        bound_bed1=DiffTFBinding/${TF}/beds/${TF}_${g1}_bound.bed
        bound_bed2=DiffTFBinding/${TF}/beds/${TF}_${g2}_bound.bed
        echo $bound_bed1
	echo "Plot aggregate footprints for ${TF} across conditions \${g1} and \${g2}"
        TOBIAS PlotAggregate \
            --TFBS ${bound_bed1} ${bound_bed2} \
            --signals ${signal1} ${signal2} \
            --output ${TF}_footprint_plot_${g1}vs${g2}.png \
            --share_y both \
            --plot_boundaries \
            --signal_labels ${g1} ${g2} \
            --smooth 5

        echo "Footprint plot between ${g1} and ${g2} generated for ${TF}!"
    done
done
