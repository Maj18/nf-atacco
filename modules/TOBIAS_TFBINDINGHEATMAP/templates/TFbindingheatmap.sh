#!/usr/bin/env bash

# Read env vars
TF="$1"
groups="$2"
# ftscores="$3"
corrected_dirs="$3"
ls ${corrected_dirs}

echo "Generate binding heatmap for ${TF}..."
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
        signal2=${corrected_dirs[$j]}/${g2}_merged_filtered_corrected.bw
        bound_bed1=DiffTFBinding/${TF}/beds/${TF}_${g1}_bound.bed
        unbound_bed1=DiffTFBinding/${TF}/beds/${TF}_${g1}_unbound.bed
        bound_bed2=DiffTFBinding/${TF}/beds/${TF}_${g2}_bound.bed
        unbound_bed2=DiffTFBinding/${TF}/beds/${TF}_${g2}_unbound.bed
	echo $unbound_bed2
        echo "Create heatmaps comparing bound and unbound sites between conditions ${g1} and ${g2} for ${TF}"
        TOBIAS PlotHeatmap \
            --TFBS ${bound_bed1} ${unbound_bed1} \
            --TFBS ${bound_bed2} ${unbound_bed2} \
            --signals ${signal1} ${signal2} \
            --output ${TF}_bindingHeatmap_${g1}vs${g2}.png \
            --signal_labels ${g1} ${g2} \
            --share_colorbar \
            --sort_by -1

        echo "Heatmaps comparing bound and unbound sites between conditions ${g1} and ${g2} generated for ${TF}!"
    done
done
