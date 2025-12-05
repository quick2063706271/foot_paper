hm_cell_lines=( "WM4324" "B16_foot_vs_skin" "all_Acral_vs_SunExposed" "primary_all_Acral_vs_SunExposed" "metastatic_all_Acral_vs_SunExposed" "acral_cutaneous_weiss" "primary_acral_cutaneous_weiss" "CU_AM_CM" "primary_CU_AM_CM" )

mus_cell_lines=( "YUMM1-7" "YUMM4-1" "B16_foot_vs_skin" "new_foot_ear" )

##### change dir of and also dir defined in run_gsea_cli.sh

for mus_cell_line in "${mus_cell_lines[@]}"; do
    bash /Users/kyu/Documents/rebecca_lab/acral_melanoma/scripts/run_gsea_cli.sh $mus_cell_line mouse
done


hm_cell_lines=( "acral_cutaneous_weiss" )
cell_lines=( "weiss_AM_primary_vs_metastatic" "weiss_CM_primary_vs_metastatic" "farshidfar_AM_primary_vs_metastatic" "farshidfar_CM_primary_vs_metastatic" )
for hm_cell_line in "${hm_cell_lines[@]}"; do
    bash /Users/kyu/Documents/rebecca_lab/acral_melanoma/scripts/run_gsea_cli.sh $hm_cell_line human
done

