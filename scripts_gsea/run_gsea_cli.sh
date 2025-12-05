name=$1
specie=$2

###### change below
gsea_software_dir=/Users/kyu/Documents/app/GSEA_4.3.3/gsea-cli.sh
rnk_file=/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/gsea/${name}_DE_gene_names.rnk
out_dir=/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/gsea/result/

#####


# gmt_dir=/Users/kyu/Documents/app/msigdb_v2024.1.Hs_files_to_download_locally/msigdb_v2024.1.Hs_GMTs
gmt_dir=ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets
# chip_dir=/Users/kyu/Documents/app/msigdb_v2024.1.Hs_files_to_download_locally/chip
chip_dir=ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations
h_gmt_file=${gmt_dir}/h.all.v2024.1.Hs.symbols.gmt
c8_gmt_file=${gmt_dir}/c8.all.v2024.1.Hs.symbols.gmt
c5_gmt_file=${gmt_dir}/c5.go.bp.v2024.1.Hs.symbols.gmt
c2_gmt_file=${gmt_dir}/c2.cp.v2024.1.Hs.symbols.gmt

gmt_files=("$h_gmt_file" "$c8_gmt_file" "$c5_gmt_file" "$c2_gmt_file")
gmt_names=(h c8 c5gobp c2)
length=${#gmt_files[@]}

if [ "$specie" = "human" ]; then
    collapse_option=No_Collapse
    # chip=ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_HGNC_ID_MSigDB.v2024.1.Hs.chip
    chip=${chip_dir}/Human_HGNC_ID_MSigDB.v2024.1.Hs.chip
else
	collapse_option=Collapse
    chip=${chip_dir}/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v2024.1.Hs.chip
fi



# /Users/kyu/Documents/app/GSEA_4.3.3/gsea-cli.sh GSEAPreranked -gmx $gmt_file -collapse No_Collapse -mode Abs_max_of_probes -norm meandiv -nperm 1000 -rnd_seed timestamp -rnk $rnk_file -scoring_scheme weighted -rpt_label my_analysis -chip $chip -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 100 -set_max 800 -set_min 15 -zip_report false -out /Users/kyu/gsea_home/output/oct07

for i in "${!gmt_files[@]}"; do
    out_name=${name}_${gmt_names[$i]}
    echo $out_name
    echo $rnk_file
    echo ${gmt_names[$i]}
    $gsea_software_dir GSEAPreranked -gmx ${gmt_files[$i]} -collapse $collapse_option -mode Abs_max_of_probes -norm meandiv -nperm 1000 -rnd_seed timestamp -rnk $rnk_file -scoring_scheme weighted -rpt_label $out_name -chip $chip -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 100 -set_max 800 -set_min 15 -zip_report false -out $out_dir
done
