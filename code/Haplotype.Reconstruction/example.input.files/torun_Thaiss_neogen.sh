marker_grid_file="/projects/GeDI/resource/haplo_reconst/support_files/marker_grid_0.02cM_plus.txt"

working_dir="/projects/compsci/corneb_cs"

snps_file="/projects/GeDI/resource/haplo_reconst/support_files/gm_uwisc_v1.csv"

inventory_file="/projects/compsci/corneb_cs/data/Univ_of_Penn_Thaiss/Neogen/Univ_of_Penn_Thaiss_neogen_inventory.file.csv"

founder_geno_file="/projects/GeDI/resource/haplo_reconst/support_files/GM_founders.csv"

output="/projects/compsci/corneb_cs/haplo.reconst/Univ_of_Penn_Thaiss"

tag="GigaMUGA"
 
project_file="/projects/compsci/corneb_cs/data/Univ_of_Penn_Thaiss/Neogen/Univ_of_Penn_Thaiss_neogen.projectfile.csv"

project=`basename $project_file`
name="$(cut -d . -f1 <<< $project)" 
output_dir=${output}/$name  

echo "bash launch_haplotype_reconstruction.sh -o ${output_dir} -i ${inventory_file} \
 -p ${project_file} -s ${snps_file} -f ${founder_geno_file} -g ${grid_file} -t ${tag}"

bash launch_haplotype_reconstruction.sh -o ${output_dir} -i ${inventory_file} \
 -p ${project_file} -s ${snps_file} -f ${founder_geno_file} -g ${marker_grid_file} -t ${tag}

