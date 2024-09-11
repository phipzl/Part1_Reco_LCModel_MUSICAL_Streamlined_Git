#!/bin/bash
# This script extracts the B1 map from a DICOM directory and saves it so it can be used in Part2 for B1 correction. 
# 
# Requires exporting the variables B1_path, tmp_dir and out_path in a parent script.
# 
# If you want to run this script as a standalone, simply define the appropriate B1_path in the shell using:
# export B1_path=/path/to/B1/dicoms
# export tmp_dir=/path/to/tmp
# export out_path=/path/to/out
# Note that the tmp directory must contain csi_template.mnc

#B1_path=/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/3DMRSIMAP_Volunteers/3DMRSIMAP_Vol_12_A/B1/
#tmp_dir=/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/plazen/WIP/StandaloneB1Corr/tmp
#out_path=/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/plazen/WIP/StandaloneB1Corr/out


echo -e "\n\n##### B1 preparations #####"
date
echo -e "\nB1 path: $B1_path"
echo -e "out path: $out_path"
sleep 0.5

echo -e "\nremoving stuff..."
[ -f "${tmp_dir}/B1map*" ] && rm -f "${tmp_dir}/B1map*"
[ -f "${out_path}/maps/Extra/B1map*" ] && rm -f "${out_path}/maps/Extra/B1map*"
sleep 0.5

##### Create B1 minc
mkdir -p ${tmp_dir}/mnc
mkdir -p ${out_path}/maps/Extra

echo -e "\ndcm2mnc..."

dcm2mnc $B1_path ${tmp_dir}/mnc 
mv ${tmp_dir}/mnc/*/* ${tmp_dir}/B1map_orig.mnc
[ -d "${tmp_dir}/mnc" ] && rm -rf "${tmp_dir}/mnc"

# cp ${tmp_dir}/B1.mnc ${out_path}/maps/Extra/B1.mnc


# echo -e "\nmincstats..."
# mincstats ${tmp_dir}/csi_template.mnc
# mincstats ${tmp_dir}/B1map_orig.mnc 

# echo -e "\nmnc2raw..."
# minctoraw ${tmp_dir}/B1map_orig.mnc -nonormalize -float > ${tmp_dir}/B1map.raw 

# echo -e "\nraw2mnc..."
# rawtominc -float -clobber -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/B1map.raw ${tmp_dir}/B1map_orig.mnc

echo -e "\nmincresample..."
mincresample -clobber -nearest_neighbour -like ${tmp_dir}/csi_template.mnc ${tmp_dir}/B1map_orig.mnc ${tmp_dir}/B1map.mnc 

echo -e "\nmnc2raw..."
minctoraw ${tmp_dir}/B1map.mnc -nonormalize -float > ${tmp_dir}/B1map.raw 

echo -e "\ncp to out_dir..."
cp ${tmp_dir}/B1map.raw ${out_path}/maps/Extra/B1map.raw 
cp ${tmp_dir}/B1map.mnc ${out_path}/maps/Extra/B1map.mnc

echo -e "\nB1 preparations completed! B1map.raw and B1map.mnc are saved in:\n$out_path/maps/Extra\n\n"

