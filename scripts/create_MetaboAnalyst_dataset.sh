#!/bin/sh

#mkdir /e/Metabolomics_Metaboanalyst_2024_Sene/Inputs
#mkdir /e/Metabolomics_Metaboanalyst_2024_Sene/metadata
#mkdir /e/Metabolomics_Metaboanalyst_2024_Sene/scripts

create_dataset () {
col_start=$1
input=$2
output=$3
col_end=$(awk 'END {print NF}' ${input})

head -n5 "${input}" |\
 tail -n 1 | cut -f 36-52 > tmp_header

cat tmp_header | sed -E 's/([A-Z0-9\.]+)([\._])([I12]+)\t/\1\t/g' |\
 sed -E 's/([A-Za-z]+\.blank[\.a-z]*)/BLANK/g' > tmp_row2

tail -n +6 "${input}" |\
 cut -f ${col_start}-${col_end} > tmp_Raw_height

cat tmp_header tmp_row2 tmp_Raw_height > tmp_file

awk 'BEGIN {FS="\t";OFS="/"; print "Sample\nLabel"}\
 NR > 5 {print $2,$3}'  "${input}"\
 > tmp_RowNames

paste tmp_RowNames tmp_file | sed 's/E30/pPTGE30/g' > ${output}

rm tmp_file tmp_RowNames tmp_Raw_height tmp_row2 tmp_header
}

#Positive 

create_dataset 32 RAW_DATA/LCMSMS_RP_Positive_MSDial.txt Inputs/LCMSMS_RP_Positive_rawHeight.txt

##Negative 
create_dataset 32 RAW_DATA/LCMSMS_RP_Negative_MSDial.txt Inputs/LCMSMS_RP_Negative_rawHeight.txt

##HILIC 
create_dataset 32 RAW_DATA/LCMSMS_HILIC_Positive_MSDial.txt Inputs/LCMSMS_HILIC_Positive_rawHeight.txt
