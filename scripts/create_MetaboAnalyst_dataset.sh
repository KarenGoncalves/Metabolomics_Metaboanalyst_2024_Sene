#!/bin/sh

#mkdir /e/Metabolomics_Metaboanalyst_2024_Sene/Inputs
#mkdir /e/Metabolomics_Metaboanalyst_2024_Sene/metadata
#mkdir /e/Metabolomics_Metaboanalyst_2024_Sene/scripts

#Positive has 17 samples: 12 replicates, 2 QCs, 3 blanks
input=/e/Samples\ and\ QCs_Positives/Height_0_2024_05_23_10_09_51.txt
output="/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/LCMSMS_RP_Positive_rawHeight.txt"
head -n5 "${input}" |\
 tail -n 1 | cut -f 36-52 > tmp_header

cat tmp_header | sed -E 's/([A-Z0-9\.]+)([\._])([I12]+)\t/\1\t/g' |\
 sed -E 's/([A-Za-z]+\.blank[\.a-z]*)/BLANK/g' > tmp_row2

tail -n +6 "${input}" |\
 cut -f 36-52 > tmp_Raw_height

cat tmp_header tmp_row2 tmp_Raw_height > tmp_file

awk 'BEGIN {FS="\t";OFS="/"; print "Sample\nLabel"}\
 NR > 5 {print $2,$3}'  "${input}"\
 > tmp_RowNames

paste tmp_RowNames tmp_file > ${output}

rm tmp_file tmp_RowNames tmp_Raw_height tmp_row2 tmp_header

##Negative has 17 samples: 12 replicates, 2 QCs, 3 blanks
input=/e/Samples\ and\ QCs_Negative/Height_0_2024_05_23_09_56_49.txt

head -n5 "${input}" |\
 tail -n 1 | cut -f 36-52 > tmp_header

cat tmp_header | sed -E 's/([A-Z0-9\.]+)([\._])([I12]+)\t/\1\t/g' |\
 sed -E 's/([A-Za-z]+\.blank[\.a-z]*)/BLANK/g' > tmp_row2

tail -n +6 "${input}" |\
 cut -f 36-52 > tmp_Raw_height

cat tmp_header tmp_row2 tmp_Raw_height > tmp_file

awk 'BEGIN {FS="\t";OFS="/"; print "Sample\nLabel"}\
 NR > 5 {print $2,$3}' "${input}"\
 > tmp_RowNames

paste tmp_RowNames tmp_file > /e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/LCMSMS_RP_Negative_rawHeight.txt
rm tmp_file tmp_RowNames tmp_Raw_height tmp_row2 tmp_header

##HILIC has 16 samples: 12 replicates, 2 QCs, 2 blanks
input=/e/Samples\ and\ QCs\ recalibrated_HILIC/Height_0_2024_05_04_11_47_27.txt

head -n5 "${input}" |\
 tail -n 1 | cut -f 36-51 > tmp_header

cat tmp_header | sed -E 's/([A-Z0-9\.]+)([\._])([I12]+)\t/\1\t/g' |\
 sed -E 's/([A-Za-z]+\.blank[\.a-z]*)/BLANK/g' > tmp_row2

tail -n +6 "${input}" |\
 cut -f 36-51 > tmp_Raw_height

cat tmp_header tmp_row2 tmp_Raw_height > tmp_file

awk 'BEGIN {FS="\t";OFS="/"; print "Sample\nLabel"}\
 NR > 5 {print $2,$3}' "${input}"\
 > tmp_RowNames

paste tmp_RowNames tmp_file > /e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/LCMSMS_HILIC_Positive_rawHeight.txt

rm tmp_file tmp_RowNames tmp_Raw_height tmp_row2 tmp_header
