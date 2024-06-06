#!/bin/sh

#Positive 
input=/e/Samples\ and\ QCs_Positives/Height_0_2024_05_23_10_09_51.txt
output="/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/LCMSMS_RP_Positive_identification.txt"
tail -n +5 "${input}" | cut -f 1-35 > ${output}

##Negative
input=/e/Samples\ and\ QCs_Negative/Height_0_2024_05_23_09_56_49.txt
output="/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/LCMSMS_RP_Negative_identification.txt"
tail -n +5 "${input}" | cut -f 1-35 > ${output}


##HILIC
input=/e/Samples\ and\ QCs\ recalibrated_HILIC/Height_0_2024_05_04_11_47_27.txt
output="/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/LCMSMS_HILIC_Positive_identification.txt"
tail -n +5 "${input}" | cut -f 1-35 > ${output}
