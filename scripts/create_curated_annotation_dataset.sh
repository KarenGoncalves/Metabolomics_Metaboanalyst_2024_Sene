#!/bin/sh

get_annotation (){
 tail -n +5 1 | cut -f 1-35 > $2
}

#Positive
input="/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/LCMSMS_RP_Positive_Identification_MSDial_format.txt"
output="/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/Corrected_LCMSMS_RP_Positive_identification.txt"


##Negative
input"/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/LCMSMS_RP_Negativee_Identification_MSDial_format.txt"
output="/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/Corrected_LCMSMS_RP_Negative_identification.txt"
get_annotation "${input}" ${output}


##HILIC
input="/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/LCMSMS_HILIC_Positive_Identification_MSDial_format.txt"
output="/e/Metabolomics_Metaboanalyst_2024_Sene/Inputs/Corrected_LCMSMS_HILIC_Positive_identification.txt"
get_annotation "${input}" ${output}
