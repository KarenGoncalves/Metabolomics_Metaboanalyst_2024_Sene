#!/bin/sh

#Positive 
input=RAW_DATA/LCMSMS_RP_Positive_MSDial.txt
output=Inputs/LCMSMS_RP_Positive_identification.txt
tail -n +5 "${input}" | cut -f 1-35 > ${output}

##Negative
input=RAW_DATA/LCMSMS_RP_Negative_MSDial.txt
output=Inputs/LCMSMS_RP_Negative_identification.txt
tail -n +5 "${input}" | cut -f 1-35 > ${output}


##HILIC
input=RAW_DATA/LCMSMS_HILIC_Positive_MSDial.txt
output=Inputs/LCMSMS_HILIC_Positive_identification.txt
tail -n +5 "${input}" | cut -f 1-35 > ${output}
