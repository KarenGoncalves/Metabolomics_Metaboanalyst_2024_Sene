library(tidyverse)

input = read_delim("Inputs/LCMSMS_RP_Negative_rawHeight.txt")
metadata = data.frame(Replicates = names(input),
											Groups = input[27,]
)
