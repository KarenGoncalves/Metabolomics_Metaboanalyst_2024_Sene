# Annotate deregulated analytes
library(tidyverse)
library(MetaboAnalystR)

#### Variables ####
Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
FDR_threshold = 0.01
FC_threshold = 1

#### Load data ####
load("Results/Siggenes_DAAs.RData")
annotation <- sapply(Analysis_modes, simplify = F, \(x) {
    paste0("Inputs/Corrected_LCMSMS_", x, "_identification.txt") %>%
        read_delim(delim = "\t", na = "null") %>%
        mutate(AnalysisMode = x)
}) %>% list_rbind() 

#### Keep only significantly DAAs ####
differential_abundance_sig <- 
    differential_abundance_all %>%
    list_rbind %>% filter(!is.na(pValue),
                          pValue < FDR_threshold,
                          abs(FoldChange) > FC_threshold
    ) %>%
    separate(col = Metabolite, 
             into = c("Rt", "Mz"), 
             sep = "/",  convert = T)

DAM_ids <- differential_abundance_sig %>%
    select(Rt, Mz, AnalysisMode) %>%
    unique

#### Get info on DAAs ####
annotation_DAAs <- 
    inner_join(annotation, DAM_ids,
               by = join_by("Average Rt(min)" == "Rt",
                            "Average Mz" == "Mz",
                            "AnalysisMode" == "AnalysisMode")) %>%
    mutate(Annotation = 
               ifelse(
                   `MS/MS matched` == T & !grepl("w/o MS2", `Metabolite name`),
                   `Metabolite name`, "Unknown"
                   ),
           Clean_name = gsub("; (LC-|CE\\d).+$", "",
                             Annotation) %>%
               gsub(pattern="\\(*[Nn]ot validated.*",
                    replacement="") %>%
               gsub(pattern="^DDAO$",
                    replacement="Decyl dimethyl amine oxide") %>%
               gsub(pattern="^l-", replacement="L-") %>%
               gsub(pattern="^(NCGC\\d+-\\d+).+$",
                    replacement="\\1") %>%
               gsub(pattern="andrographolide",
                    replacement="Andrographolide", fixed=T) %>%
               gsub(pattern="(2R,3S,4S,5R,6R)-6-(((4S,5aS,7S,11aR,12aS)-4,7-dihydroxy-3-((2R,5S)-5-(2-hydroxypropan-2-yl)-2-methyltetrahydrofuran-2-yl)-2a,5a,8,8-tetramethylhexadecahydrocyclopenta[a]cyclopropa[e]phenanthren-9-yl)oxy)-5-(((2S,3R,4R)-3,4-dihydroxy-4-(hydroxymethyl)tetrahydrofuran-2-yl)oxy)-2-(hydroxymethyl)tetrahydro-2H-pyran-3,4-diol",
                    replacement="Cucurbitacin glycoside [2]", fixed=T) %>%
               gsub(pattern="(2S,3S,4S,5S,6R)-2-(((2aR,4S,5aS,7S,11aR,12aS)-4-hydroxy-3-((2R,5S)-5-(2-hydroxypropan-2-yl)-2-methyltetrahydrofuran-2-yl)-2a,5a,8,8-tetramethyl-9-(((2S,3R,4S,5R)-3,4,5-trihydroxytetrahydro-2H-pyran-2-yl)oxy)hexadecahydrocyclopenta[a]cyclopropa[e]phen...",
                    replacement="Cucurbitacin glycoside [1]", fixed=T) %>%
               gsub(pattern="(2S,3S,4S,5S,6R)-2-(((2aR,4S,5aS,7S,11aR,12aS)-4-hydroxy-3.+",
                    replacement="Cyclocarposide", fixed=T) %>%
               gsub(pattern="(4aR,5aS,9R)-9-ethynyl-9a,11b-dimethylhexadecahydrocyclopenta[1,2]phenanthro[8a,9-b]oxirene-3,9-diol",
                    replacement="Estrane steroid", fixed=T) %>%
               gsub(pattern="N-((octahydro-1H-quinolizin-1-yl)methyl)-2,4,5,6-tetrahydrocyclopenta[c]pyrazole-3-carboxamide",
                    replacement="Quinolizine", fixed=T) %>%
               gsub(pattern="2-(3,4-dimethoxyphenyl)-7-methoxy-4H-chromen-4-one",
                    replacement="7-O-methylated flavonoid [1]", fixed=T) %>%
               gsub(pattern="(R)-((2R,3S,4S,5R,6S)-6-((3-(2,3-dihydrobenzo[b][1,4]dioxin-6-yl)-4-oxo-4H-chromen-7-yl)oxy)-3,4,5-trihydroxytetrahydro-2H-pyran-2-yl)methyl 2-((tert-butoxycarbonyl)amino)-3-phenylpropanoate",
                    replacement="Isoflavonoid O-glycoside [1]", fixed=T) %>%
               gsub(pattern="methyl 2-((4-methyl-2-oxo-2H-chromen-7-yl)oxy)propanoate",
                    replacement="Coumarin derivative [1]", fixed=T) %>%
               gsub(pattern="4-((1R,3S,5r,7r)-5,7-dimethyl-1,3-diazaadamantan-2-yl)-2-methoxyphenol",
                    replacement="Methoxyphenol-type compound [1]", fixed=T) %>%
               gsub(pattern="7-benzyl-11,14-dimethyl-16-(2-methylpropyl)-10,13-di(propan-2-yl)-17-oxa-1,5,8,11,14-pentazabicyclo[17.3.0]docosane-2,6,9,12,15,18-hexone",
                    replacement="Cyclodepsipeptide", fixed=T) %>%
               gsub(pattern="7-hydroxy-3-(4-hydroxyphenyl)-8-((2S,3R,4R,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)tetrahydro-2H-pyran-2-yl)-4H-chromen-4-one",
                    replacement="Isoflavonoid C-glycoside [1]", fixed=T)           
           
    )

# Replace with sentence case names
names_to_correct <- annotation_DAAs$Clean_name %in% 
    c("THIAMINE", "ABIETIC ACID", "THIAMINE PYROPHOSPHATE", "lappaconitine")
annotation_DAAs$Clean_name[names_to_correct] <- 
    str_to_sentence(annotation_DAAs$Clean_name[names_to_correct])

write_delim(annotation_DAAs,
            file="Results/annotation_DAAs.txt",
            delim="\t", na="NA")
