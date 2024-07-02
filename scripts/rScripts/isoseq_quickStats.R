## Iso-seq quick stats

workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/2024_05_29_mskiLab_MelanomaIsoSeq/v1.compareAMtoCM/'
setwd(workingDir)

allInfo_iq <- data.table::fread('combined_allInfo_Incomplete_wIsoquantAssignments.gz',sep = "\t")
colnames(allInfo_iq) <- c("Read","Gene","CellLine","BC","Clust","IntronChain","TSS","PolyA",
                          "ExonChain","knownT","NumIntrons","Transcript","Type","Class")
allInfo_iq <- allInfo_iq %>% mutate(Class = gsub("Classification=|;","",Class))

summary_stats <- allInfo_iq %>% group_by(CellLine,Class) %>% dplyr::select(CellLine, Class) %>% 
  add_count(name = "NumReads") %>% ungroup() %>% distinct() %>% filter(Class!="")

ggplot(summary_stats, aes(x = Class, y = NumReads, fill = CellLine)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90))

taStatus <- allInfo_iq %>% dplyr::select(CellLine, TSS,PolyA) %>% mutate(TSS = case_when(TSS != "NoTSS" ~ "Present", TRUE ~ "NoTSS"), 
                                                        PolyA = case_when(PolyA != "NoPolyA" ~ "Present", TRUE ~ "NoPolyA")) %>% 
  group_by_all() %>% add_count() %>% distinct()


ggplot(taStatus, aes(x = interaction(TSS,PolyA), y = n, fill = CellLine)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90))
