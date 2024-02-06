library(chromunity)
library(tidyr)
library(gUtils)
library(gTrack)


simplicies = readRDS('/gpfs/commons/groups/imielinski_lab/home/jorvis/Projects/testing/db/random_simplex_snippet/simplicial_snippet.rds')
view_range = GRanges('chr8:126000000-129000000')
plot(gTrack(split(simplicies, simplicies$cid) %>% unname, height = 10, name = 'simplicies'), view_range)