setwd("~/multiomic_pseudotimes/")


metpaths_obj <- 'syn47742117'
synapser::synGet(metpaths_obj, downloadLocation = "files/")
metpaths <- load(file="files/mets_DEanova_stats_path.rds")

class(path_res)
state <- c("2","3","4","5","6","7")
#library(dplyr)
path_res2 <- map2(path_res, state, ~cbind(.x, state = .y))
metpaths2 <- Reduce(full_join,path_res2)



metpaths2$leadingEdge <- as.character(metpaths2$leadingEdge)
metpaths2 <- metpaths2 %>% relocate(leadingEdge, .after = last_col())
write.table(metpaths2, file="mets_GSEAanalysis.csv", sep="\t", row.names=FALSE)
