library(grofit)
library(tidyr)
library(dplyr)

growth <-
    gather(grofit.data, timepoint, od, -V1, -V2, -V3, -strain) %>%
    mutate(timepoint=as.numeric(gsub("V", "", timepoint)) - 3) %>%
    mutate(medium=V2) %>%
    select(-V1, -V2) %>%
    write.table(file="yeast-growth.csv", sep=",", quote=FALSE, row.names=FALSE)

uni <- read.table("ecoli.txt", sep=",", header=TRUE, stringsAsFactors=FALSE)

ecoli <- group_by(uni, bnumber) %>%
    do({
        data.frame(bnumber=.$bnumber, genes=unlist(strsplit(.$symbol, ";")),
                   stringsAsFactors=FALSE)
    })
write.table(ecoli, file="ecoli.csv", quote=FALSE, row.names=FALSE, sep=",")
