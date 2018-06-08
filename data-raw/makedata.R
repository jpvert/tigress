ecoli <- list(exp=read.csv('E_coli_v3_Build_1_chips445probes4345.tab', header=TRUE, sep="\t", row.names=1), reg=read.csv('regulation.txt', header=TRUE, sep="\t"))
devtools::use_data(ecoli)
