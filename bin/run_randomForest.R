library(randomForest)

argsIn <- commandArgs(trailingOnly = TRUE)
load(argsIn[1])

#argsIn[1] <- "/data/hpc/cog_bioinf/kloosterman/users/jvalleinclan/nanopore/lowpassSV/emt026/features_table.txt"
#argsIn[2] <- 2

ftab <- read.table(argsIn[2], sep="\t", header=T)
mean_cov <- as.integer(argsIn[3])

total_cov <- (ftab$dr1 + ftab$dr2 + ftab$dv1 + ftab$dv2)
ftab$total_cov_norm <- round(total_cov/mean_cov)
ftab$vaf <- (ftab$dv1 + ftab$dv2)/total_cov
ftab <- ftab[,-which(colnames(ftab) %in% c("dr1", "dr2", "dv1", "dv2"))]

test.forest <- predict(output.forest, newdata=ftab)
ftab$pred <- predict(output.forest, newdata=ftab)
table(test.forest)

write.table(ftab[,c("id","pred")], paste(dirname(argsIn[1]),"/predict_labels.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
