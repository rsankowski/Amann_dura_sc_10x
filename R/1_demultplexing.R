library(tidyverse)
library(deMULTIplex)

pth <- file.path("data","counts")
fls <-paste("GL",1:10,sep = "-")

## load cell IDs
cellid_lst <- map(fls, function(x) {
  .df <- read.table(file.path(pth,x,"outs","filtered_feature_bc_matrix","barcodes.tsv.gz"))
  .df$V1 <- gsub("-1","",.df$V1)
  .df[[1]]
})

names(cellid_lst) <- fls

## load hashes
hashes <- read.csv(file.path("data","hashes.csv"), header = F, row.names = 1) #, col.names = c("sequence")

barcodes_lst <- lst()
barcodes_lst[[1]] <- hashes[[1]]#[1:3]
barcodes_lst[[2]] <- hashes[[1]]#[4:6]
barcodes_lst[[3]] <- hashes[[1]]#[7:9]
barcodes_lst[[4]] <- hashes[[1]]#[1:3]
barcodes_lst[[5]] <- hashes[[1]]#[4:6]
barcodes_lst[[6]] <- hashes[[1]]#[7:9]
barcodes_lst[[7]] <- hashes[[1]]#[10:12]
barcodes_lst[[8]] <- hashes[[1]]#[1:3]
barcodes_lst[[9]] <- hashes[[1]]#[10:12]
barcodes_lst[[10]] <- hashes[[1]]#[1:3]

names(barcodes_lst) <- fls

## Pre-process MULTI-seq sample barcode FASTQs
readTable <- map(1:10, function(x) {
  df <- MULTIseq.preProcess(R1 = file.path("data","fastq",paste0("CMO_",fls[x],"_S",x,"_R1_001.fastq.gz")), R2 = file.path("data","fastq",paste0("CMO_",fls[x],"_S",x,"_R2_001.fastq.gz")), cellIDs = cellid_lst[[x]], cell=c(1,16), umi=c(17,28), tag=c(1,8))
  })
names(readTable) <- fls

bar.table <- map(fls, function(x) {
  MULTIseq.align(readTable[[x]],cellid_lst[[x]],barcodes_lst[[x]])
})

bar.tsne <- barTSNE(bar.table[[9]][,1:12]) 

## plot
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

readTable2 <- bind_rows(readTable)
