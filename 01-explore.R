
library(dplyr)

getwd()
setwd("/Users/koe1/GitHub/crisprWTcas9")
data <- read.csv(file="WT_Cas9_x_IgG_cd47.csv", sep=",", header=TRUE)
expr <- "(\\w)*(?=_[-+]_[\\d\\w.-]+)"
grep(expr, head(data$sgRNA), perl=TRUE, value=TRUE)

expr <- "(_[-+]_[\\d\\w.-]+)"
data$gene <- gsub(expr,replacement="",data$sgRNA,perl=TRUE)


data <- data[data$sum > 0, ]
hist(data$sum[data$sum < 20])

nonzero_data_backup <- data

data[data$sum > 0, "sum"] %>% summary()
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 1.0       2.0       5.0    1959.3     759.8 1104308.0

data[data$sum > 20, "sum"] %>% summary()
# Min. 1st Qu.  Median    Mean  3rd Qu.    Max. 
# 21     414    1058      4703  2050       1104308 




do_counts <-  function(data, sgRNA_ct_thresh) {
  gene_counts <- as.data.frame((table(data$gene)))
  gene_counts <- gene_counts[order(gene_counts$Freq, decreasing=TRUE),]

  print(paste0("how many genes: ", length(unique(data$gene))))
  print(paste0("how many w/ > ",sgRNA_ct_thresh," sgRNA detected: ", sum(gene_counts$Freq > sgRNA_ct_thresh)))
  table(gene_counts$Freq)
}

data[data$sum > 50, ] %>% do_counts(sgRNA_ct_thresh=1)
# [1] "how many genes: 1447"            
# [1] "how many w/ >1 sgRNA detected: 339"                            
# 
#     1    2    3    4    5 
#  1048  335   56    7    1 

data[data$sum > 20, ] %>% do_counts(sgRNA_ct_thresh=1)
# [1] "how many genes: 1475"
# [1] "how many w/ >1 sgRNA detected: 412"
# 
#     1    2    3    4    5 
#  1063  345   59    6    2 

data[data$sum > 5, ] %>% do_counts(sgRNA_ct_thresh=1)
# [1] "how many genes: 1653"
# [1] "how many w/ >1 sgRNA detected: 530"
# 
#     1    2    3    4    5    7 
#  1123  417  101    9    2    1 

data[data$sum > 50, ] %>% do_counts(sgRNA_ct_thresh=3)
# [1] "how many genes: 1447"
# [1] "how many w/ > 3 sgRNA detected: 8"
# 
#    1    2    3    4    5 
# 1048  335   56    7    1 

data[data$sum > 20, ] %>% do_counts(sgRNA_ct_thresh=3)
# [1] "how many genes: 1475"
# [1] "how many w/ > 3 sgRNA detected: 8"
# 
#     1    2    3    4    5 
#  1063  345   59    6    2 

data[data$sum > 5, ] %>% do_counts(sgRNA_ct_thresh=3)
# [1] "how many genes: 1653"
# [1] "how many w/ > 3 sgRNA detected: 12"
# 
#     1    2    3    4    5    7 
#  1123  417  101    9    2    1 

library(knitr) 
stitch("01-explore.R")
