library(readxl)
library(data.table)
CONPenh <- read_excel("aat6720_CONPPutEnhs.xlsx", sheet = 2,header=TRUE)
> head(CONPenh)
# A tibble: 6 x 16
  Chrom  Start    End CONP_ID All_OP_No. iPSC_OP_No. TD0_OP_No. TD11_OP_No.
  <chr>  <dbl>  <dbl>   <dbl>      <dbl>       <dbl>      <dbl>       <dbl>
1 chr1  817175 817393       2          2           0          0           0
2 chr1  818953 819108       3          2           0          1           1
3 chr1  907577 908523       8          8           0          0           3
4 chr1  913787 914227      11          2           0          0           0
5 chr1  915886 918310      12         27           0          7          10
6 chr1  994848 995028      47          3           0          1           0
# ... with 8 more variables: TD30_OP_No. <dbl>, CTX1_OP_No. <dbl>,
#   CTX2_OP_No. <dbl>, TD0_annotation <chr>, TD11_annotation <chr>,
#   TD30_annotation <chr>, CTX1_annotation <chr>, CTX2_annotation <chr>

write.table(CONPenh,"AmiriEtAl_PutEnh.bed",quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")