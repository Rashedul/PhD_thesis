---
title: "ighv-status"
author: "Rashedul"
date: "3/23/2020"
output: 
  html_document: 
    keep_md: yes
---



# confusion matrix


```r
library(caret)
library(gdata)
library(tidyverse)

#test
numLlvs <- 2
confusionMatrix(
   factor(sample(rep(letters[1:numLlvs], 200), 50)),
   factor(sample(rep(letters[1:numLlvs], 200), 50)))
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction  a  b
##          a 13 16
##          b 12  9
##                                           
##                Accuracy : 0.44            
##                  95% CI : (0.2999, 0.5875)
##     No Information Rate : 0.5             
##     P-Value [Acc > NIR] : 0.8389          
##                                           
##                   Kappa : -0.12           
##                                           
##  Mcnemar's Test P-Value : 0.5708          
##                                           
##             Sensitivity : 0.5200          
##             Specificity : 0.3600          
##          Pos Pred Value : 0.4483          
##          Neg Pred Value : 0.4286          
##              Prevalence : 0.5000          
##          Detection Rate : 0.2600          
##    Detection Prevalence : 0.5800          
##       Balanced Accuracy : 0.4400          
##                                           
##        'Positive' Class : a               
## 
```

```r
x = read.xls("~/Documents/CRIS/comparisons/dbgap_CRIS.xlsx", sheet = 1)
x2 = data.frame(CRIS = x$dbGAP, Sanger = x$CRIS) %>% na.omit()

confusionMatrix(factor(x2$Sanger), factor(x2$CRIS)) #sanger as referene
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction  0  1
##          0 22  1
##          1  0 36
##                                           
##                Accuracy : 0.9831          
##                  95% CI : (0.9091, 0.9996)
##     No Information Rate : 0.6271          
##     P-Value [Acc > NIR] : 3.989e-11       
##                                           
##                   Kappa : 0.9641          
##                                           
##  Mcnemar's Test P-Value : 1               
##                                           
##             Sensitivity : 1.0000          
##             Specificity : 0.9730          
##          Pos Pred Value : 0.9565          
##          Neg Pred Value : 1.0000          
##              Prevalence : 0.3729          
##          Detection Rate : 0.3729          
##    Detection Prevalence : 0.3898          
##       Balanced Accuracy : 0.9865          
##                                           
##        'Positive' Class : 0               
## 
```

```r
x = read.xls("~/Documents/CRIS/comparisons/EGA_CRIS.xlsx", sheet = 1)
x2 = data.frame(CRIS = x$CRIS, Sanger = x$VDJER) %>% na.omit()

confusionMatrix(factor(x2$Sanger), factor(x2$CRIS)) #sanger as referene
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction 0 1
##          0 2 0
##          1 0 5
##                                      
##                Accuracy : 1          
##                  95% CI : (0.5904, 1)
##     No Information Rate : 0.7143     
##     P-Value [Acc > NIR] : 0.09486    
##                                      
##                   Kappa : 1          
##                                      
##  Mcnemar's Test P-Value : NA         
##                                      
##             Sensitivity : 1.0000     
##             Specificity : 1.0000     
##          Pos Pred Value : 1.0000     
##          Neg Pred Value : 1.0000     
##              Prevalence : 0.2857     
##          Detection Rate : 0.2857     
##    Detection Prevalence : 0.2857     
##       Balanced Accuracy : 1.0000     
##                                      
##        'Positive' Class : 0          
## 
```

```r
#PNAS
x = read.xls("~/Documents/CRIS/comparisons/PNAS_CRIS.xlsx", sheet = 1)
x2 = data.frame(CRIS = x$BWA_align, Sanger = x$STAR_align) %>% na.omit()

confusionMatrix(factor(x2$CRIS), factor(x2$Sanger))
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction  0  1
##          0 11  0
##          1  0  6
##                                      
##                Accuracy : 1          
##                  95% CI : (0.8049, 1)
##     No Information Rate : 0.6471     
##     P-Value [Acc > NIR] : 0.000611   
##                                      
##                   Kappa : 1          
##                                      
##  Mcnemar's Test P-Value : NA         
##                                      
##             Sensitivity : 1.0000     
##             Specificity : 1.0000     
##          Pos Pred Value : 1.0000     
##          Neg Pred Value : 1.0000     
##              Prevalence : 0.6471     
##          Detection Rate : 0.6471     
##    Detection Prevalence : 0.6471     
##       Balanced Accuracy : 1.0000     
##                                      
##        'Positive' Class : 0          
## 
```

# VDJ-region assembly completeness CEMT samples


```r
library(gdata)
library(dplyr)
library(tidyverse)

g = read.table("~/Documents/research/cll/Blueprint_dbGAP/Request2/dbGAP/igBLAST_out/ighvdj_IMGT_db.fa_seqId_len.txt", fill = T)

colnames(g) = c("gene", "length_germline")

#v-gene
x = read.xls("~/Documents/research/IGHV-status/IGHV-status_manuscript/data/IGHV_PNAS_comparison.xlsx", sheet = 2)
#
x2 = left_join(g, x, by = c( "gene" = "IGHV.gene" ))
v = tibble(id = x2$CEMT_ID, gene = x2$gene, fraction = as.numeric(as.character(x2$IGHV.len)) / x2$length_germline) %>% na.omit()

#d-gene
x2 = left_join(g, x, by = c( "gene" = "IGHD.gene" ))
d = tibble(id = x2$CEMT_ID, gene = x2$gene, fraction = as.numeric(as.character(x2$IGHD.len)) / x2$length_germline) %>% na.omit()

#j-gene
x2 = left_join(g, x, by = c( "gene" = "IGHJ.gene" ))
j = tibble(id = x2$CEMT_ID, gene = x2$gene, fraction = as.numeric(as.character(x2$IGHJ.len)) / x2$length_germline) %>% na.omit()

y = rbind(v,d,j) %>% mutate(Ig = ifelse(grepl("IGHV", gene, ignore.case = T), "IGHV",
                                        ifelse(grepl("IGHD", gene, ignore.case = T), "IGHD", 
                                        "IGHJ")))

y %>%
    ggplot(aes(Ig, fraction)) +
    geom_boxplot() +
    geom_jitter(width = 0.10) 
```

![](chapter-2_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
v_cemt = v

# V-region assembly completeness pnas
library(gdata)
library(dplyr)
library(tidyverse)

g = read.table("~/Documents/research/cll/Blueprint_dbGAP/Request2/dbGAP/igBLAST_out/ighvdj_IMGT_db.fa_seqId_len.txt", fill = T)
colnames(g) = c("gene", "length_germline")

#v-gene

x = read.xls("~/Documents/research/IGHV-status/IGHV-status_manuscript/data/IGHV_PNAS_comparison.xlsx", sheet = 1)

#
x2 = left_join(g, x, by = c( "gene" = "IGHV.gene" ))
v = tibble(id = x2$RNA.ID, gene = x2$gene, fraction = as.numeric(as.character(x2$IGHV.len)) / x2$length_germline) %>% na.omit()

#d-gene
x2 = left_join(g, x, by = c( "gene" = "IGHD.gene" )) 
d = tibble(id = x2$RNA.ID, gene = x2$gene, fraction = as.numeric(as.character(x2$IGHD.len)) / x2$length_germline) %>% na.omit()

#j-gene
x2 = left_join(g, x, by = c( "gene" = "IGHJ.gene" ))
j = tibble(id = x2$RNA.ID, gene = x2$gene, fraction = as.numeric(as.character(x2$IGHJ.len)) / x2$length_germline) %>% na.omit()

y = rbind(v,d,j) %>% mutate(Ig = ifelse(grepl("IGHV", gene, ignore.case = T), "IGHV",
                                        ifelse(grepl("IGHD", gene, ignore.case = T), "IGHD", 
                                        "IGHJ")))

y %>%
    ggplot(aes(Ig, fraction)) +
    geom_boxplot() +
    geom_jitter(width = 0.10)
```

![](chapter-2_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

```r
v_pnas = v %>% filter(id != "US-1422294")

vy = rbind(v_cemt, v_pnas) %>% mutate(Sample = ifelse(grepl("CEMT", id, ignore.case = T), "CEMT", "GSE66228"))

vy %>%
    ggplot(aes(Sample, fraction)) +
    geom_boxplot(color = "gray", outlier.colour = "white") +
    geom_jitter(width = 0.10, shape = 2) +
    ylim(0,1.1) +
    theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.text = element_text( size = 10),
        strip.background =element_rect(fill="white"),
        text = element_text(size=10),
        axis.text = element_text(color = "black"))
```

![](chapter-2_files/figure-html/unnamed-chunk-2-3.png)<!-- -->

```r
ggsave("~/Documents/research/IGHV-status/IGHV-status_manuscript/Illustrator_IGHV/IGHV_plots/IGHV_len.pdf", width = 8, height = 10, units = "cm")

#NS
t.test(v_cemt$fraction, v_pnas$fraction)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  v_cemt$fraction and v_pnas$fraction
## t = -1.5061, df = 15.08, p-value = 0.1527
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.06276599  0.01077692
## sample estimates:
## mean of x mean of y 
## 0.9714906 0.9974851
```

# Benchmark: VDJ-region assembly completeness 


```r
library(gdata)
library(dplyr)
library(tidyverse)
library(reshape2)

g = read.table("~/Documents/research/cll/Blueprint_dbGAP/Request2/dbGAP/igBLAST_out/ighvdj_IMGT_db.fa_seqId_len.txt", fill = T)

colnames(g) = c("gene", "length_germline")

#v-gene
x = read.xls("~/Documents/CRIS/comparisons/PNAS_CRIS.xlsx", sheet = 1)
x2 = data.frame(sample = x$RNA.ID, gene =x$IGHV.gene , CRIS = x$IGHV.len, VDJer = x$IGHV_len_VDJER) %>% na.omit() 
#
x3 = left_join(x2, g) %>% 
    mutate(CRIS_f = CRIS/length_germline) %>%
    mutate(VDJer_f = VDJer/length_germline) %>% 
    select(sample, CRIS_f, VDJer_f) %>%
    melt()

#mean
cris = 0.9974343
vdjer = 0.7543688

x3 %>% filter(sample != "US-1422294") %>%
    ggplot(aes(sample, 100*(value), fill = variable)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_hline(yintercept = c(100), color = "gray", linetype = "dashed") +
    #geom_hline(yintercept = c(99.74343), color = "#F8766D", linetype = "dashed") +
    #geom_hline(yintercept = c( 75.43688), color = "#00BFC4", linetype = "dashed") + 
    #scale_fill_brewer(palette="Set1") +
     theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.text = element_text( size = 10),
        strip.background =element_rect(fill="white"),
        text = element_text(size=10),
        axis.text = element_text(color = "black", angle = 90))
```

![](chapter-2_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
ggsave("~/Documents/CRIS/plot/PNAS_VDJER_CRIS_IGHV_len.pdf", width = 30, height = 14, units = "cm")

#decided not to use it

y = read.xls("~/Documents/research/IGHV-status/script/comparisons/PNAS_CRIS.xlsx", sheet = 1) %>% na.omit()
colnames(y)
```

```
##  [1] "RNA.ID"           "SRR"              "IGHV.gene.Sanger" "Sanger_IGHV"     
##  [5] "CRIS_V_gene"      "CRIS_IGHV"        "CRIS_len"         "clonotypes"      
##  [9] "transcript"       "multiclonal"      "BWA_align"        "VDJER"           
## [13] "VDJER_len"        "VDJER...identity" "TRUST"            "TRUST_len"       
## [17] "TRUST_status"
```

```r
#
y2 = left_join(y, g, by = c("CRIS_V_gene" = "gene")) %>% 
    mutate(CRIS_f = CRIS_len/length_germline) %>%
    mutate(VDJer_f = VDJER_len/length_germline) %>% 
    mutate(TRUST_f = TRUST_len/length_germline) %>% 
    select(RNA.ID, CRIS_f, VDJer_f, TRUST_f) %>%
    melt()

#trust mean = 0.4152447

y2 %>% 
    ggplot(aes(RNA.ID, 100*(value), fill = variable)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_hline(yintercept = c(99.74343), color = "#F8766D", linetype = "dashed") +
    geom_hline(yintercept = c(75.43688), color = "green", linetype = "dashed") + 
    geom_hline(yintercept = c(41.52447), color = "blue", linetype = "dashed") + 
    #scale_fill_brewer(palette="Set1") +
    xlab("") +
    ylab("") +
     theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.text = element_text( size = 10),
        strip.background =element_rect(fill="white"),
        text = element_text(size=10),
        axis.text = element_text(color = "black", angle = 90))
```

![](chapter-2_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
ggsave("~/Documents/research/IGHV-status/plot/Fraction_IGHV_len.pdf", width = 30, height = 15, units = "cm")
```

# ighv % correlation and deviation


```r
library(gdata)
library(dplyr)
library(tidyverse)

x = read.xls("~/Documents/research/IGHV-status/IGHV-status_manuscript/data/IGHV_PNAS_comparison.xlsx", sheet = 1)


x %>% 
    ggplot(aes(Sanger_IGHV, CRIS_IGHV)) +
    geom_abline(color = "gray") +
    geom_point(shape = 2) +
    xlim(-1,11) +
    ylim(-1,11) +
    theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.text = element_text( size = 10),
        strip.background =element_rect(fill="white"),
        text = element_text(size=10),
        axis.text = element_text(color = "black"))
```

![](chapter-2_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
ggsave("~/Documents/research/IGHV-status/IGHV-status_manuscript/Illustrator_IGHV/IGHV_plots/Corr.pdf", width = 8, height = 8, units = "cm")


#cor
x2 = x %>% select(Sanger_IGHV, CRIS_IGHV) %>% na.omit()
cor(x2$Sanger_IGHV, x2$CRIS_IGHV, method = "pearson")
```

```
## [1] 0.9479579
```

```r
#dev
sd(x2$Sanger_IGHV - x2$CRIS_IGHV)
```

```
## [1] 1.050976
```

```r
#confidence interval
library(Rmisc)
CI(x2$Sanger_IGHV - x2$CRIS_IGHV)
```

```
##      upper       mean      lower 
##  0.7603619  0.2200000 -0.3203619
```

```r
cor.test(x2$Sanger_IGHV, x2$CRIS_IGHV)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  x2$Sanger_IGHV and x2$CRIS_IGHV
## t = 11.531, df = 15, p-value = 7.439e-09
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.8584492 0.9814319
## sample estimates:
##       cor 
## 0.9479579
```

```r
#EGA

x = read.xls("~/Documents/CRIS/comparisons/EGA_CRIS.xlsx", sheet = 1)


x %>% 
    ggplot(aes(100-Sanger_IGHV, 100-CRIS_IGHV)) +
    geom_abline(color = "gray") +
    geom_point(shape = 2) +
    xlim(-1,11) +
    ylim(-1,11) +
    geom_text(aes(label=Sanger),hjust=0, vjust=0) +
    theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.text = element_text( size = 8),
        strip.background =element_rect(fill="white"),
        text = element_text(size=8),
        axis.text = element_text(color = "black"))
```

![](chapter-2_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```r
ggsave("~/Documents/research/IGHV-status/IGHV-status_manuscript/Illustrator_IGHV/IGHV_plots/Corr_EGA_v2.pdf", width = 14, height = 14, units = "cm")

#cor
x2 = x %>% select(Sanger_IGHV, CRIS_IGHV) %>% na.omit()
cor(100-x2$Sanger_IGHV, 100-x2$CRIS_IGHV, method = "pearson")
```

```
## [1] 0.7985379
```

```r
#dev
sd(x2$Sanger_IGHV - x2$CRIS_IGHV)
```

```
## [1] 2.036049
```

```r
#confidence interval
library(Rmisc)
CI(x2$Sanger_IGHV - x2$CRIS_IGHV)
```

```
##      upper       mean      lower 
##  0.9173169 -0.9657143 -2.8487455
```

```r
cor.test(x2$Sanger_IGHV, x2$CRIS_IGHV)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  x2$Sanger_IGHV and x2$CRIS_IGHV
## t = 2.9664, df = 5, p-value = 0.03129
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.1140833 0.9689327
## sample estimates:
##       cor 
## 0.7985379
```
