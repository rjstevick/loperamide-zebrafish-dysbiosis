---
title: "Water Loperamide CFUs"
author: "Rebecca Stevick"
output:
  html_document:
    toc: true
    keep_md: TRUE
    theme: "cerulean"
    toc_float:
      collapsed: false
      smooth_scroll: false
---

## About the Data

Individual strains survival in water.

20 mL of volvic water with each strain at final OD=0.001 (~10^6^ CFU/mL)\
DMSO at 1:10000, Loperamide at 10 mg/L

### Timepoints

Add loperamide and strains to water.
Sample water at 0 h, 3 h, 6 hours on day 1.
Then, daily until 3 days post inoculation. 24 h, 48 h, 72 h.


### Sample collection & plating

At each timepoint:

Make 0 to -3 dilutions in 96-well plates, in triplicate\
Plate 10 µL microdrops on big square plates.\
2 square plates total per timepoint

Put plates at 28C for 2 days, then count colonies.

### Conditions

**Treatment (Strain)**: Bc1, Bc2, Bc3, Bc4, Bc10, W6, W8, Mz1, Mz8, Fjohn

**Loperamide Treatment**: Control (water only), DMSO, Loperamide 10 mg/L



## Setup

### Load libraries and settings


```r
library(tidyverse)
library(scales)
library(ggpubr)
library(ungeviz)
library(ggtext)
library(rmdformats)

# set global theme
theme_set(theme_minimal()+
theme(panel.grid.major.y = element_line(color="grey80"), strip.text=element_text(size=16),
        strip.text.y = element_text(angle=0), plot.caption = element_text(size=10),
        panel.grid.major.x = element_blank(),legend.position="top",
        plot.background = element_rect(fill="transparent", color="transparent"),
        axis.ticks = element_line(inherit.blank = FALSE),
        panel.background = element_rect(color="grey50", size=2),
        legend.title = element_text(size=18),
        axis.text = element_text(size=15), axis.title = element_text(size=18),
        legend.text = element_text(size=16), plot.title = element_text(hjust=0.5)))

# set global options for code output
knitr::opts_chunk$set(echo=TRUE, warning=FALSE,message=FALSE)
```

### Import data


```r
# import individual strain data from exp 1
datacfus1 <-
   readxl::read_xlsx("LoperamideWaterSurvivalCFUs_sub.xlsx", sheet="Round1") %>%
   drop_na(DF) %>%
   mutate(LoperamideTreatment=factor(LoperamideTreatment,
                                    levels=c("None", "DMSO", "Loperamide 10 mg/L"),
                                    labels=c("Control","DMSO", "Loperamide")),
         Treatment = factor(Treatment, levels=c("Bc1","Bc2","Bc3","Bc4","Bc10"))) %>%
   add_column(Assay=1)

# import individual strain data from exp 2
datacfus2 <-
   readxl::read_xlsx("LoperamideWaterSurvivalCFUs_sub.xlsx", sheet="Round2") %>%
   drop_na(DF) %>%
   mutate(LoperamideTreatment=factor(LoperamideTreatment,
                                    levels=c("None", "DMSO", "Loperamide 10 mg/L"),
                                    labels=c("Control","DMSO", "Loperamide")),
         Treatment = factor(Treatment, levels=c("Bc1","Bc2","Bc3","Bc4","Bc10"))) %>%
   add_column(Assay=2)

# import individual strain data from exp 3
datacfus3 <-
   readxl::read_xlsx("LoperamideWaterSurvivalCFUs_sub.xlsx", sheet="Round3") %>%
   drop_na(DF) %>%
   mutate(LoperamideTreatment=factor(LoperamideTreatment,
                                    levels=c("None", "DMSO", "Loperamide 10 mg/L"),
                                    labels=c("Control","DMSO", "Loperamide")),
         Treatment = factor(Treatment, levels=c("Bc1","Bc2","Bc3","Bc4","Bc10"))) %>%
   add_column(Assay=3)


# import individual strain data from exp 4
datacfus4 <-
   readxl::read_xlsx("LoperamideWaterSurvivalCFUs_sub.xlsx", sheet="Round4") %>%
   drop_na(DF) %>%
   mutate(LoperamideTreatment=factor(LoperamideTreatment,
                                    levels=c("None", "DMSO", "Loperamide 10 mg/L"),
                                    labels=c("Control","DMSO", "Loperamide")),
         Treatment = factor(Treatment, levels=c("W6","W8","Mz8")))

# import individual strain data from exp 5
datacfus5 <-
   readxl::read_xlsx("LoperamideWaterSurvivalCFUs_sub.xlsx", sheet="Round5") %>%
   drop_na(DF) %>%
   mutate(LoperamideTreatment=factor(LoperamideTreatment,
                                    levels=c("None", "DMSO", "Loperamide 10 mg/L"),
                                    labels=c("Control","DMSO", "Loperamide")),
         Treatment = factor(Treatment, levels=c("W6","W8","Mz8","Fjohn","Mz1")))


straininfo <- readxl::read_xlsx("../../LoperamideStrainInfo.xlsx") %>%
   mutate(Strain=recode(Strain, "W6t"="W6"))


dataallCFUs <- full_join(datacfus1, datacfus2) %>%
   full_join(datacfus3) %>%
   full_join(datacfus4) %>% full_join(datacfus5) %>%
   left_join(straininfo, by=c("Treatment"="Strain")) %>%
   mutate(CodeName=factor(CodeName, levels=unique(straininfo$CodeName)))
```

------------------------------------------------------------------------


## Stats of all significant comparisons


```r
longdata <- dataallCFUs %>%
  pivot_longer(Rep1:Rep3) %>%
  mutate(RepCFUs=(1000/VolPlated_ul)*DF*value) %>%
  select(Treatment, LoperamideTreatment, Timepoint_hrs, RepCFUs)

statsbydayCFUs <- compare_means(data=longdata,
                            RepCFUs~LoperamideTreatment,
                            group.by = c("Treatment","Timepoint_hrs")) %>%
   left_join(straininfo, by=c("Treatment"="Strain")) %>%
   mutate(CodeName=factor(CodeName, levels=unique(straininfo$CodeName)))
```


## All together now

### Timeline for control conditions


```r
dataallCFUs %>% filter(LoperamideTreatment=="Control") %>%
  pivot_longer(Rep1:Rep3) %>%
  mutate(RepCFUs=(1000/VolPlated_ul)*DF*value) %>%
  ggplot(aes(x=Timepoint_hrs, y=RepCFUs,
             color=PaperCode, fill=PaperCode, shape=PaperCode))+
  stat_summary(geom="errorbar", fun.data="mean_sd", width=0.5, show.legend = FALSE)+
  stat_summary(geom="point", fun="mean", size=3) +
  stat_summary(aes(lty=PaperCode), geom="line", fun="mean", lwd=0.8) +
  scale_y_continuous(trans = 'log10', labels = trans_format('log10', math_format(10^.x)))+
  scale_shape_manual(values=c(21,22,23,24,25,21,22,23,24,25))+
  scale_fill_manual(values=c(nationalparkcolors::park_palette("Saguaro", n=6),
                             "grey","navy","lightcoral","deepskyblue"))+
  scale_color_manual(values=c(nationalparkcolors::park_palette("Saguaro", n=6),
                              "grey","navy","lightcoral","deepskyblue"))+
  labs(y="CFUs per mL", x="Time (hours)", color=NULL, shape=NULL, fill=NULL, lty=NULL,
       title="Survival in water, control conditions")+
  theme(legend.key.width = unit(1.33,"cm"), legend.position = "right")
```

![](WaterSurvivalCFUs_loperamide_figure3_files/figure-html/timelineplotcontrolall-1.png)<!-- -->


### Timeline for each strain with mean


```r
statsbydayCFUsformat <- statsbydayCFUs %>% filter(p.signif!="ns")

dataallCFUs %>% pivot_longer(Rep1:Rep3) %>%
  mutate(RepCFUs=(1000/VolPlated_ul)*DF*value) %>%
  ggplot(aes(x = Timepoint_hrs, y=RepCFUs,
             color=LoperamideTreatment, shape=LoperamideTreatment))+
  facet_wrap(.~CodeName, ncol=5)+
  stat_summary(geom="errorbar", fun.data="mean_sd", width=1, show.legend = FALSE)+
  stat_summary(geom="point", fun="mean", size=4) +
  stat_summary(aes(lty=LoperamideTreatment), geom="line", fun="mean", lwd=1) +
  geom_text(data=statsbydayCFUsformat, aes(label="*", y=5e6, color=NA, shape=NA),
            size=10, color="#0fc08e",
            show.legend=FALSE)+
  scale_y_continuous(trans = 'log10', limits=c(NA,1e7),
                     labels = trans_format('log10', math_format(10^.x)))+
   scale_x_continuous(breaks=c(0,6,24,48,72))+
  scale_linetype_manual(values=c("solid","dotted","dashed"))+
  scale_color_manual(values=c('#000000', '#1c5580', '#0fc08e'))+
  labs(y="CFUs per mL", x="Time since treatment (hours)",
       color=NULL, shape=NULL, lty=NULL)+
  theme(strip.text = element_markdown(size = 16, face="bold"),
        panel.grid.major.x = element_line(size=0.5),
        panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),
        legend.text = element_text(size=18), legend.key.width = unit(2,"cm"))
```

![](WaterSurvivalCFUs_loperamide_figure3_files/figure-html/timelineplotall-1.png)<!-- -->

```r
ggsave("Figure3_WaterInVitro_LoperamideCFUs.png", bg="transparent", width = 19, height = 9.5)
ggsave("Figure3_WaterInVitro_LoperamideCFUs.tiff", bg="transparent", width = 19, height = 9.5)
ggsave("Figure3_WaterInVitro_LoperamideCFUs.pdf", bg="transparent", width = 19, height = 9.5)
```



```r
sessionInfo()
```

```
## R version 4.1.3 (2022-03-10)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur/Monterey 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] rmdformats_1.0.3 ggtext_0.1.1     ungeviz_0.1.0    ggpubr_0.4.0    
##  [5] scales_1.1.1     forcats_0.5.1    stringr_1.4.0    dplyr_1.0.8     
##  [9] purrr_0.3.4      readr_2.1.2      tidyr_1.2.0      tibble_3.1.6    
## [13] ggplot2_3.3.5    tidyverse_1.3.1 
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.8.3             lubridate_1.8.0          assertthat_0.2.1        
##  [4] digest_0.6.29            utf8_1.2.2               R6_2.5.1                
##  [7] cellranger_1.1.0         backports_1.4.1          reprex_2.0.1            
## [10] evaluate_0.15            highr_0.9                httr_1.4.2              
## [13] pillar_1.7.0             rlang_1.0.2              readxl_1.4.0            
## [16] rstudioapi_0.13          car_3.0-12               jquerylib_0.1.4         
## [19] nationalparkcolors_0.1.0 rmarkdown_2.13           labeling_0.4.2          
## [22] gridtext_0.1.4           munsell_0.5.0            broom_0.7.12            
## [25] compiler_4.1.3           modelr_0.1.8             xfun_0.30               
## [28] pkgconfig_2.0.3          htmltools_0.5.2          tidyselect_1.1.2        
## [31] bookdown_0.25            fansi_1.0.3              crayon_1.5.1            
## [34] tzdb_0.3.0               dbplyr_2.1.1             withr_2.5.0             
## [37] grid_4.1.3               jsonlite_1.8.0           gtable_0.3.0            
## [40] lifecycle_1.0.1          DBI_1.1.2                magrittr_2.0.3          
## [43] cli_3.2.0                stringi_1.7.6            carData_3.0-5           
## [46] farver_2.1.0             ggsignif_0.6.3           fs_1.5.2                
## [49] xml2_1.3.3               bslib_0.3.1              ellipsis_0.3.2          
## [52] generics_0.1.2           vctrs_0.4.0              tools_4.1.3             
## [55] strapgod_0.0.4.9000      glue_1.6.2               markdown_1.1            
## [58] hms_1.1.1                abind_1.4-5              fastmap_1.1.0           
## [61] yaml_2.3.5               colorspace_2.0-3         rstatix_0.7.0           
## [64] rvest_1.0.2              knitr_1.38               haven_2.4.3             
## [67] sass_0.4.1
```
