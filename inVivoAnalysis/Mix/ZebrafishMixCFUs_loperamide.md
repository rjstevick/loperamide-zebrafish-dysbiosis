---
title: "Gnotobiotic Mix-Reconv Loperamide CFUs"
author: "Rebecca Stevick/David Perez Pascual"
output:
  html_document:
    toc: true
    keep_md: TRUE
    theme: "cerulean"
    toc_float:
      collapsed: false
      smooth_scroll: false
---

# About the Data

## Timepoints

Treat with loperamide at 5 dpf for 24 hours.

-   Sample timepoint 1 at 6 dpf (24 hour treatment)
-   Sample timepoint 2 at 7 dpf (24 hour treatment + 24 hour water)
-   Sample timepoint 3 at 11 dpf (24 hour treatment + 5 days water)

### Sample collection & plating

At each timepoint:

1.  Wash all fish twice by transferring into sterile volvic in a 6-well plate
2.  Add fish with 500 µL sterile volvic water into a fastprep tube
3.  Homogenize sample at 6.5 for 45 seconds

Make 3 dilutions in 1.5 mL tubes: 100 µL into 900 µL water --\> 10^0^, 10^-1^, 10^-2^, 10^-3^\
Spread 3 x 100 µL aliquots of each dilution on LB plates = 12 plates per fish\
144 plates total per timepoint

Put plates at 28C for 2 days, then count colonies.

### Conditions

|  Treatment (Strain)  | Loperamide Treatment |
|:--------------------:|:--------------------:|
| Bc1/Bc2/Bc3/Bc4/Bc10 |                      |
| Bc1/Bc2/Bc3/Bc4/Bc10 |         DMSO         |
| Bc1/Bc2/Bc3/Bc4/Bc10 |  Loperamide 10 mg/L  |

# Setup

## Load libraries



## Import data


```r
# import mix data
datacfusmix <- 
   readxl::read_xlsx("MixCFUs_LoperamideZebrafish.xlsx", sheet="FishMix") %>%
   drop_na(DF) %>%
   mutate(LoperamideTreatment=factor(LoperamideTreatment, 
                                     levels=c("None", "DMSO", "Loperamide 10 mg/L"),
                                     labels=c("Control","DMSO", "Loperamide"))) %>% 
  group_by(Fish,Treatment,LoperamideTreatment,FishNum,TrialDay, Timepoint) %>% 
  summarise_all(.funs="mean", na.rm=TRUE) %>% 
   unite("LoperamideTimepoint", LoperamideTreatment,Timepoint, remove=FALSE) %>% 
   unite("FishID", LoperamideTreatment,Timepoint,FishNum, remove=FALSE)


# import strain metadata
straininfo <- readxl::read_xlsx("../../LoperamideStrainInfo.xlsx")
```

------------------------------------------------------------------------


# Mix fish colonization

## Total CFUs


```r
mixtotalCFUplot <- datacfusmix %>% 
  ggplot(aes(x=Timepoint, y=CFUs_perFish, 
             fill=LoperamideTreatment, color=LoperamideTreatment, shape=LoperamideTreatment))+
  geom_boxplot(color="black", alpha=0.8, show.legend = FALSE)+
  geom_point(size=2, position=position_jitterdodge(jitter.width=0.3)) +
   scale_color_manual(values=c('#000000', '#1c5580', '#0fc08e'))+
   scale_fill_manual(values=c('#000000', '#1c5580', '#0fc08e'))+
  scale_y_continuous(trans = 'log10', limits=c(5e2,5e5), 
                     labels = trans_format('log10', math_format(10^.x)))+
  theme(legend.position = "none")+
  labs(y="Total CFUs per fish", x=NULL, fill=NULL, color=NULL, shape=NULL)+
   theme(legend.position = "right",legend.direction = "vertical", 
         panel.background = element_rect(size=1),
         legend.text = element_markdown(size=22))+
   guides(color=guide_legend(override.aes = list(size=4)))
mixtotalCFUplot
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/mixtotalcfus-1.png)<!-- -->


## Abundance of each strain per fish with stats


```r
library(ggh4x)

statsbyday <- compare_means(data=datacfusmix %>% 
   pivot_longer(Bc1perFish:Bc10perFish) %>% 
   mutate(name=factor(name,
                      levels=c("Bc2perFish","Bc10perFish","Bc3perFish","Bc4perFish","Bc1perFish"),
                      labels=c("S1. *Pseudomonas mossellii*",
                               "S3. *Pseudomonas nitroreducens*",
                               "S5. *Stenotrophomonas maltophila*",
                               "S6. *Aeromonas caviae*",
                               "S7. *Aeromonas veronii*")),
          percent=value/CFUs_perFish), 
   percent~LoperamideTreatment, 
   group.by = c("Treatment","Timepoint", "name"))
statsdatastrains <- statsbyday %>% filter(p.format<0.05 & group1=="DMSO") %>% 
   mutate(CFUs_perFish=1, LoperamideTreatment="Loperamide")

abundplot <- datacfusmix %>% 
   pivot_longer(Bc1perFish:Bc10perFish) %>% 
   mutate(name=factor(name,
                      levels=c("Bc2perFish","Bc10perFish","Bc3perFish","Bc4perFish","Bc1perFish"),
                      labels=c("S1. *Pseudomonas mossellii*",
                               "S3. *Pseudomonas nitroreducens*",
                               "S5. *Stenotrophomonas maltophila*",
                               "S6. *Aeromonas caviae*",
                               "S7. *Aeromonas veronii*")),
          percent=value/CFUs_perFish) %>% 
   ggplot(aes(x = LoperamideTreatment, y=percent, 
             fill=name, shape=LoperamideTreatment))+
   facet_wrap(.~Timepoint, nrow=1, scales="free_x")+
   geom_point(size=2, position=position_jitterdodge(jitter.width=0.2)) +
   geom_boxplot(alpha=0.7, show.legend = FALSE)+
   geom_text(data=statsdatastrains, aes(label=p.signif, fill=name), size=12,
             color=c("#847CA3", "#F4A65E", "#80792B"),
             y=c(1,1,1), face="bold",
             nudge_x=c(0.3,0,-0.15),
            show.legend=FALSE)+
   scale_fill_manual(values = rev(nationalparkcolors::park_palette("Saguaro", 5)))+
   scale_color_manual(values = rev(nationalparkcolors::park_palette("Saguaro", 5)))+
   scale_shape_manual(values=c(21,24,22))+
   scale_y_continuous(labels=scales::label_percent(sale=1), limits=c(0,1.05))+
   theme(strip.text=element_text(size=15), 
         strip.background = element_rect(fill="grey85", color="white"),
         legend.position="none")+
   labs(y="% CFUs per strain per fish", x=NULL, fill="Treatment", color="Treatment", shape="Treatment")

abundplot
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/strainpercent-1.png)<!-- -->

## Percent abundance of each fish


```r
library(ggh4x)

percentplot <- datacfusmix %>% 
   pivot_longer(Bc1perFish:Bc10perFish) %>% 
   mutate(name=factor(name,
                      levels=c("Bc2perFish","Bc10perFish","Bc3perFish","Bc4perFish","Bc1perFish"),
                      labels=c("S1. *Pseudomonas mossellii*",
                               "S3. *Pseudomonas nitroreducens*",
                               "S5. *Stenotrophomonas maltophila*",
                               "S6. *Aeromonas caviae*",
                               "S7. *Aeromonas veronii*"))) %>% 
   ggplot(aes(x=as.factor(FishNum), y=value, fill=name))+
   geom_col(position="fill")+
   facet_nested(.~Timepoint+LoperamideTreatment, scales="free")+
   scale_y_continuous(labels=scales::label_percent())+
   scale_fill_manual(values = rev(nationalparkcolors::park_palette("Saguaro", 5)))+
   labs(x=NULL, y="% CFUs per fish", fill=NULL)+
   theme(strip.text=element_text(size=15), axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         strip.background = element_rect(fill="grey85", color="white"),
         legend.position = "bottom", legend.direction = "vertical",
         panel.background = element_blank(),
         legend.text = element_markdown(size=22))
percentplot
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/straintotal-1.png)<!-- -->

## Summary figure


```r
library(cowplot)
percentlegend <- get_legend(percentplot)
cfulegend <- get_legend(mixtotalCFUplot)

legends <- plot_grid(cfulegend, percentlegend, ncol=1, axis = "l")

((mixtotalCFUplot+legends + plot_layout(widths=c(2,4)))  /
      percentplot) /
   abundplot+plot_annotation(tag_levels = list(c('A','','B','C')))  &
   theme(legend.position="none", plot.tag = element_text(face = "bold", size=20))
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/strainpercentsummary-1.png)<!-- -->

```r
ggsave("Figure5_MixReconvPlots.png", width=12, height=11)
ggsave("Figure5_MixReconvPlots.pdf", width=12, height=11)
ggsave("Figure5_MixReconvPlots.tiff", width=12, height=11)

plot_grid(percentlegend) / abundplot  + plot_layout(heights=c(1,2))
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/strainpercentsummary-2.png)<!-- -->

```r
ggsave("Figure5C.png", width=12, height=5.5)
```


# Alpha diversity

## Simpsons


```r
matrixcfusmix<- datacfusmix %>% ungroup() %>% 
   unite("FishID", LoperamideTreatment,Timepoint,FishNum) %>% 
   select(FishID, Bc1perFish:Bc10perFish) %>% 
   column_to_rownames(var="FishID") %>% as.matrix()

datacfusmix$Simpsons<-diversity(matrixcfusmix, index="simpson")

statsSimpsonsMix <-  compare_means(data=datacfusmix, 
                            Simpsons~LoperamideTreatment, 
                            group.by = c("Timepoint")) %>% 
   filter(p.format<0.05 & group1=="DMSO")

mixsimpson <- datacfusmix %>% 
   ggplot(aes(x=Timepoint, y=Simpsons, 
             fill=LoperamideTreatment, color=LoperamideTreatment, shape=LoperamideTreatment))+
   geom_boxplot(color="black", alpha=0.8, show.legend = FALSE)+
   geom_point(size=2, position=position_jitterdodge(jitter.width=0.3)) +
   geom_text(data=statsSimpsonsMix, aes(label="*", x=Timepoint, 
                                        y=0.9, fill=NA, shape=NA), 
             size=9, color="#0fc08e", show.legend=FALSE, nudge_x=.25)+
   scale_color_manual(values=c('#000000', '#1c5580', '#0fc08e'))+
   scale_fill_manual(values=c('#000000', '#1c5580', '#0fc08e'))+
   scale_y_continuous(limits=c(0,1.0), expand=c(0,0))+
   theme(legend.position = "none")+
   labs(x=NULL, y="Simpson's Index",fill=NULL, shape=NULL, color=NULL,
        title="Mix5 gnotobiotic fish")+
   theme(legend.position = "bottom", plot.title = element_text(size=20),
         panel.background = element_rect(color="grey80", size=0.8))+
   guides(color=guide_legend(override.aes = list(size=4)))
mixsimpson
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/alphasimp-1.png)<!-- -->


## Richness


```r
# calculate richness
Smix5 <- specnumber(matrixcfusmix) 
# add into metadata
datacfusmix$richness <- Smix5

# calculate stats
compare_means(data=datacfusmix, 
                            richness~LoperamideTreatment, 
                            group.by = c("Timepoint")) %>% 
   filter(p.format<0.05 & group1=="DMSO")
```

```
## # A tibble: 0 × 9
## # … with 9 variables: Timepoint <chr>, .y. <chr>, group1 <chr>, group2 <chr>,
## #   p <dbl>, p.adj <dbl>, p.format <chr>, p.signif <chr>, method <chr>
```

```r
mixrichness <- datacfusmix %>% 
   ggplot(aes(x=Timepoint, y=richness, 
             fill=LoperamideTreatment, color=LoperamideTreatment, shape=LoperamideTreatment))+
   geom_boxplot(color="black", alpha=0.8, show.legend = FALSE)+
   geom_point(size=2, position=position_jitterdodge(jitter.width=0.3)) +
   scale_color_manual(values=c('#000000', '#1c5580', '#0fc08e'))+
   scale_fill_manual(values=c('#000000', '#1c5580', '#0fc08e'))+
   scale_y_continuous(limits=c(0,6), expand=c(0,0))+
   theme(legend.position = "none")+
   labs(x=NULL, y="Observed richness",fill=NULL, shape=NULL, color=NULL,
        title="Mix5 gnotobiotic fish")+
   theme(legend.position = "bottom", plot.title = element_text(size=20),
         panel.background = element_rect(color="grey80", size=0.8))+
   guides(color=guide_legend(override.aes = list(size=4)))
mixrichness
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/alpharich-1.png)<!-- -->

## Evenness


```r
# calculate evenness
Jmix5 <- microbiome::evenness(t(matrixcfusmix))

# add into metadata dataframe
datacfusmix <- Jmix5 %>% rownames_to_column("FishID") %>% 
   left_join(datacfusmix)

# stats
compare_means(data=datacfusmix, 
                            pielou~LoperamideTreatment, 
                            group.by = c("Timepoint")) %>% 
   filter(p.format<0.05 & group1=="DMSO")
```

```
## # A tibble: 0 × 9
## # … with 9 variables: Timepoint <chr>, .y. <chr>, group1 <chr>, group2 <chr>,
## #   p <dbl>, p.adj <dbl>, p.format <chr>, p.signif <chr>, method <chr>
```

```r
mixeven <- datacfusmix %>% 
   ggplot(aes(x=Timepoint, y=pielou, 
             fill=LoperamideTreatment, color=LoperamideTreatment, shape=LoperamideTreatment))+
   geom_boxplot(color="black", alpha=0.8, show.legend = FALSE)+
   geom_point(size=2, position=position_jitterdodge(jitter.width=0.3)) +
   scale_color_manual(values=c('#000000', '#1c5580', '#0fc08e'))+
   scale_fill_manual(values=c('#000000', '#1c5580', '#0fc08e'))+
   scale_y_continuous(limits=c(0,1.0), expand=c(0,0))+
   theme(legend.position = "none")+
   labs(x=NULL, y="Pielou's evenness",fill=NULL, shape=NULL, color=NULL,
        title="Mix5 gnotobiotic fish")+
   theme(legend.position = "bottom", plot.title = element_text(size=20),
         panel.background = element_rect(color="grey80", size=0.8))+
   guides(color=guide_legend(override.aes = list(size=4)))
mixeven
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/alphaeven-1.png)<!-- -->

## Summary figure

Get conventional panel from 16S script


```r
convalpha + mixalpha + labs(y=NULL) + 
   plot_layout(guides="collect")+plot_annotation(tag_levels = "A") & 
   theme(legend.position="bottom",plot.tag = element_text(face = "bold", size=20))
```


Get conventional panel from 16S script


```r
(convsimpson+convrichness+conveven) / 
   (mixsimpson+mixrichness+mixeven)  + labs(y=NULL) + 
   plot_layout(guides="collect")+plot_annotation(tag_levels = "A") & 
   theme(legend.position="bottom",plot.tag = element_text(face = "bold", size=20))

ggsave("Figure6_AlphaDiversity_Conv_Mix5.png", bg="transparent", width = 11, height = 8)
ggsave("Figure6_AlphaDiversity_Conv_Mix5.pdf", bg="transparent", width = 8, height = 6)
ggsave("Figure6_AlphaDiversity_Conv_Mix5.tiff", bg="transparent", width = 8, height = 5)
```

------------------------------------------------------------------------

# Compare with Mono as a mix


## Import Mono data for these strains


```r
datamonocfus <-
   readxl::read_xlsx("../Mono/MonoCFUs_LoperamideZebrafish.xlsx", sheet="Trial49") %>%
   drop_na(DF) %>%
   mutate(LoperamideTreatment=factor(LoperamideTreatment, 
                                    levels=c("None", "DMSO", "Loperamide 10 mg/L"),
                                    labels=c("Control","DMSO", "Loperamide")),
         Treatment = factor(Treatment,
                            levels=c("Bc1","Bc2","Bc3","Bc4","Bc10")))
```


## Mono plots


```r
percentmono <- datamonocfus %>% left_join(straininfo, by=c("Treatment" = "Strain")) %>% 
   group_by(Timepoint, LoperamideTreatment, CodeName) %>% 
   summarise(meanCFUs=mean(CFUs_perFish), sdCFUs=sd(CFUs_perFish)) %>% 
   ggplot(aes(x=LoperamideTreatment, y=meanCFUs, fill=CodeName))+
   geom_col(position="fill")+
   facet_grid(.~Timepoint, space="free", scales="free")+
   theme(legend.text = element_markdown(), legend.position = "right")+
   scale_y_continuous(labels=scales::label_percent())+ #limits=c(0,1e6))+ #
   scale_fill_manual(values = rev(nationalparkcolors::park_palette("Saguaro", 5)))+
   labs(x=NULL, y="mean CFUs per strain \n(% of total CFUs per condition)", fill="Strain",
        title="Mono means as a mix")

abundmono <- datamonocfus  %>% left_join(straininfo, by=c("Treatment" = "Strain")) %>% 
   group_by(Timepoint, LoperamideTreatment, Treatment, CodeName) %>% 
   summarise(meanCFUs=mean(CFUs_perFish), sdCFUs=sd(CFUs_perFish)) %>% 
   mutate(Treatment=factor(Treatment, levels=c("Bc1","Bc2","Bc3","Bc4","Bc10"))) %>% 
   ggplot(aes(x=LoperamideTreatment, y=meanCFUs, fill=CodeName))+
   geom_col(position=position_dodge(width=0.8), width=0.8)+
   geom_errorbar(aes(ymin=pmax(meanCFUs-sdCFUs,0), ymax=meanCFUs+sdCFUs),
                 position=position_dodge(width=0.8),width=0.2)+
   facet_grid(.~Timepoint, space="free", scales="free")+
    theme(legend.text = element_markdown(), legend.position = "right")+
   scale_y_continuous(trans = 'log10', expand=c(0,0),limits=c(NA,2e6),
                      labels = trans_format('log10', math_format(10^.x)))+
   scale_fill_manual(values = rev(nationalparkcolors::park_palette("Saguaro", 5)))+
   labs(x=NULL, y="CFUs per strain per fish", fill="Strain",
        title="Mono means as a mix")

percentmono/abundmono+plot_layout(guides="collect") &
  theme(legend.position='bottom')
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/monomixcfus-1.png)<!-- -->

## Mix plots


```r
percentmix<- datacfusmix  %>% 
   pivot_longer(Bc1perFish:Bc10perFish) %>% 
   mutate(name=factor(name,
                      levels=c("Bc1perFish","Bc2perFish","Bc3perFish","Bc4perFish","Bc10perFish"),
                      labels=c("Bc1","Bc2","Bc3","Bc4","Bc10"))) %>% 
   left_join(straininfo, by=c("name" = "Strain")) %>% 
   group_by(LoperamideTreatment, Timepoint, CodeName) %>% 
   summarise(meanCFUs=mean(value)) %>% 
   ggplot(aes(x=LoperamideTreatment, y=meanCFUs, fill=CodeName))+
   geom_col(position="fill")+
   facet_grid(.~Timepoint, space="free", scales="free")+
   scale_y_continuous(labels=scales::label_percent())+ #limits=c(0,1e6))+ #
   scale_fill_manual(values = rev(nationalparkcolors::park_palette("Saguaro", 5)))+
   theme(legend.text = element_markdown(), legend.position = "right")+
   labs(x=NULL, y="CFUs per strain \n(% of total CFUs per fish)", fill="Strain",
        title="Mix5 means per condition/timepoint")

abundmix <- datacfusmix %>% 
   pivot_longer(Bc1perFish:Bc10perFish) %>% 
   mutate(name=factor(name,
                      levels=c("Bc1perFish","Bc2perFish","Bc3perFish","Bc4perFish","Bc10perFish"),
                      labels=c("Bc1","Bc2","Bc3","Bc4","Bc10"))) %>% 
   left_join(straininfo, by=c("name" = "Strain")) %>% 
   group_by(CodeName, LoperamideTreatment, Timepoint) %>% 
   summarise(meanCFUs=mean(value), sdCFUs=sd(value)) %>% 
   ggplot(aes(x=LoperamideTreatment, y=meanCFUs, fill=CodeName))+
   geom_col(position=position_dodge(width=0.8), width=0.8)+
   geom_errorbar(aes(ymin=pmax(meanCFUs-sdCFUs,0), ymax=meanCFUs+sdCFUs),
                 position=position_dodge(width=0.8),width=0.2)+
   facet_grid(.~Timepoint, space="free", scales="free")+
   scale_y_continuous(trans = 'log10', expand=c(0,0),limits=c(NA,2e6),
                      labels = trans_format('log10', math_format(10^.x)))+
   theme(legend.text = element_markdown(), legend.position = "right")+
   scale_fill_manual(values = rev(nationalparkcolors::park_palette("Saguaro", 5)))+
   labs(x=NULL, y="mean CFUs per strain per fish", fill="Strain",
        title="Mix5 means per condition/timepoint")

percentmix/abundmix+plot_layout(guides="collect") &
  theme(legend.position='bottom')
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/straintotalmean-1.png)<!-- -->


## Summary figure


```r
percentmix+percentmono+
   (abundmix+labs(title=NULL))+(abundmono+labs(title=NULL)) +
   plot_layout(guides="collect", nrow=2) + plot_annotation(tag_levels = "A") &
   theme(legend.position='bottom',plot.title = element_text(size=34),
         axis.text.x = element_markdown(size = 19),
         legend.text = element_markdown(size = 30),
         axis.text.y = element_markdown(size = 22),
         axis.title = element_markdown(size = 30),
         legend.title = element_markdown(size = 30),
         plot.tag = element_text(face = "bold", size=30),
         strip.text = element_text(size=30)) & 
   guides(fill = guide_legend(ncol = 3))
```

![](ZebrafishMixCFUs_loperamide_files/figure-html/mixmonosummary-1.png)<!-- -->

```r
ggsave("FigureS9_CompareMonoMixReconv.png", width=26, height=16, dpi=400)
ggsave("FigureS9_CompareMonoMixReconv.pdf", width=26, height=16)
ggsave("FigureS9_CompareMonoMixReconv.tiff", width=26, height=16)
```
Sum of the mono-reconv does not equal the mix-reconv.  
Also, increased colonization for each strain in mono-reconv than when part of a mix.


