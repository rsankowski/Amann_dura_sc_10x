library(ggpubr)
library(tidyquant)
library(tidyverse)

df <- read.csv("/Users/romansankowski/Downloads/Quantification_CoQ2_S96N_vs_WT - Tabellenblatt1.csv")

df <- df %>% 
  mutate(density=case_when(Staining %in% c("NeuN", "Iba1","APP") ~ Quant/Area, 
                               T ~ Quant)) %>% 
  mutate(Staining = factor(Staining, levels=c("HE","LFB-PAS","NeuN", "Iba1","APP")))

df %>% 
  filter(Anatomical_Region=="mesenc") %>% 
  ggplot(aes(x=Condition, y=density, fill=Condition))  +
  stat_summary(geom="bar", fun=mean, color="black", ldw =.1) +
  stat_summary(fun.data = mean_se) +
  geom_point(size=5, pch = 21) +
  scale_fill_tq() +
  facet_wrap(~Staining, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test") +
  theme_linedraw() +
  expand_limits(y=0) +
  theme(legend.position = "none")

ggsave(file.path("plots","barplot_mesencephalon.pdf"))

df %>% 
  filter(Anatomical_Region=="med") %>% 
  ggplot(aes(x=Condition, y=density, fill=Condition))  +
  stat_summary(geom="bar", fun=mean, color="black", lwd=.1) +
  stat_summary(fun.data = mean_se) +
  geom_point(size=5, pch = 21) +
  scale_fill_tq() +
  facet_wrap(~Staining, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test") +
  theme_linedraw() +
  expand_limits(y=0) +
  theme(legend.position = "none")

ggsave(file.path("plots","barplot_medulla.pdf"))

