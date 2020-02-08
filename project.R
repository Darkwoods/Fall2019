setwd("E:/Маша/Bioinf/Genoms")
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)


data <- as_data_frame(fread('overlap_1.csv'))
data <- as_data_frame(fread('snp.csv'))
data <- as_data_frame(fread('stat.csv'))
data %>% mutate(n=row_number()) %>% mutate(color = case_when(V1 %% 2==0 ~ 0,
                                                             V1 %% 2!=0 ~ 1)) %>%
  ggplot(aes(x=n, y=V6, fill=as.factor(color), label = V1))+
  geom_bar(stat='identity', show.legend = F)+
  theme_bw()+
  labs(y = "Number of variants in each window")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.margin = unit(0, "lines"))+
  scale_x_continuous(expand = c(0, 0))

data

set.seed(16)
rand_data <- sample_n(data, 20)


variants <- as_data_frame(fread('variants.csv'))
str(variants)

var <- filter(variants, V1 == rand_data$V1 & V4 == rand_data$V2)

so <- var %>% group_by(V1) %>% unite(group,c(V1,V4),sep='_',remove = F) %>% ungroup()

var %>% ggplot(aes(x=V4, y=V3, fill=V4))+
  geom_violin()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45))+
  scale_y_log10()

so %>% group_by(V1,group) %>% slice(1)
so %>% ggplot(aes(x=group, y=V3, fill=group))+
  geom_violin()+
  geom_boxplot(width=0.15, fill="white")+
  theme_bw()+
  scale_y_log10(name = "Frequency")+
  scale_x_discrete(name = "Interval")+
  scale_fill_discrete(name = "Interval")+
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
  

data %>% mutate(n=row_number()) %>% mutate(color = case_when(V1 %% 2==0 ~ 0,
                                                             V1 %% 2!=0 ~ 1)) %>%
  ggplot(aes(x=n, y = V4))+
  geom_line()+
  theme_bw()+
  labs(y = "Mean")


data %>% mutate(n=row_number()) %>%
  ggplot(aes(x=n, y = V3))+
  geom_line()+
  theme_bw()+
  labs(y = "Freq")



