#clean memory
# rm(list=ls(all=TRUE))

library(tidyverse)
library(gridExtra)
library(ggpmisc)
library(ggplot2)
my.formula <- y ~ x

#### Ecologycal Selection - (Asexual Growth and Top/Bottom Sel.)

competition_data<-read.table("ecol_sel_measurements.txt", T) %>% 
  mutate(counts=red+wt, 
         red=red+1, 
         wt=wt+1) %>% 
  filter(counts>500) %>% 
  select(-counts)


inference_m.evol<-competition_data %>% filter(Time=="AG") %>% mutate(Sat_ref_total=30000000, 
                            N_ref_BG=((((Sat_ref_total/300)*100)/350)*5)*30/505 ,
                            N_evo_BG=(freq.before*N_ref_BG)/(1-freq.before) , 
                            tAG=6 ,
                            N_ref_AG= N_ref_BG*exp(tAG) ,
                            r.evol.AG= wt/(wt+red) ,
                            N_evol_AG=(r.evol.AG*N_ref_AG)/(1-r.evol.AG), 
                            m.evol = log(N_evol_AG/N_evo_BG)/tAG) 



# General distribution of all replicates 
inference_m.evol %>% 
  mutate(group_ID=paste(Group, evolved, sep="__")) %>% 
  filter(Ade=="P") %>% 
  mutate(population=ifelse(Group %in% c("P3P", "P3R"), "P",as.character(Pop))) %>% 
  ggplot(aes(as.factor(population), as.numeric(m.evol), fill=Group)) +
  geom_boxplot() + 
  coord_flip() +
  geom_point(alpha=0.4, size=0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_jitterdodge()) +
  ggtitle(label = "Growth rate - pink pop." ) +
  xlab(label = "Population") +
  ylab(label = "Growth rate (m)") + 
  facet_grid(group_ID ~ . , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        legend.position="none",
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("plots/01_growth_m_perPop_pink.png"), width = 6, height = 14, dpi = 400)



# distribution of mean per population (mean between technical replicates), including parental strains
# colours_top_bottom<-c("red3","#D95F02","cornflowerblue")
colours_top_bottom<-c("red3","orange", "cornflowerblue") ; inference_m.evol %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
  #filter(Group_sample!="Parental") %>%
  #filter(Ade=="R") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  #filter(Group_sample=="Parental") %>% 
  # mutate(group_ID=paste(Group, evolved, sep="_")) %>% 
  # group_by(group_ID, Pop, Ade, Group) %>% 
  # dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
  # ungroup() %>% 
  # mutate(group_ID=ifelse(Group %in% c("P3P", "P3R"), as.character(Group),as.character(group_ID))) %>%
  # ggplot(aes(as.factor(group_ID), median_m.evol, fill=Group)) +
  ggplot(aes(Group, median_m.evol, fill=Group_sample)) + 
  geom_boxplot(alpha=0.6) +
  # geom_point(alpha=0.4, size=1) +
  geom_point(aes(fill=Group_sample), size = 1.5, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=colours_top_bottom) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  xlab(label = "Treatment") +
  ylab(label = "Relative Growth rate (m)") + 
  # facet_grid(Ade ~ . , scale="free", space="free") +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.x = element_blank() , 
        strip.background = element_blank(), 
        #legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
   ggsave(paste0("new_analyses_plots/02_0_growth_m_median_perPop_all_sim.png"), width = 8, height = 5, dpi = 400)
  # ggsave(paste0("plots/02_growth_m_median_perPop_all_sim.png"), width = 3, height = 3, dpi = 400)
  # ggsave(paste0("plots/02_growth_m_median_perPop_red_sim.png"), width = 3, height = 3, dpi = 400)


   

# Change in growth rate relative to parental strain:
colours_top_bottom<-c("red3","cornflowerblue") ; inference_m.evol %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
  #filter(Group_sample!="Parental") %>%
  #filter(Ade=="R") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  merge(inference_m.evol %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
          #filter(Group_sample!="Parental") %>%
          #filter(Ade=="R") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                              levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade) %>% 
          dplyr::summarise(parental_median_m=median(median_m.evol), 
                    parental_mean_m=mean(mean_m.evol)), 
        by="Ade", all.x = T) %>% 
  # mutate(dif_median_m.evol=median_m.evol-parental_median_m, 
  #        dif_mean_m.evol=mean_m.evol-parental_mean_m) %>% 
  # mutate(dif_median_m.evol=exp(log(median_m.evol)-log(parental_median_m)),
  #        dif_mean_m.evol=exp(log(mean_m.evol)-log(parental_mean_m))) %>%
  mutate(dif_median_m.evol=exp(median_m.evol-parental_median_m),
         dif_mean_m.evol=exp(mean_m.evol-parental_mean_m)) %>%
  # ggplot(aes(dif_median_m.evol, dif_mean_m.evol)) +
  # geom_point()
  # filter(Ade=="R") %>%
  ggplot(aes(Group, dif_median_m.evol, fill=Group_sample)) + 
  geom_boxplot(alpha=0.6) +
  # geom_point(alpha=0.4, size=1) +
  geom_point(aes(fill=Group_sample), size = 1.5, shape = 21, position = position_jitterdodge()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  ylab(label = "Relative Growth rate (m)") + 
  # facet_grid(Ade ~ . , scale="free", space="free") +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        strip.text.x = element_blank() , 
        strip.background = element_blank(), 
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("new_analyses_plots/02_1_relative_growth_m_median_pink_perPop_all_sim.png"), width = 3, height = 3, dpi = 400)
  # ggsave(paste0("new_analyses_plots/02_1_relative_growth_m_median_red_perPop_all_sim.png"), width = 3, height = 3, dpi = 400)
  ggsave(paste0("new_analyses_plots/02_1_relative_growth_m_median_perPop_all_sim.png"), width = 8, height = 5, dpi = 400)
  ggsave(paste0("new_analyses_plots/02_1_growth_m_median_perPop_all_sim.png"), width = 8, height = 5, dpi = 400)
  
  


# Normalised values per treatment using variance per treatment:
inference_m.evol %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
  #filter(Group_sample!="Parental") %>%
  #filter(Ade=="R") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  merge(inference_m.evol %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
          #filter(Group_sample!="Parental") %>%
          #filter(Ade=="R") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade) %>% 
          dplyr::summarise(parental_median_m=median(median_m.evol), 
                    parental_mean_m=mean(mean_m.evol)), 
        by="Ade", all.x = T) %>% 
  # mutate(dif_median_m.evol=median_m.evol-parental_median_m, 
  #        dif_mean_m.evol=mean_m.evol-parental_mean_m) %>% 
  mutate(dif_median_m.evol=exp(log(median_m.evol)-log(parental_median_m)), 
         dif_mean_m.evol=exp(log(mean_m.evol)-log(parental_mean_m))) %>% 
  # ggplot(aes(dif_median_m.evol, dif_mean_m.evol)) +
  # geom_point()
  merge(inference_m.evol %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
          #filter(Group_sample!="Parental") %>%
          #filter(Ade=="R") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                              levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group!="Parental") %>% 
          merge(inference_m.evol %>% 
                  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental", 
                                             as.character(evolved))) %>% 
                  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
                  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
                  #filter(Group_sample!="Parental") %>%
                  #filter(Ade=="R") %>%
                  ungroup() %>% 
                  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                                      levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
                  filter(Group=="Parental") %>% 
                  group_by(Ade) %>% 
                  summarise(parental_median_m=median(median_m.evol), 
                            parental_mean_m=mean(mean_m.evol)), 
                by="Ade", all.x = T) %>% 
          # mutate(dif_median_m.evol=median_m.evol-parental_median_m, 
          #        dif_mean_m.evol=mean_m.evol-parental_mean_m) %>% 
          mutate(dif_median_m.evol=exp(log(median_m.evol)-log(parental_median_m)), 
                 dif_mean_m.evol=exp(log(mean_m.evol)-log(parental_mean_m))) %>% 
          group_by(Ade,Group, Group_sample) %>% 
          dplyr::summarize(sd_diff_median_m.evol=sd(dif_median_m.evol), 
                           sd_diff_mean_m.evol=sd(dif_mean_m.evol)), 
        by=c("Ade", "Group", "Group_sample"), all.x = T) %>% 
  group_by(Ade, Group, Group_sample) %>% 
  dplyr::summarise(median_dif_median_m.evol=median(dif_median_m.evol), 
            median_sd_diff_median_m.evol=median(sd_diff_median_m.evol)) %>% 
  mutate(N_median_dif_median_m.evol=median_dif_median_m.evol/median_sd_diff_median_m.evol) %>% 
  # filter(Ade=="P")
  ggplot(aes(Group, N_median_dif_median_m.evol, fill=Group_sample)) + 
  # geom_point(alpha=0.4, size=1) +
  geom_point(aes(fill=Group_sample), size = 2, shape = 21, alpha=0.6) +
  scale_fill_manual(values=colours_top_bottom) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  xlab(label = "Treatment") +
  ylab(label = "N. Change Growth rate (m)") + 
  # facet_grid(Ade ~ . , scale="free", space="free") +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.x = element_blank() , 
        # strip.background = element_blank(), 
        # legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("new_analyses_plots/02_3_N_growth_m_median_perPop_all_sim.png"), width = 8, height = 5, dpi = 400)
  

# Normalised values per treatment using total variance:
inference_m.evol %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
  #filter(Group_sample!="Parental") %>%
  #filter(Ade=="R") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  # this will add the parental parameters:
  merge(inference_m.evol %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
          #filter(Group_sample!="Parental") %>%
          #filter(Ade=="R") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade) %>% 
          summarise(parental_median_m=median(median_m.evol), 
                    parental_mean_m=mean(mean_m.evol)), 
        by="Ade", all.x = T) %>% 
  # calculate difference between evolved populations and parental strains:
  # mutate(dif_median_m.evol=median_m.evol-parental_median_m, 
  #        dif_mean_m.evol=mean_m.evol-parental_mean_m) %>% 
  mutate(dif_median_m.evol=exp(log(median_m.evol)-log(parental_median_m)), 
         dif_mean_m.evol=exp(log(mean_m.evol)-log(parental_mean_m))) %>% 
  # ggplot(aes(dif_median_m.evol, dif_mean_m.evol)) +
  # geom_point()
  # this will add the variance component. In this one, I used the total variance rather than the variance per treatment.
  merge(inference_m.evol %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
          #filter(Group_sample!="Parental") %>%
          #filter(Ade=="R") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                              levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group!="Parental") %>% 
          merge(inference_m.evol %>% 
                  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental", 
                                             as.character(evolved))) %>% 
                  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
                  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
                  #filter(Group_sample!="Parental") %>%
                  #filter(Ade=="R") %>%
                  ungroup() %>% 
                  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                                      levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
                  filter(Group=="Parental") %>% 
                  group_by(Ade) %>% 
                  summarise(parental_median_m=median(median_m.evol), 
                            parental_mean_m=mean(mean_m.evol)), 
                by="Ade", all.x = T) %>% 
          # mutate(dif_median_m.evol=median_m.evol-parental_median_m, 
          #        dif_mean_m.evol=mean_m.evol-parental_mean_m) %>% 
          mutate(dif_median_m.evol=exp(log(median_m.evol)-log(parental_median_m)), 
                 dif_mean_m.evol=exp(log(mean_m.evol)-log(parental_mean_m))) %>% 
          group_by(Ade) %>% 
          dplyr::summarize(sd_diff_median_m.evol=sd(dif_median_m.evol), 
                           sd_diff_mean_m.evol=sd(dif_mean_m.evol)), 
        by=c("Ade"), all.x = T) %>% 
  group_by(Ade, Group, Group_sample) %>% 
  summarise(median_dif_median_m.evol=median(dif_median_m.evol), 
            median_sd_diff_median_m.evol=median(sd_diff_median_m.evol)) %>% 
  mutate(N_median_dif_median_m.evol=median_dif_median_m.evol/median_sd_diff_median_m.evol) %>% 
  # filter(Ade=="P")
  ggplot(aes(Group, N_median_dif_median_m.evol, fill=Group_sample)) + 
  # geom_point(alpha=0.4, size=1) +
  geom_point(aes(fill=Group_sample), size = 2, shape = 21, alpha=0.6) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  ylab(label = "N. Median Change Growth rate (m)") + 
  # facet_grid(Ade ~ . , scale="free", space="free") +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.x = element_blank() , 
        # strip.background = element_blank(), 
        # legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("new_analyses_plots/02_4_totalVarN_growth_m_median_perPop_all_sim.png"), width = 8, height = 5, dpi = 400)


# Change in growth rate without variance normalisation:
inference_m.evol %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
  #filter(Group_sample!="Parental") %>%
  #filter(Ade=="R") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  merge(inference_m.evol %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
          #filter(Group_sample!="Parental") %>%
          #filter(Ade=="R") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade) %>% 
          summarise(parental_median_m=median(median_m.evol), 
                    parental_mean_m=mean(mean_m.evol)), 
        by="Ade", all.x = T) %>% 
  # mutate(dif_median_m.evol=median_m.evol-parental_median_m, 
  #        dif_mean_m.evol=mean_m.evol-parental_mean_m) %>% 
  mutate(dif_median_m.evol=exp(log(median_m.evol)-log(parental_median_m)), 
         dif_mean_m.evol=exp(log(mean_m.evol)-log(parental_mean_m))) %>% 
  group_by(Ade, Group, Group_sample) %>% 
  summarise(median_dif_median_m.evol=median(dif_median_m.evol)) %>% 
  #mutate(N_mean_dif_median_m.evol=mean_dif_median_m.evol/meam_sd_diff_median_m.evol) %>% 
  # filter(Ade=="P")
  ggplot(aes(Group, median_dif_median_m.evol, fill=Group_sample)) + 
  # geom_point(alpha=0.4, size=1) +
  geom_point(aes(fill=Group_sample), size = 2, shape = 21, alpha=0.6) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  ylab(label = "Mean Change Growth rate (m)") + 
  # facet_grid(Ade ~ . , scale="free", space="free") +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.x = element_blank() , 
        # strip.background = element_blank(), 
        # legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("new_analyses_plots/02_5_change_growth_m_median_perPop_all_sim.png"), width = 8, height = 5, dpi = 400)


# median value per population: 

median_inference_m.evol_cor<-inference_m.evol %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
  #filter(Group_sample!="Parental") %>%
  #filter(Ade=="R") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  merge(inference_m.evol %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
          #filter(Group_sample!="Parental") %>%
          #filter(Ade=="R") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade) %>% 
          dplyr::summarise(parental_median_m=median(median_m.evol), 
                    parental_mean_m=mean(mean_m.evol)), 
        by="Ade", all.x = T) %>% 
  # mutate(dif_median_m.evol=median_m.evol-parental_median_m) %>% 
  # mutate(dif_median_m.evol=median_m.evol/parental_median_m) %>% 
  mutate(dif_median_m.evol=exp(log(median_m.evol)-log(parental_median_m)), 
         dif_mean_m.evol=exp(log(mean_m.evol)-log(parental_mean_m))) %>% 
  # mutate(test=ifelse(Group_sample=="BOTTOM", "B", ifelse(Group_sample=="TOP", "T", "Other"))) %>% 
  # mutate(test2=test==Fraction) %>% 
  # filter(test2==T) %>% 
  # head()
  select(Group, Fraction, Pop, Ade, dif_median_m.evol)


median_inference_m.evol<-inference_m.evol %>%
  mutate(group_ID=paste(Group, evolved, sep="__")) %>% 
  group_by(group_ID, Pop, Ade, Group) %>% 
  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
  select(-N_rep, -Group) %>% 
  separate(group_ID, c("Group", "evolved"), sep="__")
  

inference_m.evol %>% 
  filter(Ade=="R") %>% 
  group_by(Structure, Group, evolved, Pop, Ade) %>% 
  dplyr::summarize(N_rep=n(), mean_m.evol=mean(m.evol)) %>% 
  ggplot(aes(as.factor(evolved), as.numeric(mean_m.evol), fill=evolved)) +
  geom_boxplot() + 
  geom_point() +
  facet_grid(. ~ Group, scale="free") 

inference_m.evol %>% 
  filter(Ade=="R") %>% 
  group_by(Structure, Group, evolved, Pop, Ade) %>% 
  dplyr::summarize(N_rep=n(), mean_m.evol=mean(m.evol)) %>% 
  filter(Group=="Para") %>% 
  spread(evolved, mean_m.evol) %>%
  ggplot(aes(BOTTOM, TOP)) +
  geom_point() 
  facet_grid(. ~ Group, scale="free") 
  
ggplot(inference_m.evol %>% filter(Ade=="P"), aes(as.factor(Pop), as.numeric(m.evol), fill=evolved)) +
  geom_boxplot() + 
  facet_grid(. ~ Group, scale="free") #+
  ggsave("01_raw_growthRate_replicates_pink.png", width = 10, height = 10)
 
# inference_m.evol<-inference_m.evol %>% group_by(Structure, evolved, Pop, Ade, replicate) %>% 
#   select(Structure, Group, evolved, Pop, replicate, Ade, m.evol) 

  
  
  
  
  
  
# Ecological Selection: Top and bottom selection.
head(median_inference_m.evol)

ecological_sel<-
  competition_data %>% filter(Time!="AG") %>% 
  group_by(Structure, Group, evolved, Pop, Ade, replicate) %>% 
  merge(median_inference_m.evol, by=c("Group", "evolved", "Pop", "Ade"), all.x=TRUE) %>% 
  mutate(Sat_ref_total=30000000, 
         tAG=6 ,
         Expected_ref_BS=Sat_ref_total/6, 
         selection_co_0=0.01, 
         Expected_ref_After_Sel=Sat_ref_total*selection_co_0, 
         N_ref_BS=(Sat_ref_total*100/300)*50/350, 
         N_evo_BS=(freq.before*N_ref_BS)/(1-freq.before), 
         N_ref_AS=N_ref_BS*Expected_ref_After_Sel/Expected_ref_BS, 
         N_ref_S=N_ref_AS*exp(tAG) , 
         r.evol.S= wt/(wt+red) , 
         N_evo_S=(r.evol.S*N_ref_S)/(1-r.evol.S), 
         # N_evo_AS=N_evo_S/(exp(tAG*m.evol)),
         N_evo_AS=N_evo_S/(exp(tAG*median_m.evol)), 
         v_ref=N_ref_AS/N_ref_BS, 
         v_evo=N_evo_AS/N_evo_BS, 
         fq_BS_evo=N_evo_BS/(N_evo_BS+N_ref_BS),
         fq_AS_evo=N_evo_AS/(N_evo_AS+N_ref_AS), 
         fq_BS_ref=N_ref_BS/(N_evo_BS+N_ref_BS), 
         fq_AS_ref=N_ref_AS/(N_evo_AS+N_ref_AS), 
         s_estimate_fq=1-(((fq_BS_ref/fq_AS_ref)-fq_BS_ref)/fq_BS_evo), 
         relative_W_evo=v_evo/v_ref, 
         s_evo=1-relative_W_evo, 
         s_evol_difflog=log((fq_AS_evo/(1-fq_AS_evo)))-log((fq_BS_evo/(1-fq_BS_evo)))
         )

head(ecological_sel)

# general distribution of technical replicates within populations:
ecological_sel %>%
  mutate(group_ID=paste(Group, evolved, sep="__")) %>% 
  filter(Ade=="R") %>% 
  mutate(population=ifelse(Group %in% c("P3P", "P3R"), "P",as.character(Pop))) %>% 
  ggplot(aes(as.factor(population), as.numeric(log(relative_W_evo)), fill=Group)) +
  geom_boxplot() + 
  xlab(label = "Population") +
  ylab(label = "Log Relative W") +
  facet_grid(group_ID ~ Time, scale="free", space="free") +
  coord_flip() +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("plots/03_log_relative_W_ecologicalSel_perPop_red.png"), width = 8, height = 14, dpi = 400)


# Median values per population. Including parentals:
# colours_top_bottom<-c("red3","#D95F02","cornflowerblue")
library(scales)
asinh_trans <- function(){
    trans_new(name = 'asinh', transform = function(x) asinh(x), 
              inverse = function(x) sinh(x))
}
  
colours_top_bottom<-c("red3","orange", "cornflowerblue") ; ecological_sel %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
  filter(Ade=="R") %>%
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
            median_relative_W_evo = median(relative_W_evo)) %>%
  #filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  mutate(Time=factor(Time, levels = c("TS", "BS"))) %>% 
  # ggplot(aes(Group, log(mean_relative_W_evo), fill=Group_sample)) +
  ggplot(aes(Group, log(median_relative_W_evo), fill=Group_sample)) +
  # ggplot(aes(Time, relative_W_evo, fill=Group_sample)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  ylab(label = "Relative W") + 
  # ylab(label = "Log (Relative W)") + 
  # ylim(c(0,1000)) +
  geom_point(aes(fill=Group_sample), alpha=0.8, size = 1.5, shape = 21, position = position_jitterdodge()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
  # scale_y_continuous(breaks = seq((-3),8,1), 
  #                    labels = c(exp(seq((-3),8,1)))) +
  # scale_y_continuous(trans = 'asinh',breaks=c(-20,-10,-5,-2,-1,0,1,5,10,20,50,100))+
  scale_y_continuous(breaks = log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)),
                     labels = exp(log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)))) +
  facet_grid(Time ~ Ade, scale="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.y = element_blank() , 
        strip.background = element_blank(), 
        legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("new_analyses_plots/04_0_relative_W_ecologicalSel_all_sim.png", width = 3, height = 5.7, dpi = 400)
  # ggsave("new_analyses_plots/04_0_relative_W_ecologicalSel_pink_sim.png", width = 7, height = 5.5, dpi = 400)
  ggsave("new_analyses_plots/04_0_relative_W_ecologicalSel_red_sim.png", width = 7, height = 5.5, dpi = 400)



# Change in relative fitness per population
# colours_top_bottom<-c("red3","#D95F02","cornflowerblue")
library(scales)
asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}


colours_top_bottom<-c("red3","cornflowerblue") ; ecological_sel %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
  filter(Ade=="R") %>%
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
            median_relative_W_evo = median(relative_W_evo)) %>%
  #filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  merge(ecological_sel %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
          #filter(Ade=="R") %>%
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
          summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
                    median_relative_W_evo = median(relative_W_evo)) %>%
          #filter(Group_sample!="Parental") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                              levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade, Time) %>% 
          summarise(parental_median_relative_W_evo=median(median_relative_W_evo), 
                    parental_mean_relative_W_evo=mean(median_relative_W_evo)), 
        by=c("Ade","Time"), all.x = T) %>% 
  # mutate(dif_median_relative_W.evol=median_relative_W_evo-parental_median_relative_W_evo, 
  #        dif_mean_relative_W.evol=mean_relative_W_evo-parental_mean_relative_W_evo) %>% 
  # using the ratio intead?
  mutate(dif_median_relative_W.evol=median_relative_W_evo/parental_median_relative_W_evo, 
         dif_mean_relative_W.evol=mean_relative_W_evo/parental_mean_relative_W_evo) %>% 
  mutate(Time=factor(Time, levels = c("TS", "BS"))) %>% 
  ggplot(aes(Group, log(dif_median_relative_W.evol), fill=Group_sample)) +
  # ggplot(aes(Time, relative_W_evo, fill=Group_sample)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  ylab(label = "Relative W") + 
  # ylab(label = "Log (Relative W)") + 
  # ylim(c(0,1000)) +
  geom_point(aes(fill=Group_sample), alpha=0.8, size = 1.5, shape = 21, position = position_jitterdodge()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  # scale_y_continuous(breaks = seq((-3),8,1),
  #                    labels = c(exp(seq((-3),8,1)))) +
  scale_y_continuous(breaks = log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)),
                     labels = exp(log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)))) +
  # scale_y_continuous(trans = 'asinh',breaks=c(-20,-10,-5,-2,-1,0,1,5,10,20,50,100))+
  facet_grid(Time ~ Ade, scale="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        legend.position="none",
        strip.text = element_blank() ,
        strip.background = element_blank(),
        # legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("new_analyses_plots/04_1_change_relative_W_ecologicalSel_pink_sim_fig.png"), width = 3, height = 5.3, dpi = 400)
  # ggsave(paste0("new_analyses_plots/04_1_change_relative_W_ecologicalSel_red_sim_fig.png"), width = 3, height = 5.3, dpi = 400)
  # ggsave("new_analyses_plots/04_1_relative_W_ecologicalSel_all_sim.png", width = 3, height = 5.7, dpi = 400)
  # ggsave("new_analyses_plots/04_1_change_relative_W_ecologicalSel_pink_sim.png", width = 7, height = 5.5, dpi = 400)
  ggsave("new_analyses_plots/04_1_change_relative_W_ecologicalSel_red_sim.png", width = 7, height = 5.5, dpi = 400)


# Normalised values per treatment using variance using total variance:
colours_top_bottom<-c("red3","cornflowerblue") ; ecological_sel %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
  # filter(Ade=="R") %>%
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
            median_relative_W_evo = median(relative_W_evo)) %>%
  #filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  # add parental values and caluclate the change in median from parental to evolved population:
  merge(ecological_sel %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
          #filter(Ade=="R") %>%
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
          summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
                    median_relative_W_evo = median(relative_W_evo)) %>%
          #filter(Group_sample!="Parental") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                              levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade, Time) %>% 
          summarise(parental_median_relative_W_evo=median(median_relative_W_evo), 
                    parental_mean_relative_W_evo=mean(median_relative_W_evo)), 
        by=c("Ade","Time"), all.x = T) %>% 
  mutate(dif_median_relative_W.evol=median_relative_W_evo/parental_median_relative_W_evo, 
         dif_mean_relative_W.evol=mean_relative_W_evo/parental_mean_relative_W_evo) %>% 
  mutate(Time=factor(Time, levels = c("TS", "BS"))) %>% 
  # add variance parameter:
  merge(ecological_sel %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
          # filter(Ade=="R") %>%
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
          summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
                    median_relative_W_evo = median(relative_W_evo)) %>%
          #filter(Group_sample!="Parental") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                              levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group!="Parental") %>% 
          # add parental values and caluclate the change in median from parental to evolved population:
          merge(ecological_sel %>% 
                  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",
                                             as.character(evolved))) %>%
                  #filter(Ade=="R") %>%
                  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
                  summarise(mean_relative_W_evo = mean(relative_W_evo), 
                            sd_relative_W_evo = sd(relative_W_evo), 
                            median_relative_W_evo = median(relative_W_evo)) %>%
                  #filter(Group_sample!="Parental") %>%
                  ungroup() %>% 
                  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                                      levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
                  filter(Group=="Parental") %>% 
                  group_by(Ade, Time) %>% 
                  summarise(parental_median_relative_W_evo=median(median_relative_W_evo), 
                            parental_mean_relative_W_evo=mean(median_relative_W_evo)), 
                by=c("Ade","Time"), all.x = T) %>% 
          mutate(dif_median_relative_W.evol=median_relative_W_evo/parental_median_relative_W_evo, 
                 dif_mean_relative_W.evol=mean_relative_W_evo/parental_mean_relative_W_evo) %>% 
          mutate(Time=factor(Time, levels = c("TS", "BS"))) %>% 
          group_by(Ade, Time) %>% 
          dplyr::summarize(sd_dif_median_relative_W.evol=sd(dif_median_relative_W.evol), 
                           sd_dif_mean_relative_W.evol=sd(dif_mean_relative_W.evol)), 
        by=c("Ade", "Time"), all.x = T) %>% 
  group_by(Ade, Time, Group, Group_sample) %>% 
  summarise(median_dif_median_relative_W.evol=median(dif_median_relative_W.evol), 
            median_sd_dif_median_relative_W.evol=median(sd_dif_median_relative_W.evol)) %>% 
  mutate(N_median_dif_median_relative_W.evol=
           median_dif_median_relative_W.evol/median_sd_dif_median_relative_W.evol) %>% 
  ggplot(aes(Group, (N_median_dif_median_relative_W.evol), fill=Group_sample)) +
  # ggplot(aes(Time, relative_W_evo, fill=Group_sample)) + 
  geom_point(aes(fill=Group_sample), size = 2, shape = 21, alpha=0.6) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  ylab(label = "N. Ratio Change Relative W") + 
  # ylab(label = "Log (Relative W)") + 
  # ylim(c(0,1000)) +
  # geom_point(aes(fill=Group_sample), alpha=0.8, size = 1.5, shape = 21, position = position_jitterdodge()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  # scale_y_continuous(breaks = seq((-3),8,1), 
  #                    labels = c(exp(seq((-3),8,1)))) +
  # scale_y_continuous(breaks = log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)),
  #                    labels = exp(log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)))) +
  # scale_y_continuous(trans = 'asinh',breaks=c(-20,-10,-5,-2,-1,0,1,5,10,20,50,100))+
  facet_grid(Time ~ Ade, scale="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.y = element_blank() , 
        strip.background = element_blank(), 
        legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("new_analyses_plots/04_2_relative_W_ecologicalSel_all_sim.png", width = 3, height = 5.7, dpi = 400)
  # ggsave("new_analyses_plots/04_2_relative_W_ecologicalSel_pink_sim.png", width = 7, height = 5.5, dpi = 400)
  ggsave("new_analyses_plots/04_4_N_totalVar_change_relative_W_ecologicalSel_sim.png", width = 7, height = 5.5, dpi = 400)



colours_top_bottom<-c("red3","cornflowerblue") ; ecological_sel %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
  # filter(Ade=="R") %>%
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
            median_relative_W_evo = median(relative_W_evo)) %>%
  #filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  # add parental values and caluclate the change in median from parental to evolved population:
  merge(ecological_sel %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
          #filter(Ade=="R") %>%
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
          summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
                    median_relative_W_evo = median(relative_W_evo)) %>%
          #filter(Group_sample!="Parental") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                              levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade, Time) %>% 
          summarise(parental_median_relative_W_evo=median(median_relative_W_evo), 
                    parental_mean_relative_W_evo=mean(median_relative_W_evo)), 
        by=c("Ade","Time"), all.x = T) %>% 
  mutate(dif_median_relative_W.evol=median_relative_W_evo/parental_median_relative_W_evo, 
         dif_mean_relative_W.evol=mean_relative_W_evo/parental_mean_relative_W_evo) %>% 
  mutate(Time=factor(Time, levels = c("TS", "BS"))) %>% 
  group_by(Ade, Time, Group, Group_sample) %>% 
  summarise(median_dif_median_relative_W.evol=median(dif_median_relative_W.evol)) %>% 
  ggplot(aes(Group, (median_dif_median_relative_W.evol), fill=Group_sample)) +
  # ggplot(aes(Time, relative_W_evo, fill=Group_sample)) + 
  geom_point(aes(fill=Group_sample), size = 2, shape = 21, alpha=0.6) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  ylab(label = "Median Ratio Change Relative W") + 
  # ylab(label = "Log (Relative W)") + 
  # ylim(c(0,1000)) +
  # geom_point(aes(fill=Group_sample), alpha=0.8, size = 1.5, shape = 21, position = position_jitterdodge()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  # scale_y_continuous(breaks = seq((-3),8,1), 
  #                    labels = c(exp(seq((-3),8,1)))) +
  # scale_y_continuous(breaks = log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)),
  #                    labels = exp(log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)))) +
  # scale_y_continuous(trans = 'asinh',breaks=c(-20,-10,-5,-2,-1,0,1,5,10,20,50,100))+
  facet_grid(Time ~ Ade, scale="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.y = element_blank() , 
        strip.background = element_blank(), 
        legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("new_analyses_plots/04_2_relative_W_ecologicalSel_all_sim.png", width = 3, height = 5.7, dpi = 400)
  # ggsave("new_analyses_plots/04_2_relative_W_ecologicalSel_pink_sim.png", width = 7, height = 5.5, dpi = 400)
  ggsave("new_analyses_plots/04_5_change_relative_W_ecologicalSel_sim.png", width = 7, height = 5.5, dpi = 400)




#######

ecological_sel %>% filter(Group=="Allo" & Ade=="R") %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_W_evo), fill=evolved)) + 
  geom_boxplot() + 
  geom_point(aes(fill=evolved), size = 1.5, shape = 21, position = position_jitterdodge()) +
  facet_grid(Time ~ evolved, scale="free") #+
  ggsave("02_log_relative_W_allo_red.png", width = 10, height = 8)

ecological_sel %>% filter(Group=="Allo" & Ade=="P") %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_W_evo), fill=evolved)) + 
  geom_boxplot() + 
  geom_point(aes(fill=evolved), size = 1.5, shape = 21, position = position_jitterdodge()) +
  facet_grid(Time ~ evolved, scale="free") #+
  ggsave("02_log_relative_W_allo_pink.png", width = 10, height = 8)

ecological_sel %>% filter(Group=="Para") %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_W_evo), fill=evolved)) + 
  geom_boxplot() + 
  geom_point(aes(fill=evolved), size = 1.5, shape = 21, position = position_jitterdodge()) +
  facet_grid(Time ~ Ade, scale="free") #+
  ggsave("02_log_relative_W_para.png", width = 10, height = 8)

ecological_sel %>% filter(Group=="LM") %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_W_evo), fill=evolved)) + 
  geom_boxplot() + 
  geom_point(aes(fill=evolved), size = 1.5, shape = 21, position = position_jitterdodge()) +
  facet_grid(Time ~ Ade, scale="free") #+
  ggsave("02_log_relative_W_LM.png", width = 10, height = 8)

ecological_sel %>% filter(Group=="Sym") %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_W_evo), fill=evolved)) + 
  geom_boxplot() + 
  geom_point(aes(fill=evolved), size = 1.5, shape = 21, position = position_jitterdodge()) +
  facet_grid(Time ~ Ade, scale="free") +
  ggsave("02_log_relative_W_SYM.png", width = 10, height = 8)

ecological_sel %>% group_by(Structure, evolved, Pop, Time) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo)) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Structure) %>% 
  select(-sd_relative_W_evo) %>% 
  spread(Time, mean_relative_W_evo) %>% 
  ggplot(aes(log(TS), log(BS), colour=evolved)) + 
  geom_point() + geom_smooth(method='lm',formula=y~x) +
  facet_grid(Ade ~ Group )
  


head(inference_m.evol)

sem<-function(x){sd(x)/sqrt(length(x))}



# mean_inference_m.evol <- 
inference_m.evol %>% # filter(!(Group %in% c("P3P", "P3R"))) %>% 
  group_by(Structure, Group, evolved, Ade, Pop) %>% 
  summarise(mean_m.evol=mean(m.evol), low_lim=mean_m.evol-2*sem(m.evol), 
            high_lim=mean_m.evol+2*sem(m.evol), sd.evol=sd(m.evol)) %>% 
  ggplot(aes(as.factor(Group), as.numeric(mean_m.evol), colour=evolved)) +
  geom_boxplot() + 
  geom_jitter(width = 0.4) + 
  facet_grid(Ade ~ ., scale="free") #+
  ggsave("01_growth_rate.png", width = 8, height = 8)

head(mean_inference_m.evol)

mean_inference_m.evol <- inference_m.evol %>% filter(!(Group %in% c("P3P", "P3R"))) %>% 
  group_by(Structure, Group, evolved, Ade, Pop) %>% 
  summarise(mean_m.evol=mean(m.evol), low_lim=mean_m.evol-2*sem(m.evol), 
            high_lim=mean_m.evol+2*sem(m.evol), sd.evol=sd(m.evol))

mean_inference_m.evol_parentals<-inference_m.evol %>% filter(Group %in% c("P3P", "P3R")) %>% 
  group_by(Group) %>% 
  summarise(mean_m.evol=mean(m.evol), low_lim=mean_m.evol-2*sem(m.evol), 
            high_lim=mean_m.evol+2*sem(m.evol), sd.evol=sd(m.evol))

head(mean_inference_m.evol)

####  Mating measurements 
mating_data<-read.table("mating_data.txt", T) %>% 
  mutate(counts=red+wt, 
         red=red+1, 
         wt=wt+1) %>% 
  filter(counts>500) %>% 
  select(-counts)

mating_data_evolved_pop<- mating_data %>% 
  # filter(!(Group %in% c("P3P", "P3R"))) %>% 
  # group_by(Structure, Group, evolved, Ade, Pop) %>% 
  # merge(mean_inference_m.evol, by=c("Structure", "Group", "evolved", "Pop", "Ade")) 
  # group_by(Structure, Group, evolved, Pop, Ade, replicate) %>% 
  merge(median_inference_m.evol, by=c("Group", "evolved", "Pop", "Ade"), all.x=TRUE) 

# mating_data_parental_pop<- mating_data %>% filter(Group %in% c("P3P", "P3R")) %>% 
#   group_by(Group) %>% 
#   merge(mean_inference_m.evol_parentals, by="Group") 


# mating_data_ed<-rbind(mating_data_evolved_pop, mating_data_parental_pop) %>% 
mating_data_ed<-mating_data_evolved_pop %>% 
  mutate(saturate_ref=30000000, 
          N_ref_BM=(((saturate_ref/300)*20)*10/100)*25/1000, 
          N_evol_BM=(freq.before*N_ref_BM)/(1-freq.before), 
          N_cell_AM=100000, #it does not matter which number is used here
          N_evo_AMG=freq.after*N_cell_AM, 
          N_ref_AMG=N_cell_AM-N_evo_AMG, 
          growth_time_generations=(21/3.5)-1, 
          N_evol_AM=N_evo_AMG/(exp(median_m.evol*growth_time_generations)), 
          N_ref_AM=N_ref_AMG/(exp(growth_time_generations)), 
          r.AM=N_evol_AM/(N_evol_AM+N_ref_AM), 
          relative_v=(r.AM*N_ref_BM)/(N_evol_BM*(1-r.AM)))


mating_data_ed %>% 
  mutate(group_ID=paste(Group, evolved, sep="__")) %>% 
  filter(Ade=="R") %>% 
  mutate(population=ifelse(Group %in% c("P3P", "P3R"), "R",as.character(Pop))) %>% 
  ggplot(aes(as.factor(population), as.numeric(log(relative_v)), fill=Group)) + 
  geom_boxplot() + 
  xlab(label = "Population") +
  ylab(label = "Log (Relative W)") +
  facet_grid(group_ID ~ ., scale="free", space="free") +
  coord_flip() +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("plots/05_mating_perPop_red.png"), width = 6, height = 14, dpi = 400)


# colours_top_bottom<-c("red3","#D95F02","cornflowerblue")
colours_top_bottom<-c("red3","orange", "cornflowerblue") ; mating_data_ed %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
            sd_relative_v = sd(relative_v), 
            median_relative_v = median(relative_v)) %>% 
  #filter(Group_sample!="Parental") %>% 
  ungroup() %>% 
  #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  # filter(Ade=="R") %>%
  ggplot(aes(Group, log(median_relative_v), fill=Group_sample)) + 
  # ggplot(aes(Time, relative_W_evo, fill=Group_sample)) + 
  geom_boxplot(alpha=0.6) + 
  geom_point(aes(fill=Group_sample), alpha=0.8, size = 1.5, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  # ylab(label = "Log (Relative v)") + 
  ylab(label = "Relative viability (v)") + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  scale_y_continuous(breaks = log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)),
                     labels = exp(log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)))) +
  # ylim(c(0,1000)) +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.x = element_blank() , 
        strip.background = element_blank(), 
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("new_analyses_plots/06_0_relative_v_mating_all.png", width = 8, height = 5, dpi = 400)
  # ggsave("new_analyses_plots/06_0_relative_v_mating_red.png", width = 3, height = 3, dpi = 400)
  

colours_top_bottom<-c("red3","cornflowerblue") ; mating_data_ed %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
            sd_relative_v = sd(relative_v), 
            median_relative_v = median(relative_v)) %>% 
  filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  merge(mating_data_ed %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
                    sd_relative_v = sd(relative_v), 
                    median_relative_v = median(relative_v)) %>% 
          #filter(Group_sample!="Parental") %>% 
          ungroup() %>% 
          filter(Group_sample=="Parental") %>% 
          #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          group_by(Ade) %>% 
          summarise(parental_median_relative_v=median(median_relative_v), 
                    parental_mean_relative_v=mean(median_relative_v)), 
        by="Ade", all.x = T) %>% 
  mutate(dif_median_relative_v=median_relative_v/parental_median_relative_v, 
         dif_mean_relative_v=mean_relative_v/parental_mean_relative_v) %>% 
  filter(Ade=="R") %>%
  ggplot(aes(Group, log(dif_median_relative_v), fill=Group_sample)) + 
  # ggplot(aes(Time, relative_W_evo, fill=Group_sample)) + 
  geom_boxplot(alpha=0.6) + 
  geom_point(aes(fill=Group_sample), alpha=0.8, size = 1.5, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  # ylab(label = "Log (Relative v)") + 
  ylab(label = "Relative Mating Viability (v)") + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  scale_y_continuous(breaks = log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)),
                     labels = exp(log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)))) +
  # ylim(c(0,1000)) +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        strip.text.x = element_blank() ,
        strip.background = element_blank(),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("new_analyses_plots/06_1_change_relative_v_mating_pink.png"), width = 3, height = 3, dpi = 400)
  # ggsave(paste0("new_analyses_plots/06_1_change_relative_v_mating_red.png"), width = 3, height = 3, dpi = 400)
  # ggsave("new_analyses_plots/06_1_change_relative_v_mating_all.png", width = 8, height = 5, dpi = 400)
# ggsave("new_analyses_plots/06_1_change_relative_v_mating_red.png", width = 3, height = 3, dpi = 400)




colours_top_bottom<-c("red3","cornflowerblue") ; mating_data_ed %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
            sd_relative_v = sd(relative_v), 
            median_relative_v = median(relative_v)) %>% 
  filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  # this will add the parental parameters:
  merge(mating_data_ed %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
                    sd_relative_v = sd(relative_v), 
                    median_relative_v = median(relative_v)) %>% 
          #filter(Group_sample!="Parental") %>% 
          ungroup() %>% 
          filter(Group_sample=="Parental") %>% 
          #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          group_by(Ade) %>% 
          summarise(parental_median_relative_v=median(median_relative_v), 
                    parental_mean_relative_v=mean(median_relative_v)), 
        by="Ade", all.x = T) %>% 
  mutate(dif_median_relative_v=median_relative_v/parental_median_relative_v, 
         dif_mean_relative_v=mean_relative_v/parental_mean_relative_v) %>% 
  # this will add the variance component. In this one, I used the total variance rather than the variance per treatment.
  merge(mating_data_ed %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
                    sd_relative_v = sd(relative_v), 
                    median_relative_v = median(relative_v)) %>% 
          filter(Group_sample!="Parental") %>%
          ungroup() %>% 
          #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          # this will add the parental parameters:
          merge(mating_data_ed %>% 
                  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
                  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
                  summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
                            sd_relative_v = sd(relative_v), 
                            median_relative_v = median(relative_v)) %>% 
                  #filter(Group_sample!="Parental") %>% 
                  ungroup() %>% 
                  filter(Group_sample=="Parental") %>% 
                  #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
                  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
                  group_by(Ade) %>% 
                  summarise(parental_median_relative_v=median(median_relative_v), 
                            parental_mean_relative_v=mean(median_relative_v)), 
                by="Ade", all.x = T) %>% 
          mutate(dif_median_relative_v=median_relative_v/parental_median_relative_v, 
                 dif_mean_relative_v=mean_relative_v/parental_mean_relative_v) %>% 
          group_by(Ade) %>% 
          dplyr::summarize(sd_dif_median_relative_v=sd(dif_median_relative_v), 
                           sd_dif_mean_relative_v=sd(dif_mean_relative_v)), 
        by="Ade", all.x = T) %>% 
  group_by(Ade, Group, Group_sample) %>% 
  summarise(median_dif_median_relative_v=median(dif_median_relative_v), 
            median_sd_dif_median_relative_v=median(sd_dif_median_relative_v)) %>% 
  mutate(N_median_dif_median_relative_v=median_dif_median_relative_v/median_sd_dif_median_relative_v) %>% 
  ggplot(aes(Group, (N_median_dif_median_relative_v), fill=Group_sample)) + 
  geom_point(aes(fill=Group_sample), size = 2, shape = 21, alpha=0.6) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  # ylab(label = "Log (Relative v)") + 
  ylab(label = "N. Change Ratio Relative Viability (v)") + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  # scale_y_continuous(breaks = log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)),
  #                    labels = exp(log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)))) +
  # ylim(c(0,1000)) +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.x = element_blank() , 
        strip.background = element_blank(), 
        # legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("new_analyses_plots/06_4_totalVar_N_change_relative_v_mating_all.png", width = 8, height = 5, dpi = 400)
# ggsave("new_analyses_plots/06_4_totalVar_N_change_relative_v_mating_red.png", width = 3, height = 3, dpi = 400)


colours_top_bottom<-c("red3","cornflowerblue") ; mating_data_ed %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
            sd_relative_v = sd(relative_v), 
            median_relative_v = median(relative_v)) %>% 
  filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  # this will add the parental parameters:
  merge(mating_data_ed %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
                    sd_relative_v = sd(relative_v), 
                    median_relative_v = median(relative_v)) %>% 
          #filter(Group_sample!="Parental") %>% 
          ungroup() %>% 
          filter(Group_sample=="Parental") %>% 
          #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          group_by(Ade) %>% 
          summarise(parental_median_relative_v=median(median_relative_v), 
                    parental_mean_relative_v=mean(median_relative_v)), 
        by="Ade", all.x = T) %>% 
  mutate(dif_median_relative_v=median_relative_v/parental_median_relative_v, 
         dif_mean_relative_v=mean_relative_v/parental_mean_relative_v) %>% 
  # this will add the variance component. In this one, I used the total variance rather than the variance per treatment.
  merge(mating_data_ed %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
                    sd_relative_v = sd(relative_v), 
                    median_relative_v = median(relative_v)) %>% 
          filter(Group_sample!="Parental") %>%
          ungroup() %>% 
          #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          # this will add the parental parameters:
          merge(mating_data_ed %>% 
                  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
                  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
                  summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
                            sd_relative_v = sd(relative_v), 
                            median_relative_v = median(relative_v)) %>% 
                  #filter(Group_sample!="Parental") %>% 
                  ungroup() %>% 
                  filter(Group_sample=="Parental") %>% 
                  #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
                  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
                  group_by(Ade) %>% 
                  summarise(parental_median_relative_v=median(median_relative_v), 
                            parental_mean_relative_v=mean(median_relative_v)), 
                by="Ade", all.x = T) %>% 
          mutate(dif_median_relative_v=median_relative_v/parental_median_relative_v, 
                 dif_mean_relative_v=mean_relative_v/parental_mean_relative_v) %>% 
          group_by(Ade) %>% 
          dplyr::summarize(sd_dif_median_relative_v=sd(dif_median_relative_v), 
                           sd_dif_mean_relative_v=sd(dif_mean_relative_v)), 
        by="Ade", all.x = T) %>% 
  group_by(Ade, Group, Group_sample) %>% 
  summarise(median_dif_median_relative_v=median(dif_median_relative_v)) %>% 
  #mutate(N_mean_dif_median_relative_v=mean_dif_median_relative_v/meam_sd_dif_median_relative_v) %>% 
  ggplot(aes(Group, (median_dif_median_relative_v), fill=Group_sample)) + 
  geom_point(aes(fill=Group_sample), size = 2, shape = 21, alpha=0.6) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  # ylab(label = "Log (Relative v)") + 
  ylab(label = "Change Ratio Relative Viability (v)") + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  # scale_y_continuous(breaks = log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)),
  #                    labels = exp(log(c(0.1,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)))) +
  # ylim(c(0,1000)) +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        # strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # strip.text.x = element_blank() , 
        strip.background = element_blank(), 
        # legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("new_analyses_plots/06_5_change_relative_v_mating_all.png", width = 8, height = 5, dpi = 400)
# ggsave("new_analyses_plots/06_5_change_relative_v_mating_red.png", width = 3, height = 3, dpi = 400)



####




mating_data_ed %>% mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  summarise(mean_relative_v = mean(relative_v), sd_relative_v = sd(relative_v)) %>% 
  ggplot(aes(Group, log(mean_relative_v), fill=Group_sample)) + 
  geom_boxplot() + 
  geom_point(aes(fill=Group_sample), size = 1.5, shape = 21, position = position_jitterdodge()) +
  facet_grid(. ~ Ade , scale="free", space="free") #+
  ggsave("03_log_relative_v_mating.png", width = 10, height = 8)



mating_data_ed %>%
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_v), fill=evolved)) + 
  geom_boxplot() + 
  facet_grid(Ade ~ Group, scale="free_x" )


mating_data_ed %>% filter(Group=="Allo") %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_v), fill=evolved)) + 
  geom_boxplot() + 
  facet_grid(Ade ~ evolved, scale="free_x" ) #+
  ggsave("03_log_relative_v_mating_allo.png", width = 10, height = 8)

mating_data_ed %>% filter(Group=="Para") %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_v), fill=evolved)) + 
  geom_boxplot() + 
  facet_grid(. ~ Ade, scale="free_x" ) +
  ggsave("03_log_relative_v_mating_para.png", width = 10, height = 8)

mating_data_ed %>% filter(Group=="LM") %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_v), fill=evolved)) + 
  geom_boxplot() + 
  facet_grid(. ~ Ade, scale="free_x" ) +
  ggsave("03_log_relative_v_mating_LM.png", width = 10, height = 8)

mating_data_ed %>% filter(Group=="Sym") %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  ggplot(aes(Pop, log(relative_v), fill=evolved)) + 
  geom_boxplot() + 
  facet_grid(. ~ Ade, scale="free_x" ) + 
  ggsave("03_log_relative_v_mating_SYM.png", width = 10, height = 8)

mating_data_ed %>% filter(!(Group %in% c("Allo", "P3P", "P3R"))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(mean_relative_v = mean(relative_v), sd_relative_v = sd(relative_v)) %>% 
  group_by(Pop, Group, Ade) %>% select(-sd_relative_v, -sample, -Fraction) %>% 
  spread(evolved, mean_relative_v) %>% mutate(dif_relative_v_bottom_top=log(BOTTOM)-log(TOP)) %>% 
  ggplot(aes(Group, dif_relative_v_bottom_top)) + 
  geom_boxplot() + 
  geom_jitter(size = 1.5, shape = 21) +
  facet_grid(. ~ Ade, scale="free_x" ) +
  ggsave("03_log_relative_v_mating_SYM.png", width = 10, height = 8)



#### ecologycal selection Vs Mating

eco_data_summary<-ecological_sel %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
  # filter(Ade=="R") %>%
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
            median_relative_W_evo = median(relative_W_evo)) %>%
  #filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  merge(ecological_sel %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
          #filter(Ade=="R") %>%
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
          summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
                    median_relative_W_evo = median(relative_W_evo)) %>%
          #filter(Group_sample!="Parental") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                              levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade, Time) %>% 
          summarise(parental_median_relative_W_evo=median(median_relative_W_evo), 
                    parental_mean_relative_W_evo=mean(median_relative_W_evo)), 
        by=c("Ade","Time"), all.x = T) %>% 
  # using the difference between parental and evolved:
  # mutate(dif_median_relative_W.evol=median_relative_W_evo-parental_median_relative_W_evo, 
  #        dif_mean_relative_W.evol=mean_relative_W_evo-parental_mean_relative_W_evo) %>% 
  # using the ratio intead:
  mutate(dif_median_relative_W.evol=median_relative_W_evo/parental_median_relative_W_evo, 
         dif_mean_relative_W.evol=mean_relative_W_evo/parental_mean_relative_W_evo) %>% 
  mutate(Time=factor(Time, levels = c("TS", "BS"))) %>% 
  select(Fraction, Pop, sample, Group, Ade, evolved, 
         Structure, Group_sample, Time, dif_median_relative_W.evol) 


mating_data_summary<-mating_data_ed %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
            sd_relative_v = sd(relative_v), 
            median_relative_v = median(relative_v)) %>% 
  filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  merge(mating_data_ed %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
                    sd_relative_v = sd(relative_v), 
                    median_relative_v = median(relative_v)) %>% 
          #filter(Group_sample!="Parental") %>% 
          ungroup() %>% 
          filter(Group_sample=="Parental") %>% 
          #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          group_by(Ade) %>% 
          summarise(parental_median_relative_v=median(median_relative_v), 
                    parental_mean_relative_v=mean(median_relative_v)), 
        by="Ade", all.x = T) %>% 
  mutate(dif_median_relative_v=median_relative_v/parental_median_relative_v, 
         dif_mean_relative_v=mean_relative_v/parental_mean_relative_v) %>% 
  select(Fraction, Pop, sample, Group, Ade, evolved, 
         Structure, Group_sample, Time, dif_median_relative_v) 



# eco_mating_data<-
mating_data_ed %>%
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, 
            evolved, Structure, Group_sample) %>% 
  summarise(N_replicates=n(), median_relative_v = median(relative_v), 
            sd_relative_v = sd(relative_v)) %>% 
  select(-N_replicates) %>% 
  merge((ecological_sel %>% 
           mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), 
                                      "Parental",as.character(evolved))) %>%
           group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
           summarise(median_relative_W_evo = median(relative_W_evo), 
                     sd_relative_W_evo = sd(relative_W_evo))), 
        by=c("Fraction", "Pop", "sample", "Group", "Ade", "evolved", 
             "Structure", "Group_sample"), all=TRUE) %>% head()


colours_top_bottom<-c("red3","cornflowerblue")

treatment_used<-"Allo"
eco_mating_data %>% 
  filter(Group==treatment_used) %>% 
  ggplot(aes(x=log(median_relative_W_evo), y=log(median_relative_v))) +
  geom_point(aes(colour=evolved)) +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  scale_colour_manual(values=colours_top_bottom) +
  xlab(label = "Log (Relative W)") +
  ylab(label = "Log (Relative V)") + 
  ggtitle(treatment_used) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"))+ 
  facet_grid(Ade ~ Time, scale="free") #+
  ggsave(paste0("plots/07_EcologicalSel_vs_mating_", treatment_used, ".png"), width = 8, height = 6)
  


eco_mating_data %>% 
  filter(Group=="Para") %>% 
  group_by(Pop, Ade, Time) %>% 
  summarise(N_pop=n(), relative_v = mean(median_relative_v), 
            relative_w = mean(median_relative_W_evo)) %>% 
  ggplot(aes(x=log(relative_w), y=log(relative_v))) +
  geom_point() +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  # scale_colour_manual(values=colours_top_bottom) +
  xlab(label = "Log (Relative W)") +
  ylab(label = "Log (Relative V)") + 
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(Ade ~ Time, scale="free") 


colours_top_bottom<-c("red3","#D95F02","cornflowerblue")
mating_data_ed %>% mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  summarise(N_replicates=n(), mean_relative_v = median(relative_v), 
            sd_relative_v = sd(relative_v)) %>% 
  ggplot(aes(Group, log(mean_relative_v), fill=Group_sample)) + 
  # ggplot(aes(Time, relative_W_evo, fill=Group_sample)) + 
  geom_boxplot(alpha=0.6) + 
  geom_point(aes(fill=Group_sample), alpha=0.8, size = 1.5, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values=colours_top_bottom) +
  xlab(label = "Treatment") +
  ylab(label = "Log (Relative W)") + 
  # ylim(c(0,1000)) +
  facet_grid(. ~ Ade , scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  ggsave("plots/06_log_relative_W_mating.png", width = 7, height = 5)



ecological_sel %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo)) %>%
  names()
  head()
  
eco_mating_data<-merge((mating_data_ed %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(mean_relative_v = mean(relative_v), mean_m.evol=mean_m.evol[1], sd_relative_v = sd(relative_v)) %>% 
    ungroup()), 
  (ecological_sel %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo))%>% 
    ungroup()), 
  by=c("Fraction", "Pop", "sample", "Group", "Ade", "evolved", "Structure"))


head(eco_mating_data)

eco_mating_data %>% 
  filter(Group=="Allo") %>% 
  ggplot(aes(x=log(mean_relative_W_evo), y=log(mean_relative_v))) +
  geom_point(aes(colour=evolved)) +
  geom_smooth(method='lm',formula=y~x) + 
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  facet_grid(Ade ~ Time.y, scale="free") #+
  ggsave("04_EcologicalSel_vs_mating_allo.png", width = 8, height = 8)


eco_mating_data %>% 
  filter(Group=="Para") %>% 
  ggplot(aes(x=log(mean_relative_W_evo), y=log(mean_relative_v))) +
  geom_point(aes(shape=evolved)) +
  #geom_smooth(method='lm',formula=y~x) + 
  geom_line(aes(colour=Pop)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  facet_grid(Ade ~ Time.y, scale="free") + 
  ggsave("04_EcologicalSel_vs_mating_Para.png", width = 8, height = 8)

eco_mating_data %>% 
  filter(Group=="LM") %>% 
  ggplot(aes(x=log(mean_relative_W_evo), y=log(mean_relative_v))) +
  geom_point(aes(shape=evolved)) +
  #geom_smooth(method='lm',formula=y~x) + 
  geom_line(aes(colour=Pop)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  facet_grid(Ade ~ Time.y, scale="free") +
  ggsave("04_EcologicalSel_vs_mating_LM.png", width = 8, height = 8)
  
eco_mating_data %>% 
  filter(Group=="Sym") %>% 
  ggplot(aes(x=log(mean_relative_W_evo), y=log(mean_relative_v))) +
  geom_point(aes(shape=evolved)) +
  #geom_smooth(method='lm',formula=y~x) + 
  geom_line(aes(colour=Pop)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  facet_grid(Ade ~ Time.y, scale="free") +
  ggsave("04_EcologicalSel_vs_mating_SYM.png", width = 8, height = 8)







library(ggpmisc)
my.formula <- y ~ x
eco_mating_data %>% 
  filter(Group=="Allo") %>% 
  ggplot(aes(x=log(mean_relative_W_evo), y=log(mean_relative_v), colour=evolved)) +
  geom_point() +
  geom_smooth(method='lm',formula=y~x) + 
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  facet_grid(Ade ~ Time.y)


eco_mating_data %>% 
  filter(Group=="Allo") %>% 
  ggplot(aes(x=log(mean_relative_W_evo), y=(mean_m.evol), colour=evolved)) +
  geom_point() +
  geom_smooth(method='lm',formula=y~x) + 
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  facet_grid(Ade ~ Time.y)






#######

# Genotype Phenotype correlation:
library("RColorBrewer")
library(ggfortify)
library(factoextra)

# Allopatric treatments:
# red populations:
phaseIII_table_total_PCA_red<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>%
  filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))


plot_PCA_red<-autoplot(PCA_data, 
                       data=phaseIII_table_total_PCA_red, 
                       label = F , 
                       colour = 'treatment', 
                       fill =  'treatment', 
                       shape=T, label.size = 5, 
                       frame=T)


plot_PCA_red

phaseIII_table_total_Genetic_PCA_PCvalues_red <- data.frame(treatment=plot_PCA_red$data$treatment, time=plot_PCA_red$data$time,  pop_ID=plot_PCA_red$data$pop_ID, pc1=as.numeric(plot_PCA_red$data$PC1), pc2=as.numeric(plot_PCA_red$data$PC2)) 

# pink populations:
phaseIII_table_total_PCA_pink<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P", "SYM", "LM", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>%
  filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("pink3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_pink %>% select(-treatment, -time, -pop_ID))


plot_PCA_pink<-autoplot(PCA_data, 
                        data=phaseIII_table_total_PCA_pink, 
                        label = F , 
                        colour = 'treatment', 
                        fill =  'treatment', 
                        shape=T, label.size = 5, 
                        frame=T)

phaseIII_table_total_Genetic_PCA_PCvalues_pink <- data.frame(treatment=plot_PCA_pink$data$treatment, time=plot_PCA_pink$data$time,  pop_ID=plot_PCA_pink$data$pop_ID, pc1=as.numeric(plot_PCA_pink$data$PC1), pc2=as.numeric(plot_PCA_pink$data$PC2)) 



## Variance covariance matrix:

head(median_inference_m.evol_cor)
head(mating_data_summary)
head(eco_data_summary)

str(eco_data_summary)

sim_growth<-median_inference_m.evol_cor %>% 
  mutate(pop_ID=paste(Pop, Group, Fraction, sep="_")) %>% 
  ungroup() %>% 
  # mutate(log_median_m.evol=log(dif_median_m.evol)) %>% 
  # select(pop_ID, Ade, log_median_m.evol)
  mutate(log_median_m.evol=log(dif_median_m.evol)) %>%
  select(pop_ID, Ade, log_median_m.evol) 

head(eco_data_summary)
sim_ecological_sel<-eco_data_summary %>% 
  select(-Group_sample) %>% 
  ungroup() %>% 
  spread(Time, dif_median_relative_W.evol) %>% 
  mutate(pop_ID=paste(Pop, Group, Fraction, sep="_")) %>% 
  mutate(relativeWBS=BS, relativeWTS=TS) %>% 
  mutate(log_relativeWBS=log(relativeWBS), log_relativeWTS=log(relativeWTS)) %>% 
  select(pop_ID, Ade, log_relativeWBS, log_relativeWTS)

sim_phenotype_completed<-mating_data_summary %>% 
  mutate(pop_ID=paste(Pop, Group, Fraction, sep="_"), 
         log_median_relative_v=log(dif_median_relative_v)) %>% 
  #select(pop_ID, Ade, log_median_relative_v) %>% 
  merge(sim_growth, 
        by=c("Ade", "pop_ID"), all = T) %>% 
  merge(sim_ecological_sel, 
        by=c("Ade", "pop_ID"), all = T) 

head(sim_phenotype_completed)
Ade_used<-"P"
Cov_M <- list()
Cor_M <- list()
for (treatment in c("Allo","Para", "LM", "Sym")) {
  sim_phenotype <- sim_phenotype_completed %>%
    filter(Ade==Ade_used, Group==treatment)
  # sim_phenotype <- sim_phenotype_completed %>% 
  #   filter(Ade==Ade_used, Group!="Parental")
  # N_v<-c((sim_phenotype$log_median_relative_v-mean(sim_phenotype$log_median_relative_v))/sd(sim_phenotype$log_median_relative_v))
  # N_m<-c((sim_phenotype$median_m.evol-mean(sim_phenotype$median_m.evol))/sd(sim_phenotype$median_m.evol))
  # N_WBS<-c((sim_phenotype$log_relativeWBS-mean(sim_phenotype$log_relativeWBS))/sd(sim_phenotype$log_relativeWBS))
  # N_WTS<-c((sim_phenotype$log_relativeWTS-mean(sim_phenotype$log_relativeWTS))/sd(sim_phenotype$log_relativeWTS))
  N_v<-sim_phenotype$log_median_relative_v
  N_m<-sim_phenotype$log_median_m.evol
  N_WBS<-sim_phenotype$log_relativeWBS
  N_WTS<-sim_phenotype$log_relativeWTS
  table_perTreatment<-sim_phenotype %>% 
    cbind(N_v, N_m, N_WBS, N_WTS) %>% 
    filter(Group==treatment) %>% 
    select(-dif_median_relative_v, -log_median_m.evol, -log_median_relative_v, 
           -log_relativeWBS, -log_relativeWTS)
  table_perTreatment_sim<-table_perTreatment %>% select(N_v, N_m, N_WBS, N_WTS)
  Cov_M[[treatment]]<-cov(table_perTreatment_sim)
  Cor_M[[treatment]]<-cov2cor(cov(table_perTreatment_sim))
}


# for red populations only
Cov_M_R<-Cov_M
Cor_M_R<-Cor_M

# for pink populations only
Cov_M_P<-Cov_M
Cor_M_P<-Cor_M

####

Ade_used<-"R"
Cov_M <- data.frame()
Cor_M <- data.frame()
for (treatment in c("Allo","Para", "LM", "Sym")) {
treatment<-"Allo"
  sim_phenotype <- sim_phenotype_completed %>%
    filter(Ade==Ade_used, Group==treatment)
  N_mating<-sim_phenotype$log_median_relative_v
  N_growth<-sim_phenotype$log_median_m.evol
  N_WBS<-sim_phenotype$log_relativeWBS
  N_WTS<-sim_phenotype$log_relativeWTS
  table_perTreatment<-sim_phenotype %>% 
    cbind(N_mating, N_growth, N_WBS, N_WTS) %>% 
    filter(Group==treatment) %>% 
    select(-dif_median_relative_v, -log_median_m.evol, -log_median_relative_v, 
           -log_relativeWBS, -log_relativeWTS)
  table_perTreatment_sim<-table_perTreatment %>% select(N_mating, N_growth, N_WBS, N_WTS)
  Cov_M<-rbind(Cov_M, data.frame(cov(table_perTreatment_sim)) %>% 
                 mutate(treatment=treatment, 
                        comp_var=row.names(cov(table_perTreatment_sim))))
  Cor_M<-rbind(Cor_M, data.frame(cov2cor(cov(table_perTreatment_sim))) %>% 
                 mutate(treatment=treatment, 
                        comp_var=row.names(cov2cor(cov(table_perTreatment_sim)))))
}



Cor_M %>% 
  gather("comp_var2", "value" ,c(N_mating, N_growth, N_WBS, N_WTS)) %>% 
  mutate(treatment=factor(treatment, levels=c("Allo","Para", "LM", "Sym")), 
         comp_var=factor(comp_var, levels=c("N_growth", "N_mating", "N_WTS", "N_WBS")), 
         comp_var2=factor(comp_var2, levels=c("N_growth", "N_mating", "N_WTS", "N_WBS"))) %>%
  filter(comp_var!=comp_var2) %>%
  filter((comp_var2=="N_growth" & comp_var %in% c("N_mating", "N_WTS", "N_WBS")) |
           (comp_var2=="N_mating" & comp_var %in% c("N_WTS", "N_WBS")) |
           (comp_var2=="N_WTS" & comp_var %in% c("N_WBS"))) %>%
  #Using Bars and colours:
  ggplot(aes(treatment, value, fill=value)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(-1,1)) +
  scale_fill_gradient2(low = "#3F007D", mid = "white", high = "#D95F02", 
                       midpoint = 0, limits = c(-1,1)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  facet_grid(comp_var ~ comp_var2) +
  theme_classic() +
  theme(axis.line = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    plot.title = element_blank(), 
    strip.background = element_blank(), 
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0, size=10, colour="black"), 
    axis.text.y = element_text(size=8, colour="black")) #+
  ggsave(paste0("new_analyses_plots/10_correlation_phenotypes_", Ade_used,".png"), 
         width = 4, height = 2.5, dpi = 400)
#


  
data.frame(matrix(unlist(Cor_M_R), nrow=length(Cor_M_R), byrow=F))

Cor_M
for ( i in Cor_M_R){
  print("A")
  print(i)
}

cor(Cor_M_R$LM, Cor_M_R$Para)
eigen(Cor_M_R$Allo)

Cor_M_R$Allo
Cor_M$Allo
Cov_M_R$Allo
Cov_M

eigen(Cov_M_R$Allo)
eigen(Cov_M_R$Sym)

KrzCor(Cov_M_R)

MatrixDistance(Cov_M_R, distance = "RiemannDist")
MatrixDistance(Cov_M_R, distance = "OverlapDist")

MatrixDistance(Cov_M_P, distance = "RiemannDist")
MatrixDistance(Cov_M_P, distance = "OverlapDist")

MatrixDistance(Cov_M_R, Cov_M_P[["Allo"]])

library(evolqg)
KrzCor(Cov_M_P, ret.dim = 1)
KrzCor(Cov_M_R[["Allo"]], Cov_M_P[["Allo"]])

KrzCor(c(Cov_M_R, Cov_M_P))
c(Cov_M_R, Cov_M_P)

library(car)
lower.tri(Cor_M_R$Allo)
MantelCor(lower.tri(Cor_M_R$Allo), Cor_M_R$Sym, permutations = 1000)
KrzCor(Cov_M_R[["Allo"]], Cov_M_R[["Sym"]])
sqrt(0.528)

laply(Cov_M_R, MonteCarloRep, 30, KrzCor, correlation = TRUE)




library(evolqg)
c1 <- RandomMatrix(10, 1, 1, 10)
c2 <- RandomMatrix(10, 1, 1, 10)
c3 <- RandomMatrix(10, 1, 1, 10)
KrzCor(c1, c2)

KrzCor(Cor_M_P)

reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRep, 10, KrzCor, iterations = 10))
KrzCor(list(c1, c2, c3), repeat.vector = reps)

c4 <- RandomMatrix(10)
KrzCor(list(c1, c2, c3), c4)



library("RColorBrewer")
library(ggfortify)
library(factoextra)

Ade_used<-"R"
treatment<-"Allo" # "Allo","Para", "LM", "Sym"
for (treatment in c("Allo","Para", "LM", "Sym")) {
sim_phenotype <- sim_phenotype_completed %>%
    filter(Ade==Ade_used, Group==treatment) 
  # sim_phenotype <- sim_phenotype_completed %>% 
  #   filter(Ade==Ade_used, Group!="Parental")
N_v<-c((sim_phenotype$log_median_relative_v-mean(sim_phenotype$log_median_relative_v))/sd(sim_phenotype$log_median_relative_v))
N_m<-c((sim_phenotype$log_median_m.evol-mean(sim_phenotype$log_median_m.evol))/sd(sim_phenotype$log_median_m.evol))
N_WBS<-c((sim_phenotype$log_relativeWBS-mean(sim_phenotype$log_relativeWBS))/sd(sim_phenotype$log_relativeWBS))
N_WTS<-c((sim_phenotype$log_relativeWTS-mean(sim_phenotype$log_relativeWTS))/sd(sim_phenotype$log_relativeWTS))

sim_phenotype_PCAtable<-sim_phenotype %>% 
  cbind(N_v, N_m, N_WBS, N_WTS) %>% 
  #filter(Group=="Allo") %>% 
  select(-dif_median_relative_v, -log_median_m.evol, -log_median_relative_v, 
         -log_relativeWBS, -log_relativeWTS) 
PCA_data_phenotype<-prcomp(sim_phenotype_PCAtable %>% select(N_v, N_m, N_WBS, N_WTS))

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")
autoplot(PCA_data_phenotype, 
         data=sim_phenotype_PCAtable, 
         label = F , 
         colour = 'evolved', 
         fill =  'evolved', 
         shape=T, label.size = 5, 
         frame=T) +
  scale_colour_manual(values=colours_PCA) +
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+

fviz_pca_var(PCA_data_phenotype,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = T ,
             repel = F ) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_blank(), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
# ggsave("plots/08_PCA_loads_allo_red_nl.png", width = 4, height = 3, dpi = 400)

# To save the final plots:
PCA_point<-autoplot(PCA_data_phenotype, 
         data=sim_phenotype_PCAtable, 
         label = F , 
         colour = 'evolved', 
         fill =  'evolved', 
         shape=T, label.size = 5, 
         frame=F) +
  scale_colour_manual(values=colours_PCA) 

pPCA_point <- ggplot_build(PCA_point)

PCA_point_data<-pPCA_point$data[[1]] %>% 
  mutate(PC1=as.numeric(x), PC2=as.numeric(y))

head(pPCA_point$data[[1]])
head(PCA_point_data)

fviz_pca_var(PCA_data_phenotype,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label = F ,
             repel = T ) +
  geom_point(aes(PC1,PC2), 
             colour=as.factor(PCA_point_data$colour), 
             data=PCA_point_data, size=2, alpha=0.7) +
  # scale_colour_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_blank(), 
        # legend.position="none",
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) +
  ggsave(paste0("new_analyses_plots/08_PCA_loads_", treatment, "_", Ade_used, ".png"), width = 4, height = 3, dpi = 400)


fviz_pca_var(PCA_data_phenotype,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #label = F ,
             repel = T ) +
  geom_point(aes(PC1,PC2), colour=as.factor(PCA_point_data$colour), data=PCA_point_data, size=2, alpha=0.7) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_blank(), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) +
  ggsave(paste0("new_analyses_plots/08_PCA_loads_", treatment, "_", Ade_used, "_loadIDs.png"), width = 4, height = 3, dpi = 400)
}




###### Final Genotype - Phenotype correlation :

library("RColorBrewer")
library(ggfortify)
library(factoextra)

# Genotype tables:

# Allopatric treatments:
# red populations:
phaseIII_table_total_PCA_red<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>%
  filter(treatment %in% c("Allo_T", "Allo_B", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  # filter(treatment %in% c("Para_T", "Para_B")) %>%
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))

plot_PCA_red<-autoplot(PCA_data, 
                       data=phaseIII_table_total_PCA_red, 
                       label = F , 
                       colour = 'treatment', 
                       fill =  'treatment', 
                       shape=T, label.size = 5, 
                       frame=T)


phaseIII_table_total_Genetic_PCA_PCvalues_red <- data.frame(treatment=plot_PCA_red$data$treatment, time=plot_PCA_red$data$time,  pop_ID=plot_PCA_red$data$pop_ID, pc1=as.numeric(plot_PCA_red$data$PC1), pc2=as.numeric(plot_PCA_red$data$PC2)) 


# pink populations:
phaseIII_table_total_PCA_pink<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P", "SYM", "LM", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  filter(treatment %in% c("Para_T", "Para_B")) %>%
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("pink3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_pink %>% select(-treatment, -time, -pop_ID))


plot_PCA_pink<-autoplot(PCA_data, 
                        data=phaseIII_table_total_PCA_pink, 
                        label = F , 
                        colour = 'treatment', 
                        fill =  'treatment', 
                        shape=T, label.size = 5, 
                        frame=T)

plot_PCA_pink

phaseIII_table_total_Genetic_PCA_PCvalues_pink <- data.frame(treatment=plot_PCA_pink$data$treatment, time=plot_PCA_pink$data$time,  pop_ID=plot_PCA_pink$data$pop_ID, pc1=as.numeric(plot_PCA_pink$data$PC1), pc2=as.numeric(plot_PCA_pink$data$PC2)) 

# phenotype table:
Ade_used<-"R"
# treatment<-"Para" # "Allo","Para", "LM", "Sym"
# for (treatment in c("Allo","Para", "LM", "Sym")) {

# sim_phenotype <- sim_phenotype_completed %>%
#   filter(Ade==Ade_used, Group==treatment) 

sim_phenotype <- sim_phenotype_completed %>%
  filter(Ade==Ade_used, Group!="Parental")

# sim_phenotype <- sim_phenotype_completed %>%
# filter(Ade==Ade_used, Group %in% c("Allo","Sym"))

N_v<-c((sim_phenotype$log_median_relative_v-mean(sim_phenotype$log_median_relative_v))/sd(sim_phenotype$log_median_relative_v))
N_m<-c((sim_phenotype$log_median_m.evol-mean(sim_phenotype$log_median_m.evol))/sd(sim_phenotype$log_median_m.evol))
N_WBS<-c((sim_phenotype$log_relativeWBS-mean(sim_phenotype$log_relativeWBS))/sd(sim_phenotype$log_relativeWBS))
N_WTS<-c((sim_phenotype$log_relativeWTS-mean(sim_phenotype$log_relativeWTS))/sd(sim_phenotype$log_relativeWTS))

sim_phenotype_PCAtable<-sim_phenotype %>% 
  cbind(N_v, N_m, N_WBS, N_WTS) %>% 
  #filter(Group=="Allo") %>% 
  select(-dif_median_relative_v, -log_median_m.evol, -log_median_relative_v, 
         -log_relativeWBS, -log_relativeWTS) 

PCA_data_phenotype<-prcomp(sim_phenotype_PCAtable %>% select(N_v, N_m, N_WBS, N_WTS))

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")
plot_PCA_phenotype<-autoplot(PCA_data_phenotype, 
                             data=sim_phenotype_PCAtable, 
                             label = F , 
                             colour = 'evolved', 
                             fill =  'evolved', 
                             shape=T, label.size = 5, 
                             frame=T) 

plot_PCA_phenotype

table_total_phenotype_PCA_PCvalues<-data.frame(pop_ID=plot_PCA_phenotype$data$pop_ID, 
                                               pc1_phenotype=as.numeric(plot_PCA_phenotype$data$PC1), 
                                               pc2_phenotype=as.numeric(plot_PCA_phenotype$data$PC2))


table_total_phenotype_PCA_PCvalues %>% 
  separate(pop_ID,c("pop","treatment2", "fraction"), sep="_") %>% 
  mutate(treatment2=factor(treatment2, levels=c("Allo","Para", "LM", "Sym"))) %>% 
  mutate(phe_dis=sqrt((pc1_phenotype^2)+(pc2_phenotype^2))) %>% 
  # ggplot(aes(treatment2, phe_dis)) +
  # geom_boxplot() +
  ggplot(aes(pc1_phenotype,pc2_phenotype)) +
  geom_point(aes(colour=fraction), size=2, alpha=0.8)+
  geom_line(aes(group=pop), data=table_total_phenotype_PCA_PCvalues %>% 
              separate(pop_ID,c("pop","treatment2", "fraction"), sep="_") %>% 
              mutate(treatment2=factor(treatment2, levels=c("Allo","Para", "LM", "Sym"))) %>% 
              filter(treatment2!="Allo")) +
  scale_colour_manual(values=c("red3", "cornflowerblue"))+
  facet_grid(. ~ treatment2) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_blank(), 
        #legend.position="none", 
        axis.text.x=element_text(angle = 45, hjust = 1,size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
  # ggsave("new_analyses_plots/12_PCA_phenotype_pairPop_red.png", width = 10, height = 3, dpi = 400)
  # ggsave("new_analyses_plots/12_PCA_phenotype_pairPop_pink.png", width = 10, height = 3, dpi = 400)


table_total_phenotype_PCA_PCvalues_para_random<-c()
for (i in seq(1,100)) {
  table_total_phenotype_PCA_PCvalues_para_random[i]<-table_total_phenotype_PCA_PCvalues %>% 
  separate(pop_ID,c("pop","treatment2", "fraction"), sep="_") %>% 
  mutate(treatment2=factor(treatment2, levels=c("Allo","Para", "LM", "Sym"))) %>% 
  filter(treatment2=="Para") %>% 
  filter(fraction=="T") %>%
  rowwise() %>% 
  mutate(samplepop=runif(1, 0, 1)) %>% 
  arrange(samplepop) %>% 
  mutate(fraction.x= fraction, pc1_phenotype.x=pc1_phenotype, pc2_phenotype.x=pc2_phenotype) %>% 
  mutate(treatment2="Para*") %>% 
  select(treatment2, pop, fraction.x, pc1_phenotype.x, pc2_phenotype.x) %>% 
  cbind(table_total_phenotype_PCA_PCvalues %>% 
          separate(pop_ID,c("pop","treatment2", "fraction"), sep="_") %>% 
          mutate(treatment2=factor(treatment2, levels=c("Allo","Para", "LM", "Sym"))) %>% 
          filter(treatment2=="Para") %>% 
          filter(fraction=="B") %>%
          rowwise() %>% 
          mutate(samplepop=runif(1, 0, 1)) %>% 
          arrange(samplepop) %>% 
          mutate(fraction.y= fraction, 
                 pc1_phenotype.y=pc1_phenotype, 
                 pc2_phenotype.y=pc2_phenotype) %>% 
          select(fraction.y, pc1_phenotype.y, pc2_phenotype.y)) %>% 
  mutate(phe_dis=sqrt(((pc1_phenotype.x-pc1_phenotype.y)^2)+((pc2_phenotype.x-pc2_phenotype.y)^2))) %>% 
  ungroup() %>% 
  select(phe_dis) %>% 
  unlist() %>% 
  as.vector() %>% 
  mean()
}

table_total_phenotype_PCA_PCvalues_allo_random<-c()
for (i in seq(1,100)) {
table_total_phenotype_PCA_PCvalues_allo_random[i]<-table_total_phenotype_PCA_PCvalues %>%
  separate(pop_ID,c("pop","treatment2", "fraction"), sep="_") %>% 
  mutate(treatment2=factor(treatment2, levels=c("Allo","Para", "LM", "Sym"))) %>% 
  filter(treatment2=="Allo") %>% 
  filter(fraction=="T") %>%
  rowwise() %>% 
  mutate(samplepop=runif(1, 0, 1)) %>% 
  arrange(samplepop) %>% 
  mutate(fraction.x= fraction, pc1_phenotype.x=pc1_phenotype, pc2_phenotype.x=pc2_phenotype) %>% 
  mutate(treatment2="Allo*") %>% 
  select(treatment2, pop, fraction.x, pc1_phenotype.x, pc2_phenotype.x) %>% 
  cbind(table_total_phenotype_PCA_PCvalues %>% 
          separate(pop_ID,c("pop","treatment2", "fraction"), sep="_") %>% 
          mutate(treatment2=factor(treatment2, levels=c("Allo","Para", "LM", "Sym"))) %>% 
          filter(treatment2=="Allo") %>% 
          filter(fraction=="B") %>%
          rowwise() %>% 
          mutate(samplepop=runif(1, 0, 1)) %>% 
          arrange(samplepop) %>% 
          mutate(fraction.y= fraction, 
                 pc1_phenotype.y=pc1_phenotype, 
                 pc2_phenotype.y=pc2_phenotype) %>% 
          select(fraction.y, pc1_phenotype.y, pc2_phenotype.y)) %>% 
  mutate(phe_dis=sqrt(((pc1_phenotype.x-pc1_phenotype.y)^2)+((pc2_phenotype.x-pc2_phenotype.y)^2))) %>% 
  ungroup() %>% 
  select(phe_dis) %>% 
  unlist() %>% 
  as.vector() %>% 
  mean()
}

table_total_phenotype_PCA_PCvalues %>% 
  separate(pop_ID,c("pop","treatment2", "fraction"), sep="_") %>% 
  mutate(treatment2=factor(treatment2, levels=c("Allo","Para", "LM", "Sym"))) %>% 
  filter(fraction=="T") %>%
  merge(table_total_phenotype_PCA_PCvalues %>% 
          separate(pop_ID,c("pop","treatment2", "fraction"), sep="_") %>% 
          mutate(treatment2=factor(treatment2, levels=c("Allo","Para", "LM", "Sym"))) %>% 
          filter(fraction=="B"), 
        by=c("treatment2","pop"), all.x = T) %>% 
  filter(treatment2!="Allo") %>%
  mutate(phe_dis=sqrt(((pc1_phenotype.x-pc1_phenotype.y)^2)+((pc2_phenotype.x-pc2_phenotype.y)^2))) %>% 
  select(treatment2, phe_dis) %>% 
  rbind(data.frame(treatment2="Allo*", 
                   phe_dis=table_total_phenotype_PCA_PCvalues_allo_random)) %>% 
  rbind(data.frame(treatment2="Para*", 
                   phe_dis=table_total_phenotype_PCA_PCvalues_para_random)) %>% 
  mutate(treatment2=factor(treatment2, levels=c("Allo","Allo*","Para","Para*","LM", "Sym"))) %>%
  ggplot(aes(treatment2, phe_dis)) +
  geom_boxplot() +
  # geom_point()+
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_blank(), 
        #legend.position="none", 
        axis.text.x=element_text(angle = 45, hjust = 1,size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) +
ggsave("new_analyses_plots/13_PCA_phenotype_dis_pairPop_red.png", width = 5, height = 4, dpi = 400)
# ggsave("new_analyses_plots/13_PCA_phenotype_dis_pairPop_pink.png", width = 5, height = 4, dpi = 400)



  
  
  
  
  

# for Red populations:
phaseIII_table_total_Genetic_phenotype_PCA_PCvalues <- 
  phaseIII_table_total_Genetic_PCA_PCvalues_red %>% 
  merge(table_total_phenotype_PCA_PCvalues, by=c("pop_ID"), all=TRUE) 

# for pink populations:
phaseIII_table_total_Genetic_phenotype_PCA_PCvalues <- 
  phaseIII_table_total_Genetic_PCA_PCvalues_pink %>% 
  merge(table_total_phenotype_PCA_PCvalues, by=c("pop_ID"), all=TRUE) 


mod1 <- glm(pc1_phenotype ~ pc1, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues, 
           family=gaussian())

mod2 <- glm(pc2_phenotype ~ pc2, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues, 
           family=gaussian())

summary(mod1)
summary(mod2)

library(nlme)
library(geiger)
fit<-lm(pc1_phenotype ~ pc1, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues)
summary(fit)

fit<-lm(pc2_phenotype ~ pc2, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues)
summary(fit)

# for Para populations:
fit<-lm(pc1_phenotype ~ pc1 + treatment, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues)
summary(fit)



cor(phaseIII_table_total_Genetic_phenotype_PCA_PCvalues$pc1, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues$pc1_phenotype)
cor(phaseIII_table_total_Genetic_phenotype_PCA_PCvalues$pc2, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues$pc2_phenotype)


phaseIII_table_total_Genetic_phenotype_PCA_PCvalues %>% 
  # filter(pc1<0.4) %>%
  ggplot(aes(pc1, pc1_phenotype)) +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  geom_point(size=2, alpha=0.8) +
  ylab(label = "PC1 Phenotype Var.") +
  xlab(label = "PC1 Genotype Var.") + 
  # ylim(c(0,5)) +
  # stat_poly_eq(formula = my.formula,
  #              eq.with.lhs = "italic(hat(y))~`=`~",
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_blank(), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+ 
  ggsave(paste0("new_analyses_plots/11_PCA_Genotype_Phenotype_", treatment, "_", Ade_used, "_PC1.png"), width = 5, height = 4, dpi = 400)


phaseIII_table_total_Genetic_phenotype_PCA_PCvalues %>% 
  ggplot(aes(pc2*(-1), pc2_phenotype)) +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  geom_point(size=2, alpha=0.8) +
  ylab(label = "PC2 Phenotype Var.") +
  xlab(label = "PC2 Genotype Var.") + 
  # ylim(c(0,5)) +
  # stat_poly_eq(formula = my.formula,
  #              eq.with.lhs = "italic(hat(y))~`=`~",
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_blank(), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
  # ggsave(paste0("new_analyses_plots/11_PCA_Genotype_Phenotype_", treatment, "_", Ade_used, "_PC2.png"), width = 5, height = 4, dpi = 400)

####

### script for parapatric populations:


library("RColorBrewer")
library(ggfortify)
library(factoextra)

# Genotype tables:

# Allopatric treatments:
# red populations:
phaseIII_table_total_PCA_red<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  filter(treatment %in% c("Para_T", "Para_B")) %>%
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))

plot_PCA_red<-autoplot(PCA_data, 
                       data=phaseIII_table_total_PCA_red, 
                       label = F , 
                       colour = 'treatment', 
                       fill =  'treatment', 
                       shape=T, label.size = 5, 
                       frame=T)


phaseIII_table_total_Genetic_PCA_PCvalues_red <- data.frame(treatment=plot_PCA_red$data$treatment, time=plot_PCA_red$data$time,  pop_ID=plot_PCA_red$data$pop_ID, pc1=as.numeric(plot_PCA_red$data$PC1), pc2=as.numeric(plot_PCA_red$data$PC2)) 

# pink populations:
phaseIII_table_total_PCA_pink<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P", "SYM", "LM", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  filter(treatment %in% c("Para_T", "Para_B")) %>%
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("pink3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_pink %>% select(-treatment, -time, -pop_ID))


plot_PCA_pink<-autoplot(PCA_data, 
                        data=phaseIII_table_total_PCA_pink, 
                        label = F , 
                        colour = 'treatment', 
                        fill =  'treatment', 
                        shape=T, label.size = 5, 
                        frame=T)

plot_PCA_pink

phaseIII_table_total_Genetic_PCA_PCvalues_pink <- data.frame(treatment=plot_PCA_pink$data$treatment, time=plot_PCA_pink$data$time,  pop_ID=plot_PCA_pink$data$pop_ID, pc1=as.numeric(plot_PCA_pink$data$PC1), pc2=as.numeric(plot_PCA_pink$data$PC2)) 

# phenotype table:
Ade_used<-"R"
treatment<-"Para" # "Allo","Para", "LM", "Sym"
# for (treatment in c("Allo","Para", "LM", "Sym")) {
sim_phenotype <- sim_phenotype_completed %>%
  # filter(Ade==Ade_used, Group %in% c("Allo", "Para"))
  filter(Ade==Ade_used, Group==treatment)
# sim_phenotype <- sim_phenotype_completed %>% 
#   filter(Ade==Ade_used, Group!="Parental")
N_v<-c((sim_phenotype$log_median_relative_v-mean(sim_phenotype$log_median_relative_v))/sd(sim_phenotype$log_median_relative_v))
N_m<-c((sim_phenotype$log_median_m.evol-mean(sim_phenotype$log_median_m.evol))/sd(sim_phenotype$log_median_m.evol))
N_WBS<-c((sim_phenotype$log_relativeWBS-mean(sim_phenotype$log_relativeWBS))/sd(sim_phenotype$log_relativeWBS))
N_WTS<-c((sim_phenotype$log_relativeWTS-mean(sim_phenotype$log_relativeWTS))/sd(sim_phenotype$log_relativeWTS))

sim_phenotype_PCAtable<-sim_phenotype %>% 
  cbind(N_v, N_m, N_WBS, N_WTS) %>% 
  #filter(Group=="Allo") %>% 
  select(-dif_median_relative_v, -log_median_m.evol, -log_median_relative_v, 
         -log_relativeWBS, -log_relativeWTS) 

PCA_data_phenotype<-prcomp(sim_phenotype_PCAtable %>% select(N_v, N_m, N_WBS, N_WTS))

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")
plot_PCA_phenotype<-autoplot(PCA_data_phenotype, 
                             data=sim_phenotype_PCAtable, 
                             label = F , 
                             colour = 'evolved', 
                             fill =  'evolved', 
                             shape=T, label.size = 5, 
                             frame=T) 

plot_PCA_phenotype

table_total_phenotype_PCA_PCvalues<-data.frame(pop_ID=plot_PCA_phenotype$data$pop_ID, 
                                               pc1_phenotype=as.numeric(plot_PCA_phenotype$data$PC1), 
                                               pc2_phenotype=as.numeric(plot_PCA_phenotype$data$PC2))


# for Red populations:
phaseIII_table_total_Genetic_phenotype_PCA_PCvalues <- 
  phaseIII_table_total_Genetic_PCA_PCvalues_red %>% 
  merge(table_total_phenotype_PCA_PCvalues, by=c("pop_ID"), all=TRUE) 

# for pink populations:
phaseIII_table_total_Genetic_phenotype_PCA_PCvalues <- 
  phaseIII_table_total_Genetic_PCA_PCvalues_pink %>% 
  merge(table_total_phenotype_PCA_PCvalues, by=c("pop_ID"), all=TRUE) 

phaseIII_table_total_Genetic_phenotype_PCA_PCvalues

mod1 <- glm(pc1_phenotype ~ pc1, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues, 
            family=gaussian())

mod2 <- glm(pc2_phenotype ~ pc2, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues, 
            family=gaussian())

summary(mod1)
summary(mod2)

library(nlme)
library(geiger)
head(phaseIII_table_total_Genetic_phenotype_PCA_PCvalues)

fit<-lm(pc1_phenotype ~ pc1 + treatment, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues)
summary(fit)
fit<-lm(pc2_phenotype ~ pc2 + treatment, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues)
summary(fit)

fit<-lm(pc1_phenotype ~ pc1, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues)
summary(fit)

fit<-lm(pc2_phenotype ~ pc2, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues)
summary(fit)




### phenotype with other treatments:

sim_growth<- median_inference_m.evol %>% 
  mutate(pop_ID=paste(Pop, Group, evolved, sep="_")) %>% 
  ungroup() %>% 
  mutate(log_median_m.evol=log(median_m.evol)) %>% 
  select(pop_ID, Ade, log_median_m.evol) 


sim_ecological_sel<-
  ecological_sel %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo)) %>% 
  ungroup() %>% 
  select(-sd_relative_W_evo) %>% 
  spread(Time, mean_relative_W_evo) %>% 
  mutate(pop_ID=paste(Pop, Group, evolved, sep="_")) %>% 
  ungroup() %>% 
  mutate(relativeWBS=BS, relativeWTS=TS) %>% 
  mutate(log_relativeWBS=log(relativeWBS), log_relativeWTS=log(relativeWTS)) %>% 
  select(pop_ID, Ade, log_relativeWBS, log_relativeWTS)


sim_phenotype<-
  mating_data_ed %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(log_median_relative_v = log(median(relative_v))) %>% 
  ungroup() %>% 
  mutate(pop_ID=paste(Pop, Group, evolved, sep="_")) %>% 
  merge(sim_growth, by=c("pop_ID", "Ade"), all=TRUE) %>% 
  merge(sim_ecological_sel, by=c("pop_ID", "Ade"), all=TRUE) %>% 
  filter(Ade=="R") %>%
  # filter(Group %in% c("Allo")) %>%
  # filter(Group %in% c("Allo", "Para")) %>%
  # filter(Group %in% c("LM", "Sym")) %>%
  # select(Structure, Fraction, Pop, log_median_relative_v, log_median_m.evol, log_relativeWBS, log_relativeWTS)
  select(Structure, Fraction, Pop, log_median_m.evol, log_median_relative_v, log_relativeWBS, log_relativeWTS) %>% 
  mutate(treatment=paste(Structure, Fraction, sep="_")) %>% 
  mutate(pop_ID=paste(Pop, Structure, Fraction, sep="_"))

head(sim_phenotype)

N_median_relative_v<-c((sim_phenotype$log_median_relative_v-mean(sim_phenotype$log_median_relative_v))/sd(sim_phenotype$log_median_relative_v))
N_median_m.evol<-c((sim_phenotype$log_median_m.evol-mean(sim_phenotype$log_median_m.evol))/sd(sim_phenotype$log_median_m.evol))
N_relativeWBS<-c((sim_phenotype$log_relativeWBS-mean(sim_phenotype$log_relativeWBS))/sd(sim_phenotype$log_relativeWBS))
N_relativeWTS<-c((sim_phenotype$log_relativeWTS-mean(sim_phenotype$log_relativeWTS))/sd(sim_phenotype$log_relativeWTS))


N_median_relative_v<-ppls::normalize.vector(sim_phenotype$log_median_relative_v)
N_median_m.evol<-ppls::normalize.vector(sim_phenotype$log_median_m.evol)
N_relativeWBS<-ppls::normalize.vector(sim_phenotype$log_relativeWBS)
N_relativeWTS<-ppls::normalize.vector(sim_phenotype$log_relativeWTS)


# sim_phenotype_PCAtable<-sim_phenotype %>%
#   cbind(N_median_relative_v, N_median_m.evol, N_relativeWBS, N_relativeWTS) %>%
#   select(-median_relative_v, -median_m.evol, -relativeWBS, -relativeWTS)

sim_phenotype_PCAtable<-
sim_phenotype %>%
  cbind(N_median_m.evol, N_median_relative_v, N_relativeWBS, N_relativeWTS) %>%
  select(-log_median_m.evol, -log_median_relative_v, -log_relativeWBS, -log_relativeWTS) %>% 
  filter(Structure=="Para")
  

PCA_data_phenotype<-prcomp(sim_phenotype_PCAtable %>% select(-Structure, -Fraction, -Pop, -pop_ID, -treatment))

library("RColorBrewer")
library(ggfortify)
library(factoextra)

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")
autoplot(PCA_data_phenotype, 
         data=sim_phenotype_PCAtable, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) +
  scale_colour_manual(values=colours_PCA) +
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
  ggsave("plots/08_PCA_allo_red.png", width = 7, height = 5, dpi = 400)


head(PCA_data_phenotype)
fviz_pca_var(PCA_data_phenotype,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = T ,
             repel = F ) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_blank(), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
  # ggsave("plots/08_PCA_loads_allo_red_nl.png", width = 4, height = 3, dpi = 400)

PCA_point<-
autoplot(PCA_data_phenotype, 
         data=sim_phenotype_PCAtable, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=F) +
  scale_colour_manual(values=colours_PCA) 

pPCA_point <- ggplot_build(PCA_point)

PCA_point_data<-pPCA_point$data[[1]] %>% 
  mutate(PC1=as.numeric(x), PC2=as.numeric(y))
  
fviz_pca_var(PCA_data_phenotype,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label = F ,
             repel = T ) +
  geom_point(aes(PC1,PC2), colour=as.factor(PCA_point_data$colour), data=PCA_point_data, size=2, alpha=0.6) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_blank(), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
  ggsave("plots/08_PCA_loads_Sym_red_nl_points.png", width = 4, height = 3, dpi = 400)




phenotype_plot_PCA_red<-autoplot(PCA_data_phenotype, 
                                 data=sim_phenotype_PCAtable, 
                                 label = F , 
                                 colour = 'treatment', 
                                 fill =  'treatment', 
                                 shape=T, label.size = 5, 
                                 frame=T) 



phaseIII_table_total_Phenotype_PCA_PCvalues_red <- 
  data.frame(treatment=phenotype_plot_PCA_red$data$treatment, 
             pop_ID=phenotype_plot_PCA_red$data$pop_ID, 
             pc1_phenotype=as.numeric(phenotype_plot_PCA_red$data$PC1), 
             pc2_phenotype=as.numeric(phenotype_plot_PCA_red$data$PC2)) 



head(phaseIII_table_total_Genetic_PCA_PCvalues_red)
head(phaseIII_table_total_Phenotype_PCA_PCvalues_red)


phaseIII_table_total_Genetic_phenotype_PCA_PCvalues_red <- 
  phaseIII_table_total_Genetic_PCA_PCvalues_red %>% 
  merge(phaseIII_table_total_Phenotype_PCA_PCvalues_red, by=c("pop_ID", "treatment"), all=TRUE) 

mod <- glm(pc2_phenotype ~ pc1+pc2, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues_red, family=gaussian())
summury_mod<-summary(mod)
summury_mod


# allopatric pink

sim_phenotype<-
  mating_data_ed %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(log_median_relative_v = log(median(relative_v))) %>% 
  ungroup() %>% 
  mutate(pop_ID=paste(Pop, Group, evolved, sep="_")) %>% 
  merge(sim_growth, by=c("pop_ID", "Ade"), all=TRUE) %>% 
  merge(sim_ecological_sel, by=c("pop_ID", "Ade"), all=TRUE) %>% 
  filter(Ade=="P") %>% 
  filter(Group %in% c("Allo")) %>%
  # filter(Group %in% c("Allo", "Para")) %>%
  # filter(Group %in% c("LM", "Sym")) %>%
  # select(Structure, Fraction, Pop, log_median_relative_v, log_median_m.evol, log_relativeWBS, log_relativeWTS)
  select(Structure, Fraction, Pop, log_median_m.evol, log_median_relative_v, log_relativeWBS, log_relativeWTS) %>% 
  mutate(treatment=paste(Structure, Fraction, sep="_")) %>% 
  mutate(pop_ID=paste(Pop, Structure, Fraction, sep="_"))



N_median_relative_v<-c((sim_phenotype$log_median_relative_v-mean(sim_phenotype$log_median_relative_v))/sd(sim_phenotype$log_median_relative_v))
N_median_m.evol<-c((sim_phenotype$log_median_m.evol-mean(sim_phenotype$log_median_m.evol))/sd(sim_phenotype$log_median_m.evol))
N_relativeWBS<-c((sim_phenotype$log_relativeWBS-mean(sim_phenotype$log_relativeWBS))/sd(sim_phenotype$log_relativeWBS))
N_relativeWTS<-c((sim_phenotype$log_relativeWTS-mean(sim_phenotype$log_relativeWTS))/sd(sim_phenotype$log_relativeWTS))


# sim_phenotype_PCAtable<-sim_phenotype %>%
#   cbind(N_median_relative_v, N_median_m.evol, N_relativeWBS, N_relativeWTS) %>%
#   select(-median_relative_v, -median_m.evol, -relativeWBS, -relativeWTS)

sim_phenotype_PCAtable<-sim_phenotype %>%
  cbind(N_median_m.evol, N_median_relative_v, N_relativeWBS, N_relativeWTS) %>%
  select(-log_median_m.evol, -log_median_relative_v, -log_relativeWBS, -log_relativeWTS)

PCA_data_phenotype<-prcomp(sim_phenotype_PCAtable %>% select(-Structure, -Fraction, -Pop, -pop_ID, -treatment))

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")
autoplot(PCA_data_phenotype, 
         data=sim_phenotype_PCAtable, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) +
  scale_colour_manual(values=colours_PCA) +
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) +
  ggsave("plots/08_PCA_allo_pink.png", width = 7, height = 5, dpi = 400)



fviz_pca_var(PCA_data_phenotype,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = F , 
             repel = T ) +
  ggsave("plots/08_PCA_loads_allo_pink.png", width = 7, height = 5, dpi = 400)


phenotype_plot_PCA_pink<-autoplot(PCA_data_phenotype, 
                                  data=sim_phenotype_PCAtable, 
                                  label = F , 
                                  colour = 'treatment', 
                                  fill =  'treatment', 
                                  shape=T, label.size = 5, 
                                  frame=T) 

phaseIII_table_total_Phenotype_PCA_PCvalues_pink <- 
  data.frame(treatment=phenotype_plot_PCA_pink$data$treatment, 
             pop_ID=phenotype_plot_PCA_pink$data$pop_ID, 
             pc1_phenotype=as.numeric(phenotype_plot_PCA_pink$data$PC1), 
             pc2_phenotype=as.numeric(phenotype_plot_PCA_pink$data$PC2)) 



head(phaseIII_table_total_Genetic_PCA_PCvalues_pink)
head(phaseIII_table_total_Phenotype_PCA_PCvalues_pink)


phaseIII_table_total_Genetic_phenotype_PCA_PCvalues_pink <- 
  phaseIII_table_total_Genetic_PCA_PCvalues_pink %>% 
  merge(phaseIII_table_total_Phenotype_PCA_PCvalues_pink, by=c("pop_ID", "treatment"), all=TRUE) 

mod <- glm(pc2_phenotype ~ pc1+pc2, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues_pink, family=gaussian())
summury_mod<-summary(mod)
summury_mod



# Other treatments:
# red populations:

phaseIII_table_total_PCA_red<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>%
  filter(treatment %in% c("Allo_T", "Allo_B", "SYM", "LM", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))


plot_PCA_red<-autoplot(PCA_data, 
                       data=phaseIII_table_total_PCA_red, 
                       label = F , 
                       colour = 'treatment', 
                       fill =  'treatment', 
                       shape=T, label.size = 5, 
                       frame=T)

phaseIII_table_total_Genetic_PCA_PCvalues_red <- data.frame(treatment=plot_PCA_red$data$treatment, time=plot_PCA_red$data$time,  pop_ID=plot_PCA_red$data$pop_ID, pc1=as.numeric(plot_PCA_red$data$PC1), pc2=as.numeric(plot_PCA_red$data$PC2)) 

# pink populations:
phaseIII_table_total_PCA_pink<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(treatment %in% c("Allo_T", "Allo_B", "SYM", "LM", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P", "SYM", "LM", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "Para_T", "Para_B")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>%
  # filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("pink3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_pink %>% select(-treatment, -time, -pop_ID))


plot_PCA_pink<-autoplot(PCA_data, 
                        data=phaseIII_table_total_PCA_pink, 
                        label = F , 
                        colour = 'treatment', 
                        fill =  'treatment', 
                        shape=T, label.size = 5, 
                        frame=T)

phaseIII_table_total_Genetic_PCA_PCvalues_pink <- data.frame(treatment=plot_PCA_pink$data$treatment, time=plot_PCA_pink$data$time,  pop_ID=plot_PCA_pink$data$pop_ID, pc1=as.numeric(plot_PCA_pink$data$PC1), pc2=as.numeric(plot_PCA_pink$data$PC2)) 




### phenotype with other treatments:

sim_growth<- median_inference_m.evol %>% 
  mutate(pop_ID=paste(Pop, Group, evolved, sep="_")) %>% 
  ungroup() %>% 
  mutate(log_median_m.evol=log(median_m.evol)) %>% 
  select(pop_ID, Ade, log_median_m.evol) 

sim_ecological_sel<-
  ecological_sel %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo)) %>% 
  ungroup() %>% 
  select(-sd_relative_W_evo) %>% 
  spread(Time, mean_relative_W_evo) %>% 
  mutate(pop_ID=paste(Pop, Group, evolved, sep="_")) %>% 
  ungroup() %>% 
  mutate(relativeWBS=BS, relativeWTS=TS) %>% 
  mutate(log_relativeWBS=log(relativeWBS), log_relativeWTS=log(relativeWTS)) %>% 
  select(pop_ID, Ade, log_relativeWBS, log_relativeWTS)

sim_phenotype<-
  mating_data_ed %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(log_median_relative_v = log(median(relative_v))) %>% 
  ungroup() %>% 
  mutate(pop_ID=paste(Pop, Group, evolved, sep="_")) %>% 
  # merge(sim_growth, by=c("pop_ID", "Ade"), all=TRUE) %>%
  merge(sim_ecological_sel, by=c("pop_ID", "Ade"), all=TRUE) %>% 
  filter(Ade=="R") %>% 
  # filter(Group %in% c("Allo")) %>%
  # filter(Group %in% c("Allo", "Para")) %>%
  # filter(Group %in% c("LM", "Sym")) %>%
  # select(Structure, Fraction, Pop, log_median_relative_v, log_median_m.evol, log_relativeWBS, log_relativeWTS)
  select(Structure, Fraction, Pop, log_median_relative_v, log_relativeWBS, log_relativeWTS) %>% 
  mutate(treatment=paste(Structure, Fraction, sep="_")) %>% 
  mutate(pop_ID=paste(Pop, Structure, Fraction, sep="_"))


N_median_relative_v<-c((sim_phenotype$log_median_relative_v-mean(sim_phenotype$log_median_relative_v))/sd(sim_phenotype$log_median_relative_v))
# N_median_m.evol<-c((sim_phenotype$log_median_m.evol-mean(sim_phenotype$log_median_m.evol))/sd(sim_phenotype$log_median_m.evol))
N_relativeWBS<-c((sim_phenotype$log_relativeWBS-mean(sim_phenotype$log_relativeWBS))/sd(sim_phenotype$log_relativeWBS))
N_relativeWTS<-c((sim_phenotype$log_relativeWTS-mean(sim_phenotype$log_relativeWTS))/sd(sim_phenotype$log_relativeWTS))


# sim_phenotype_PCAtable<-sim_phenotype %>% 
#   cbind(N_median_relative_v, N_median_m.evol, N_relativeWBS, N_relativeWTS) %>% 
#   select(-median_relative_v, -median_m.evol, -relativeWBS, -relativeWTS) 

sim_phenotype_PCAtable<-sim_phenotype %>%
  cbind(N_median_relative_v, N_relativeWBS, N_relativeWTS) %>%
  select(-log_median_relative_v, -log_relativeWBS, -log_relativeWTS)


PCA_data_phenotype<-prcomp(sim_phenotype_PCAtable %>% select(-Structure, -Fraction, -Pop, -pop_ID, -treatment))

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")
autoplot(PCA_data_phenotype, 
         data=sim_phenotype_PCAtable, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) +
  scale_colour_manual(values=colours_PCA) +
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
ggsave("plots/08_PCA_para_red.png", width = 7, height = 5, dpi = 400)



fviz_pca_var(PCA_data_phenotype,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = F , 
             repel = T ) +
  ggsave("plots/08_PCA_loads_para_red.png", width = 7, height = 5, dpi = 400)


phenotype_plot_PCA_red<-autoplot(PCA_data_phenotype, 
                                 data=sim_phenotype_PCAtable, 
                                 label = F , 
                                 colour = 'treatment', 
                                 fill =  'treatment', 
                                 shape=T, label.size = 5, 
                                 frame=T) 

phaseIII_table_total_Phenotype_PCA_PCvalues_red <- 
  data.frame(treatment=phenotype_plot_PCA_red$data$treatment, 
             pop_ID=phenotype_plot_PCA_red$data$pop_ID, 
             pc1_phenotype=as.numeric(phenotype_plot_PCA_red$data$PC1), 
             pc2_phenotype=as.numeric(phenotype_plot_PCA_red$data$PC2)) 



head(phaseIII_table_total_Genetic_PCA_PCvalues_red)
head(phaseIII_table_total_Phenotype_PCA_PCvalues_red)


phaseIII_table_total_Genetic_phenotype_PCA_PCvalues_red <- 
  phaseIII_table_total_Genetic_PCA_PCvalues_red %>% 
  merge(phaseIII_table_total_Phenotype_PCA_PCvalues_red, by=c("pop_ID", "treatment"), all=TRUE) 

mod <- glm(pc2_phenotype ~ pc1+pc2, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues_red, family=gaussian())
summury_mod<-summary(mod)
summury_mod


# pink

sim_phenotype<-
  mating_data_ed %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure) %>% 
  summarise(log_median_relative_v = log(median(relative_v))) %>% 
  ungroup() %>% 
  mutate(pop_ID=paste(Pop, Group, evolved, sep="_")) %>% 
  # merge(sim_growth, by=c("pop_ID", "Ade"), all=TRUE) %>%
  merge(sim_ecological_sel, by=c("pop_ID", "Ade"), all=TRUE) %>% 
  filter(Ade=="P") %>% 
  # filter(Group %in% c("Allo")) %>%
  filter(Group %in% c("Allo", "Para")) %>%
  # filter(Group %in% c("LM", "Sym")) %>%
  # select(Structure, Fraction, Pop, log_median_relative_v, log_median_m.evol, log_relativeWBS, log_relativeWTS)
  select(Structure, Fraction, Pop, log_median_relative_v, log_relativeWBS, log_relativeWTS) %>% 
  mutate(treatment=paste(Structure, Fraction, sep="_")) %>% 
  mutate(pop_ID=paste(Pop, Structure, Fraction, sep="_"))

N_median_relative_v<-c((sim_phenotype$log_median_relative_v-mean(sim_phenotype$log_median_relative_v))/sd(sim_phenotype$log_median_relative_v))
# N_median_m.evol<-c((sim_phenotype$log_median_m.evol-mean(sim_phenotype$log_median_m.evol))/sd(sim_phenotype$log_median_m.evol))
N_relativeWBS<-c((sim_phenotype$log_relativeWBS-mean(sim_phenotype$log_relativeWBS))/sd(sim_phenotype$log_relativeWBS))
N_relativeWTS<-c((sim_phenotype$log_relativeWTS-mean(sim_phenotype$log_relativeWTS))/sd(sim_phenotype$log_relativeWTS))


# sim_phenotype_PCAtable<-sim_phenotype %>% 
#   cbind(N_median_relative_v, N_median_m.evol, N_relativeWBS, N_relativeWTS) %>% 
#   select(-median_relative_v, -median_m.evol, -relativeWBS, -relativeWTS) 

sim_phenotype_PCAtable<-sim_phenotype %>%
  cbind(N_median_relative_v, N_relativeWBS, N_relativeWTS) %>%
  select(-log_median_relative_v, -log_relativeWBS, -log_relativeWTS)


PCA_data_phenotype<-prcomp(sim_phenotype_PCAtable %>% select(-Structure, -Fraction, -Pop, -pop_ID, -treatment))

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")
autoplot(PCA_data_phenotype, 
         data=sim_phenotype_PCAtable, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) +
  scale_colour_manual(values=colours_PCA) +
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) +
  ggsave("plots/08_PCA_para_pink.png", width = 7, height = 5, dpi = 400)



fviz_pca_var(PCA_data_phenotype,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = F , 
             repel = T ) +
  ggsave("plots/08_PCA_loads_para_pink.png", width = 7, height = 5, dpi = 400)


phenotype_plot_PCA_pink<-autoplot(PCA_data_phenotype, 
                                  data=sim_phenotype_PCAtable, 
                                  label = F , 
                                  colour = 'treatment', 
                                  fill =  'treatment', 
                                  shape=T, label.size = 5, 
                                  frame=T) 

phaseIII_table_total_Phenotype_PCA_PCvalues_pink <- 
  data.frame(treatment=phenotype_plot_PCA_pink$data$treatment, 
             pop_ID=phenotype_plot_PCA_pink$data$pop_ID, 
             pc1_phenotype=as.numeric(phenotype_plot_PCA_pink$data$PC1), 
             pc2_phenotype=as.numeric(phenotype_plot_PCA_pink$data$PC2)) 



head(phaseIII_table_total_Genetic_PCA_PCvalues_pink)
head(phaseIII_table_total_Phenotype_PCA_PCvalues_pink)


phaseIII_table_total_Genetic_phenotype_PCA_PCvalues_pink <- 
  phaseIII_table_total_Genetic_PCA_PCvalues_pink %>% 
  merge(phaseIII_table_total_Phenotype_PCA_PCvalues_pink, by=c("pop_ID", "treatment"), all=TRUE) 

mod <- glm(pc2_phenotype ~ pc1+pc2, phaseIII_table_total_Genetic_phenotype_PCA_PCvalues_pink, family=gaussian())
summury_mod<-summary(mod)
summury_mod




library("lme4")
library("languageR")
library("afex")
library("sjstats")

# GLM for difference in growth rate:
growth_table<-inference_m.evol %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
  dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>%
  filter(Group_sample!="Parental") %>%
  filter(Ade=="P") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  merge(inference_m.evol %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          dplyr::summarize(N_rep=n(), median_m.evol=median(m.evol), mean_m.evol=mean(m.evol)) %>% 
          #filter(Group_sample!="Parental") %>%
          #filter(Ade=="R") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade) %>% 
          summarise(parental_median_m=median(median_m.evol), 
                    parental_mean_m=mean(mean_m.evol)), 
        by="Ade", all.x = T) %>% 
  # mutate(dif_median_m.evol=median_m.evol-parental_median_m, 
  #        dif_mean_m.evol=mean_m.evol-parental_mean_m) %>% 
  mutate(dif_median_m.evol=exp(log(median_m.evol)-log(parental_median_m))) 

head(growth_table)
mod_growth<-lmer(dif_median_m.evol ~ 1 + Structure/Fraction + (1|Pop), data=growth_table)
res_mod<-summary(mod_growth)
#icc(modtest2)
growth_table$dif_median_m.evol_pre<-predict(mod_growth)

growth_table %>% 
  mutate(Structure=factor(Structure, levels=c("Allo", "Para", "LM", "Sym"))) %>% 
  ggplot(aes(Structure, (dif_median_m.evol_pre), colour=Fraction)) +
  geom_point()

# write.table(res_mod$coefficients, "model_test.growth_red.txt", sep = "\t", quote = F)
# write.table(res_mod$coefficients, "model_test.growth_pink.txt", sep = "\t", quote = F)
res_mod
######


# GLM for difference in mating viability:

mating_table <-mating_data_ed %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
  summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
            sd_relative_v = sd(relative_v), 
            median_relative_v = median(relative_v)) %>% 
  filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  merge(mating_data_ed %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>% 
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>% 
          summarise(N_replicates=n(), mean_relative_v = mean(relative_v), 
                    sd_relative_v = sd(relative_v), 
                    median_relative_v = median(relative_v)) %>% 
          #filter(Group_sample!="Parental") %>% 
          ungroup() %>% 
          filter(Group_sample=="Parental") %>% 
          #mutate(Group=factor(Group, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          group_by(Ade) %>% 
          summarise(parental_median_relative_v=median(median_relative_v), 
                    parental_mean_relative_v=mean(median_relative_v)), 
        by="Ade", all.x = T) %>% 
  mutate(dif_median_relative_v=log(median_relative_v/parental_median_relative_v), 
         dif_mean_relative_v=mean_relative_v/parental_mean_relative_v) %>% 
  filter(Ade=="P") 

head(mating_table)
mod_mating<-lmer(dif_median_relative_v ~ 1 + Structure/Fraction + (1|Pop), data=mating_table)
res_mod<-summary(mod_mating)
#icc(modtest2)
mating_table$log_dif_median_relative_v_pre<-predict(mod_mating)

mating_table %>% 
  mutate(Structure=factor(Structure, levels=c("Allo", "Para", "LM", "Sym"))) %>% 
  ggplot(aes(Structure, exp(log_dif_median_relative_v_pre), colour=Fraction)) +
  geom_point()

# write.table(res_mod$coefficients, "model_test.mating_red.txt", sep = "\t", quote = F)
# write.table(res_mod$coefficients, "model_test.mating_pink.txt", sep = "\t", quote = F)
res_mod

#####


# GLM for difference in ecological fitness:

eco_table <- ecological_sel %>% 
  mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
  filter(Ade=="P") %>%
  group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
  summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
            median_relative_W_evo = median(relative_W_evo)) %>%
  #filter(Group_sample!="Parental") %>%
  ungroup() %>% 
  mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
  filter(Group!="Parental") %>% 
  merge(ecological_sel %>% 
          mutate(Group_sample=ifelse(Group %in% c("P3P", "P3R"), "Parental",as.character(evolved))) %>%
          #filter(Ade=="R") %>%
          group_by( Fraction, Pop, sample, Group, Ade, evolved, Time, Structure, Group_sample) %>%
          summarise(mean_relative_W_evo = mean(relative_W_evo), sd_relative_W_evo = sd(relative_W_evo), 
                    median_relative_W_evo = median(relative_W_evo)) %>%
          #filter(Group_sample!="Parental") %>%
          ungroup() %>% 
          mutate(Group=factor(ifelse(Group_sample!="Parental",as.vector(Group),"Parental"), 
                              levels = c("Parental","Allo", "Para", "LM", "Sym"))) %>% 
          filter(Group=="Parental") %>% 
          group_by(Ade, Time) %>% 
          summarise(parental_median_relative_W_evo=median(median_relative_W_evo), 
                    parental_mean_relative_W_evo=mean(median_relative_W_evo)), 
        by=c("Ade","Time"), all.x = T) %>% 
  # mutate(dif_median_relative_W.evol=median_relative_W_evo-parental_median_relative_W_evo, 
  #        dif_mean_relative_W.evol=mean_relative_W_evo-parental_mean_relative_W_evo) %>% 
  # using the ratio intead?
  mutate(dif_median_relative_W.evol=log(median_relative_W_evo/parental_median_relative_W_evo)) %>% 
  mutate(Time=factor(Time, levels = c("TS", "BS"))) %>% 
  filter(Time=="TS")

head(eco_table)
mod_eco<-lmer(dif_median_relative_W.evol ~ 1 + Structure/Fraction + (1|Pop), data=eco_table)
res_mod<-summary(mod_eco)
#icc(modtest2)
eco_table$dif_median_relative_W.evol_pre<-predict(mod_mating)

eco_table %>% 
  mutate(Structure=factor(Structure, levels=c("Allo", "Para", "LM", "Sym"))) %>% 
  ggplot(aes(Structure, (dif_median_relative_W.evol_pre), colour=Fraction)) +
  geom_point()

# write.table(res_mod$coefficients, "model_test.eco_TS_red.txt", sep = "\t", quote = F)
# write.table(res_mod$coefficients, "model_test.eco_BS_red.txt", sep = "\t", quote = F)

# write.table(res_mod$coefficients, "model_test.eco_TS_pink.txt", sep = "\t", quote = F)
# write.table(res_mod$coefficients, "model_test.eco_BS_pink.txt", sep = "\t", quote = F)
res_mod

#####



#### Correlations between phenotypes:

sim_phenotype_PCAtable %>% 
  ungroup() %>% 
  mutate(Structure=factor(Structure, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  # filter(Structure=="Allo") %>%
  ggplot(aes((log_median_m.evol), (log_median_relative_v), colour=Ade)) + #, colour=Structure)) + 
  geom_point(alpha=0.6) +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  # scale_colour_manual(values=colours_top_bottom) +
  xlab(label = "Log (Growth rate (m))") +
  ylab(label = "Log (Relative v)") + 
  # ylim(c(0,5)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  facet_grid(Structure ~ .) +
  theme_classic() + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        # strip.placement = "outside",
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  # ggsave(paste0("plots/07_Growth_vs_mating_bytreatment.png"), width = 4, height = 8, dpi = 400)
  ggsave(paste0("plots/07_Growth_vs_mating.png"), width = 5, height = 4, dpi = 400)



sim_phenotype_PCAtable %>% 
  ungroup() %>% 
  mutate(Structure=factor(Structure, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  # filter(Structure=="Allo") %>%
  ggplot(aes((log_relativeWTS), (log_median_relative_v), colour=Ade)) + #, colour=Structure)) + 
  geom_point(alpha=0.6) +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  # scale_colour_manual(values=colours_top_bottom) +
  xlab(label = "Log (Relative W Top Sel.)") +
  ylab(label = "Log (Relative v)") + 
  # ylim(c(0,5)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  facet_grid(Structure ~ .) +
  theme_classic() + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        # strip.placement = "outside",
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  ggsave(paste0("plots/07_WTS_vs_mating_bytreatment.png"), width = 4, height = 8, dpi = 400)
# ggsave(paste0("plots/07_WTS_vs_mating.png"), width = 5, height = 4, dpi = 400)



sim_phenotype_PCAtable %>% 
  ungroup() %>% 
  mutate(Structure=factor(Structure, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  # filter(Structure=="Allo") %>%
  ggplot(aes((log_relativeWBS), (log_median_relative_v), colour=Ade)) + #, colour=Structure)) + 
  geom_point(alpha=0.6) +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  # scale_colour_manual(values=colours_top_bottom) +
  xlab(label = "Log (Relative W Bottom Sel.)") +
  ylab(label = "Log (Relative v)") + 
  # ylim(c(0,5)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  facet_grid(Structure ~ .) +
  theme_classic() + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        # strip.placement = "outside",
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  # ggsave(paste0("plots/07_WBS_vs_mating_bytreatment.png"), width = 4, height = 8, dpi = 400)
  ggsave(paste0("plots/07_WBS_vs_mating.png"), width = 5, height = 4, dpi = 400)




sim_phenotype_PCAtable %>% 
  ungroup() %>% 
  mutate(Structure=factor(Structure, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  # filter(Structure=="Allo") %>%
  ggplot(aes((log_relativeWTS), (log_relativeWBS), colour=Ade)) + #, colour=Structure)) + 
  geom_point(alpha=0.6) +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  # scale_colour_manual(values=colours_top_bottom) +
  xlab(label = "Log (Relative W Top Sel.)") +
  ylab(label = "Log (Relative W Bottom Sel.)") + 
  # ylim(c(0,5)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  facet_grid(Structure ~ .) +
  theme_classic() + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        # strip.placement = "outside",
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  # ggsave(paste0("plots/07_WTS_vs_WBS_bytreatment.png"), width = 4, height = 8, dpi = 400)
  ggsave(paste0("plots/07_WTS_vs_WBS.png"), width = 5, height = 4, dpi = 400)




sim_phenotype_PCAtable %>% 
  ungroup() %>% 
  mutate(Structure=factor(Structure, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  # filter(Structure=="Allo") %>%
  ggplot(aes((log_median_m.evol), (log_relativeWTS), colour=Ade)) + #, colour=Structure)) + 
  geom_point(alpha=0.6) +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  # scale_colour_manual(values=colours_top_bottom) +
  xlab(label = "Log (Growth rate (m))") +
  ylab(label = "Log (Relative W Top Sel.)") + 
  # ylim(c(0,5)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  # facet_grid(Structure ~ .) +
  theme_classic() + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        # strip.placement = "outside",
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  # ggsave(paste0("plots/07_Growth_vs_WTS_bytreatment.png"), width = 4, height = 8, dpi = 400)
  ggsave(paste0("plots/07_Growth_vs_WTS.png"), width = 5, height = 4, dpi = 400)



sim_phenotype_PCAtable %>% 
  ungroup() %>% 
  mutate(Structure=factor(Structure, levels = c("Allo", "Para", "LM", "Sym"))) %>% 
  # filter(Structure=="Allo") %>%
  ggplot(aes((log_median_m.evol), (log_relativeWBS), colour=Ade)) + #, colour=Structure)) + 
  geom_point(alpha=0.6) +
  geom_smooth(method='lm',formula=y~x, alpha=0.3) +
  # scale_colour_manual(values=colours_top_bottom) +
  xlab(label = "Log (Growth rate (m))") +
  ylab(label = "Log (Relative W Bottom Sel.)") + 
  # ylim(c(0,5)) +
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  facet_grid(Structure ~ .) +
  theme_classic() + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"), 
        axis.line = element_line(colour = "black"), 
        strip.text.y = element_text(angle = 360),
        # legend.position="bottom",
        # legend.position="none",
        # strip.placement = "outside",
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  ggsave(paste0("plots/07_Growth_vs_WBS_bytreatment.png"), width = 4, height = 8, dpi = 400)
# ggsave(paste0("plots/07_Growth_vs_WBS.png"), width = 5, height = 4, dpi = 400) 



mod1<-lm(log_median_relative_v ~ log_median_m.evol * Ade * Structure * Fraction, data=sim_phenotype_PCAtable)
summary(mod1)
mod1<-lm(log_median_relative_v ~ log_relativeWTS * Ade * Structure * Fraction, data=sim_phenotype_PCAtable)
mod1<-lm(log_median_relative_v ~ log_relativeWTS * Ade *  Fraction, data=sim_phenotype_PCAtable %>% filter(Structure=="Allo"))
summary(mod1)
mod1<-lm(log_median_relative_v ~ log_relativeWBS * Ade * Structure * Fraction, data=sim_phenotype_PCAtable)
summary(mod1)
mod1<-lm(log_relativeWTS ~ log_relativeWBS * Ade * Structure * Fraction, data=sim_phenotype_PCAtable)
summary(mod1)
mod1<-lm(log_median_m.evol ~ log_relativeWTS * Ade * Structure * Fraction, data=sim_phenotype_PCAtable)
summary(mod1)
mod1<-lm(log_median_m.evol ~ log_relativeWBS * Ade * Structure * Fraction, data=sim_phenotype_PCAtable)
summary(mod1)


