#=================================================#
# Plots for the publication, after reviewer updates
# Created by EK 2/24/2023
#=================================================#
library(tidyverse)
#### FIGURE 2 DISEASE DYNAMICS ####
#Figure 2A
room_fill_xtra_0.A <- data.frame(t2 = c(13.99, 27.99, 41.99),
                               med_S = c(996, 1988, 2977),
                               med_E = c(0, 0, 0),
                               med_I = c(0, 0, 0),
                               med_R = c(0, 0, 0),
                               LB_S = c(992, 1980, 2964),
                               UB_S = c(998, 1994, 2988),
                               LB_E = c(0, 0, 0),
                               UB_E = c(0, 0, 0),
                               LB_I = c(0, 0, 0),
                               UB_I = c(0, 0, 0),
                               LB_R = c(0, 0, 0),
                               UB_R = c(0, 0, 0))
results <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_0VaccEff.rds")
colors <- c("Susceptible" = "blue", "Exposed" = "mediumorchid1", "Infected" = "red", "Recovered" = "black")
A <- results %>%
  mutate(t2 = round(t, digits = 0),
         tot_S = S1 + S2 + S3 + S4,
         tot_E = E1 + E2 + E3 + E4,
         tot_I = I1 + I2 + I3 + I4,
         tot_R = R1 + R2 + R3 + R4) %>% 
  group_by(t2) %>%
  summarize(med_S = median(tot_S),
            LB_S = quantile(tot_S, probs = 0.05),
            UB_S = quantile(tot_S, probs = 0.95),
            med_E = median(tot_E),
            LB_E = quantile(tot_E, probs = 0.05),
            UB_E = quantile(tot_E, probs = 0.95),
            med_I = median(tot_I),
            LB_I = quantile(tot_I, probs = 0.05),
            UB_I = quantile(tot_I, probs = 0.95),
            med_R = median(tot_R),
            LB_R = quantile(tot_R, probs = 0.05),
            UB_R = quantile(tot_R, probs = 0.95)) %>%
  ungroup() %>%
  rbind(room_fill_xtra_0.A) %>% 
  arrange(t2) %>% 
  ggplot(aes(x = t2))+ 
  geom_vline(xintercept = 14, lty = 2, size = 1.1, col = "gray")+
  geom_text(aes(x = 14, y = 3600, label = "Room 2 Fill", fontface = "italic"), angle=90, vjust=-1)+
  geom_vline(xintercept = 28, lty = 2, size = 1.1, col = "gray")+
  geom_text(aes(x = 28, y = 3600, label = "Room 3 Fill", fontface = "italic"), angle=90, vjust=-1)+
  geom_vline(xintercept = 42, lty = 2, size = 1.1, col = "gray")+
  geom_text(aes(x = 42, y = 3600, label = "Room 4 Fill", fontface = "italic"), angle=90, vjust=-1)+
  geom_ribbon(aes(ymin = LB_S, ymax = UB_S), fill = "blue", alpha = 0.3)+
  geom_line(aes(y = med_S, color = "Susceptible"), size = 1.1)+
  # geom_line(aes(y = LB_S, color = "Susceptible"), lty = 2, alpha = .5)+
  # geom_line(aes(y = UB_S, color = "Susceptible"), lty = 2, alpha = .5)+
  geom_ribbon(aes(ymin = LB_R, ymax = UB_R), fill = "black", alpha = 0.3)+
  geom_line(aes(y = med_R, color = "Recovered"), size = 1.1)+
  # geom_line(aes(y = LB_R, color = "Recovered"), lty = 2, alpha = .5)+
  # geom_line(aes(y = UB_R, color = "Recovered"), lty = 2, alpha = .5)+
  geom_ribbon(aes(ymin = LB_I, ymax = UB_I), fill = "red", alpha = 0.3)+
  geom_line(aes(y = med_I, color = "Infected"), size = 1.1)+
  # geom_line(aes(y = LB_I, color = "Infected"), lty = 2, alpha = .5)+
  # geom_line(aes(y = UB_I, color = "Infected"), lty = 2, alpha = .5)+
  geom_ribbon(aes(ymin = LB_E, ymax = UB_E), fill = "mediumorchid1", alpha = 0.3)+
  geom_line(aes(y = med_E, color = "Exposed"), size = 1.1)+
  # geom_line(aes(y = LB_E, color = "Exposed"), lty = 2, alpha = .5)+
  # geom_line(aes(y = UB_E, color = "Exposed"), lty = 2, alpha = .5)+
  scale_x_continuous(name = "Time(Days)", limits = c(0,91),
                     breaks = c(1, 30, 60, 90))+
  scale_y_continuous(name = "Total Pigs", limits = c(0, 4000))+
  labs(x = "Time(Days)",
       y = "Total Pigs",
       color = "")+
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(#legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12),
        # legend.title.align = 0.5,
        legend.position = "bottom",
        strip.background = element_blank(),
        # strip.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12, color = "black"))
#Figure 2B
room_fill_xtra_0.B <- data.frame(t2 = c(13.99, 27.99, 41.99),
                                 med_S = c(996, 1987, 2976),
                                 med_E = c(0, 0, 0),
                                 med_I = c(0, 0, 0),
                                 med_R = c(0, 0, 0),
                                 LB_S = c(993, 1976, 2958),
                                 UB_S = c(999, 1994, 2991),
                                 LB_E = c(0, 0, 0),
                                 UB_E = c(0, 0, 0),
                                 LB_I = c(0, 0, 0),
                                 UB_I = c(0, 0, 0),
                                 LB_R = c(0, 0, 0),
                                 UB_R = c(0, 0, 0))
results <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_0VaccEff.rds")
colors <- c("Susceptible" = "blue", "Exposed" = "mediumorchid1", "Infected" = "red", "Recovered" = "black")
B <- results %>%
  mutate(t2 = round(t, digits = 0),
         tot_S = S1 + S2 + S3 + S4,
         tot_E = E1 + E2 + E3 + E4,
         tot_I = I1 + I2 + I3 + I4,
         tot_R = R1 + R2 + R3 + R4) %>% 
  group_by(t2) %>%
  summarize(med_S = median(tot_S),
            LB_S = quantile(tot_S, probs = 0.05),
            UB_S = quantile(tot_S, probs = 0.95),
            med_E = median(tot_E),
            LB_E = quantile(tot_E, probs = 0.05),
            UB_E = quantile(tot_E, probs = 0.95),
            med_I = median(tot_I),
            LB_I = quantile(tot_I, probs = 0.05),
            UB_I = quantile(tot_I, probs = 0.95),
            med_R = median(tot_R),
            LB_R = quantile(tot_R, probs = 0.05),
            UB_R = quantile(tot_R, probs = 0.95)) %>%
  ungroup() %>%
  rbind(room_fill_xtra_0.B) %>% 
  arrange(t2) %>% 
  ggplot(aes(x = t2))+ 
  geom_vline(xintercept = 14, lty = 2, size = 1.1, col = "gray")+
  geom_text(aes(x = 14, y = 3600, label = "Room 2 Fill", fontface = "italic"), angle=90, vjust=-1)+
  geom_vline(xintercept = 28, lty = 2, size = 1.1, col = "gray")+
  geom_text(aes(x = 28, y = 3600, label = "Room 3 Fill", fontface = "italic"), angle=90, vjust=-1)+
  geom_vline(xintercept = 42, lty = 2, size = 1.1, col = "gray")+
  geom_text(aes(x = 42, y = 3600, label = "Room 4 Fill", fontface = "italic"), angle=90, vjust=-1)+
  geom_ribbon(aes(ymin = LB_S, ymax = UB_S), fill = "blue", alpha = 0.3)+
  geom_line(aes(y = med_S, color = "Susceptible"), size = 1.1)+
  # geom_line(aes(y = LB_S, color = "Susceptible"), lty = 2, alpha = .5)+
  # geom_line(aes(y = UB_S, color = "Susceptible"), lty = 2, alpha = .5)+
  geom_ribbon(aes(ymin = LB_R, ymax = UB_R), fill = "black", alpha = 0.3)+
  geom_line(aes(y = med_R, color = "Recovered"), size = 1.1)+
  # geom_line(aes(y = LB_R, color = "Recovered"), lty = 2, alpha = .5)+
  # geom_line(aes(y = UB_R, color = "Recovered"), lty = 2, alpha = .5)+
  geom_ribbon(aes(ymin = LB_I, ymax = UB_I), fill = "red", alpha = 0.3)+
  geom_line(aes(y = med_I, color = "Infected"), size = 1.1)+
  # geom_line(aes(y = LB_I, color = "Infected"), lty = 2, alpha = .5)+
  # geom_line(aes(y = UB_I, color = "Infected"), lty = 2, alpha = .5)+
  geom_ribbon(aes(ymin = LB_E, ymax = UB_E), fill = "mediumorchid1", alpha = 0.3)+
  geom_line(aes(y = med_E, color = "Exposed"), size = 1.1)+
  # geom_line(aes(y = LB_E, color = "Exposed"), lty = 2, alpha = .5)+
  # geom_line(aes(y = UB_E, color = "Exposed"), lty = 2, alpha = .5)+
  scale_x_continuous(name = "Time(Days)", limits = c(0,110),
                     breaks = c(1, 25, 50, 75, 100))+
  scale_y_continuous(name = "Total Pigs", limits = c(0, 4000))+
  labs(x = "Time(Days)",
       y = "Total Pigs",
       color = "")+
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(#legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12),
        # legend.title.align = 0.5,
        legend.position = "bottom",
        strip.background = element_blank(),
        # strip.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12, color = "black"))


tiff(filename = "../Writing/Final Documents/Jan2023 Reviewer/RU Plots/Feb 28 Update/Figure2AB.tiff",
     width = 3500,
     height = 1600,
     units = 'px',
     #pointsize = 24,
     res = 300)
cowplot::plot_grid(A, B, labels = "AUTO")
dev.off()

#### FIGURE 3 VACCINATION No MDAs ####
#A. Total number infected pigs
R10_0PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Vacc_Eff = as.factor(0)) %>%  
  select(iter, Vacc_Eff, t, Tot_inf)
R10_20PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_20VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Vacc_Eff = as.factor(20)) %>%  
  select(iter, Vacc_Eff, t, Tot_inf)
R10_40PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_40VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Vacc_Eff = as.factor(40)) %>%  
  select(iter, Vacc_Eff, t, Tot_inf)
R10_60PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_60VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Vacc_Eff = as.factor(60)) %>%  
  select(iter, Vacc_Eff, t, Tot_inf)

VaccEff <- rbind(R10_0PerVacc, R10_20PerVacc, R10_40PerVacc, R10_60PerVacc)
rm(R10_0PerVacc, R10_20PerVacc, R10_40PerVacc, R10_60PerVacc)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
A <- VaccEff %>% 
  group_by(Vacc_Eff, iter) %>% 
  summarise(Tot_inf2 = max(Tot_inf)) %>% 
  ungroup() %>% 
  group_by(Vacc_Eff) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0", "20", "40", "60"), labels = c("0%", "20%", "40%", "60%"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Vacc_Eff, y = Tot_inf2, fill = Vacc_Eff))+
  # geom_bar(stat = "identity", position = "dodge", color = "black")+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Vaccine Efficacy",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("0%", "20%", "40%", "60%"))+ 
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 4000),
  #                    breaks = c(0, 1000, 2000, 3000, 4000),
  #                    labels = c("0" = "0", "1000" = "1000", "2000" = "2000", "3000" = "3000", "4000" = "4000"))+
  labs(x = "Vaccine Efficacy",
       y = "Total Infected Pigs")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
# #HISTOGRAM OF TOTAL INFECTED
# # tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
# #      width = 988,
# #      height = 800,
# #      units = 'px',
# #      #pointsize = 24,
# #      res = 300)
# VaccEff %>%
#   group_by(Vacc_Eff, iter) %>%
#   summarise(Tot_inf2 = max(Tot_inf)) %>%
#   ungroup() %>%
#   group_by(Vacc_Eff) %>%
#   mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0", "20", "40", "60"), labels = c("0%", "20%", "40%", "60%"))) %>%
#   ungroup() %>%
#   ggplot(aes(x = Tot_inf2, fill = Vacc_Eff))+
#   geom_histogram(bins = 8)+
#   # geom_histogram(binwidth = 500)+
#   facet_grid(cols = vars(Vacc_Eff))+
#   scale_fill_manual(name = "Vaccine Efficacy",
#                     values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
#                     labels = c("0%", "20%", "40%", "60%"))+ 
#   labs(x = "Total Pigs Infect",
#        y = "Count")+
#   theme_classic()+
#   theme(panel.spacing.x = unit(1.25, "lines"),
#         # legend.title = element_text(face = "bold", size = 28),
#         # legend.text = element_text(face = "bold", size = 26),
#         # legend.title.align = 0.5,
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 6),
#         axis.title = element_text(face = "bold", size = 10),
#         axis.text = element_text(face = "bold", size = 8, color = "black"))
# # dev.off()
# #988 x 800
#B. Probability of WF infection
## Prob of WF infection
R10_0PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_0VaccEff.rds") %>% 
  mutate(Vacc_Eff = as.factor(0)) %>% 
  select(iter, Vacc_Eff, t, HI)
R10_20PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_20VaccEff.rds")%>% 
  mutate(Vacc_Eff = as.factor(20)) %>% 
  select(iter, Vacc_Eff, t, HI)
R10_40PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_40VaccEff.rds")%>% 
  mutate(Vacc_Eff = as.factor(40)) %>% 
  select(iter, Vacc_Eff, t, HI)
R10_60PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_60VaccEff.rds")%>% 
  mutate(Vacc_Eff = as.factor(60)) %>% 
  select(iter, Vacc_Eff, t, HI)
VaccEff <- rbind(R10_0PerVacc, R10_20PerVacc, R10_40PerVacc, R10_60PerVacc)
rm(R10_0PerVacc, R10_20PerVacc, R10_40PerVacc, R10_60PerVacc)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
B <- VaccEff %>% 
  group_by(Vacc_Eff, iter) %>% 
  summarise(total_I = sum(HI)) %>%
  filter(total_I != 0) %>% 
  ungroup() %>% 
  group_by(Vacc_Eff) %>% 
  summarise(countwithWFInf = n(),
            PerwithWFInf = (n()/5000)*100,
            PerWOWFinf = (1-(n()/5000))*100) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0", "20", "40", "60"), labels = c("0%", "20%", "40%", "60%"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Vacc_Eff, y = PerwithWFInf, fill = Vacc_Eff))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_fill_manual(name = "Vaccine Efficacy",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("0%", "20%", "40%", "60%"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = ".25", "50" = ".5", "75" = ".75", "100" = "1"))+
  labs(x = "Vaccine Efficacy",
       y = "Probability of Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
#C. Days until first WF infection
# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3B.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
C <- VaccEff %>% 
  arrange(Vacc_Eff, iter, t) %>%
  filter(HI == 1) %>%
  group_by(Vacc_Eff, iter) %>%
  filter(row_number() == 1) %>%
  mutate(TimeSineInfIntro = t - 42) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0", "20", "40", "60"), labels = c("0%", "20%", "40%", "60%"))) %>%
  ggplot(aes(x = Vacc_Eff, y = TimeSineInfIntro, fill = Vacc_Eff))+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Vaccine Efficacy",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("0%", "20%", "40%", "60%"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = "25", "50" = "50", "75" = "75", "100" = "100"))+
  labs(x = "Vaccine Efficacy",
       y = "Days to 1st Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800

# tiff(filename = "../Writing/Final Documents/Jan2023 Reviewer/RU Plots/Figure3ABC.tiff",
#      width = 1850,
#      height = 1950,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
# cowplot::plot_grid(A, B, C, labels = "AUTO")
# dev.off()

#### FIGURE 4 VACCINATION with MDAs ####
#A. Total number infected pigs
R10_0PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Vacc_Eff = as.factor(0)) %>%  
  select(iter, Vacc_Eff, t, Tot_inf)
R10_20PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_20VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Vacc_Eff = as.factor(20)) %>%  
  select(iter, Vacc_Eff, t, Tot_inf)
R10_40PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_40VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Vacc_Eff = as.factor(40)) %>%  
  select(iter, Vacc_Eff, t, Tot_inf)
R10_60PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_60VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Vacc_Eff = as.factor(60)) %>%  
  select(iter, Vacc_Eff, t, Tot_inf)

VaccEff <- rbind(R10_0PerVacc, R10_20PerVacc, R10_40PerVacc, R10_60PerVacc)
rm(R10_0PerVacc, R10_20PerVacc, R10_40PerVacc, R10_60PerVacc)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
D <- VaccEff %>% 
  group_by(Vacc_Eff, iter) %>% 
  summarise(Tot_inf2 = max(Tot_inf)) %>% 
  ungroup() %>% 
  group_by(Vacc_Eff) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0", "20", "40", "60"), labels = c("0%", "20%", "40%", "60%"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Vacc_Eff, y = Tot_inf2, fill = Vacc_Eff))+
  # geom_bar(stat = "identity", position = "dodge", color = "black")+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Vaccine Efficacy",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("0%", "20%", "40%", "60%"))+ 
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 4000),
  #                    breaks = c(0, 1000, 2000, 3000, 4000),
  #                    labels = c("0" = "0", "1000" = "1000", "2000" = "2000", "3000" = "3000", "4000" = "4000"))+
  labs(x = "Vaccine Efficacy",
       y = "Total Infected Pigs")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
# #HISTOGRAM OF TOTAL INFECTED
# # tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
# #      width = 988,
# #      height = 800,
# #      units = 'px',
# #      #pointsize = 24,
# #      res = 300)
# VaccEff %>%
#   group_by(Vacc_Eff, iter) %>%
#   summarise(Tot_inf2 = max(Tot_inf)) %>%
#   ungroup() %>%
#   group_by(Vacc_Eff) %>%
#   mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0", "20", "40", "60"), labels = c("0%", "20%", "40%", "60%"))) %>%
#   ungroup() %>%
#   ggplot(aes(x = Tot_inf2, fill = Vacc_Eff))+
#   geom_histogram(bins = 8)+
#   # geom_histogram(binwidth = 500)+
#   facet_grid(cols = vars(Vacc_Eff))+
#   scale_fill_manual(name = "Vaccine Efficacy",
#                     values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
#                     labels = c("0%", "20%", "40%", "60%"))+ 
#   labs(x = "Total Pigs Infect",
#        y = "Count")+
#   theme_classic()+
#   theme(panel.spacing.x = unit(1.25, "lines"),
#         # legend.title = element_text(face = "bold", size = 28),
#         # legend.text = element_text(face = "bold", size = 26),
#         # legend.title.align = 0.5,
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 6),
#         axis.title = element_text(face = "bold", size = 10),
#         axis.text = element_text(face = "bold", size = 8, color = "black"))
# # dev.off()
# #988 x 800
#B. Probability of WF infection
## Prob of WF infection
R10_0PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_0VaccEff.rds") %>% 
  mutate(Vacc_Eff = as.factor(0)) %>% 
  select(iter, Vacc_Eff, t, HI)
R10_20PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_20VaccEff.rds")%>% 
  mutate(Vacc_Eff = as.factor(20)) %>% 
  select(iter, Vacc_Eff, t, HI)
R10_40PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_40VaccEff.rds")%>% 
  mutate(Vacc_Eff = as.factor(40)) %>% 
  select(iter, Vacc_Eff, t, HI)
R10_60PerVacc <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_60VaccEff.rds")%>% 
  mutate(Vacc_Eff = as.factor(60)) %>% 
  select(iter, Vacc_Eff, t, HI)
VaccEff <- rbind(R10_0PerVacc, R10_20PerVacc, R10_40PerVacc, R10_60PerVacc)
rm(R10_0PerVacc, R10_20PerVacc, R10_40PerVacc, R10_60PerVacc)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
E <- VaccEff %>% 
  group_by(Vacc_Eff, iter) %>% 
  summarise(total_I = sum(HI)) %>%
  filter(total_I != 0) %>% 
  ungroup() %>% 
  group_by(Vacc_Eff) %>% 
  summarise(countwithWFInf = n(),
            PerwithWFInf = (n()/5000)*100,
            PerWOWFinf = (1-(n()/5000))*100) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0", "20", "40", "60"), labels = c("0%", "20%", "40%", "60%"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Vacc_Eff, y = PerwithWFInf, fill = Vacc_Eff))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_fill_manual(name = "Vaccine Efficacy",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("0%", "20%", "40%", "60%"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = ".25", "50" = ".5", "75" = ".75", "100" = "1"))+
  labs(x = "Vaccine Efficacy",
       y = "Probability of Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
#C. Days until first WF infection
# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3B.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
F <- VaccEff %>% 
  arrange(Vacc_Eff, iter, t) %>%
  filter(HI == 1) %>%
  group_by(Vacc_Eff, iter) %>%
  filter(row_number() == 1) %>%
  mutate(TimeSineInfIntro = t - 42) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0", "20", "40", "60"), labels = c("0%", "20%", "40%", "60%"))) %>%
  ggplot(aes(x = Vacc_Eff, y = TimeSineInfIntro, fill = Vacc_Eff))+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Vaccine Efficacy",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("0%", "20%", "40%", "60%"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = "25", "50" = "50", "75" = "75", "100" = "100"))+
  labs(x = "Vaccine Efficacy",
       y = "Days to 1st Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800

tiff(filename = "../Writing/Final Documents/Jan2023 Reviewer/RU Plots/Feb 28 Update/Figure3ABCDEF.tiff",
     width = 3000,
     height = 2600,
     units = 'px',
     #pointsize = 24,
     res = 300)
cowplot::plot_grid(A, B, C, D, E, F, labels = "AUTO")
dev.off()

#### FIGURE 5 QUARATINE NO MDAs####
#A. total infected pigs
R10_0Qrate <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         QRate = as.factor(0)) %>%  
  select(iter, QRate, t, Tot_inf)
R10_.33Qrate <- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/No MDAs/R10_.33Qrate_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         QRate = as.factor(0.33)) %>%  
  select(iter, QRate, t, Tot_inf)
R10_.5Qrate <- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/No MDAs/R10_.5Qrate_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         QRate = as.factor(0.5)) %>%  
  select(iter, QRate, t, Tot_inf)
R10_1Qrate<- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/No MDAs/R10_1Qrate_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         QRate = as.factor(1)) %>%  
  select(iter, QRate, t, Tot_inf)

QRate <- rbind(R10_0Qrate, R10_.33Qrate, R10_.5Qrate, R10_1Qrate)
rm(R10_0Qrate, R10_.33Qrate, R10_.5Qrate, R10_1Qrate)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
A <- QRate %>% 
  group_by(QRate, iter) %>% 
  summarise(Tot_inf2 = max(Tot_inf)) %>% 
  ungroup() %>% 
  group_by(QRate) %>% 
  mutate(QRate = factor(QRate, levels = c("0", "0.33", "0.5", "1"), labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))) %>% 
  ungroup() %>% 
  ggplot(aes(x = QRate, y = Tot_inf2, fill = QRate))+
  # geom_bar(stat = "identity", position = "dodge", color = "black")+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Quarantine Rate",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c("0" = "0", "1000" = "1000", "2000" = "2000", "3000" = "3000", "4000" = "4000"))+
  labs(x = "Quarantine Rate",
       y = "Total Infected Pigs")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
# #HISTOGRAM OF TOTAL INFECTED
# # tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
# #      width = 988,
# #      height = 800,
# #      units = 'px',
# #      #pointsize = 24,
# #      res = 300)
# QRate %>%
#   group_by(QRate, iter) %>%
#   summarise(Tot_inf2 = max(Tot_inf)) %>%
#   ungroup() %>%
#   group_by(QRate) %>%
#   mutate(QRate = factor(QRate, levels = c("0", "0.33", "0.5", "1"), labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))) %>%
#   ungroup() %>%
#   ggplot(aes(x = Tot_inf2, fill = QRate))+
#   geom_histogram(bins = 8)+
#   # geom_histogram(binwidth = 500)+
#   facet_grid(cols = vars(QRate))+
#   scale_fill_manual(name = "Quarantine Rate",
#                     values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
#                     labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))+
#   labs(x = "Total Pigs Infect",
#        y = "Count")+
#   theme_classic()+
#   theme(panel.spacing.x = unit(1.25, "lines"),
#         # legend.title = element_text(face = "bold", size = 28),
#         # legend.text = element_text(face = "bold", size = 26),
#         # legend.title.align = 0.5,
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 6),
#         axis.title = element_text(face = "bold", size = 10),
#         axis.text = element_text(face = "bold", size = 8, color = "black"))
# # dev.off()
# #988 x 800
#B. Probability of WF infection
R10_0Qrate <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/No MDAs/R10_0VaccEff.rds") %>% 
  mutate(QRate = as.factor(0))%>% 
  select(iter, QRate, t, HI)
R10_.33Qrate <- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/No MDAs/R10_.33Qrate_0VaccEff.rds") %>% 
  mutate(QRate = as.factor(0.33))%>% 
  select(iter, QRate, t, HI)
R10_.5Qrate <- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/No MDAs/R10_.5Qrate_0VaccEff.rds") %>% 
  mutate(QRate = as.factor(0.5))%>% 
  select(iter, QRate, t, HI)
R10_1Qrate<- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/No MDAs/R10_1Qrate_0VaccEff.rds") %>% 
  mutate(QRate = as.factor(1))%>% 
  select(iter, QRate, t, HI)

QRate <- rbind(R10_0Qrate, R10_.33Qrate, R10_.5Qrate, R10_1Qrate)
rm(R10_0Qrate, R10_.33Qrate, R10_.5Qrate, R10_1Qrate)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure4A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
B <- QRate %>% 
  group_by(QRate, iter) %>% 
  summarise(total_I = sum(HI)) %>%
  filter(total_I != 0) %>% 
  ungroup() %>% 
  group_by(QRate) %>% 
  summarise(countwithWFInf = n(),
            PerwithWFInf = (n()/5000)*100,
            PerWOWFinf = (1-(n()/5000))*100) %>% 
  mutate(QRate = factor(QRate, levels = c("0", "0.33", "0.5", "1"), labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))) %>% 
  ungroup() %>% 
  ggplot(aes(x = QRate, y = PerwithWFInf, fill = QRate))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_fill_manual(name = "Quarantine Rate",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = ".25", "50" = ".5", "75" = ".75", "100" = "1"))+
  labs(x = "Mean time from onset of infectiousness",
       y = "Probability of Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
#C. Days until first WF infection
# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure4B.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
C <- QRate %>% 
  arrange(QRate, iter, t) %>%
  filter(HI == 1) %>%
  group_by(QRate, iter) %>%
  filter(row_number() == 1) %>%
  mutate(TimeSineInfIntro = t - 42) %>% 
  mutate(QRate = factor(QRate, levels = c("0", "0.33", "0.5", "1"), labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))) %>% 
  ggplot(aes(x = QRate, y = TimeSineInfIntro, fill = QRate))+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Quarantine Rate",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = "25", "50" = "50", "75" = "75", "100" = "100"))+
  labs(x = "Mean time from onset of infectiousness",
       y = "Days to 1st Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
# tiff(filename = "../Writing/Final Documents/Jan2023 Reviewer/RU Plots/Figure5ABC.tiff",
#      width = 1850,
#      height = 1950,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
# cowplot::plot_grid(A, B, C, labels = "AUTO")
# dev.off()

#### FIGURE 6 QUARATINE WITH MDAs####
#A. total infected pigs
R10_0Qrate <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         QRate = as.factor(0)) %>%  
  select(iter, QRate, t, Tot_inf)
R10_.33Qrate <- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/with MDAs/R10_.33Qrate_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         QRate = as.factor(0.33)) %>%  
  select(iter, QRate, t, Tot_inf)
R10_.5Qrate <- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/with MDAs/R10_.5Qrate_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         QRate = as.factor(0.5)) %>%  
  select(iter, QRate, t, Tot_inf)
R10_1Qrate<- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/with MDAs/R10_1Qrate_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         QRate = as.factor(1)) %>%  
  select(iter, QRate, t, Tot_inf)

QRate <- rbind(R10_0Qrate, R10_.33Qrate, R10_.5Qrate, R10_1Qrate)
rm(R10_0Qrate, R10_.33Qrate, R10_.5Qrate, R10_1Qrate)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
D <- QRate %>% 
  group_by(QRate, iter) %>% 
  summarise(Tot_inf2 = max(Tot_inf)) %>% 
  ungroup() %>% 
  group_by(QRate) %>% 
  mutate(QRate = factor(QRate, levels = c("0", "0.33", "0.5", "1"), labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))) %>% 
  ungroup() %>% 
  ggplot(aes(x = QRate, y = Tot_inf2, fill = QRate))+
  # geom_bar(stat = "identity", position = "dodge", color = "black")+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Quarantine Rate",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c("0" = "0", "1000" = "1000", "2000" = "2000", "3000" = "3000", "4000" = "4000"))+
  labs(x = "Quarantine Rate",
       y = "Total Pigs Infected")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
# #HISTOGRAM OF TOTAL INFECTED
# # tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
# #      width = 988,
# #      height = 800,
# #      units = 'px',
# #      #pointsize = 24,
# #      res = 300)
# QRate %>%
#   group_by(QRate, iter) %>%
#   summarise(Tot_inf2 = max(Tot_inf)) %>%
#   ungroup() %>%
#   group_by(QRate) %>%
#   mutate(QRate = factor(QRate, levels = c("0", "0.33", "0.5", "1"), labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))) %>%
#   ungroup() %>%
#   ggplot(aes(x = Tot_inf2, fill = QRate))+
#   geom_histogram(bins = 8)+
#   # geom_histogram(binwidth = 500)+
#   facet_grid(cols = vars(QRate))+
#   scale_fill_manual(name = "Quarantine Rate",
#                     values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
#                     labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))+
#   labs(x = "Total Pigs Infect",
#        y = "Count")+
#   theme_classic()+
#   theme(panel.spacing.x = unit(1.25, "lines"),
#         # legend.title = element_text(face = "bold", size = 28),
#         # legend.text = element_text(face = "bold", size = 26),
#         # legend.title.align = 0.5,
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 6),
#         axis.title = element_text(face = "bold", size = 10),
#         axis.text = element_text(face = "bold", size = 8, color = "black"))
# # dev.off()
# #988 x 800
#B. Probability of WF infection
R10_0Qrate <- readRDS(file = "Current code/Reviewer Update Code and Models/Vaccination/with MDAs/R10_0VaccEff.rds") %>% 
  mutate(QRate = as.factor(0))%>% 
  select(iter, QRate, t, HI)
R10_.33Qrate <- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/with MDAs/R10_.33Qrate_0VaccEff.rds") %>% 
  mutate(QRate = as.factor(0.33))%>% 
  select(iter, QRate, t, HI)
R10_.5Qrate <- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/with MDAs/R10_.5Qrate_0VaccEff.rds") %>% 
  mutate(QRate = as.factor(0.5))%>% 
  select(iter, QRate, t, HI)
R10_1Qrate<- readRDS(file = "Current Code/Reviewer Update Code and Models/Quarantine/with MDAs/R10_1Qrate_0VaccEff.rds") %>% 
  mutate(QRate = as.factor(1))%>% 
  select(iter, QRate, t, HI)

QRate <- rbind(R10_0Qrate, R10_.33Qrate, R10_.5Qrate, R10_1Qrate)
rm(R10_0Qrate, R10_.33Qrate, R10_.5Qrate, R10_1Qrate)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure4A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
E <- QRate %>% 
  group_by(QRate, iter) %>% 
  summarise(total_I = sum(HI)) %>%
  filter(total_I != 0) %>% 
  ungroup() %>% 
  group_by(QRate) %>% 
  summarise(countwithWFInf = n(),
            PerwithWFInf = (n()/5000)*100,
            PerWOWFinf = (1-(n()/5000))*100) %>% 
  mutate(QRate = factor(QRate, levels = c("0", "0.33", "0.5", "1"), labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))) %>% 
  ungroup() %>% 
  ggplot(aes(x = QRate, y = PerwithWFInf, fill = QRate))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_fill_manual(name = "Quarantine Rate",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = ".25", "50" = ".5", "75" = ".75", "100" = "1"))+
  labs(x = "Mean time from onset of infectiousness",
       y = "Probability of Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
#C. Days until first WF infection
# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure4B.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
F <- QRate %>% 
  arrange(QRate, iter, t) %>%
  filter(HI == 1) %>%
  group_by(QRate, iter) %>%
  filter(row_number() == 1) %>%
  mutate(TimeSineInfIntro = t - 42) %>% 
  mutate(QRate = factor(QRate, levels = c("0", "0.33", "0.5", "1"), labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))) %>% 
  ggplot(aes(x = QRate, y = TimeSineInfIntro, fill = QRate))+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Quarantine Rate",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("No Quarantine", "3 Days", "2 Days", "1 Day"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = "25", "50" = "50", "75" = "75", "100" = "100"))+
  labs(x = "Mean time from onset of infectiousness",
       y = "Days to 1st Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800
tiff(filename = "../Writing/Final Documents/Jan2023 Reviewer/RU Plots/Feb 28 Update/Figure4ABCDEF.tiff",
     width = 3000,
     height = 2600,
     units = 'px',
     #pointsize = 24,
     res = 300)
cowplot::plot_grid(A, B, C, D, E, F, labels = "AUTO")
dev.off()

#### FIGURE 7 WF FLOW NO MDAs####
#A. boxplots of total infected pigs
R10_room4 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/No MDAs/R10_Room4_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Room = as.factor(4)) %>%  
  select(iter, Room, t, Tot_inf)
R10_room3 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/No MDAs/R10_Room3_0VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Room = as.factor(3)) %>%  
  select(iter, Room, t, Tot_inf)
R10_room2 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/No MDAs/R10_Room2_0VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Room = as.factor(2)) %>%  
  select(iter, Room, t, Tot_inf)
R10_room1 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/No MDAs/R10_Room1_0VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Room = as.factor(1)) %>%  
  select(iter, Room, t, Tot_inf)
Room <- rbind(R10_room4, R10_room3, R10_room2, R10_room1)
rm(R10_room4, R10_room3, R10_room2, R10_room1)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
A <- Room %>% 
  group_by(Room, iter) %>% 
  summarise(Tot_inf2 = max(Tot_inf)) %>% 
  ungroup() %>% 
  group_by(Room) %>% 
  mutate(Room = factor(Room, levels = c("4", "3", "2", "1"), labels = c("4", "3", "2", "1"))) %>% 
  ungroup() %>% 
  ggplot(aes(x = Room, y = Tot_inf2, fill = Room))+
  # geom_bar(stat = "identity", position = "dodge", color = "black")+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Infection Roomn/Introduction",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("Room 4", "Room 3", "Room 2", "Room 1"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c("0" = "0", "1000" = "1000", "2000" = "2000", "3000" = "3000", "4000" = "4000"))+
  labs(x = "Infected Pig Introduction Room",
       y = "Probability of Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
# #HISTOGRAM OF TOTAL INFECTED
# # tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
# #      width = 988,
# #      height = 800,
# #      units = 'px',
# #      #pointsize = 24,
# #      res = 300)
# Room %>%
#   group_by(Room, iter) %>%
#   summarise(Tot_inf2 = max(Tot_inf)) %>%
#   ungroup() %>%
#   group_by(Room) %>%
#   mutate(Room = factor(Room, levels = c("4", "3", "2", "1"), labels = c("4", "3", "2", "1"))) %>% 
#   ungroup() %>%
#   ggplot(aes(x = Tot_inf2, fill = Room))+
#   geom_histogram(bins = 8)+
#   # geom_histogram(binwidth = 500)+
#   facet_grid(cols = vars(Room))+
#   scale_fill_manual(name = "Infection Roomn/Introduction",
#                     values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
#                     labels = c("Room 4", "Room 3", "Room 2", "Room 1"))+ 
#   labs(x = "Total Pigs Infect",
#        y = "Count")+
#   theme_classic()+
#   theme(panel.spacing.x = unit(1.25, "lines"),
#         # legend.title = element_text(face = "bold", size = 28),
#         # legend.text = element_text(face = "bold", size = 26),
#         # legend.title.align = 0.5,
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 6),
#         axis.title = element_text(face = "bold", size = 10),
#         axis.text = element_text(face = "bold", size = 8, color = "black"))
# # dev.off()
# #988 x 800

#B. Probability of WF Infection
R10_room4 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/No MDAs/R10_Room4_0VaccEff.rds") %>% 
  mutate(Room = as.factor(4)) %>%  
  select(iter, Room, t, HI)
R10_room3 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/No MDAs/R10_Room3_0VaccEff.rds")%>% 
  mutate(Room = as.factor(3)) %>%  
  select(iter, Room, t, HI)
R10_room2 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/No MDAs/R10_Room2_0VaccEff.rds")%>% 
  mutate(Room = as.factor(2)) %>%  
  select(iter, Room, t, HI)
R10_room1 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/No MDAs/R10_Room1_0VaccEff.rds")%>% 
  mutate(Room = as.factor(1)) %>%  
  select(iter, Room, t, HI)
Room <- rbind(R10_room4, R10_room3, R10_room2, R10_room1)
rm(R10_room4, R10_room3, R10_room2, R10_room1)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure5A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
B <- Room %>% 
  group_by(Room, iter) %>% 
  summarise(total_I = sum(HI)) %>%
  filter(total_I != 0) %>% 
  ungroup() %>% 
  group_by(Room) %>% 
  summarise(countwithWFInf = n(),
            PerwithWFInf = (n()/5000)*100,
            PerWOWFinf = (1-(n()/5000))*100) %>% 
  mutate(Room = factor(Room, levels = c("4", "3", "2", "1"), labels = c("4", "3", "2", "1"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Room, y = PerwithWFInf, fill = Room))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_fill_manual(name = "Infection Roomn/Introduction",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("Room 4", "Room 3", "Room 2", "Room 1"))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = ".25", "50" = ".5", "75" = ".75", "100" = "1"))+
  labs(x = "Infected Pig Introduction Room",
       y = "Probability of Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800

## Time to WF Infection
# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure5B.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
C <- Room %>%  
  arrange(Room, iter, t) %>%
  filter(HI == 1) %>%
  group_by(Room, iter) %>%
  filter(row_number() == 1) %>%
  mutate(TimeSineInfIntro = ifelse(Room == 4, t - 42,
                                   ifelse(Room == 3, t - 28,
                                          ifelse(Room == 2, t - 14, t)))) %>% 
  mutate(Room = factor(Room, levels = c("4", "3", "2", "1"), labels = c("4", "3", "2", "1"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Room, y = TimeSineInfIntro, fill = Room))+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Infection Roomn/Introduction",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("Room 4", "Room 3", "Room 2", "Room 1"))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = "25", "50" = "50", "75" = "75", "100" = "100"))+
  labs(x = "Infected Pig Introduction Room",
       y = "Days to 1st Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
#dev.off()
#988 x 800
# tiff(filename = "../Writing/Final Documents/Jan2023 Reviewer/RU Plots/Figure7ABC.tiff",
#      width = 1850,
#      height = 1950,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
# cowplot::plot_grid(A, B, C, labels = "AUTO")
# dev.off()

#### FIGURE 8 WF FLOW WITH MDAs####
#A. boxplots of total infected pigs
R10_room4 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/with MDAs/R10_Room4_0VaccEff.rds") %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Room = as.factor(4)) %>%  
  select(iter, Room, t, Tot_inf)
R10_room3 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/with MDAs/R10_Room3_0VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Room = as.factor(3)) %>%  
  select(iter, Room, t, Tot_inf)
R10_room2 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/with MDAs/R10_Room2_0VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Room = as.factor(2)) %>%  
  select(iter, Room, t, Tot_inf)
R10_room1 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/with MDAs/R10_Room1_0VaccEff.rds")%>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Room = as.factor(1)) %>%  
  select(iter, Room, t, Tot_inf)
Room <- rbind(R10_room4, R10_room3, R10_room2, R10_room1)
rm(R10_room4, R10_room3, R10_room2, R10_room1)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
D <- Room %>% 
  group_by(Room, iter) %>% 
  summarise(Tot_inf2 = max(Tot_inf)) %>% 
  ungroup() %>% 
  group_by(Room) %>% 
  mutate(Room = factor(Room, levels = c("4", "3", "2", "1"), labels = c("4", "3", "2", "1"))) %>% 
  ungroup() %>% 
  ggplot(aes(x = Room, y = Tot_inf2, fill = Room))+
  # geom_bar(stat = "identity", position = "dodge", color = "black")+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Infection Roomn/Introduction",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("Room 4", "Room 3", "Room 2", "Room 1"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c("0" = "0", "1000" = "1000", "2000" = "2000", "3000" = "3000", "4000" = "4000"))+
  labs(x = "Infected Pig Introduction Room",
       y = "Probability of Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
# #HISTOGRAM OF TOTAL INFECTED
# # tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure3A.tiff",
# #      width = 988,
# #      height = 800,
# #      units = 'px',
# #      #pointsize = 24,
# #      res = 300)
# Room %>%
#   group_by(Room, iter) %>%
#   summarise(Tot_inf2 = max(Tot_inf)) %>%
#   ungroup() %>%
#   group_by(Room) %>%
#   mutate(Room = factor(Room, levels = c("4", "3", "2", "1"), labels = c("4", "3", "2", "1"))) %>% 
#   ungroup() %>%
#   ggplot(aes(x = Tot_inf2, fill = Room))+
#   geom_histogram(bins = 8)+
#   # geom_histogram(binwidth = 500)+
#   facet_grid(cols = vars(Room))+
#   scale_fill_manual(name = "Infection Roomn/Introduction",
#                     values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
#                     labels = c("Room 4", "Room 3", "Room 2", "Room 1"))+ 
#   labs(x = "Total Pigs Infect",
#        y = "Count")+
#   theme_classic()+
#   theme(panel.spacing.x = unit(1.25, "lines"),
#         # legend.title = element_text(face = "bold", size = 28),
#         # legend.text = element_text(face = "bold", size = 26),
#         # legend.title.align = 0.5,
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 6),
#         axis.title = element_text(face = "bold", size = 10),
#         axis.text = element_text(face = "bold", size = 8, color = "black"))
# # dev.off()
# #988 x 800

#B. Probability of WF Infection
R10_room4 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/with MDAs/R10_Room4_0VaccEff.rds") %>% 
  mutate(Room = as.factor(4)) %>%  
  select(iter, Room, t, HI)
R10_room3 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/with MDAs/R10_Room3_0VaccEff.rds")%>% 
  mutate(Room = as.factor(3)) %>%  
  select(iter, Room, t, HI)
R10_room2 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/with MDAs/R10_Room2_0VaccEff.rds")%>% 
  mutate(Room = as.factor(2)) %>%  
  select(iter, Room, t, HI)
R10_room1 <- readRDS(file = "Current code/Reviewer Update Code and Models/WF Flow/with MDAs/R10_Room1_0VaccEff.rds")%>% 
  mutate(Room = as.factor(1)) %>%  
  select(iter, Room, t, HI)
Room <- rbind(R10_room4, R10_room3, R10_room2, R10_room1)
rm(R10_room4, R10_room3, R10_room2, R10_room1)

# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure5A.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
E <- Room %>% 
  group_by(Room, iter) %>% 
  summarise(total_I = sum(HI)) %>%
  filter(total_I != 0) %>% 
  ungroup() %>% 
  group_by(Room) %>% 
  summarise(countwithWFInf = n(),
            PerwithWFInf = (n()/5000)*100,
            PerWOWFinf = (1-(n()/5000))*100) %>% 
  mutate(Room = factor(Room, levels = c("4", "3", "2", "1"), labels = c("4", "3", "2", "1"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Room, y = PerwithWFInf, fill = Room))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_fill_manual(name = "Infection Roomn/Introduction",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("Room 4", "Room 3", "Room 2", "Room 1"))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = ".25", "50" = ".5", "75" = ".75", "100" = "1"))+
  labs(x = "Infected Pig Introduction Room",
       y = "Probability of Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
# dev.off()
#988 x 800

## Time to WF Infection
# tiff(filename = "../Writing/Manuscript Plots/Final figures/Figure5B.tiff",
#      width = 988,
#      height = 800,
#      units = 'px',
#      #pointsize = 24,
#      res = 300)
F <- Room %>%  
  arrange(Room, iter, t) %>%
  filter(HI == 1) %>%
  group_by(Room, iter) %>%
  filter(row_number() == 1) %>%
  mutate(TimeSineInfIntro = ifelse(Room == 4, t - 42,
                                   ifelse(Room == 3, t - 28,
                                          ifelse(Room == 2, t - 14, t)))) %>% 
  mutate(Room = factor(Room, levels = c("4", "3", "2", "1"), labels = c("4", "3", "2", "1"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Room, y = TimeSineInfIntro, fill = Room))+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Infection Roomn/Introduction",
                    values = c("#882255", "#DDCC77", "#999933", "#6699CC"),
                    labels = c("Room 4", "Room 3", "Room 2", "Room 1"))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0", "25" = "25", "50" = "50", "75" = "75", "100" = "100"))+
  labs(x = "Infected Pig Introduction Room",
       y = "Days to 1st Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        # legend.title = element_text(face = "bold", size = 28),
        # legend.text = element_text(face = "bold", size = 26),
        # legend.title.align = 0.5,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8, color = "black"))
#dev.off()
#988 x 800
tiff(filename = "../Writing/Final Documents/Jan2023 Reviewer/RU Plots/Feb 28 Update/Figure5ABCDEF.tiff",
     width = 3000,
     height = 2600,
     units = 'px',
     #pointsize = 24,
     res = 300)
cowplot::plot_grid(A, B, C, D, E, F, labels = "AUTO")
dev.off()



