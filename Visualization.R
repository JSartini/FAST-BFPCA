library(ggplot2)
library(ggridges)
library(ggforce)
library(scales)
library(tidyverse)
library(ggpubr)

read_directory <- function(directory){
  files = list.files(directory, full.names = T)
  max_sample = 0
  all_results = map(files, function(x){
    ind_data = read.csv(x)
    ind_data$Sample = ind_data$Sample + max_sample
    max_sample <<- max(ind_data$Sample)
    return(ind_data)
  }) %>% list_rbind()
  return(all_results)
}

# Create output directory if not present
if(!dir.exists("Figures")){
  dir.create("Figures")
}

# Set colors 
levord = c("FAST", "POLAR", "GFSR", "VMP")
labord = c("FAST", "POLAR", "GFSR", "VMP")
colors = pal_hue()(5)[1:4]
names(colors) = labord

# Read in cluster results
EF = rbind(
  read_directory("Results/EF/Can_N50_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical", N = 50, M = 50),
  read_directory("Results/EF/CGM_N50_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 50, M = 50)
) %>%
  mutate(Simulation = factor(Simulation, levels = c("CGM", "Cannonical"), 
                             labels = c("S1", "S2")), 
         Method = factor(Method, levels = levord, labels = labord))

FE = rbind(
  read_directory("Results/FE/Can_N50_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical", N = 50, M = 50),
  read_directory("Results/FE/CGM_N50_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 50, M = 50)
) %>%
  mutate(Simulation = factor(Simulation, levels = c("CGM", "Cannonical"), 
                             labels = c("S1", "S2")), 
         Method = factor(Method, levels = levord, labels = labord))

Score = rbind(
  read_directory("Results/Score/Can_N50_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical", N = 50, M = 50),
  read_directory("Results/Score/CGM_N50_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 50, M = 50)
) %>%
  mutate(Simulation = factor(Simulation, levels = c("CGM", "Cannonical"), 
                             labels = c("S1", "S2")), 
         Method = factor(Method, levels = levord, labels = labord))

# Functional accuracy
efs = EF %>%
  mutate(FName = paste0("phi[", substring(FPC_Num, 5), "](t)")) %>%
  select(-c(FPC_Num))

fes = FE %>%
  mutate(FName = "mu(t)")

efs %>%
  ggplot() + 
  geom_boxplot(aes(x = Method, y = MSE, fill = Method)) +
  theme_bw() + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_grid(Simulation~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  scale_fill_manual(values = colors) + 
  labs(y = "ISE")
ggsave("Figures/FPC_ISE.png", height = 5, width = 9)

fes %>%
  ggplot() + 
  geom_boxplot(aes(x = Method, y = MSE, fill = Method)) +
  theme_bw() + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_grid(Simulation~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  scale_fill_manual(values = colors) + 
  labs(y = "ISE")
ggsave("Figures/Mu_ISE.png", height = 5, width = 9)

# Functional coverage
efs %>%
  drop_na(Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Method, fill = Method, x = Cov, color = Method), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage") + 
  scale_fill_manual(values = colors) + 
  facet_grid(Simulation~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  xlim(0, 1)
ggsave("Figures/FPC_Cov.png", height = 5, width = 9)

fes %>%
  drop_na(Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Method, fill = Method, x = Cov, color = Method), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage") + 
  scale_fill_manual(values = colors) + 
  facet_grid(Simulation~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  xlim(0, 1)
ggsave("Figures/Mu_Cov.png", height = 5, width = 9)

# Visualize score accuracy
Score_plot_df = Score %>%
  mutate(FName = paste0("phi[", substring(FPC, 5), "](t)")) %>%
  select(-c(FPC))

Score_plot_df %>%
  drop_na(Cov) %>%
  ggplot() + 
  geom_boxplot(aes(x = Method, y = MSE, fill = Method)) +
  theme_bw() + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_grid(Simulation~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  scale_fill_manual(values = colors) + 
  labs(y = "ISE")
ggsave("Figures/Score_ISE.png", height = 5, width = 9)

# Visualize score coverage
Score_plot_df %>%
  drop_na(Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Method, fill = Method, x = Cov, color = Method), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage") + 
  scale_fill_manual(values = colors) + 
  facet_grid(Simulation~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  xlim(0, 1)
ggsave("Figures/Score_Cov.png", height = 5, width = 9)

#---------------------------------#
# Computation times               #
#---------------------------------#

M_timing = rbind(
  read_directory("Timing_Results/Time/CGM_N50_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 50, M = 50),
  read_directory("Timing_Results/Time/CGM_N50_M100") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 50, M = 100),
  read_directory("Timing_Results/Time/CGM_N50_M250") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 50, M = 250) ,
  read_directory("Timing_Results/Time/CGM_N50_M500") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 50, M = 500)
) %>%
  mutate(Minutes = Seconds/60)

N_timing = rbind(
  read_directory("Timing_Results/Time/CGM_N50_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 50, M = 50),
  read_directory("Timing_Results/Time/CGM_N100_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 100, M = 50),
  read_directory("Timing_Results/Time/CGM_N250_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 250, M = 50),
  read_directory("Timing_Results/Time/CGM_N500_M50") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", N = 500, M = 50)
) %>%
  mutate(Minutes = Seconds/60)

# Timing table
M_timing %>%
  group_by(N, M, Method) %>%
  summarize(est = median(Minutes), 
            LB = min(Minutes), 
            UB = max(Minutes),
            n_obs = n())

N_timing %>%
  group_by(N, M, Method) %>%
  summarize(est = median(Minutes), 
            LB = min(Minutes), 
            UB = max(Minutes), 
            n_obs = n())

#---------------------------------#
# Multilevel extension            #
#---------------------------------#

# Reset colors
levord = c("FAST", "POLAR", "GFSR", "VMP")
labord = c("FAST", "POLAR", "GFSR", "VMP")
colors = pal_hue()(5)[1:4]
names(colors) = labord

EF = read_directory("Results/EF/MLev") %>%
  select(-c(X)) %>%
  mutate(Method = factor(Method, levels = levord, labels = labord), 
         FName = paste0("phi[", substring(FPC_Num, 5), "](t)"))

Score = read_directory("Results/Score/MLev") %>%
  select(-c(X)) %>%
  mutate(Method = factor(Method, levels = levord, labels = labord),
         FName = paste0("phi[", substring(FPC, 5), "](t)"))

EF %>%
  ggplot() + 
  geom_boxplot(aes(x = Method, y = MSE, fill = Method)) +
  theme_bw() + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_grid(Level~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  scale_fill_manual(values = colors) + 
  labs(y = "ISE")
ggsave("Figures/ML_FPC_ISE.png", height = 5, width = 9)

EF %>%
  drop_na(Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Method, fill = Method, x = Cov, color = Method), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage") + 
  scale_fill_manual(values = colors) + 
  facet_grid(Level~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  xlim(0, 1)
ggsave("Figures/ML_FPC_Cov.png", height = 5, width = 9)

Score %>%
  drop_na(Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Method, fill = Method, x = Cov, color = Method), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage") + 
  scale_fill_manual(values = colors) + 
  facet_grid(Level~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  xlim(0, 1)
ggsave("Figures/ML_Score_Cov.png", height = 5, width = 9)
