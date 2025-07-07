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

#---------------------------------#
# Visualize Q simulations         #
#---------------------------------#

EF_Q = rbind(
  read_directory("Results/EF/Can_Q5") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical", Q = 5),
  read_directory("Results/EF/Can_Q10") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical", Q = 10),
  read_directory("Results/EF/Can_Q20") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical", Q = 20),
  read_directory("Results/EF/Can_Q30") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical", Q = 30),
  read_directory("Results/EF/Can_Q40") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical", Q = 40),
  read_directory("Results/EF/CGM_Q5") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", Q = 5),
  read_directory("Results/EF/CGM_Q10") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", Q = 10),
  read_directory("Results/EF/CGM_Q20") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", Q = 20),
  read_directory("Results/EF/CGM_Q30") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", Q = 30),
  read_directory("Results/EF/CGM_Q40") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM", Q = 40)
)

EF_Q %>%
  filter(Method == "FAST") %>%
  mutate(FName = paste0("phi[", substring(FPC_Num, 5), "](t)")) %>%
  select(-c(FPC_Num, Method)) %>%
  mutate(Simulation = factor(Simulation, levels = c("CGM", "Cannonical"), 
                             labels = c("S1", "S2"))) %>%
  ggplot() + 
  geom_boxplot(aes(x = Q, y = MSE, group = Q)) +
  theme_bw() + 
  theme(legend.position = "none",  
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_grid(Simulation~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  labs(y = "ISE", x = "Spline Basis Dimension Q")
ggsave("Figures/Q_ISE.png", height = 5, width = 9)

EF_Q %>%
  filter(Method == "FAST") %>%
  mutate(FName = paste0("phi[", substring(FPC_Num, 5), "](t)")) %>%
  select(-c(FPC_Num, Method)) %>%
  mutate(Simulation = factor(Simulation, levels = c("CGM", "Cannonical"), 
                             labels = c("S1", "S2"))) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Q, x = Cov, group = Q), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...) mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage", y = "Spline Basis Dimension Q") + 
  facet_grid(Simulation~FName, labeller = labeller(FName = label_parsed)) + 
  xlim(0, 1)
ggsave("Figures/Q_Cov.png", height = 5, width = 9)

#---------------------------------#
# Visualize K simulations         #
#---------------------------------#

EF_K = rbind(
  read_directory("Results/EF/Can_K2") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/Can_K3") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/Can_K4") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/Can_K5") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/Can_K6") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/CGM_K2") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM"),
  read_directory("Results/EF/CGM_K3") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM"),
  read_directory("Results/EF/CGM_K4") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM"),
  read_directory("Results/EF/CGM_K5") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM"),
  read_directory("Results/EF/CGM_K6") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM")
)

EF_K %>%
  filter(Method == "FAST") %>%
  mutate(FName = paste0("phi[", substring(FPC_Num, 5), "](t)")) %>%
  select(-c(FPC_Num, Method)) %>%
  mutate(Simulation = factor(Simulation, levels = c("CGM", "Cannonical"), 
                             labels = c("S1", "S2"))) %>%
  ggplot() + 
  geom_boxplot(aes(x = K, y = MSE, group = K)) +
  theme_bw() + 
  theme(legend.position = "none",  
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_grid(Simulation~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  labs(y = "ISE", x = "FPC Dimension K")
ggsave("Figures/K_ISE.png", height = 5, width = 9)

EF_K %>%
  filter(Method == "FAST") %>%
  mutate(FName = paste0("phi[", substring(FPC_Num, 5), "](t)")) %>%
  select(-c(FPC_Num, Method)) %>%
  mutate(Simulation = factor(Simulation, levels = c("CGM", "Cannonical"), 
                             labels = c("S1", "S2"))) %>%
  ggplot() + 
  geom_density_ridges(aes(y = K, x = Cov, group = K), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...) mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage", y = "FPC Dimension K") + 
  facet_grid(Simulation~FName, labeller = labeller(FName = label_parsed)) + 
  xlim(0, 1)
ggsave("Figures/K_Cov.png", height = 5, width = 9)
  
#---------------------------------#
# Visualize alpha simulations     #
#---------------------------------#

EF_alpha = rbind(
  read_directory("Results/EF/Can_Alpha0.01") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/Can_Alpha0.05") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/Can_Alpha0.1") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/Can_Alpha0.2") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/Can_Alpha0.3") %>%
    select(-c(X)) %>%
    mutate(Simulation = "Cannonical"),
  read_directory("Results/EF/CGM_Alpha0.01") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM"),
  read_directory("Results/EF/CGM_Alpha0.05") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM"),
  read_directory("Results/EF/CGM_Alpha0.1") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM"),
  read_directory("Results/EF/CGM_Alpha0.2") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM"),
  read_directory("Results/EF/CGM_Alpha0.3") %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM")
)

EF_alpha %>%
  filter(Method == "FAST") %>%
  mutate(FName = paste0("phi[", substring(FPC_Num, 5), "](t)")) %>%
  select(-c(FPC_Num, Method)) %>%
  mutate(Simulation = factor(Simulation, levels = c("CGM", "Cannonical"), 
                             labels = c("S1", "S2"))) %>%
  ggplot() + 
  geom_boxplot(aes(x = Alpha, y = MSE, group = Alpha), width = 0.025) +
  theme_bw() + 
  theme(legend.position = "none",  
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_grid(Simulation~FName, scales = "free_y", 
             labeller = labeller(FName = label_parsed)) + 
  labs(y = "ISE", x = parse(text = "Penalty~Parameter~alpha"))
ggsave("Figures/Alpha_ISE.png", height = 5, width = 9)

EF_alpha %>%
  filter(Method == "FAST") %>%
  mutate(FName = paste0("phi[", substring(FPC_Num, 5), "](t)")) %>%
  select(-c(FPC_Num, Method)) %>%
  mutate(Simulation = factor(Simulation, levels = c("CGM", "Cannonical"), 
                             labels = c("S1", "S2"))) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Alpha, x = Cov, group = Alpha), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...) mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage", y = parse(text = "Penalty~Parameter~alpha")) + 
  facet_grid(Simulation~FName, labeller = labeller(FName = label_parsed)) + 
  xlim(0, 1)
ggsave("Figures/Alpha_Cov.png", height = 5, width = 9)

#---------------------------------#
# Timing table                    #
#---------------------------------#

Q_vals = c(20, 30, 40)
K_vals = c(3, 4, 5) 
grid_vals = expand.grid(Q_vals, K_vals)
colnames(grid_vals) = c("Q", "K")
time_df = pmap(grid_vals, function(Q, K){
  dir_name = paste0("Timing_Results/Time/CGM_Q", Q, "_K", K)
  output = read_directory(dir_name) %>%
    select(-c(X)) %>%
    mutate(Simulation = "CGM")
  return(output)
}) %>% list_rbind()

time_df %>%
  group_by(K, Q) %>%
  summarize(Time = max(Seconds)/60)
