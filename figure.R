
rm(list=ls())
setwd("C:/Users/Administrator.DESKTOP-36N49P8/Desktop/Simulation studies/testn500")

################################################################
### boxplots ###
################################################################
library(ggplot2)
library(latex2exp)

fts_data1 <- read.csv("AR_result_fts_n400.csv", header = T)
Orts_data1 <- read.csv("AR_result_Orts_n400.csv", header = T)
pls_data1 <- read.csv("AR_result_pls_n400.csv", header = T)
Orls_data1 <- read.csv("AR_result_Orls_n400.csv", header = T)

rho <- rep(pls_data1$rho, 4)
methods <- rep(c("FPLS", "FPLS-Oracle", "FPWLS", "FPWLS-Oracle"), each =length(pls_data1$rho))
MSEf_value1 <- c(pls_data1$ISE_mean, Orls_data1$ISE_mean, 
                fts_data1$ISE_mean, Orts_data1$ISE_mean)
MSEs_value1 <- c(pls_data1$MSE_mean, Orls_data1$MSE_mean, 
                fts_data1$MSE_mean, Orts_data1$MSE_mean)
PE_value1 <- c(pls_data1$PE_mean, Orls_data1$PE_mean, 
              fts_data1$PE_mean, Orts_data1$PE_mean)

AR_value <- cbind(rho, methods, MSEf_value1, MSEs_value1, PE_value1)
colnames(AR_value) <- c("rho", "methods", "MSEf_value", "MSEs_value", "PE_value")
AR_data <- data.frame(AR_value)
AR_data$MSEf_value <- as.numeric(AR_data$MSEf_value)
AR_data$MSEs_value <- as.numeric(AR_data$MSEs_value)
AR_data$PE_value <- as.numeric(AR_data$PE_value)

fts_data2 <- read.csv("EX_result_fts_n400.csv", header = T)
Orts_data2 <- read.csv("EX_result_Orts_n400.csv", header = T)
pls_data2 <- read.csv("EX_result_pls_n400.csv", header = T)
Orls_data2 <- read.csv("EX_result_Orls_n400.csv", header = T)

MSEf_value2 <- c(pls_data2$ISE_mean, Orls_data2$ISE_mean, 
                 fts_data2$ISE_mean, Orts_data2$ISE_mean)
MSEs_value2 <- c(pls_data2$MSE_mean, Orls_data2$MSE_mean, 
                 fts_data2$MSE_mean, Orts_data2$MSE_mean)
PE_value2 <- c(pls_data2$PE_mean, Orls_data2$PE_mean, 
               fts_data2$PE_mean, Orts_data2$PE_mean)

EX_value <- cbind(rho, methods, MSEf_value2, MSEs_value2, PE_value2)
colnames(EX_value) <- c("rho", "methods", "MSEf_value", "MSEs_value", "PE_value")
EX_data <- data.frame(EX_value)
EX_data$MSEf_value <- as.numeric(EX_data$MSEf_value)
EX_data$MSEs_value <- as.numeric(EX_data$MSEs_value)
EX_data$PE_value <- as.numeric(EX_data$PE_value)


#########################################################
## plot of AR Sample ##
#########################################################
##plot of MSEf ##
AR_MSEf_plot <- ggplot(AR_data, aes(x = factor(rho), y = MSEf_value, 
                          group = methods, color = methods,
                          linetype = methods, shape = methods))+
  geom_line(size = 0.6)+ 
  geom_point(size = 1.7)+
  scale_x_discrete(expand = c(0.02, 0.02))+
  scale_y_continuous(expand = c(0.002, 0.002))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.box.spacing = unit(-1, "pt"))+
  labs(#title = "",
       x = TeX(r'($rho$)'),
       y = TeX(r'(MSE$_{f}$)'))

## plot of MSEs ##
AR_MSEs_plot <- ggplot(AR_data, aes(x = factor(rho), y = MSEs_value, 
                              group = methods, color = methods,
                            linetype = methods, shape = methods))+
  geom_line(size = 0.6)+ 
  geom_point(size = 1.7)+
  scale_x_discrete(expand = c(0.02, 0.02))+
  scale_y_continuous(expand = c(0.003, 0.003))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.box.spacing = unit(-1, "pt"))+
  labs(x = TeX(r'($rho$)'), y = TeX(r'(MSE$_{s}$)'))

## plot of PE ##
AR_PE_plot <- ggplot(AR_data, aes(x = factor(rho), y = PE_value, 
                              group = methods, linetype = methods, 
                              shape = methods, color = methods))+
  geom_line(size = 0.6)+ 
  geom_point(size = 1.7)+
  scale_x_discrete(expand = c(0.02, 0.02))+
  scale_y_continuous(expand = c(0.005, 0.005))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.box.spacing = unit(-1, "pt"))+
  labs(x = TeX(r'($rho$)'), y = TeX(r'(Prediction Error)'))

#########################################################
## plot of EX Sample ##
#########################################################
##plot of MSEf ##
EX_MSEf_plot <- ggplot(EX_data, aes(x = factor(rho), y = MSEf_value, 
                          group = methods, color = methods,
                          linetype = methods, shape = methods))+
  geom_line(size = 0.6)+ 
  geom_point(size = 1.7)+
  scale_x_discrete(expand = c(0.02, 0.02))+
  scale_y_continuous(expand = c(0.002, 0.002))+
  # theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.box.spacing = unit(-1, "pt"))+
  labs(#title = "",
    x = TeX(r'($rho$)'),
    y = TeX(r'(MSE$_{f}$)'))

## plot of MSEs ##
EX_MSEs_plot <- ggplot(EX_data, aes(x = factor(rho), y = MSEs_value, 
                          group = methods, color = methods,
                          linetype = methods, shape = methods))+
  geom_line(size = 0.6)+ 
  geom_point(size = 1.7)+
  scale_x_discrete(expand = c(0.02, 0.02))+
  scale_y_continuous(expand = c(0.003, 0.003))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.box.spacing = unit(-1, "pt"))+
  labs(x = TeX(r'($rho$)'), y = TeX(r'(MSE$_{s}$)'))

## plot of PE ##
EX_PE_plot <- ggplot(AR_data, aes(x = factor(rho), y = PE_value, 
                          group = methods, linetype = methods, 
                          shape = methods, color = methods))+
  geom_line(size = 0.6)+ 
  geom_point(size = 1.7)+
  scale_x_discrete(expand = c(0.02, 0.02))+
  scale_y_continuous(expand = c(0.005, 0.005))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.box.spacing = unit(-1, "pt"))+
  labs(x = TeX(r'($rho$)'), y = TeX(r'(Prediction Error)'))

library(patchwork)
plot_AR <- AR_MSEf_plot + labs(title = "AR error with sample size n=400")+ 
           AR_MSEs_plot + AR_PE_plot
plot_EX <- EX_MSEf_plot + labs(title = "EX error with sample size n=400")+ 
           EX_MSEs_plot + EX_PE_plot
p <- cowplot::plot_grid(plot_AR, plot_EX, nrow = 2) 
ggsave(p, filename = "plotn400.pdf", width = 14, height = 9)
