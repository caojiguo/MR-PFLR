
rm(list=ls())
setwd("C:/~")

Y <- read.csv("Y.csv")
Z <- read.csv("Z.csv")
Y <- as.matrix(Y)[, -1]
Z <- as.matrix(Z)[, -1]

source("functionalTS.R")
source("PenalizedLS.R")
#Load ROI data set
load("ROI_data.rda")

W <- NULL
W.time <- NULL
W$left_thalamus <- ROI_data$left_thalamus$Trajectories
W$right_lateral_ventricle <- ROI_data$right_lateral_ventricle$Trajectories
W.time$left_thalamus <- ROI_data$left_thalamus$Time
W.time$right_lateral_ventricle <- ROI_data$right_lateral_ventricle$Time

## centered functional and scalar predictors ##
Y <- scale(Y)
Z[, -c(1, 3)] <- scale(Z[, -c(1, 3)])
W$left_thalamus <- apply(W$left_thalamus, 2, function(x) scale(x, scale = FALSE))
W$right_lateral_ventricle <- apply(W$right_lateral_ventricle, 2, function(x) scale(x, scale = FALSE))

##########################################################
###### local linear fit for functional observations ######
# {
# local.fit <- function(W, W.time, tgrid.len){
#   x.fit <- NULL
#   x.tgrid <- NULL
#   for(j in 1 : length(W)){
#     Wj <- W[[j]]
#     wj.time <- W.time[[j]]
#     sample_size <- nrow(Wj)
#     xj.tgrid <- seq(min(wj.time), max(wj.time), len = tgrid.len)
#     xj.fit <- matrix(NA, sample_size, tgrid.len)
#     for (i in 1 : sample_size){
#       wij <- as.vector(Wj[i, ])
#       Wij.data <- data.frame(wij)
#       Wij.data$t <- wj.time
#       xij.fit <- locpol::locpol(wij ~ t, Wij.data, xeval = xj.tgrid, bw = 0.05)
#       xj.fit[i, ] <- xij.fit$lpFit$wij
#     }
#     x.fit[[j]] <- xj.fit
#     x.tgrid[[j]] <- xj.tgrid
#   }
#   return(list(x = x.fit, tgrid = x.tgrid))
# }
# x.fit <- local.fit(W, W.time, tgrid.len = 100)
# x.hat <- x.fit$x
# x.tgrid <- x.fit$tgrid
# }


##################################################################
####### Estimation for all samples #######
{
case_result <- NULL
result.fts <- FTSfun(Y, W, Z, W.time, c(5, 5))
result.pls <- PLSfun(Y, W, Z, W.time, c(5, 5))
Yfit.fts <- result.fts$fit
Yfit.pls <- result.pls$fit
Yfit.ipls <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
for (i in 1 : ncol(Y)){
  result.ipls <- PLSfun(Y[, i], W, Z, W.time, c(5, 5))
  Yfit.ipls[, i] <- result.ipls$fit
}
res0sum.fts <- apply(abs(Yfit.fts - Y), 2, sum)
res0sum.pls <- apply(abs(Yfit.pls - Y), 2, sum)
res0sum.ipls <- apply(abs(Yfit.ipls - Y), 2, sum)
beta_fts <- result.fts$beta
case_result$beta1.grid <- W.time[[1]]
case_result$beta2.grid <- W.time[[2]]
case_result$alpha <- result.fts$alpha
case_result$beta1 <- beta_fts[[1]]
case_result$beta2 <- beta_fts[[2]]
}
save(case_result, file = "case_result.rda")

snps_label <- read.csv("snps_select_label_combine.csv")
gene_label <- paste0("(", snps_label$Gene, ")")
snps_gene_label <- paste(snps_label$SNPs, gene_label, sep = " ")
z_label <- c(colnames(Z)[1:5], snps_gene_label)

load("case_result.rda")

alpha_est <- case_result$alpha
z_select <- unlist(apply(alpha_est, 2, function(x) which(x != 0)))
z_select_ind <- names(unlist(table(z_select))[unlist(table(z_select)) >= 6])
z_select_ind <- c(2, as.numeric(z_select_ind))
alpha_select <- alpha_est[z_select_ind, ] * 10^2
alpha_select <- round(alpha_select, 3)
z_label_select <- z_label[z_select_ind]
alpha_select_result <- cbind(z_label_select, alpha_select)
colnames(alpha_select_result) <- c("covariate", colnames(Y))
write.csv(alpha_select_result, "alpha_select_result.csv", row.names = FALSE)
alpha_clinic <- round(alpha_est[1:4, ] * 10, 4)
colnames(alpha_clinic) <- colnames(Y)
rownames(alpha_clinic) <- colnames(Z)[1:4]
write.csv(alpha_clinic, "alpha_clinic_result.csv", row.names = TRUE)

beta1 <- case_result$beta1
beta2 <- case_result$beta2
beta1_grid <- case_result$beta1.grid
beta2_grid <- case_result$beta2.grid

# create data
beta_data <- data.frame(
  beta1_MMSE = beta1[1, ],
  beta1_RAVLT.lea = beta1[4, ],
  beta1_ADAS11 = beta1[2, ],
  beta1_ADASQ4 = beta1[3, ],
  beta1_RAVLT.for = beta1[5, ],
  beta1_FAQ = beta1[6, ],
  beta2_MMSE <- beta2[1, ],
  beta2_RAVLT.lea <- beta2[4, ],
  beta2_ADAS11 <- beta2[2, ],
  beta2_ADASQ4 <- beta2[3, ],
  beta2_RAVLT.for <- beta2[5, ],
  beta2_FAQ <- beta2[6, ]
)

library(ggplot2)
library(latex2exp)

# plot curves using ggplot2 
beta1_MR <- ggplot(beta_data, aes(beta1_grid)) +
  geom_line(linewidth = 0.8, aes(y = beta1_MMSE, color = "beta1_MMSE")) +
  geom_line(linewidth = 0.8, aes(y = beta1_RAVLT.lea, color = "beta1_RAVLT.lea")) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))+
  scale_color_discrete(name = "Legend", labels = c("MMSE", "RAVLT.learning")) +
  theme(legend.position = "bottom")+
  labs(#title = "",
    x = "Brain volume",
    y = TeX(r'($\hat{beta}_{1}(\cdot)$)'))

beta1_AD <- ggplot(beta_data, aes(beta1_grid)) +
  geom_line(linewidth = 0.8, aes(y = beta1_ADAS11, color = "beta1_ADAS11")) +
  geom_line(linewidth = 0.8, aes(y = beta1_ADASQ4, color = "beta1_ADASQ4")) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))+
  scale_color_discrete(name = "Legend", labels = c("ADAS11", "ADASQ4")) +
  #ggtitle('Plot 1')+
  theme(legend.position = "bottom")+
  labs(#title = "",
    x = "Brain volume",
    y = TeX(r'($\hat{beta}_{1}(\cdot)$)'))

beta1_RF <- ggplot(beta_data, aes(beta1_grid)) +
  geom_line(linewidth = 0.8, aes(y = beta1_RAVLT.for, color = "beta1_RAVLT.for")) +
  geom_line(linewidth = 0.8, aes(y = beta1_FAQ, color = "beta1_FAQ")) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))+
  scale_color_discrete(name = "Legend", labels = c("RAVLT.forgetting", "FAQ")) +
  theme(legend.position = "bottom")+
  labs(#title = "",
    x = "Brain volume",
    y = TeX(r'($\hat{beta}_{1}(\cdot)$)'))

beta2_MR <- ggplot(beta_data, aes(beta2_grid)) +
  geom_line(size = 0.8, aes(y = beta2_MMSE, color = "beta2_MMSE")) +
  geom_line(size = 0.8, aes(y = beta2_RAVLT.lea, color = "beta2_RAVLT.lea")) +
  scale_x_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5, 1.8))+
  scale_color_discrete(name = "Legend", labels = c("MMSE", "RAVLT.learning")) +
  theme(legend.position = "bottom")+
  labs(#title = "",
    x = "Brain volume",
    y = TeX(r'($\hat{beta}_{2}(\cdot)$)'))

beta2_AD <- ggplot(beta_data, aes(beta2_grid)) +
  geom_line(size = 0.8, aes(y = beta2_ADAS11, color = "beta2_ADAS11")) +
  geom_line(size = 0.8, aes(y = beta2_ADASQ4, color = "beta2_ADASQ4")) +
  scale_x_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5, 1.8))+
  scale_color_discrete(name = "Legend", labels = c("ADAS11", "ADASQ4")) +
  theme(legend.position = "bottom")+
  labs(#title = "",
    x = "Brain volume",
    y = TeX(r'($\hat{beta}_{2}(\cdot)$)'))

beta2_RF <- ggplot(beta_data, aes(beta2_grid)) +
  geom_line(size = 0.8, aes(y = beta2_RAVLT.for, color = "beta2_RAVLT.for")) +
  geom_line(size = 0.8, aes(y = beta2_FAQ, color = "beta2_FAQ")) +
  scale_x_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5, 1.8))+
  scale_color_discrete(name = "Legend", labels = c("RAVLT.forgetting", "FAQ")) +
  theme(legend.position = "bottom")+
  labs(#title = "",
    x = "Brain volume",
    y = TeX(r'($\hat{beta}_{2}(\cdot)$)'))

library(patchwork)
plot_beta1 <- beta1_MR + beta1_AD + beta1_RF
plot_beta1 <- plot_beta1 + plot_annotation(title = "Estimated functional coefficients for the ROI:left thalamus")

plot_beta2 <- beta2_MR + beta2_AD + beta2_RF
plot_beta2 <- plot_beta2 + plot_annotation(title = "Estimated functional coefficients for the ROI: right lateral ventricle")
p <- cowplot::plot_grid(plot_beta1, plot_beta2, nrow = 2) 

ggsave(p, filename = "BetaEst.pdf", width = 12, height = 8)


#######################################################
##### Estimation by 500 Bootstrap #####
{
result.BS = NULL
B.n <- 500
for (B in 1 : B.n) {
  BS <- sample(1 : nrow(Y), nrow(Y), replace = TRUE)
  Y.BS <- Y[BS, ]
  Z.BS <- Z[BS, ]
  # FD data ↓
  xx <- NULL
  xx.time  <- NULL
  for (j in 1 : length(W)) {
    xx[[j]] <- W[[j]][BS, ]
    xx.time[[j]] <- W.time[[j]]
  }
  fit.fts <- FTSfun(Y.BS, xx, Z.BS, xx.time, c(5, 5))
  result.BS[[B]] <- fit.fts$alpha
  print(B)
}
}
#### save the results

save(result.BS, file = "boot_result.rda")

load("boot_result.rda")

cogn_names <- colnames(Y)
alpha_snps <- NULL
alpha_freq_select <- NULL
snps_freq_select <- NULL
snps_label_select <- NULL
z_label_select <- NULL
snps_gene_cogn <- NULL
z_cogn <- NULL
alpha_array <- array(unlist(result.BS$alpha), dim = c(dim(result.BS$alpha[[1]]), length(result.BS$alpha)))
for(l in 1 : ncol(Y)){
  alpha_freq <- table(unlist(apply(alpha_array[, l, ], 2, function(x) which(x != 0))))
  alpha_freq_sort <- sort(unlist(alpha_freq), decreasing = TRUE)
  alpha_freq_select[[l]] <- alpha_freq_sort[alpha_freq_sort >= 250]
  z_label_select[[l]] <- z_label[as.numeric(names(alpha_freq_select[[l]]))]
  z_cogn[[l]] <- rep(cogn_names[l], length(z_label_select[[l]]))
}
z_freq_data <- cbind(unlist(z_cogn), unlist(z_label_select), as.vector(unlist(alpha_freq_select)))
colnames(z_freq_data) <- c("cogn", "covariates", "freq")
z_important <- table(z_freq_data[, 2])[table(z_freq_data[, 2]) >= 6]


for(l in c(1, 2, 4)){
  alpha_snps[[l]] <- alpha_array[-c(1:5), l, ]
  snps_freq <- table(unlist(apply(alpha_snps[[l]], 2, function(x) which(x != 0))))
  snps_freq_sort <- sort(unlist(snps_freq), decreasing = TRUE)
  snps_freq_select[[l]] <- snps_freq_sort[snps_freq_sort >= 250]
  snps_label_select[[l]] <- snps_gene_label[as.numeric(names(snps_freq_select[[l]]))]
  snps_gene_cogn[[l]] <- rep(cogn_names[l], length(snps_label_select[[l]]))
}
snps_freq_data <- cbind(unlist(snps_gene_cogn), unlist(snps_label_select), as.vector(unlist(snps_freq_select)))
colnames(snps_freq_data) <- c("cogn", "snps", "freq")
snps_freq_data <- data.frame(snps_freq_data)
snps_freq_data$freq <- as.numeric((snps_freq_data$freq))
snps_freq_data$snps <- factor((snps_freq_data$snps))

snps_plot1 <- ggplot(snps_freq_data, aes(x = freq, y = reorder(snps, freq), fill = cogn)) + 
  geom_bar(stat="identity", color = "black", position="dodge",
           size = 0.55, width = 0.55) +
  geom_text(aes_(x =~ freq + 5, label =~ abs(freq)), color="black", size = 2) + 
  scale_x_continuous(label = abs, 
                     expand = expansion(mult = c(.05, .05))) + 
  scale_y_discrete(expand = c(0.02, 0.02))+
  coord_cartesian(xlim = c(255, 410)) +
  theme_classic() + 
  ylab("") + 
  labs(x="frequency", y = "SNPs(Gene)") +
  theme(axis.text.y = element_text(size = 6, color = "black", family = "serif"),
        axis.title.x = element_text(face = "bold.italic", size = 10), 
        axis.title.y = element_text(face = "bold", size = 10),
        legend.position = "none") +      
  facet_grid(cogn ~ ., scales="free", space="free") +
  theme(strip.text = element_text(face = "bold", size = 10), 
        strip.background = element_rect(fill = "gray", 
                                        colour = "black", size = 0.8)) +  
  scale_fill_brewer(palette="Set1")

alpha_snps <- NULL
snps_freq_select <- NULL
snps_label_select <- NULL
snps_gene_cogn <- NULL
for(l in c(3, 5, 6)){
  alpha_snps[[l]] <- alpha_array[-c(1:5), l, ]
  snps_freq <- table(unlist(apply(alpha_snps[[l]], 2, function(x) which(x != 0))))
  snps_freq_sort <- sort(unlist(snps_freq), decreasing = TRUE)
  snps_freq_select[[l]] <- snps_freq_sort[snps_freq_sort >= 250]
  snps_label_select[[l]] <- snps_gene_label[as.numeric(names(snps_freq_select[[l]]))]
  snps_gene_cogn[[l]] <- rep(cogn_names[l], length(snps_label_select[[l]]))
}
snps_freq_data <- cbind(unlist(snps_gene_cogn), unlist(snps_label_select), as.vector(unlist(snps_freq_select)))
colnames(snps_freq_data) <- c("cogn", "snps", "freq")
snps_freq_data <- data.frame(snps_freq_data)
snps_freq_data$freq <- as.numeric((snps_freq_data$freq))
snps_freq_data$snps <- factor((snps_freq_data$snps))

snps_plot2 <- ggplot(snps_freq_data, aes(x = freq, y = reorder(snps, freq), fill = cogn)) + 
  geom_bar(stat="identity", color = "black", position="dodge",
           size = 0.55, width = 0.55) +
  #geom_col(width = 0.7) +       
  geom_text(aes_(x =~ freq + 5, label =~ abs(freq)), color="black", size = 2) + 
  scale_x_continuous(label = abs, 
                     expand = expansion(mult = c(.05, .05))) + 
  scale_y_discrete(expand = c(0.02, 0.02))+
  coord_cartesian(xlim = c(255, 390)) +
  theme_classic() + 
  ylab("") + 
  labs(x="frequency", y = "SNPs(Gene)") +
  theme(axis.text.y = element_text(size = 6, color = "black", family = "serif"),
        axis.title.x = element_text(face = "bold.italic", size = 10), 
        axis.title.y = element_text(face = "bold", size = 10),
        legend.position = "none") +      
  facet_grid(cogn ~ ., scales="free", space="free") +
  theme(strip.text = element_text(face = "bold", size = 10), 
        strip.background = element_rect(fill = "gray", 
                                        colour = "black", size = 0.8)) +  
  #facet_wrap(~cogn, nrow = 2, strip.position = "right")
  scale_fill_brewer(palette="Set2")

library(patchwork)
snps_plot <- snps_plot1 + snps_plot2
ggsave(snps_plot, filename = "snps_plot.pdf", width = 12, height = 14)


#### plot of functional predictors ####
library(fda.usc)
pdf("logROI_plot.pdf", width = 10, height = 4.5)
par(mfrow =c(1,2))
lt_plot1 <- plot.fdata(fdata(W$left_thalamus, W.time$left_thalamus),
                       ylab = "Log-transformed density curve", xlab = "Brain volume", 
                       main = "ROI: Left thalamus", adj = 0.5, mgp = c(1.7, 0.5, 0), 
                       xaxs = "i", cex.main = 1.2, cex.lab = 1.2, cex.axis = 0.9)

rlv_plot1 <- plot.fdata(fdata(W$right_lateral_ventricle, W.time$right_lateral_ventricle),
                        ylab = "Log-transformed density curve", xlab = "Brain volume", 
                        main = "ROI: Right lateral ventricle", adj = 0.5, mgp = c(1.7, 0.5, 0), 
                        xaxs = "i", cex.main = 1.2, cex.lab = 1.2, cex.axis = 0.9)
dev.off()


pdf("corrplot.pdf")
corrplot::corrplot(cor(Y), method = "square", tl.col="black", 
                   tl.srt=45, tl.cex = 0.8, order = "FPC")
dev.off()
