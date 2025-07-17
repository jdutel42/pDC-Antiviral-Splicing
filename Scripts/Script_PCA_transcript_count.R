ctrl_1_f <- "~/LocalPedago/Ctrl/7010/abundance.tsv"
ctrl_2_f <- "~/LocalPedago/Ctrl/7011//abundance.tsv"
ctrl_3_f <- "~/LocalPedago/Ctrl/7023//abundance.tsv"
inf_1_f <- "~/LocalPedago/Infected/7005/abundance.tsv"
inf_2_f <- "~/LocalPedago/Infected/7014/abundance.tsv"
inf_3_f <- "~/LocalPedago/Infected/7017/abundance.tsv"
inf_4_f <- "~/LocalPedago/Infected/7020/abundance.tsv"
noninf_1_f <- "~/LocalPedago/NonInfected/7001/abundance.tsv"
noninf_2_f <- "~/LocalPedago/NonInfected/7004/abundance.tsv"
noninf_3_f <- "~/LocalPedago/NonInfected/7007/abundance.tsv"

counts_ctrl_1 <- read.table(ctrl_1_f,header = T, sep = "\t")
counts_ctrl_2 <- read.table(ctrl_2_f,header = T, sep = "\t")
counts_ctrl_3 <- read.table(ctrl_3_f,header = T, sep = "\t")
counts_inf_1 <- read.table(inf_1_f,header = T, sep = "\t")
counts_inf_2 <- read.table(inf_2_f,header = T, sep = "\t")
counts_inf_3 <- read.table(inf_3_f,header = T, sep = "\t")
counts_inf_4 <- read.table(inf_4_f,header = T, sep = "\t")
counts_noninf_1 <- read.table(noninf_1_f,header = T, sep = "\t")
counts_noninf_2 <- read.table(noninf_2_f,header = T, sep = "\t")
counts_noninf_3 <- read.table(noninf_3_f,header = T, sep = "\t")

counts_ctrl_1 <- counts_ctrl_1[, c("target_id", "tpm")]
counts_ctrl_2 <- counts_ctrl_2[, c("target_id", "tpm")]
counts_ctrl_3 <- counts_ctrl_3[, c("target_id", "tpm")]
counts_inf_1 <- counts_inf_1[, c("target_id", "tpm")]
counts_inf_2 <- counts_inf_2[, c("target_id", "tpm")]
counts_inf_3 <- counts_inf_3[, c("target_id", "tpm")]
counts_inf_4 <- counts_inf_4[, c("target_id", "tpm")]
counts_noninf_1 <- counts_noninf_1[, c("target_id", "tpm")]
counts_noninf_2 <- counts_noninf_2[, c("target_id", "tpm")]
counts_noninf_3 <- counts_noninf_3[, c("target_id", "tpm")]

library(dplyr)
library(purrr)

counts_ctrl_1 <- counts_ctrl_1 %>% rename(tpm_ctrl_1 = tpm)
counts_ctrl_2 <- counts_ctrl_2 %>% rename(tpm_ctrl_2 = tpm)
counts_ctrl_3 <- counts_ctrl_3 %>% rename(tpm_ctrl_3 = tpm)
counts_inf_1 <- counts_inf_1 %>% rename(tpm_inf_1 = tpm)
counts_inf_2 <- counts_inf_2 %>% rename(tpm_inf_2 = tpm)
counts_inf_3 <- counts_inf_3 %>% rename(tpm_inf_3 = tpm)
counts_inf_4 <- counts_inf_4 %>% rename(tpm_inf_4 = tpm)
counts_noninf_1 <- counts_noninf_1 %>% rename(tpm_noninf_1 = tpm)
counts_noninf_2 <- counts_noninf_2 %>% rename(tpm_noninf_2 = tpm)
counts_noninf_3 <- counts_noninf_3 %>% rename(tpm_noninf_3 = tpm)

donnees_fusionnees <- reduce(list(counts_ctrl_1, counts_ctrl_2, counts_ctrl_3, 
                                  counts_inf_1, counts_inf_2, counts_inf_3, counts_inf_4,
                                  counts_noninf_1, counts_noninf_2, counts_noninf_3), full_join, by = "target_id")




library(dplyr)

donnees_filtrees <- donnees_fusionnees %>%
  # Filtrer les lignes où toutes les valeurs de TPM sont > 0.5
  filter(apply(.[, -which(names(.) == "target_id")], 1, function(x) all(x > 0.5)))



donnees_pca <- donnees_filtrees %>% select(-target_id)
pca_result <- prcomp(donnees_pca, center = T, scale. = T)

# Tracer les scores des deux premières composantes principales
plot(pca_result$x[, 1], pca_result$x[, 2], xlab = "PC1", ylab = "PC2", main = "PCA sur les valeurs TPM filtrées")
# text(pca_result$x[, 1], pca_result$x[, 2], labels = row.names(pca_result$x), cex = 0.6, pos = 4) ça va faire bugger

summary(pca_result)


conditions <- c("contrôle", "contrôle", "contrôle", "inf", "inf", "inf", "inf", "non-inf", "non-inf", "non-inf")

# Création du data frame pour la visualisation
dataframe_pca <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], Condition = conditions)

library(ggplot2)

ggplot(dataframe_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA par Condition", x = "PC1", y = "PC2") +
  scale_color_manual(values = c("contrôle" = "blue", "inf" = "red", "non-inf" = "green"))









