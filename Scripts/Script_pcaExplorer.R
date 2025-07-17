#------------------------------------------------------------------------------------------------
## traitement des réplicats condition controle
#------------------------------------------------------------------------------------------------

#chargement des tableaux de comptages et suppression des 4 premièreslignes
counts_ctrl_1 <- read.table("~/LocalPedago2/7010_ReadsPerGene.out.tab",header = F, sep = "\t")
counts_ctrl_1 = counts_ctrl_1[c(5:nrow(counts_ctrl_1)),]

counts_ctrl_2 <- read.table("~/LocalPedago2/7011_ReadsPerGene.out.tab",header = F, sep = "\t")
counts_ctrl_2 = counts_ctrl_2[c(5:nrow(counts_ctrl_2)),]

counts_ctrl_3 <- read.table("~/LocalPedago2/7023_ReadsPerGene.out.tab",header = F, sep = "\t")
counts_ctrl_3 = counts_ctrl_3[c(5:nrow(counts_ctrl_3)),]

#renommage des colonnes
colnames(counts_ctrl_1) <- c("gene_ID","unstranded","R1","R2")
colnames(counts_ctrl_2) <- c("gene_ID","unstranded","R1","R2")
colnames(counts_ctrl_3) <- c("gene_ID","unstranded","R1","R2")

#selection des 2 colonnes qui nous interessent ici gene_id et unstranded
ctrl_1 <- counts_ctrl_1[, c("gene_ID", "unstranded")]
ctrl_2 <- counts_ctrl_2[,c("gene_ID", "unstranded")]
ctrl_3 <- counts_ctrl_3[,c("gene_ID", "unstranded")]






#------------------------------------------------------------------------------------------------
## traitement des réplicats condition infecté
#------------------------------------------------------------------------------------------------

# Lecture des tableaux de comptages pour la condition "inf"
counts_inf_1 <- read.table("~/LocalPedago/Infecte/7005_S1_ReadsPerGene.out.tab", header = FALSE, sep = "\t")
counts_inf_2 <- read.table("~/LocalPedago/Infecte/7014_ReadsPerGene.out.tab", header = FALSE, sep = "\t")
counts_inf_3 <- read.table("~/LocalPedago/Infecte/7017_ReadsPerGene.out.tab", header = FALSE, sep = "\t")
counts_inf_4 <- read.table("~/LocalPedago/Infecte/7020_ReadsPerGene.out.tab", header = FALSE, sep = "\t")

# Supprimer les premières lignes inutiles
counts_inf_1 <- counts_inf_1[5:nrow(counts_inf_1), ]
counts_inf_2 <- counts_inf_2[5:nrow(counts_inf_2), ]
counts_inf_3 <- counts_inf_3[5:nrow(counts_inf_3), ]
counts_inf_4 <- counts_inf_4[5:nrow(counts_inf_4), ]

# Renommage des colonnes
colnames(counts_inf_1) <- c("gene_ID", "unstranded", "strand1", "strand2")
colnames(counts_inf_2) <- c("gene_ID", "unstranded", "strand1", "strand2")
colnames(counts_inf_3) <- c("gene_ID", "unstranded", "strand1", "strand2")
colnames(counts_inf_4) <- c("gene_ID", "unstranded", "strand1", "strand2")

#selection des 2 colonnes qui nous interessent ici gene_id et unstranded
inf_1 <- counts_inf_1[, c("gene_ID", "unstranded")]
inf_2 <- counts_inf_2[, c("gene_ID", "unstranded")]
inf_3 <- counts_inf_3[, c("gene_ID", "unstranded")]
inf_4 <- counts_inf_4[, c("gene_ID", "unstranded")]





#------------------------------------------------------------------------------------------------
## traitement des réplicats condition non infecté
#------------------------------------------------------------------------------------------------

# Lecture des tableaux de comptages 
counts_nonInf_1 <- read.table("~/LocalPedago/Non_Infecte/7001_ReadsPerGene.out.tab", header = FALSE, sep = "\t")
counts_nonInf_2 <- read.table("~/LocalPedago/Non_Infecte/7004_ReadsPerGene.out.tab", header = FALSE, sep = "\t")
counts_nonInf_3 <- read.table("~/LocalPedago/Non_Infecte/7007_ReadsPerGene.out.tab", header = FALSE, sep = "\t")

# Supprimer les premières lignes inutiles
counts_nonInf_1 <- counts_nonInf_1[5:nrow(counts_nonInf_1), ]
counts_nonInf_2 <- counts_nonInf_2[5:nrow(counts_nonInf_2), ]
counts_nonInf_3 <- counts_nonInf_3[5:nrow(counts_nonInf_3), ]

# Renommage des colonnes
colnames(counts_nonInf_1) <- c("gene_ID", "unstranded", "strand1", "strand2")
colnames(counts_nonInf_2) <- c("gene_ID", "unstranded", "strand1", "strand2")
colnames(counts_nonInf_3) <- c("gene_ID", "unstranded", "strand1", "strand2")

#selection des 2 colonnes qui nous interessent ici gene_id et unstranded
nonInf_1 <- counts_nonInf_1[, c("gene_ID", "unstranded")]
nonInf_2 <- counts_nonInf_2[, c("gene_ID", "unstranded")]
nonInf_3 <- counts_nonInf_3[, c("gene_ID", "unstranded")]





#-----------------------------------------------------------------------
##Création d'un seul tableau de comptage et analyse des données
#-----------------------------------------------------------------------


#-------------------Création d'1 seul tableau de comptage ---------------------------------------

# Listes des noms des réplicats
noms_replicats <- c("ctrl_1", "ctrl_2", "ctrl_3", "inf_1", "inf_2", "inf_3", "inf_4", "nonInf_1", "nonInf_2", "nonInf_3")

# Initialisation du tableau final avec la colonne gene_id
merged_counts_all <- unique(ctrl_1[, "gene_ID", drop = FALSE])

# Fusion des tableaux restants
for (nom_replicat in noms_replicats) {
  # Récupérer le tableau correspondant au réplicat
  tableau_replicat <- get(nom_replicat)
  
  # Renommer la colonne unstranded avec le nom du réplicat
  colname <- paste("unstranded", nom_replicat, sep = "_")
  colnames(tableau_replicat)[colnames(tableau_replicat) == 'unstranded'] <- colname
  
  # Fusion avec le tableau final
  merged_counts_all <- merge(merged_counts_all, tableau_replicat, by = "gene_ID", all = TRUE)
}


# Remplacer les valeurs NA par 0
merged_counts_all[is.na(merged_counts_all)] <- 0

#conversion des comptages en données numériques
merged_counts_all[, -1] <- lapply(merged_counts_all[, -1], as.numeric)



## ACP
# Chargement des librairies nécessaires
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pcaExplorer")
library("pcaExplorer")
browseVignettes("pcaExplorer")

install.packages("markdown")

BiocManager::install("airway")

BiocManager::install("org.Hs.eg.db")

#création d'un tableau coldata

coldata <- data.frame(
  Sample_ID = c("ctrl_1", "ctrl_2", "ctrl_3", "inf_1", "inf_2", "inf_3", "inf_4", "nonInf_1", "nonInf_2", "nonInf_3"),
  Condition = c(rep("control", 3), rep("infected", 4), rep("non-infected", 3)),
  Replicate = c(1, 2, 3, 1, 2, 3, 4, 1, 2, 3),
  StringsAsFactors = FALSE
)

counts <- merged_counts_all[,-1]


library(rtracklayer)

gtf <- import("~/Documents/gencode.v45.annotation.gtf")
gtf_df <- as.data.frame(gtf)
gtf_df_unique <- gtf_df[!duplicated(gtf_df$gene_id), ]

# ACP sur les données normalisées

pcaExplorer(countmatrix = counts, coldata = coldata, annotation = gtf_df_unique)
