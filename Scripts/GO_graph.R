ExonSkip = read.table(file = "~/SE.MATS.JC.txt", header = T)

ExonSkip_SortPvalues <- ExonSkip[order(ExonSkip$PValue), ]



# Chargement des packages nécessaires
library(readr)
library(dplyr)
if (!requireNamespace("biomaRt", quietly = TRUE))
  install.packages("biomaRt")

library(biomaRt)

# Lecture des fichiers CSV
Kallisto7010 = read.csv(file = "~/7010_output_kallisto_count/abundance.tsv", sep = "\t")
Kallisto7014 = read.csv(file = "~/7014_output_kallisto_count/abundance.tsv", sep = "\t")

# Préparation des data frames en renommant la colonne tpm pour éviter les conflits lors de la fusion
Kallisto7010 <- rename(Kallisto7010, tmp_7010 = tpm)
Kallisto7014 <- rename(Kallisto7014, tmp_7014 = tpm)

# Fusion interne des deux data frames sur target_id pour garder uniquement les target_id en commun
donnees_combinees <- inner_join(Kallisto7010, Kallisto7014, by = "target_id")

# Récupération des GeneID avec biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
transcript_ids <- donnees_combinees$target_id
attributs <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
                   filters = "ensembl_transcript_id",
                   values = transcript_ids,
                   mart = ensembl)


# Fusion de donnees_combinees avec attributs pour associer les GeneID
donnees_combinees <- left_join(donnees_combinees, attributs, by = c("target_id" = "ensembl_transcript_id"))

# Mise à jour de la colonne target_id avec les GeneID
donnees_finale <- select(donnees_combinees, ensembl_gene_id, tmp_7010, tmp_7014)
donnees_finale <- rename(donnees_finale, target_id = ensembl_gene_id)

# Sauvegarde du nouveau data frame dans un fichier TSV
write_tsv(donnees_finale, "~/transcriptome_count_gene_id.tsv")










library(biomaRt)

# Utilisation de ENSEMBL
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Supposons que transcript_ids contient vos 266 000 identifiants
transcript_ids <- c(...) # Remplacez ceci par votre liste réelle d'identifiants

# Diviser la liste des identifiants en sous-ensembles
chunks <- split(transcript_ids, ceiling(seq_along(transcript_ids)/1000))

# Fonction pour interroger biomaRt pour un sous-ensemble d'identifiants
get_gene_ids <- function(ids) {
  getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
        filters = "ensembl_transcript_id",
        values = ids,
        mart = ensembl)
}

# Appliquer la fonction à chaque sous-ensemble
results <- lapply(chunks, get_gene_ids)

# Combiner les résultats obtenus en un seul data frame
final_results <- do.call(rbind, results)

# Affichage des premières lignes des résultats
head(final_results)


SE = read.table(file = "~/LocalPedago2/SE.MATS.JCEC.txt", h = T)
which.max(SE$IJC_SAMPLE_1)
which.max(SE$SJC_SAMPLE_1)
which.max(SE$IJC_SAMPLE_2)
which.max(SE$SJC_SAMPLE_2)

SE$Sum_IJC_SJC = SE$IJC_SAMPLE_1 + SE$SJC_SAMPLE_1 + SE$IJC_SAMPLE_2 + SE$SJC_SAMPLE_2
which.max(SE$Sum_IJC_SJC)


GO = read.table(file = "~/Documents/Alligned_ReadsPerGene.out.tab", sep = "\t")
GO <- GO[-c(3, 4)]
GO <- GO[-c(1:4), ]
GO <- GO[order(-GO[,2]), ]
top_500 <- head(GO, 2999)

write.table(top_500, "~/Documents/top_500.txt", row.names=FALSE, sep="\t", col.names=FALSE, quote=FALSE)


premiere_colonne <- top_500[, 1, drop=FALSE]
write.table(premiere_colonne, "~/Documents/premiere_colonne.txt", row.names=FALSE, sep="\t", col.names=FALSE, quote=FALSE)

# Créer une copie du dataframe top_500 pour éviter de modifier les données originales
top_500_modifie <- top_500

# Enlever ce qui se trouve après le point dans la première colonne
top_500_modifie$V1 <- sub("\\..*", "", top_500_modifie$V1)

# Écrire uniquement la première colonne modifiée dans un fichier texte
write.table(top_500_modifie$V1, "~/Documents/premiere_colonne_modifiee.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


















