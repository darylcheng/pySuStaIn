# Cleaning features for pySuStaIn
# Inputs:
#  - CSV from pyradiomics
# Outputs:
# - Cleaned features for pySuStaIn

library(tidyverse)
library(ggcorrplot)
library(caret)
library(randomForest)
library(ggplot2)
library(cluster)
library(ggdendro)
library(factoextra)

source("/home/daryl/Documents/ascent-clinical/creation.R")

malignant_flag <- ascent_NoduleList |>
  select(scan_name, malignancy_diagnosis_ascent, NoduleType, tstg, nstg, mstg, finalhisttyp, lymphinvasion, sprdthrairspace)

# ------------- Feature cleaning ----------------

rad_features <- read_csv("/home/daryl/Documents/data/radiomics/features_2024-04-19_14-41.csv") |>
  select(-Image, -starts_with("diagnostics")) |>
  # extract ppatient ID from Mask column
  mutate(summit_id = str_extract(Mask, "(summit-\\d{4}-[a-z]{3})")) |>
  # mutate(lesion_id = str_extract(Mask, "(\\d+).nii.gz")) |>
  mutate(scan_name = str_extract(Mask, "summit-\\d{4}-[a-z]{3}.+Lesion\\d+.nii.gz")) |>
  left_join(malignant_flag, by = "scan_name") |>
  select(summit_id, scan_name, malignancy_diagnosis_ascent, everything()) |>
  select(-Mask)
  # scale all columns beginning with 'original' to have mean 0 and variance 1
  # mutate_at(vars(starts_with("original")), scale)

ps_rad_feat <- rad_features |>
  select(-summit_id) |>
  filter(!is.na(malignancy_diagnosis_ascent)) |>
  mutate(malignancy_diagnosis_ascent = as.factor(malignancy_diagnosis_ascent))

ps_rad_feat_only <- ps_rad_feat |>
  select(starts_with("original")) |>
  mutate_at(vars(starts_with("original")), scale)

# Correlation matrix

correlation_matrix <- cor(ps_rad_feat_only)
p.mat <- cor_pmat(ps_rad_feat_only)
ggcorrplot(correlation_matrix, hc.order = TRUE, p.mat = p.mat, insig = "blank") 
highly_correlated <- findCorrelation(correlation_matrix, cutoff = 0.8)
print(highly_correlated)

# ----------------- Hierarchical clustering ------------

m <- c("average", "single", "complete", "ward")
names(m) <- c("average", "single", "complete", "ward")

# Compute agglomerative coefficient to pick best method
ac <- function(x) {
  agnes(ps_rad_feat_only, method = x)$ac
}

dg_plot <- ggdendrogram(dg)
ggsave("/home/daryl/Documents/pysustain/pySuStaIn/plots/dendrogram.png", dg_plot,
  width = 15, height = 12, dpi = 300)

# Compute silhouette scores to find optimal clusters
silhouette_scores <- numeric(20)

for (k in 2:20) {
  cluster_assignments <- cutree(model, k = k)
  sil <- silhouette(cluster_assignments, dissimilarity_matrix)
  silhouette_scores[k] <- mean(sil[, "sil_width"])
}

# Plot the average silhouette scores
ggplot(data = data.frame(K = 2:20, Silhouette = silhouette_scores[2:20]), aes(x = K, y = Silhouette)) +
  geom_line() +
  geom_point() +
  labs(title = "Average Silhouette Score by Number of Clusters",
       x = "Number of Clusters",
       y = "Average Silhouette Score")

# Select clsuter number k
num_clusters <- 4
clust <- cutree(model, k = num_clusters)

sil_scores <- silhouette(clust, as.dist(dissimilarity_matrix))

# Visualize silhouette scores
sil_plot <- fviz_silhouette(sil_scores)
print(sil_plot)
ggsave("/home/daryl/Documents/pysustain/pySuStaIn/plots/silhouette_plot.png", sil_plot, 
       width = 8, height = 6, dpi = 300)

# Visualise clusters
cluster_plot <- fviz_cluster(
  list(data = as.dist(dissimilarity_matrix), cluster = clust), 
  geom = "point") +
  labs(title = "Hierarchical clustering of radiomic features",
       subtitle = "Ward method with 4 clusters",
       caption = "Dissimilarity matrix used for clustering")

ggsave("/home/daryl/Documents/pysustain/pySuStaIn/plots/cluster_plot.png", cluster_plot,
  width = 6, height = 4, dpi = 300)

# ------------- PCA for feature selection ---------------

pca_result <- prcomp(ps_rad_feat_only, scale. = TRUE)

# Visualize the importance of components
scree_plot <- fviz_eig(pca_result,
  addlabels = TRUE, ylim = c(0, 50),
  main = "Scree plot of PCA")

ggsave("/home/daryl/Documents/pysustain/pySuStaIn/plots/scree_plot.png", scree_plot,
  width = 6, height = 4, dpi = 300, bg = "white")

# visualise pca results for individuals
pca_ind_plot <- fviz_pca_ind(pca_result, label = "none") +
  labs(title = "PCA of radiomic features",
       subtitle = "First two principal components",
       caption = "Principal component analysis of radiomic features")

ggsave("/home/daryl/Documents/pysustain/pySuStaIn/plots/pca_ind_plot.png", pca_ind_plot,
  width = 6, height = 6, dpi = 300, bg = "white")

summary(pca_result)

# Extract loadings for the first few principal components (e.g., the first 3)
loadings <- pca_result$rotation[, 1:3]

# Find the top features based on absolute loadings across the selected components
top_features <- apply(abs(loadings), 1, max) # Get the maximum loading per feature
top_pca_feature_names <- names(sort(top_features, decreasing = TRUE)[1:5]) # Adjust number to select top N features

# Subset the original data to include only these top features
pca_reduced_data_adenos <- ps_rad_feat |>
  select(scan_name, original_shape_MeshVolume, malignancy_diagnosis_ascent, NoduleType, tstg, nstg, mstg, finalhisttyp, lymphinvasion, sprdthrairspace, all_of(top_pca_feature_names)) |>
  # only keep adenos
  filter(str_detect(finalhisttyp, "adeno"))

write.csv(pca_reduced_data_adenos, "/home/daryl/Documents/data/radiomics/pysustain_pca_reduced_features_adenos.csv")


# ------------- Reduced feature set ---------------

top_variables <- results$optVariables
reduced_features <- ps_rad_feat |> 
  select(scan_name, original_shape_MeshVolume, malignancy_diagnosis_ascent, NoduleType, tstg, nstg, mstg, finalhisttyp, lymphinvasion, sprdthrairspace, all_of(top_variables)) |>
  # drop rows in malignant diagnosis is NA
  filter(!is.na(malignancy_diagnosis_ascent))

reduced_features_adeno_only <- reduced_features |>
  filter(str_detect(finalhisttyp, "adenocarcinoma")) |>
  mutate(invasion = ifelse(finalhisttyp == 'minimal'))

# Write reduced features to CSV

write.csv(reduced_features, "/home/daryl/Documents/data/radiomics/pysustain_reduced_features.csv")
