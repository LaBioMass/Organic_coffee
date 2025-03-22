### COMEÇAR A CALCULAR OS CALCULOS DO T² ###

# Carregar os resultados do PCA
X <- pca_result_org$x  # Scores do PCA
mean_X <- colMeans(X)  # Médias dos scores
V <- pca_result_org$rotation  # Componentes principais

# Número de amostras e componentes principais
n <- nrow(X)  # Número de amostras
p <- min(23, ncol(V))  # Número de componentes principais a serem usados (até 23)

# Calcular T² usando a fórmula correta com a limitação de componentes principais
T2_values <- rowSums((X[, 1:p] - mean_X[1:p]) %*% solve(t(V[, 1:p]) %*% V[, 1:p]) * (X[, 1:p] - mean_X[1:p]))

# Exibir os valores de T²
print("Valores de T²")
print(T2_values)

# Cálculo do T² crítico usando a distribuição F
alpha <- 0.05  # Nível de significância
T2_critical <- qf(1 - alpha, p, n - p)  # Calculando o valor crítico de T²

# Exibir o T² crítico
print("Valor crítico de T²")
print(T2_critical)

# Comparar os valores de T² com o valor crítico
outliers <- T2_values > T2_critical  # Amostras que têm T² maior que o valor crítico
print("Amostras consideradas outliers (T² > T² crítico):")
print(outliers)

### COMEÇAR A CALCULAR OS CALCULOS DO Q RESIDUALS DOS MODELOS ### min(nrow(pca_result_qc$x) - 1, ncol(pca_result_qc$x))

p <- min(nrow(pca_result_qc$x) - 1, ncol(pca_result_qc$x))
X_scores <- pca_result_qc$x[, 1:p]
V <- pca_result_qc$rotation[, 1:p]
X_reconstructed <- X_scores %*% t(V) 
X_reconstructed <- scale(X_reconstructed, center = -colMeans(pca_matrix_qc), scale = FALSE)

# Calcular os resíduos quadrados (Q residual)
Q_values <- rowSums((pca_matrix_qc - X_reconstructed)^2)

# Exibir os primeiros valores de Q residual
print(Q_values)

Q_df <- data.frame(Sample = rownames(pca_matrix_qc), Q_Residual = Q_values)

# Plotar os Q residuals
ggplot(Q_df, aes(x = Sample, y = Q_Residual)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  ggtitle("Q Residuals por Amostra") +
  xlab("Amostras") +
  ylab("Q Residual") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Definir nível de confiança (exemplo: 95%)
alpha <- 0.95
z_alpha <- qnorm(alpha)  # Quantil da normal padrão

# Obter os autovalores do PCA
eigenvalues <- pca_result_qc$sdev^2  

# Definir número de PCs usados
m <- length(eigenvalues)  # Número total de variáveis

# Calcular θ1 e θ2 (baseados nos autovalores não usados)
theta1 <- sum(eigenvalues[(p+1):m])
theta2 <- sum(eigenvalues[(p+1):m]^2)

# Calcular o limite crítico de Q residual
Q_critical <- theta1 * ((z_alpha * sqrt(2 * theta2) / theta1) + 1 + (theta2 / theta1))^(theta1 / theta2)

# Mostrar o valor crítico
print(paste("Limite crítico do Q residual:", round(Q_critical, 4)))

# Comparar Q residual com o limite crítico
outliers <- which(Q_values > Q_critical)
print("Amostras com Q residual acima do limite crítico:")
print(outliers)

#####################################################

# Calcular o índice de silhueta para os clusters formados
dist_matrix <- dist(scores_df_qc[, c("PC1", "PC2")])  # Matriz de distâncias entre as amostras
sil_score <- silhouette(kmeans_result$cluster, dist_matrix)  # Índice de silhueta

# Mostrar o índice de silhueta
print(sil_score)

# Plotar o gráfico do índice de silhueta
plot(sil_score, main = "Silhouette Plot para K-means Clustering")

# Calcular a média do índice de silhueta
mean_silhouette <- mean(sil_score[, 3])
print(paste("Média do índice de silhueta:", round(mean_silhouette, 3)))