# Carregar os pacotes necessários
library(MSnbase)
library(xcms)
library(ggplot2)
library(writexl)
library(patchwork)
library(limma)
library(cluster)
library(pls)
library(caret)
library(plotly)
library(writexl)
library(MASS)
library(Rtsne)
library(dplyr)
library(Spectra)
library(CAMERA)

# Ler e processar os dados com o pacote xcms que vamos usar

data_dir <- "E:/UEM/Resultados/Cafe_31_01_25/POSITIVO/ORGANICO/mzml/R/QC"
files <- sort(list.files(data_dir, pattern = "\\.mzML", full.names = TRUE))
print(files)
raw_data <- readMSData(files, mode = "onDisk")
ms_levels <- msLevel(raw_data)
print(ms_levels)
ms1_count <- sum(ms_levels == 1)
ms2_count <- sum(ms_levels == 2)

# Mostrar os resultados
cat("Número de espectros MS1:", ms1_count, "\n")
cat("Número de espectros MS2:", ms2_count, "\n")

#Plot do cromatograma do qc

bpc <- chromatogram(raw_data, aggregationFun = "max")

plot(bpc, 
     col = "goldenrod", 
     main = "BPC - QC in Organic Positive Mode", 
     xlab = "Retetion Time (s)", 
     ylab = "Intensity", 
     xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))

#Visualizando o branco

data_dir2 <- "E:/UEM/Resultados/Cafe_31_01_25/POSITIVO/ORGANICO/mzml/R/Blancks"
files2 <- sort(list.files(data_dir2, pattern = "\\.mzML", full.names = TRUE))
print(files2)
raw_data2 <- readMSData(files2, mode = "onDisk")

bpc_blanck <- chromatogram(raw_data2, aggregationFun = "max")

plot(bpc_blanck, col = "black", add = TRUE)

legend("topright", 
       legend = c("QC", "Branco"), 
       col = c("goldenrod", "black"), 
       lty = 1, 
       lwd = 2, 
       bty = "n")  # Remove a borda da legenda


#Ajustar os Parâmetros COM TESTES 

# Definir os parâmetros de CentWave
cwp_low <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 5, fitgauss = TRUE)
cwp_mid <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, fitgauss = TRUE)
cwp_high <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 20, fitgauss = TRUE)

# Aplicar o findChromPeaks usando lapply em cada conjunto de dados filtrados
xset_low_list <- findChromPeaks(raw_data, param = cwp_low)
xset_mid_list <- findChromPeaks(raw_data, param = cwp_mid)
xset_high_list <- findChromPeaks(raw_data, param = cwp_high)

cat("Número de picos detectados com snthresh = 5:", nrow(chromPeaks(xset_low_list)), "\n")
cat("Número de picos detectados com snthresh = 10:", nrow(chromPeaks(xset_mid_list)), "\n")
cat("Número de picos detectados com snthresh = 20:", nrow(chromPeaks(xset_high_list)), "\n")

# Criar os cromatogramas
bpc_t5 <- chromatogram(xset_low_list , aggregationFun = "max")
bpc_t10 <- chromatogram(xset_mid_list, aggregationFun = "max")
bpc_t20<- chromatogram(xset_high_list, aggregationFun = "max")

par(mfrow = c(1, 3))  # Define layout: 1 linha, 3 colunas
par(mfrow = c(1, 1))  # Define layout: 1 linha, 1 colunas

plot(bpc_t5, col = "firebrick", main = "Comparação dos Parâmetros - SNTHRESH 5",
     xlab = "Retention Time (s)", ylab = "Intensity", xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))

plot(bpc_t10, col = "dodgerblue", main = "Comparação dos Parâmetros - SNTHRESH 10",
     xlab = "Retention Time (s)", ylab = "Intensity", xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))

plot(bpc_t20, col = "seagreen", main = "Comparação dos Parâmetros - SNTHRESH 20",
     xlab = "Retention Time (s)", ylab = "Intensity", xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))

# O SNTHRESH ESCOLHIDO FOI 10

# Configurações com Noise
cwp_default <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, fitgauss = TRUE)
cwp_noise1 <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, noise = 10000, fitgauss = TRUE)
cwp_noise2 <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, noise = 50000, fitgauss = TRUE)

# Aplicar Noise em cada conjunto de dados
xset_default_list <- findChromPeaks(raw_data, param = cwp_default)
xset_noise1_list <- findChromPeaks(raw_data, param = cwp_noise1)
xset_noise2_list <- findChromPeaks(raw_data, param = cwp_noise2)

# Comparar o número de picos detectados para cada configuração de Noise
cat("Número de picos detectados com cwp_default_list", nrow(chromPeaks(xset_default_list)), "\n")
cat("Número de picos detectados com cwp_noise1_list", nrow(chromPeaks(xset_noise1_list)), "\n")
cat("Número de picos detectados com cwp_noise2_list", nrow(chromPeaks(xset_noise2_list)), "\n")

# Criar os cromatogramas
bpc_def <- chromatogram(xset_default_list, aggregationFun = "max")
bpc_n1 <- chromatogram(xset_noise1_list, aggregationFun = "max")
bpc_n2 <- chromatogram(xset_noise2_list, aggregationFun = "max")

plot(bpc_def, col = "firebrick", main = "Comparação dos Parâmetros - Noise 500",
     xlab = "Retention Time (s)", ylab = "Intensity", xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))

plot(bpc_n1, col = "dodgerblue", main = "Comparação dos Parâmetros - Noise 10000",
     xlab = "Retention Time (s)", ylab = "Intensity", xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))

plot(bpc_n2, col = "seagreen", main = "Comparação dos Parâmetros - Noise 50000",
     xlab = "Retention Time (s)", ylab = "Intensity", xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))


#O VALOR DE NOISE PARA DESCARTE DE RUIDO É DE 1000, pois apesar de aperecer um picos de ruido, eu consigo usar a função adjustedRtime. os brancos serão reduzidos com a matriz.

#tempo de retenção (nao esquecer de carregar o setup escolhido)

data_dir_3 <- "E:/UEM/Resultados/Cafe_31_01_25/POSITIVO/ORGANICO/mzml"

files_amostras <- sort(list.files(data_dir_3, pattern = "\\.mzML", full.names = TRUE))

print(files_amostras)

raw_data_completed <- readMSData(files_amostras, mode = "onDisk")

n <- 3
I <- 1000
cwp_setup <- CentWaveParam(ppm = 10, mzdiff = -0.001, peakwidth = c(5, 20), snthresh = 10, prefilter = c(n, I), noise = I,  fitgauss = TRUE)

# Aplicando findChromPeaks() para cada conjunto de dados
xset_setup <- findChromPeaks(raw_data_completed, param = cwp_setup)

#AGRUPAMENTO DOS PICOS
# Criando o PeakDensityParam

pdp <- PeakDensityParam (sampleGroups = rep(1, length(files_amostras), bw = 20, minFraction = 0.01, minSamples = 1))

xset_grouped <- groupChromPeaks(xset_setup, param = pdp)
chromPeaks(xset_grouped)

#aqui é meu código para retenção (ATÉ AQUI NAO MEXER MAIS LKKKSDKJD)

# Ajustar tempo de retenção
pdp_adjust <- PeakGroupsParam(minFraction = 0.1, smooth = "linear", span = 0.9)

xset_adjusted <- adjustRtime(xset_grouped, param = pdp_adjust)

#Comparação dos tempos de retenção antes e depois do ajuste
rt_before <- rtime(xset_grouped)
head(rt_before)
rt_after <- adjustedRtime(xset_adjusted)
head(rt_after)

# Escolher o m/z para visualização
rt_range <- range(rtime(xset_grouped))
print(rt_range)

error_mz <- 0.05
mz_target <- 195.0853

#Cromatrograma antes do ajuste
chr_before_sample <- chromatogram(xset_grouped, mz = c(mz_target - error_mz, mz_target + error_mz), rt = rt_range)

# Cromatograma depois do ajuste para o m/z específico
chr_after_sample <- chromatogram(xset_adjusted, mz = c(mz_target - error_mz, mz_target + error_mz), rt = rt_range)

# Plots dos cromatogramas antes e depois do ajuste para a amostra específica
plot(chr_before_sample, col = "green", main = paste("Cromatograma antes e depois do ajuste para m/z", mz_target, "Amostra"), xlim = c(0, 1200),ylim = c(0, 2e5))

plot(chr_after_sample, col = "blue", add = T)

# Criando o PeakDensityParam dps do adjusted

pdp <- PeakDensityParam (sampleGroups = rep(1, length(files_amostras), bw = 15, minFraction = 0.02, minSamples = 1))

xset_grouped1 <- groupChromPeaks(xset_adjusted, param = pdp)
chromPeaks(xset_grouped1)

#fill peakings

std_filled <- fillChromPeaks(xset_grouped1, param = FillChromPeaksParam(ppm = 10))
matrix_data <- featureValues(std_filled)
View(matrix_data)

# Verificar a estrutura da matriz extraída
head(matrix_data)

#Estatistica já é aqui

#picos [M+] e outros aduteros 

# Criando um objeto CAMERA a partir do XCMS
xset_converted <- as(std_filled, "xcmsSet")

# Criar o objeto CAMERA
xset_annotated <- annotate(xset_converted)
print(xset_annotated)

isotope_data <- getPeaklist(xset_annotated)
head(isotope_data)
isotope_free_matrix <- as.data.frame(isotope_data)
View(isotope_free_matrix)




# Converter colunas vazias em NA (caso ainda não tenha sido feito)
isotope_free_matrix[, 58][isotope_free_matrix[, 58] == ""] <- NA

# Substituir valores NA por "[M]+"
isotope_free_matrix[, 58][is.na(isotope_free_matrix[, 58])] <- "[M]+"
View(isotope_free_matrix)
resume0 <- isotope_free_matrix[,-c(10:57)]
View(resume0)
     
# Filtrar apenas as linhas que contêm "[M]+" ou "[x][M]+"
isotope_mplus <- isotope_free_matrix[grep("(\\[M\\]\\+|\\[\\d+\\]\\[M\\]\\+)", isotope_free_matrix[, 58]), ]
View(isotope_mplus)
resume <- isotope_mplus[,-c(10:57)]
View(resume)

# Criar uma cópia da matriz para armazenar os resultados
isotope_sums <- isotope_mplus  

# Para cada linha de isotope_mplus, somar as intensidades de [M+1]+ e [M+2]+ do mesmo grupo
for (i in 1:nrow(isotope_mplus)) {
  
  # Identificar o padrão do M+ atual (exemplo: "[1][M]+")
  current_mplus <- isotope_mplus[i, 58]
  
  # Criar os padrões correspondentes a [M+1]+ e [M+2]+
  m_plus_1 <- gsub("\\[M\\]\\+", "[M+1]+", current_mplus)
  m_plus_2 <- gsub("\\[M\\]\\+", "[M+2]+", current_mplus)
  
  # Encontrar as linhas correspondentes em isotope_free_matrix
  matching_rows <- which(isotope_free_matrix[, 58] %in% c(m_plus_1, m_plus_2))
  
  # Definir apenas as colunas de 10 a 57 para a soma
  cols_to_sum <- 10:57  
  
  # Somar as intensidades das linhas correspondentes
  if (length(matching_rows) > 0) {
    isotope_sums[i, cols_to_sum] <- isotope_sums[i, cols_to_sum] + 
      colSums(isotope_free_matrix[matching_rows, cols_to_sum], na.rm = TRUE)
  }
}

View(isotope_mplus)
View(isotope_sums)
write_xlsx(isotope_sums, "isotope_sums.xlsx")

resume <- isotope_free_matrix[,-c(10:57)]
View(resume)

linhas <- c(5,9,10,11,15,24,26,28,29,33,35,36,65,66,67,68,82,92,102,97)
novamatrix <- isotope_mplus[-c(linhas), ]
View(novamatrix)

# Supondo que as 6 primeiras amostras são os brancos
blank_values <- novamatrix[,10:15]  # Extrair os brancos
View(blank_values)

# Calcular o maior valor dos brancos para cada m/z (linha)
max_blank_values <- apply(blank_values, 1, max)

# Subtrair os valores máximos dos brancos das amostras e do controle de qualidade
novamatrix[, 16:57] <- lapply(novamatrix[, 16:57], as.numeric)

# Verifique se max_blank_values tem o mesmo número de linhas que novamatrix
if (length(max_blank_values) == nrow(novamatrix)) {
  # Subtrair os valores máximos dos brancos das amostras
  for (i in 16:57) {
    novamatrix[, i] <- novamatrix[, i] - max_blank_values
  }
} else {
  stop("O número de elementos em max_blank_values não corresponde ao número de linhas em novamatrix.")
}

# Substituir NA e valores menores que zero por zero nas colunas 16 a 57
novamatrix[, 16:57][is.na(novamatrix[, 16:57]) | novamatrix[, 16:57] < 0] <- 0

# Remover as linhas que possuem apenas zeros entre as colunas 16 e 58
novamatrix <- novamatrix[apply(novamatrix[, 16:57], 1, function(x) any(x != 0)), ]

# Função para normalizar por LOESS
normalize_loess <- function(x) {
  loess_fit <- loess(x ~ seq_along(x))  # Ajusta um modelo LOESS
  normalized_values <- predict(loess_fit)  # Prediz os valores normalizados
  return(normalized_values)
}

# Aplicar a normalização LOESS nas colunas 16 a 57
novamatrix[, 16:57] <- apply(novamatrix[, 16:57], 2, normalize_loess)

# Verifique se as colunas foram normalizadas corretamente
head(novamatrix[, 16:57])

#Trocando nome das amostras
# Criando os nomes fixos para as colunas QC
colunas_qc <- c("QC1", "QC2", "QC3")

# Atribuindo os novos nomes de coluna para as colunas 16 a 18
colnames(novamatrix)[16:18] <- colunas_qc

# Criando os nomes fixos para as colunas QC
colunas_cb <- c("CBE1", "CBE2", "QBE3","CBM1", "CBM2", "QBM3","CBV1", "CBV2", "QBV3")
colnames(novamatrix)[19:27] <- colunas_cb

colunas_j <- c("J1", "J2", "J3")
colnames(novamatrix)[28:30] <- colunas_j

colunas_itf <- c("25.1", "25.2", "25.3","35.1","35.2","35.3","45.1","45.2","45.3","55.1","55.2","55.3","65.1","65.2","65.3","75.1","75.2","75.3","85.1","85.2","85.3","95.1","95.2","95.3")
colnames(novamatrix)[31:54] <- colunas_itf

colunas_sp <- c("sp1","sp2","sp3")
colnames(novamatrix)[55:57] <- colunas_sp

#REPETIR ATÉ CHEGAR NAS 57
View(novamatrix)

# Extraindo a coluna de m/z (primeira coluna)
mz_values <- novamatrix[-c(1,11,13), 1]  # valores de m/z

# Selecionando as colunas de 16 a 58 (amostras para PCA)
data_for_pca <- novamatrix[, 16:57]

data_for_pca <- data_for_pca[,-c(10:12)]
data_for_pca <- data_for_pca[-c(1,11,13),]

# Atribuindo os valores de m/z como nomes das linhas
rownames(data_for_pca) <- mz_values
View(data_for_pca)

#USAR data_for_pca considerando as triplicatas

num_cols <- ncol(data_for_pca)
grouped_matrix <- sapply(seq(1, num_cols, by = 3), function(i) {
  rowMeans(data_for_pca[, i:min(i+2, num_cols)], na.rm = TRUE)  # Evita erro caso o número de colunas não seja múltiplo de 3
})

# Transformar em data.frame, se necessário
grouped_matrix <- as.data.frame(grouped_matrix)
View(grouped_matrix)
# Exibir a nova matriz
head(grouped_matrix)

# Atribuindo os valores de m/z como nomes das linhas
rownames(grouped_matrix) <- mz_values

# Criando os nomes fixos para as colunas QC
qc_name <- c("QC")
colnames(grouped_matrix)[1] <- qc_name

cb_name <- c("C75", "C65", "CV")
colnames(grouped_matrix)[2:4] <- cb_name

j <- c("J75")
colnames(grouped_matrix)[5] <- j

itf <- c("25","35","45","55","65","75","85","95")
colnames(grouped_matrix)[6:13] <- itf

sp <- c("SP65")
colnames(grouped_matrix)[14] <- sp

media_pca <- grouped_matrix[,-c(4)]
View(media_pca)

View(grouped_matrix)


# Aplicar o PCA (com controle de qualidade)
pca_matrix_qc <- t(data_for_pca)
pca_result_qc <- prcomp(pca_matrix_qc, center = TRUE, scale. = TRUE)

# Ver o resumo dos scores para todos os PCs
summary(pca_result_qc)

# Criar DataFrame com os scores do PCA
scores_df_qc <- as.data.frame(pca_result_qc$x)
scores_df_qc$Sample <- rownames(scores_df_qc)

# Criando um grupo que muda a cada 3 amostras (PLOT DA TRIPLICATAS)
scores_df_qc$Group <- factor(1:13)

ggplot(scores_df_qc, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 5) +
  theme_minimal() + 
  ggtitle("PCA - Positive Mode Organic Fase - with Quality Control (QC)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_qc)$importance[2, 2], 1), "%)")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  # Diminui o tamanho do texto da legenda
    legend.title = element_text(size = 10), # Ajusta o tamanho do título da legenda
    legend.background = element_rect(fill = "white", color = "black") # Adiciona a caixa ao redor da legenda
  ) +
  scale_color_manual(
    name = "Sample Group",  
    values = c(
      "red", "blue", "green", "purple", "orange", "pink", "cyan", 
      "brown", "yellow", "magenta", "gray", "darkgreen", "darkblue","gold"))


#Criar o Kmeans 

set.seed(123)
wcss <- sapply(1:10, function(k) {
  kmeans(scores_df_qc[, c("PC1", "PC2")], centers = k, nstart = 25)$tot.withinss
})

print(wcss)  # Checar se os valores são todos zeros
plot(1:10, wcss, type = "b", pch = 19, col = "blue",
     xlab = "Número de clusters", ylab = "Soma dos Quadrados Dentro dos Clusters (WCSS)",
     main = "Método do Cotovelo")

set.seed(123) # Garantir reprodutibilidade
k <- 4 # escolhemos clusters com base no método do cotovelo
kmeans_result <- kmeans(scores_df_qc[, c("PC1", "PC2")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_qc$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_qc, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 4) +  # Diminui o tamanho dos nomes das amostras
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 4.0, linetype = 4) +
  theme_minimal() +
  ggtitle(label = "K-means Clustering in PCA (With Quality Control)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_qc)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("sienna2", "navyblue","red","green")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  # Reduz o tamanho do texto da legenda
    legend.title = element_text(size = 10), # Ajusta o tamanho do título da legenda
    legend.background = element_rect(fill = "white", color = "black") # Adiciona caixa ao redor da legenda
  ) 

colours()

# Calcular a matriz de distâncias (usando a matriz PCA reduzida)
dist_matrix <- dist(scores_df_qc[, c("PC1", "PC2", "PC3")])

# Aplicar o método de agrupamento hierárquico (usando "ward.D2")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plotar o dendrograma
plot(hclust_result, main = "Dendrogram with QC", 
     xlab = "Samples", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 3, border = "red")

#Loadings

loadings_df <- as.data.frame(pca_result_qc$rotation)  # Loadings dos PCs
loadings_df$Feature <- rownames(loadings_df)  # Adiciona nomes das variáveis
head(loadings_df)

loadings_pc1_pc2 <- loadings_df[, c("PC1", "PC2", "Feature")]

# Selecionar os 10 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2$Importance <- abs(loadings_pc1_pc2$PC1) + abs(loadings_pc1_pc2$PC2)
top20_loadings_qc <- loadings_pc1_pc2[order(-loadings_pc1_pc2$Importance), ][1:20, ]
View(top20_loadings_qc)

top20_loadings_qc$Feature <- sub("^X", "", top20_loadings_qc$Feature)

# Criar o gráfico corrigido
ggplot(top20_loadings_qc, aes(x = PC1, y = PC2, color = Feature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               size = 1) +
  scale_color_manual(values = rainbow(length(top20_loadings_qc$Feature))) +  # Escolhe cores para cada Feature
  theme_minimal() +
  ggtitle("Loadings - Top 10 m/z with QC") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings") +
  theme(legend.position = "right",    
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "black"),  
        axis.ticks = element_line(colour = "black"))  

top20_mz <- data.frame(mz = top20_loadings_qc$mz)

# Salvar em um arquivo Excel
write_xlsx(top20_mz, "top20_loadings_mz_R_positive_organic_cromatografia.xlsx")

#Gráfico em 3D

p <- ggplot(scores_df_qc, aes(x = PC1, y = PC2, z = PC3, color = Cluster)) +
  geom_point(size = 4) +
  labs(title = "PCA 3D (PC1, PC2, PC3)") +
  theme_minimal()

ggplotly(p)

plot_ly(scores_df_qc, 
        x = ~PC1, y = ~PC2, z = ~PC3, 
        color = ~Cluster, colors = c("red", "blue", "green", "orange"),
        text = ~Sample, hoverinfo = "text") %>%
  add_markers() %>%
  add_text(text = ~Sample, showlegend = FALSE, textposition = "top center") %>%
  layout(title = "PCA 3D (PC1, PC2, PC3)",
         scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")))


# Supondo que as componentes principais estão em scores_df_qc (excluindo a coluna de identificação da amostra)
pca_data <- scores_df_qc[, -which(names(scores_df_qc) == "Sample")]
View(pca_data)

# Aplicando t-SNE no PCA
set.seed(42)  # Definindo a semente para reprodutibilidade
tsne_result <- Rtsne(pca_data, dims = 2, perplexity = 10)  # Reduzindo para 2 dimensões

# Criando o dataframe com os resultados do t-SNE
tsne_df <- data.frame(
  tSNE1 = tsne_result$Y[, 1],
  tSNE2 = tsne_result$Y[, 2],
  Sample = scores_df_qc$Sample
)

# Plotando o gráfico
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "t-SNE a partir do PCA",
       x = "tSNE1",
       y = "tSNE2")

# Aplicando o k-means nos resultados do t-SNE
set.seed(123)  # Para reprodutibilidade
kmeans_result <- kmeans(tsne_result$Y, centers = 3)  

#Adicionando os resultados do k-means ao dataframe
tsne_df <- data.frame(
  tSNE1 = tsne_result$Y[, 1],
  tSNE2 = tsne_result$Y[, 2],
  Cluster = factor(kmeans_result$cluster))

tsne_df$Sample <- scores_df_qc$Sample


# Visualizando com ggplot2 e adicionando elipses
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  geom_point(size = 4) +  # Adiciona os pontos
  stat_ellipse(aes(group = Cluster), level = 0.95, linetype = 2) +  # Adiciona as elipses
  geom_text(aes(label = Sample), vjust = -1, size = 3) +  # Adiciona os nomes das amostras
  theme_minimal() +
  ggtitle("t-SNE com K-means e Elipses de Agrupamento") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  theme(plot.title = element_text(hjust = 0.5))

# Garantir que scores_df_qc[[mz_column]] e tsne_df$tSNE1 sejam numéricos
# Certifique-se de que as amostras estão ordenadas da mesma forma
# Verifique se o identificador de amostra é o mesmo para PCA e t-SNE
identificadores_pca <- rownames(pca_data)  # Supondo que as amostras de PCA tenham um identificador de linha
identificadores_tsne <- tsne_df$Sample  # Supondo que as amostras de t-SNE estão na coluna 'Sample'

# Alinhe as amostras de acordo com o identificador
pca_data <- pca_data[identificadores_pca %in% identificadores_tsne, ]
tsne_df <- tsne_df[identificadores_tsne %in% identificadores_pca, ]

# Verifique novamente as dimensões
dim(pca_data)  # O número de amostras deve ser o mesmo
dim(tsne_df)   # O número de amostras deve ser o mesmo

# Calcule a correlação entre as PCs do PCA e as duas componentes do t-SNE
View(scores_df_qc)

# Inicialize o dataframe para armazenar as correlações
# Obter os componentes principais (scores) do PCA
pca_scores <- pca_result$x

#ANALpca_rotation_df#ANALISE SEM O QC

# Normalizar os dados usando LOESS considerando os QC
# Aqui, estamos dizendo que as colunas 7 a 9 são as QC
# Supondo que você tenha uma lista de amostras de controle chamada `qc_samples`

abc <- data_for_pca[,4:39]
View(abc)

# Remover as linhas da matriz
feature_matrix_without_qc_normalized <- normalizeBetweenArrays(abc, method = "cyclicloess")
View(feature_matrix_without_qc_normalized)

# Aplicar o PCA (sem controle de qualidade)
pca_matrix_no_qc <- t(feature_matrix_without_qc_normalized)
pca_result_no_qc <- prcomp(pca_matrix_no_qc, center = TRUE, scale. = TRUE)

# Criar DataFrame com os scores do PCA
scores_df_no_qc <- as.data.frame(pca_result_no_qc$x)
scores_df_no_qc$Sample <- rownames(scores_df_no_qc)
scores_df_no_qc$Group <- factor(rep(1:13, each = 3, length.out = nrow(scores_df_no_qc)))

# Plotar o PCA
ggplot(scores_df_no_qc, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 0) +
  theme_minimal() + 
  ggtitle("PCA - Positive Mode Organic Fase - without Quality Control (QC)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_no_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_no_qc)$importance[2, 2], 1), "%)")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  # Diminui o tamanho do texto da legenda
    legend.title = element_text(size = 10), # Ajusta o tamanho do título da legenda
    legend.background = element_rect(fill = "white", color = "black") # Adiciona a caixa ao redor da legenda
  ) +
  scale_color_manual(
    name = "Sample Group",  
    values = c(
      "blue", "green", "purple", "orange", "pink", "cyan", 
      "brown", "yellow", "magenta", "gray", "darkgreen", "darkblue","gold"),
    labels = scores_df_no_qc$Sample[seq(1, length(scores_df_no_qc$Sample), by = 3)] 
  )

# Calcular a soma dos quadrados dentro dos clusters para diferentes k

set.seed(123)
wcss <- sapply(1:10, function(k) {
  kmeans(scores_df_no_qc[, c("PC1", "PC2")], centers = k, nstart = 25)$tot.withinss
})

print(wcss)  # Checar se os valores são todos zeros
plot(1:10, wcss, type = "b", pch = 19, col = "blue",
     xlab = "Número de clusters", ylab = "Soma dos Quadrados Dentro dos Clusters (WCSS)",
     main = "Método do Cotovelo")

set.seed(123) # Garantir reprodutibilidade
k <- 3 # escolhemos clusters com base no método do cotovelo
kmeans_result <- kmeans(scores_df_no_qc[, c("PC1", "PC2")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_no_qc$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_no_qc, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 4) +  # Diminui o tamanho dos nomes das amostras
  stat_ellipse(geom = "polygon", alpha = 0.2) +
  theme_minimal() +
  ggtitle("K-means Clustering in PCA - Organic Positive Mode (Without QC)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_no_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_no_qc)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("orange", "blue", "magenta", "firebrick")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  # Reduz o tamanho do texto da legenda
    legend.title = element_text(size = 10), # Ajusta o tamanho do título da legenda
    legend.background = element_rect(fill = "white", color = "black") # Adiciona caixa ao redor da legenda
  )

# Calcular a matriz de distâncias (usando a matriz PCA reduzida)
dist_matrix <- dist(scores_df_no_qc[, c("PC1", "PC2", "PC3")])

# Aplicar o método de agrupamento hierárquico (usando "ward.D2")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plotar o dendrograma
plot(hclust_result, main = "Dendrogram - Without QC", 
     xlab = "Samples", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 3, border = "saddlebrown")

#Loadings

loadings_df_amostral <- as.data.frame(pca_result_qc$rotation)  # Loadings dos PCs
loadings_df_amostral$Feature <- rownames(loadings_df_amostral)  # Adiciona nomes das variáveis
head(loadings_df_amostral)

loadings_pc1_pc2_amostral <- loadings_df_amostral[, c("PC1", "PC2", "Feature")]

# Selecionar os 10 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2_amostral$Importance <- abs(loadings_pc1_pc2_amostral$PC1) + abs(loadings_pc1_pc2_amostral$PC2)
top100_loadings_amostras <- loadings_pc1_pc2_amostral[order(-loadings_pc1_pc2_amostral$Importance), ][1:5, ]
top100_loadings_amostras$Feature <- sub("^X", "", top100_loadings_amostras$Feature)

# Criar o gráfico com apenas os 10 m/z mais importante eabsa jd as

ggplot(top100_loadings_amostras, aes(x = PC1, y = PC2, color = Feature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               size = 1) +
  scale_color_manual(values = rainbow(length(top100_loadings_amostras$Feature))) +  # Escolhe cores para cada Feature
  theme_minimal() +
  ggtitle("Loadings - Top 20 m/z (PC1 vs PC2) without QC") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings") +
  theme(legend.position = "right",    # Mantém a legenda
        panel.grid = element_blank(), # Remove a grade de fundo
        axis.line = element_line(colour = "black"),  # Adiciona linhas dos eixos
        axis.ticks = element_line(colour = "black"))  # Adiciona as marcações dos eixos

#Gráfico em 3D

p <- ggplot(scores_df_no_qc, aes(x = PC1, y = PC2, z = PC3, color = Cluster)) +
  geom_point(size = 4) +
  labs(title = "PCA 3D (PC1, PC2, PC3)") +
  theme_minimal()

ggplotly(p)

plot_ly(scores_df_no_qc, 
        x = ~PC1, y = ~PC2, z = ~PC3, 
        color = ~Cluster, colors = c("red", "blue", "green", "orange"),
        text = ~Sample, hoverinfo = "text") %>%
  add_markers() %>%
  add_text(text = ~Sample, showlegend = FALSE, textposition = "top center") %>%
  layout(title = "PCA 3D (PC1, PC2, PC3)",
         scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")))


# Supondo que as componentes principais estão em scores_df_qc (excluindo a coluna de identificação da amostra)
pca_data1 <- scores_df_no_qc[, -which(names(scores_df_no_qc) == "Sample")]
View(pca_data1)

# Aplicando t-SNE no PCA
set.seed(42)  # Definindo a semente para reprodutibilidade
tsne_result1 <- Rtsne(pca_data1, dims = 2, perplexity = 10)  # Reduzindo para 2 dimensões

# Criando o dataframe com os resultados do t-SNE
tsne_df1 <- data.frame(
  tSNE1 = tsne_result1$Y[, 1],
  tSNE2 = tsne_result1$Y[, 2],
  Sample = scores_df_no_qc$Sample
)

# Plotando o gráfico
ggplot(tsne_df1, aes(x = tSNE1, y = tSNE2, color = Sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "t-SNE a partir do PCA",
       x = "tSNE1",
       y = "tSNE2")

# Aplicando o k-means nos resultados do t-SNE
set.seed(123)  # Para reprodutibilidade
kmeans_result1 <- kmeans(tsne_result1$Y, centers = 3)  

#Adicionando os resultados do k-means ao dataframe
tsne_df1 <- data.frame(
  tSNE1 = tsne_result1$Y[, 1],
  tSNE2 = tsne_result1$Y[, 2],
  Cluster = factor(kmeans_result1$cluster))

tsne_df1$Sample <- scores_df_no_qc$Sample


# Visualizando com ggplot2 e adicionando elipses
ggplot(tsne_df1, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  geom_point(size = 4) +  # Adiciona os pontos
  stat_ellipse(aes(group = Cluster), level = 0.95, linetype = 2) +  # Adiciona as elipses
  geom_text(aes(label = Sample), vjust = -1, size = 3) +  # Adiciona os nomes das amostras
  theme_minimal() +
  ggtitle("t-SNE com K-means e Elipses de Agrupamento") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  theme(plot.title = element_text(hjust = 0.5))


#PCA pela torra
View(data_for_pca)
feature_matrix_without_cb <- data_for_pca[, -c(1:12)]
View(feature_matrix_without_cb)
feature_matrix_org_values_real <- feature_matrix_without_cb[,-c(25:27)]
View(feature_matrix_org_values_real)


feature_matrix_org_values_real2 <- normalizeBetweenArrays(feature_matrix_org_values_real, method = "cyclicloess")
View(feature_matrix_org_values_real2)

# Aplicar o PCA (sem controle de qualidade)
pca_org_values <- t(feature_matrix_org_values_real2)
pca_result_org <- prcomp(pca_org_values, center = TRUE, scale = TRUE)

# Criar DataFrame com os scores do PCA
scores_df_org <- as.data.frame(pca_result_org$x)
scores_df_org$Sample <- rownames(scores_df_org)
scores_df_org$Group <- factor(rep(1:8, each = 3, length.out = nrow(scores_df_org)))

# Plotar o PCA
ggplot(scores_df_org, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 2) +
  theme_minimal() + 
  ggtitle(label = "PCA") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_org)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_org)$importance[2, 2], 1), "%)"))+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  # Diminui o tamanho do texto da legenda
    legend.title = element_text(size = 10), # Ajusta o tamanho do título da legenda
    legend.background = element_rect(fill = "white", color = "black") # Adiciona a caixa ao redor da legenda
  ) +
  scale_color_manual(
    name = "Sample Group",  
    values = c(
     "pink", "cyan", "brown", "black", "gray", "darkgreen", "darkblue","gold"),
    labels = scores_df_org$Sample[seq(1, length(scores_df_org$Sample), by = 3)] 
)

# Calcular a soma dos quadrados dentro dos clusters para diferentes k

set.seed(123)
wcss <- sapply(1:10, function(k) {
  kmeans(scores_df_org[, c("PC1", "PC2")], centers = k, nstart = 25)$tot.withinss
})

print(wcss)  # Checar se os valores são todos zeros
plot(1:10, wcss, type = "b", pch = 19, col = "blue",
     xlab = "Número de clusters", ylab = "Soma dos Quadrados Dentro dos Clusters (WCSS)",
     main = "Método do Cotovelo")

set.seed(123) # Garantir reprodutibilidade
k <- 3 # Suponha que escolhemos 3 clusters com base no método do cotovelo
kmeans_result <- kmeans(scores_df_org[, c("PC1", "PC2")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_org$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_org, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 4) +  # Diminui o tamanho dos nomes das amostras
  stat_ellipse(geom = "polygon", alpha = 0.2) +
  theme_minimal() +
  ggtitle("K-means Clustering in PCA - Positive Mode Organic Fase - 8 different toasts") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_org)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_org)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c( "magenta", "firebrick","cyan")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  # Reduz o tamanho do texto da legenda
    legend.title = element_text(size = 10), # Ajusta o tamanho do título da legenda
    legend.background = element_rect(fill = "white", color = "black") # Adiciona caixa ao redor da legenda
  )

# Calcular a matriz de distâncias (usando a matriz PCA reduzida)
dist_matrix <- dist(scores_df_org[, c("PC1", "PC2", "PC3", "PC4")])

# Aplicar o método de agrupamento hierárquico (usando "ward.D2")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plotar o dendrograma
plot(hclust_result, main = "Dendrogram - Organic Samples", 
     xlab = "Samples", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 4, border = "turquoise")

#Loadings

loadings_df_org <- as.data.frame(pca_result_org$rotation)  # Loadings dos PCs
loadings_df_org$Feature <- rownames(loadings_df_org)  # Adiciona nomes das variáveis
head(loadings_df_org)

loadings_pc1_pc2_org <- loadings_df_org[, c("PC1", "PC2", "Feature")]

# Selecionar os 100 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2_org$Importance <- abs(loadings_pc1_pc2_org$PC1) + abs(loadings_pc1_pc2_org$PC2)
top100_loadings_org <- loadings_pc1_pc2_org[order(-loadings_pc1_pc2_org$Importance), ][1:5, ]
top100_loadings_org$Feature <- sub("^X", "", top100_loadings_org$Feature)

# Criar o gráfico com apenas os 10 m/z mais importantes
ggplot(top100_loadings_org, aes(x = PC1, y = PC2, color = Feature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               size = 1) +
  scale_color_manual(values = rainbow(length(top100_loadings_org$Feature))) +  # Escolhe cores para cada Feature
  theme_minimal() +
  ggtitle("Loadings - Top 20 m/z - 8 different roasted") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings") +
  theme(legend.position = "right",    
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "black"),  
        axis.ticks = element_line(colour = "black"))  

View(top100_loadings_org)
top100_loadings_org_data <- as.data.frame(top100_loadings_org)
write_xlsx(top100_loadings_org_data, "top20_loadings_mz_for_organic.xlsx")

#Gráfico em 3D

p <- ggplot(scores_df_org, aes(x = PC1, y = PC2, z = PC3, color = Cluster)) +
  geom_point(size = 4) +
  labs(title = "PCA 3D (PC1, PC2, PC3)") +
  theme_minimal()

ggplotly(p)

plot_ly(scores_df_org, 
        x = ~PC1, y = ~PC2, z = ~PC3, 
        color = ~Cluster, colors = c("red", "blue", "green", "orange"),
        text = ~Sample, hoverinfo = "text") %>%
  add_markers() %>%
  add_text(text = ~Sample, showlegend = FALSE, textposition = "top center") %>%
  layout(title = "PCA 3D (PC1, PC2, PC3)",
         scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")))


# Supondo que as componentes principais estão em scores_df_qc (excluindo a coluna de identificação da amostra)
pca_data2 <- scores_df_org[, -which(names(scores_df_org) == "Sample")]
View(pca_data2)

# Aplicando t-SNE no PCA
set.seed(42)  # Definindo a semente para reprodutibilidade
tsne_result2 <- Rtsne(pca_data2, dims = 2, perplexity = 5)  # Reduzindo para 2 dimensões

# Criando o dataframe com os resultados do t-SNE
tsne_df2 <- data.frame(
  tSNE1 = tsne_result2$Y[, 1],
  tSNE2 = tsne_result2$Y[, 2],
  Sample = scores_df_org$Sample
)

# Plotando o gráfico
ggplot(tsne_df2, aes(x = tSNE1, y = tSNE2, color = Sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "t-SNE a partir do PCA",
       x = "tSNE1",
       y = "tSNE2")

# Aplicando o k-means nos resultados do t-SNE
set.seed(123)  # Para reprodutibilidade
kmeans_result2 <- kmeans(tsne_result2$Y, centers = 2)  

#Adicionando os resultados do k-means ao dataframe
tsne_df2 <- data.frame(
  tSNE1 = tsne_result2$Y[, 1],
  tSNE2 = tsne_result2$Y[, 2],
  Cluster = factor(kmeans_result2$cluster))

tsne_df2$Sample <- scores_df_org$Sample


# Visualizando com ggplot2 e adicionando elipses
ggplot(tsne_df2, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  geom_point(size = 4) +  # Adiciona os pontos
  stat_ellipse(aes(group = Cluster), level = 0.95, linetype = 2) +  # Adiciona as elipses
  geom_text(aes(label = Sample), vjust = -1, size = 3) +  # Adiciona os nomes das amostras
  theme_minimal() +
  ggtitle("t-SNE com K-means e Elipses de Agrupamento") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  theme(plot.title = element_text(hjust = 0.5))











#PCA Torra média

View(data_for_pca)
feature_matrix_torramedia <- data_for_pca[, c(7,8,9,10,11,12,22,23,24,25,26,27,28,29,30,37,38,39)]
View(feature_matrix_torramedia)

feature_matrix_matrix_torramedia <- normalizeBetweenArrays(feature_matrix_torramedia, method = "cyclicloess")
View(feature_matrix_matrix_torramedia)

# Aplicar o PCA 
pca_torramedia <- t(feature_matrix_matrix_torramedia)
pca_torramedia <- prcomp(pca_torramedia, center = TRUE, scale = TRUE)

# Criar DataFrame com os scores do PCA
scores_df_torramedia <- as.data.frame(pca_torramedia$x)
scores_df_torramedia$Sample <- rownames(scores_df_torramedia)
scores_df_torramedia$Group <- factor(rep(1:6, each = 3, length.out = nrow(scores_df_torramedia)))

# Plotar o PCA
ggplot(scores_df_torramedia, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 2) +
  theme_minimal() + 
  ggtitle("PCA") +
  xlab(paste0("PC1 (", round(100 * summary(pca_torramedia)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_torramedia)$importance[2, 2], 1), "%)")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  
    legend.title = element_text(size = 10), 
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  scale_color_manual(
    name = "Sample Group",  
    values = c("black", "gray", "darkgreen", "darkblue", "gold", "firebrick")
  )


# Calcular a soma dos quadrados dentro dos clusters para diferentes k
set.seed(123)

# Calcular o WCSS para diferentes valores de k
wcss <- sapply(1:10, function(k) {
  kmeans(scores_df_torramedia[, c("PC1", "PC2")], centers = k, nstart = 25)$tot.withinss
})

print(wcss)  # Checar se os valores são todos zeros
plot(1:10, wcss, type = "b", pch = 19, col = "blue",
     xlab = "Número de clusters", ylab = "Soma dos Quadrados Dentro dos Clusters (WCSS)",
     main = "Método do Cotovelo")

set.seed(123) # Garantir reprodutibilidade
k <- 3 # Suponha que escolhemos 3 clusters com base no método do cotovelo
kmeans_result <- kmeans(scores_df_torramedia[, c("PC1", "PC2")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_torramedia$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_torramedia, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 4) +  # Diminui o tamanho dos nomes das amostras
  stat_ellipse(geom = "polygon", alpha = 0.2) +
  theme_minimal() +
  ggtitle(label = "K-means Clustering in PCA - Positive Mode Organic Fase - 8 different toasts") +
  xlab(paste0("PC1 (", round(100 * summary(pca_torramedia)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_torramedia)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("magenta", "firebrick","gold")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  # Reduz o tamanho do texto da legenda
    legend.title = element_text(size = 10), # Ajusta o tamanho do título da legenda
    legend.background = element_rect(fill = "white", color = "black") # Adiciona caixa ao redor da legenda
  )

# Calcular a matriz de distâncias (usando a matriz PCA reduzida)
dist_matrix <- dist(scores_df_torramedia[, c("PC1", "PC2", "PC3", "PC4")])

# Aplicar o método de agrupamento hierárquico (usando "ward.D2")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plotar o dendrograma
plot(hclust_result, main = "Dendrogram - Organic Samples", 
     xlab = "Samples", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 3, border = "turquoise")

#Loadings

loadings_df_torramedia <- as.data.frame(pca_torramedia$rotation)  # Loadings dos PCs
loadings_df_torramedia$Feature <- rownames(loadings_df_torramedia)  # Adiciona nomes das variáveis
head(loadings_df_torramedia)

loadings_pc1_pc2_torramedia <- loadings_df_torramedia[, c("PC1", "PC2", "Feature")]

# Selecionar os 100 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2_torramedia$Importance <- abs(loadings_pc1_pc2_torramedia$PC1) + abs(loadings_pc1_pc2_torramedia$PC2)
top20_loadings_torramedia <- loadings_pc1_pc2_org[order(-loadings_pc1_pc2_torramedia$Importance), ][1:5, ]
top20_loadings_torramedia$Feature <- sub("^X", "", top20_loadings_torramedia$Feature)

# Criar o gráfico com apenas os 10 m/z mais importantes
ggplot(top20_loadings_torramedia, aes(x = PC1, y = PC2, color = Feature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               size = 1) +
  scale_color_manual(values = rainbow(length(top20_loadings_torramedia$Feature))) +  # Escolhe cores para cada Feature
  theme_minimal() +
  ggtitle("Loadings - Top 20 m/z - 8 different roasted") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings") +
  theme(legend.position = "right",    
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "black"),  
        axis.ticks = element_line(colour = "black"))  

#Gráfico em 3D

p <- ggplot(scores_df_torramedia, aes(x = PC1, y = PC2, z = PC3, color = Cluster)) +
  geom_point(size = 4) +
  labs(title = "PCA 3D (PC1, PC2, PC3)") +
  theme_minimal()

ggplotly(p)

plot_ly(scores_df_torramedia, 
        x = ~PC1, y = ~PC2, z = ~PC3, 
        color = ~Cluster, colors = c("red", "blue", "green", "orange"),
        text = ~Sample, hoverinfo = "text") %>%
  add_markers() %>%
  add_text(text = ~Sample, showlegend = FALSE, textposition = "top center") %>%
  layout(title = "PCA 3D (PC1, PC2, PC3)",
         scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")))


# Supondo que as componentes principais estão em scores_df_qc (excluindo a coluna de identificação da amostra)
pca_data3 <- scores_df_torramedia[, -which(names(scores_df_torramedia) == "Sample")]
View(pca_data3)

# Aplicando t-SNE no PCA
set.seed(42)  # Definindo a semente para reprodutibilidade
tsne_result3 <- Rtsne(pca_data3, dims = 2, perplexity = 5)  # Reduzindo para 2 dimensões

# Criando o dataframe com os resultados do t-SNE
tsne_df3 <- data.frame(
  tSNE1 = tsne_result3$Y[, 1],
  tSNE2 = tsne_result3$Y[, 2],
  Sample = scores_df_torramedia$Sample
)

# Plotando o gráfico
ggplot(tsne_df3, aes(x = tSNE1, y = tSNE2, color = Sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_text(aes(label = Sample), vjust = -1, size = 3) +  # Adiciona os nomes das amostras
  labs(title = "t-SNE a partir do PCA",
       x = "tSNE1",
       y = "tSNE2")

# Aplicando o k-means nos resultados do t-SNE
set.seed(123)  # Para reprodutibilidade
kmeans_result3 <- kmeans(tsne_result3$Y, centers = 3)  

#Adicionando os resultados do k-means ao dataframe
tsne_df3 <- data.frame(
  tSNE1 = tsne_result3$Y[, 1],
  tSNE2 = tsne_result3$Y[, 2],
  Cluster = factor(kmeans_result1$cluster))

tsne_df3$Sample <- scores_df_torramedia$Sample


# Visualizando com ggplot2 e adicionando elipses
ggplot(tsne_df3, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  geom_point(size = 4) +  # Adiciona os pontos
  geom_text(aes(label = Sample), vjust = -1, size = 3) +  # Adiciona os nomes das amostras
  theme_minimal() +
  ggtitle("t-SNE com K-means e Elipses de Agrupamento") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  theme(plot.title = element_text(hjust = 0.5))









#Supondo que eu queria visualizar os espectros os top 10 loadings

sps_sciex <- Spectra(files, source = MsBackendMzR())
length(sps_sciex)
spectraVariables(sps_sciex)
View(sps_sciex)

spectra_df <- data.frame(mz = mz(sps_sciex), intensity = intensity(sps_sciex))
head(spectra_df)
View(spectra_df)

library(dplyr)
spectra_filtered <- spectra_df %>%
  group_by(mz.group) %>%
  slice_tail(n = 1) %>%  # 
  ungroup()
View(spectra_filtered)

#Plotando um espectro aleatorio
plotSpectra(sps_sciex[28246])

# Definir os m/z para filtrar (top loadings)
mz_to_plot <- c(564.4401)

# Filtrando os espectros para m/z específicos
selected_spectra <- sps_sciex[sapply(sps_sciex, function(x) any(mz(x) %in% mz_to_plot))]

#Plotando eles lado a lado
plotSpectra(selected_spectra)

# Usar plotSpectraOverlay para plotar os espectros sobrepostos
colors <- c("bisque", "mediumpurple3", "tomato1", "yellowgreen")

plotSpectraOverlay(selected_spectra,
                   main = "Picos Selecionados sobrepostos",
                   xlab = "m/z", 
                   ylab = "Intensidade",
                   col = colors)

legend("topright", 
       legend = c("mz 1", "mz 2", "mz 3", "mz 4"), 
       col = c(colors), 
       lty = 1, 
       lwd = 2, 
       bty = "n")  

# Supondo que feature_matrix_withqc é o seu conjunto de dados com QC e amostras
# Como exemplo, vou criar uma matriz fictícia com o mesmo formato
# feature_matrix_withqc <- <sua matriz de dados>

amostras <- c(rep("Organico", 9), rep("Tradicional", 3), rep("Organico", 24), rep("Tradicional", 3))

# Imprime o vetor de rótulos
print(amostras)

# Divisão dos dados (verifique se trainIndex está correto)
set.seed(123)
trainIndex <- createDataPartition(amostras, p = 0.66, list = FALSE, times = 1)
X_train <- pca_matrix_no_qc[trainIndex, ]
Y_train <- amostras[trainIndex]
X_test <- pca_matrix_no_qc[-trainIndex, ]
Y_test <- amostras[-trainIndex]

# Verifique e converta Y_train para fator
if (!is.factor(Y_train)) {
  Y_train <- as.factor(Y_train)
}

# Crie o modelo PLS-DA
pls_model <- plsda(X_train, Y_train, ncomp = 10)

# Faça previsões
Y_pred <- predict(pls_model, X_test, ncomp = 10, type = "prob")
Y_pred_class <- apply(Y_pred, 1, which.max)
Y_pred_class <- levels(Y_train)[Y_pred_class]

# Gere a tabela de confusão
conf_matrix <- table(predicted = Y_pred_class, actual = Y_test)
print(conf_matrix)

conf_matrix <- matrix(c(16, 0, 3, 3), nrow = 2, byrow = TRUE,
                      dimnames = list(c("Organico", "Tradicional"),
                                      c("Organico", "Tradicional")))

# Precisão (Organico)
precisao_organico <- conf_matrix[1, 1] / sum(conf_matrix[1, ])
print(paste("Precisão (Organico):", precisao_organico))

# Revocação (Organico)
revocacao_organico <- conf_matrix[1, 1] / sum(conf_matrix[, 1])
print(paste("Revocação (Organico):", revocacao_organico))

# F1-Score (Organico)
f1_organico <- 2 * (precisao_organico * revocacao_organico) / (precisao_organico + revocacao_organico)
print(paste("F1-Score (Organico):", f1_organico))

# Precisão (Tradicional)
precisao_tradicional <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
print(paste("Precisão (Tradicional):", precisao_tradicional))

# Revocação (Tradicional)
revocacao_tradicional <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
print(paste("Revocação (Tradicional):", revocacao_tradicional))

# F1-Score (Tradicional)
f1_tradicional <- 2 * (precisao_tradicional * revocacao_tradicional) / (precisao_tradicional + revocacao_tradicional)
print(paste("F1-Score (Tradicional):", f1_tradicional))

# Como exemplo, vou criar uma matriz fictícia com o mesmo formato
# por torras

View(pca_matrix_no_qc)
print(files)
matrix_plc_torra <- pca_matrix_no_qc[-c(7:9),]
View(matrix_plc_torra)
amostrasorg <- c(rep("Escura", 3), rep("Media", 6), rep("Escura", 9), rep("Media", 9), rep("Clara", 6), rep("Escura", 3))

# Imprime o vetor de rótulos
print(amostrasorg)

# Divisão dos dados (verifique se trainIndex está correto)
set.seed(123)
trainIndex <- createDataPartition(amostrasorg, p = 0.50, list = FALSE, times = 1)
X_train1 <- matrix_plc_torra[trainIndex, ]
Y_train1 <- amostrasorg[trainIndex]
X_test1 <- matrix_plc_torra[-trainIndex, ]
Y_test1 <- amostrasorg[-trainIndex]

# Verifique e converta Y_train para fator
if (!is.factor(Y_train1)) {
  Y_train1 <- as.factor(Y_train1)
}

# Crie o modelo PLS-DA
pls_model1 <- plsda(X_train1, Y_train1, ncomp = 10)

# Faça previsões
Y_pred1 <- predict(pls_model1, X_test1, ncomp = 10, type = "prob")
Y_pred_class1 <- apply(Y_pred1, 1, which.max)
Y_pred_class1 <- levels(Y_train1)[Y_pred_class1]

# Gere a tabela de confusão
conf_matrix_torra <- table(predicted = Y_pred_class1, actual = Y_test1)
print(conf_matrix_torra)

# Precisão (Escura)
precisao_escura <- conf_matrix_torra[1, 1] / sum(conf_matrix_torra[1, ])
print(paste("Precisão (Escura):", precisao_escura))

# Revocação (Escura)
revocacao_escura <- conf_matrix_torra[1, 1] / sum(conf_matrix_torra[1, ])
print(paste("Revocação (Escura):", revocacao_escura))

# F1-Score (Escura)
f1_escura <- 2 * (precisao_escura * revocacao_escura) / (precisao_escura + revocacao_escura)
print(paste("F1-Score (Escura):", f1_escura))

# Precisão (Média)
precisao_media <- conf_matrix_torra[2, 2] / sum(conf_matrix_torra[2, ])
print(paste("Precisão (Media):", precisao_escura))

# Revocação (Media)
revocacao_media <- conf_matrix_torra[2, 2] / sum(conf_matrix_torra[2, ])
print(paste("Revocação (Media):", revocacao_media))

# F1-Score (Media)
f1_media <- 2 * (precisao_media * revocacao_media) / (precisao_media+ revocacao_media)
print(paste("F1-Score (Media):", f1_media))

# Precisão (Clara)
precisao_clara <- conf_matrix_torra[3, 3] / sum(conf_matrix_torra[3, ])
print(paste("Precisão (Clara):", precisao_clara))

# Revocação (Clara)
revocacao_clara <- conf_matrix_torra[3, 3] / sum(conf_matrix_torra[3, ])
print(paste("Revocação (Clara):", revocacao_clara))

# F1-Score (Clara)
f1_clara <- 2 * (precisao_clara * revocacao_clara) / (precisao_clara + revocacao_clara)
print(paste("F1-Score (Clara):", f1_clara))

#MONA POSITIVO


