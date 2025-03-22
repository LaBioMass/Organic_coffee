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

data_dir <- "E:/UEM/Resultados/Cafe_31_01_25/NEGATIVO/ORGANICO/mzml/R/QC"
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
     main = "BPC - QC in Organic Phase - Negative Mode", 
     xlab = "Retetion Time (s)", 
     ylab = "Intensity", 
     xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))

#Visualizando o branco

data_dir2 <- "E:/UEM/Resultados/Cafe_31_01_25/NEGATIVO/ORGANICO/mzml/R/Blanck"
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
cwp_low <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 5, noise = 1000, fitgauss = TRUE)
cwp_mid <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, noise = 1000, fitgauss = TRUE)
cwp_high <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 20, noise = 1000, fitgauss = TRUE)

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

par(mfrow = c(1, 3))  # Define layout: 1 linha, 2 colunas
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

plot(bpc_def, col = "firebrick", main = "Comparação dos Parâmetros - Noise 1000",
     xlab = "Retention Time (s)", ylab = "Intensity", xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))

plot(bpc_n1, col = "dodgerblue", main = "Comparação dos Parâmetros - Noise 10000",
     xlab = "Retention Time (s)", ylab = "Intensity", xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))

plot(bpc_n2, col = "seagreen", main = "Comparação dos Parâmetros - Noise 50000",
     xlab = "Retention Time (s)", ylab = "Intensity", xlim = c(0, 1200), 
     ylim = c(0, 1.5e6))


#O VALOR DE NOISE PARA DESCARTE DE RUIDO É DE 500, padrão da função
#tempo de retenção (nao esquecer de carregar o setup escolhido)

data_dir_3 <- "E:/UEM/Resultados/Cafe_31_01_25/NEGATIVO/ORGANICO/mzml"

files_amostras <- sort(list.files(data_dir_3, pattern = "\\.mzML", full.names = TRUE))

print(files_amostras)

raw_data_completed <- readMSData(files_amostras, mode = "onDisk")

n <- 3
I <- 10000
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

std_filled <- fillChromPeaks(xset_grouped, param = FillChromPeaksParam(ppm = 10))
matrix_data <- featureValues(std_filled)
View(matrix_data)

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
write_xlsx(isotope_free_matrix, "isotope_free_organic_negative.xlsx")

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

#PENSAR SE DEVE SOMAR, MAS NA REAL, M/Z DOS CONTAMINANTES

# Substituir valores vazios por NA
isotope_mplus[isotope_mplus == ""] <- NA  

# Substituir NA por 0
isotope_mplus[is.na(isotope_mplus)] <- 0  

# Eliminar todas as linhas onde a coluna 59 seja diferente de 0
isotope_mplus <- isotope_mplus[isotope_mplus[, 59] == 0, ]

# Exibir a matriz modificada
print(isotope_mplus)
View(isotope_mplus)

#todos os adutos eliminados

#subir a pasta com contaminantes
# Carregar os dados
mz_isotope <- as.numeric(isotope_mplus[, 1])  # m/z do tratamento estatístico
mz_contaminantes <- Pasta2[, 1]   # m/z dos contaminantes

# Remover valores NA
mz_isotope <- mz_isotope[!is.na(mz_isotope)]
mz_contaminantes <- mz_contaminantes[!is.na(mz_contaminantes)]

# Verificar se ainda há NAs após a limpeza
if (any(is.na(mz_isotope)) || any(is.na(mz_contaminantes))) {
  stop("Existem valores NA nos dados após a limpeza.")
}

# Identificar quais m/z em isotope_mplus estão dentro do intervalo de ±0.005 dos contaminantes
linhas_para_remover <- sapply(mz_isotope, function(x) any(abs(x - mz_contaminantes) <= 0.05))

# Remover as linhas em isotope_mplus que estão dentro do intervalo
isotope_mplus_filtrado <- isotope_mplus[!linhas_para_remover, ]

# Exibir o número de linhas após a remoção
cat("Número de linhas após remoção:", nrow(isotope_mplus_filtrado), "\n")
View(isotope_mplus_filtrado)  # Para visualizar o resultado no RStudio

#ACIMA OS CONTAMINANTES FORAM REMOVIDOS
# Supondo que as 6 primeiras amostras são os brancos

# Supondo que as colunas 10 a 15 são os brancos
nova_matrix <- isotope_mplus_filtrado
blank_values <- nova_matrix[, 10:15]  # Extrair os valores dos brancos
View(blank_values)

# Calcular o maior valor dos brancos para cada m/z (linha)
max_blank_values <- apply(blank_values, 1, max, na.rm = TRUE)

# Verificar se o tamanho de max_blank_values bate com o número de linhas de nova_matrix
if (length(max_blank_values) == nrow(nova_matrix)) {
  
  # Converter as colunas 16 a 57 para numérico
  nova_matrix[, 16:57] <- lapply(nova_matrix[, 16:57], function(x) as.numeric(as.character(x)))
  
  # Subtrair max_blank_values das colunas 16 a 57
  nova_matrix[, 16:57] <- nova_matrix[, 16:57] - max_blank_values
  
  # Substituir NA e valores menores que zero por zero
  nova_matrix[, 16:57][is.na(nova_matrix[, 16:57]) | nova_matrix[, 16:57] < 0] <- 0
  
  # Remover as linhas que possuem apenas zeros entre as colunas 16 e 57
  nova_matrix <- nova_matrix[apply(nova_matrix[, 16:57], 1, function(x) any(x != 0)), ]
  
} else {
  stop("Erro: O número de valores em max_blank_values não corresponde ao número de linhas da matriz.")
}

View(nova_matrix)

#Eliminando manualmente

nova_matrix <- nova_matrix[-c(3,4,6,9,11,12,13,15,16,17,26,37,41,58,56,52,51),]

#normalizacao
normalize_pareto <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)  # Calcula a média da variável
  sd_x <- sd(x, na.rm = TRUE)  # Calcula o desvio padrão
  pareto_scaled <- (x - mean_x) / sqrt(sd_x)  # Aplica a fórmula de Pareto Scaling
  return(pareto_scaled)
}

# Aplicar a normalização Pareto nas colunas 16 a 57
nova_matrix[, 16:57] <- apply(nova_matrix[, 16:57], 2, normalize_pareto)

# Verifique se as colunas foram normalizadas corretamente
head(nova_matrix[, 16:57])

#Trocando nome das amostras
# Criando os nomes fixos para as colunas QC
colunas_qc <- c("QC1", "QC2", "QC3")

# Atribuindo os novos nomes de coluna para as colunas 16 a 18
colnames(nova_matrix)[16:18] <- colunas_qc

# Criando os nomes fixos para as colunas QC
colunas_cb <- c("CBE1", "CBE2", "CBE3","CBM1", "CBM2", "QBM3","CBV1", "CBV2", "QBV3")
colnames(nova_matrix)[19:27] <- colunas_cb

colunas_j <- c("J1", "J2", "J3")
colnames(nova_matrix)[28:30] <- colunas_j

colunas_itf <- c("25.1", "25.2", "25.3","35.1","35.2","35.3","45.1","45.2","45.3","55.1","55.2","55.3","65.1","65.2","65.3","75.1","75.2","75.3","85.1","85.2","85.3","95.1","95.2","95.3")
colnames(nova_matrix)[31:54] <- colunas_itf

colunas_sp <- c("SP1","SP2","SP3")
colnames(nova_matrix)[55:57] <- colunas_sp
View(nova_matrix)

# Extraindo a coluna de m/z (primeira coluna)
mz_values <- nova_matrix[,1]  # valores de m/z

# Selecionando as colunas de 16 a 57 (amostras para PCA)
data_for_pca <- nova_matrix[, 16:57]

# Atribuindo os valores de m/z como nomes das linhas
rownames(data_for_pca) <- mz_values
View(data_for_pca)

#USAR data_for_pca considerando as triplicatas

#APARTIR DAQUI EH A MEDIA DA TRIPLICATA
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
View(data_for_pca)
data_for_pca <- data_for_pca[, -c(10:12)]
pca_matrix_qc <- t(data_for_pca)
pca_result_qc <- prcomp(pca_matrix_qc, center = TRUE, scale. = TRUE)

# Ver o resumo dos scores para todos os PCs
summary(pca_result_qc)

# Criar DataFrame com os scores do PCA
scores_df_qc <- as.data.frame(pca_result_qc$x)
scores_df_qc$Sample <- rownames(scores_df_qc)

par(mfrow = c(1, 1))  # Define layout: 1 linha, 1 colunas

# Criar os grupos que mudam a cada 3 amostras
scores_df_qc$Group <- factor(rep(1:(nrow(scores_df_qc) / 3), each = 3, length.out = nrow(scores_df_qc)), label = group_labels)

group_labels <- c("Grupo A", "Grupo B", "Grupo C", "Grupo D", "Grupo E", 
                  "Grupo F", "Grupo G", "Grupo H", "Grupo I", "Grupo J",
                  "Grupo K", "Grupo L", "Grupo M")

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
      "firebrick", "blue", "green", "purple", "orange", "pink", "cyan", 
      "brown", "magenta", "gray", "darkgreen", "darkblue","gold"))


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
k <- 3 # escolhemos clusters com base no método do cotovelo
kmeans_result <- kmeans(scores_df_qc[, c("PC1", "PC2","PC3","PC4", "PC5","PC6","PC7", "PC8")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_qc$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_qc, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 4) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.90, linetype = 4) +
  theme_minimal() +
  ggtitle(label = "K-means Clustering in PCA (With Quality Control)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_qc)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("sienna2", "navyblue","green")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  # Reduz o tamanho do texto da legenda
    legend.title = element_text(size = 10), # Ajusta o tamanho do título da legenda
    legend.background = element_rect(fill = "white", color = "black") # Adiciona caixa ao redor da legenda
  ) 

colours()

# Calcular a matriz de distâncias (usando a matriz PCA reduzida)
dist_matrix <- dist(scores_df_qc[, c("PC1", "PC2","PC3","PC4", "PC5","PC6","PC7", "PC8")])

# Aplicar o método de agrupamento hierárquico (usando "ward.D2")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plotar o dendrograma
plot(hclust_result, main = "Dendrogram with QC", 
     xlab = "Samples", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 3, border = "steelblue4")

#Loadings

loadings_df <- as.data.frame(pca_result_qc$rotation)  # Loadings dos PCs
loadings_df$Feature <- rownames(loadings_df)  # Adiciona nomes das variáveis
head(loadings_df)

loadings_pc1_pc2 <- loadings_df[, c("PC1", "PC2", "Feature")]

# Selecionar os 10 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2$Importance <- abs(loadings_pc1_pc2$PC1) + abs(loadings_pc1_pc2$PC2)
top10_loadings_qc <- loadings_pc1_pc2[order(-loadings_pc1_pc2$Importance), ][1:10, ]
View(top10_loadings_qc)

top10_loadings_qc$Feature <- sub("^X", "", top10_loadings_qc$Feature)

# Criar o gráfico corrigido
ggplot(top10_loadings_qc, aes(x = PC1, y = PC2, color = Feature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               linewidth = 1) +
  scale_color_manual(values = rainbow(length(top10_loadings_qc$Feature))) +  # Escolhe cores para cada Feature
  theme_minimal() +
  ggtitle("Loadings - Top 10 m/z with QC") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings") +
  theme(legend.position = "right",    
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "black"),  
        axis.ticks = element_line(colour = "black"))  

top10_mz <- data.frame(mz = top10_loadings_qc$mz)

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
  ggtitle("t-SNE Grouping for Positive Mode - Org with QC") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  theme(plot.title = element_text(hjust = 0.5))

#AANALISE SEM O QC

# Normalizar os dados usando LOESS considerando os QC
# Aqui, estamos dizendo que as colunas 7 a 9 são as QC
# Supondo que você tenha uma lista de amostras de controle chamada `qc_samples`
View(data_for_pca)
abc <- data_for_pca[,c(4:39)]
View(abc)

# Remover as linhas da matriz
feature_matrix_without_qc_normalized <- apply(abc, 2, normalize_pareto)
View(feature_matrix_without_qc_normalized)
summary(pca_result_no_qc)

# Aplicar o PCA (sem controle de qualidade)
pca_matrix_no_qc <- t(feature_matrix_without_qc_normalized)
pca_result_no_qc <- prcomp(pca_matrix_no_qc, center = TRUE, scale. = TRUE)

# Criar DataFrame com os scores do PCA
# Definir rótulos dos grupos primeiro
group_labels1 <- c("Grupo B", "Grupo C", "Grupo D", "Grupo E", 
                   "Grupo F", "Grupo G", "Grupo H", "Grupo I", "Grupo J",
                   "Grupo K", "Grupo L", "Grupo M")

scores_df_no_qc <- as.data.frame(pca_result_no_qc$x)
scores_df_no_qc$Sample <- rownames(scores_df_no_qc)
scores_df_no_qc$Group <- factor(rep(1:(nrow(scores_df_no_qc) / 3), each = 3, length.out = nrow(scores_df_no_qc)), labels = group_labels1)

# Cores para os grupos
cores <- c("blue", "green", "purple", "orange", "pink", "cyan", 
           "brown", "yellow", "magenta", "gray", "darkgreen", "darkblue")

# Plotar o PCA
ggplot(scores_df_no_qc, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 4) +  # Adiciona nome dos grupos acima dos pontos
  theme_minimal() + 
  ggtitle("PCA - Positive Mode Organic Fase - without Quality Control (QC)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_no_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_no_qc)$importance[2, 2], 1), "%)")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  scale_color_manual(name = "Sample Group", values = cores)

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
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
  theme_minimal() +
  ggtitle("K-means Clustering in PCA - Organic Positive Mode (Without QC)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_no_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_no_qc)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("blue", "magenta", "firebrick")) +
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
top100_loadings_amostras <- loadings_pc1_pc2_amostral[order(-loadings_pc1_pc2_amostral$Importance), ][1:10, ]
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
  ggtitle("Loadings - Top 10 m/z (PC1 vs PC2) without QC") +
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


feature_matrix_org_values_real <- apply(feature_matrix_org_values_real, 2, normalize_pareto)
View(feature_matrix_org_values_real)

# Aplicar o PCA (sem controle de qualidade)
pca_org_values <- t(feature_matrix_org_values_real)
pca_result_org <- prcomp(pca_org_values, center = TRUE, scale = TRUE)


group_labels2 <- c("Grupo E", "Grupo F", "Grupo G", "Grupo H", "Grupo I", "Grupo J",
                   "Grupo K", "Grupo L")

# Criar DataFrame com os scores do PCA
scores_df_org <- as.data.frame(pca_result_org$x)
scores_df_org$Sample <- rownames(scores_df_org)
scores_df_org$Group <- factor(rep(1:(nrow(scores_df_org) / 3), each = 3, length.out = nrow(scores_df_org)), labels = group_labels2)

# Cores para os grupos
cores3 <- c("firebrick", "wheat2", 
           "maroon4", "yellow", "magenta", "darkkhaki", "darkgreen", "chocolate")

# Plotar o PCA
ggplot(scores_df_org, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 4) +
  theme_minimal() + 
  ggtitle(label = "PCA with Organic Samples") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_org)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_org)$importance[2, 2], 1), "%)"))+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),  # Diminui o tamanho do texto da legenda
    legend.title = element_text(size = 10), # Ajusta o tamanho do título da legenda
    legend.background = element_rect(fill = "white", color = "black") # Adiciona a caixa ao redor da legenda
  ) +
  scale_color_manual(name = "Sample Group", values = cores3)

colours()

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
k <- 2 # Suponha que escolhemos 3 clusters com base no método do cotovelo
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
  scale_color_manual(values = c( "magenta", "firebrick")) +
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

rect.hclust(hclust_result, k = 2, border = "brown")

#Loadings

loadings_df_org <- as.data.frame(pca_result_org$rotation)  # Loadings dos PCs
loadings_df_org$Feature <- rownames(loadings_df_org)  # Adiciona nomes das variáveis
head(loadings_df_org)

loadings_pc1_pc2_org <- loadings_df_org[, c("PC1", "PC2", "Feature")]

# Selecionar os 100 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2_org$Importance <- abs(loadings_pc1_pc2_org$PC1) + abs(loadings_pc1_pc2_org$PC2)
top100_loadings_org <- loadings_pc1_pc2_org[order(-loadings_pc1_pc2_org$Importance), ][1:10, ]
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
  ggtitle("Loadings - Top 10 m/z - 8 different roasted") +
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

feature_matrix_torramedia <- apply(feature_matrix_torramedia, 2, normalize_pareto)
View(feature_matrix_torramedia)

# Aplicar o PCA 
pca_torramedia <- t(feature_matrix_torramedia)
pca_torramedia <- prcomp(pca_torramedia, center = TRUE, scale = TRUE)

# Criar DataFrame com os scores do PCA
scores_df_torramedia <- as.data.frame(pca_torramedia$x)
scores_df_torramedia$Sample <- rownames(scores_df_torramedia)

#ARRUMANDO
group_labels3 <- c("Grupo C", "Grupo D", "Grupo H", "Grupo I", "Grupo K", "Grupo M")


scores_df_torramedia$Group <- factor(rep(1:(nrow(scores_df_torramedia) / 3), each = 3, length.out = nrow(scores_df_torramedia)), labels = group_labels3)

# Cores para os grupos
cores5 <- c("firebrick", "black", 
            "darkblue", "darkkhaki", "darkgreen", "chocolate")


# Plotar o PCA
ggplot(scores_df_torramedia, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 3.5) +
  theme_minimal() + 
  ggtitle("PCA for medium coffee roasts") +
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
    values = cores5)


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
kmeans_result <- kmeans(scores_df_torramedia[, c("PC1", "PC2","PC3","PC4")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_torramedia$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_torramedia, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 4) +  # Diminui o tamanho dos nomes das amostras
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
  theme_minimal() +
  ggtitle(label = "K-means Clustering in PCA - Positive Mode Organic Fase - for medium coffee roasts") +
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
plot(hclust_result, main = "Dendrogram - for medium coffee roasts", 
     xlab = "Samples", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 3, border = "chocolate")

#Loadings

loadings_df_torramedia <- as.data.frame(pca_torramedia$rotation)  # Loadings dos PCs
loadings_df_torramedia$Feature <- rownames(loadings_df_torramedia)  # Adiciona nomes das variáveis
head(loadings_df_torramedia)

loadings_pc1_pc2_torramedia <- loadings_df_torramedia[, c("PC1", "PC2", "Feature")]

# Selecionar os 100 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2_torramedia$Importance <- abs(loadings_pc1_pc2_torramedia$PC1) + abs(loadings_pc1_pc2_torramedia$PC2)
top20_loadings_torramedia <- loadings_pc1_pc2_org[order(-loadings_pc1_pc2_torramedia$Importance), ][1:10, ]
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
  ggtitle("Loadings - Top 10 m/z - for medium coffee roasts") +
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
  Cluster = factor(kmeans_result3$cluster))

tsne_df3$Sample <- scores_df_torramedia$Sample


# Visualizando com ggplot2 e adicionando elipses
ggplot(tsne_df3, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  geom_point(size = 6) +  # Adiciona os pontos
  geom_text(aes(label = Sample), vjust = -1, size = 3.5) +  # Adiciona os nomes das amostras
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95)
  theme_minimal() +
  ggtitle("t-SNE for medium coffee roasts") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  theme(plot.title = element_text(hjust = 0.5))

