rm(list = ls())

#----- Input Data -----#
data <- read.csv(file.choose(),header=T,sep=";",dec=",")
View(data)
str(data)
attach(data)
library(dplyr)
head(data)

#Matriks Data
Y <- as.matrix(data[,2:8])
Y

#Standarisasi Data
y1 <- (Y[,1]-mean(Y[,1]))/sd(Y[,1])
y2 <- (Y[,2]-mean(Y[,2]))/sd(Y[,2])
y3 <- (Y[,3]-mean(Y[,3]))/sd(Y[,3])
y4 <- (Y[,4]-mean(Y[,4]))/sd(Y[,4])
y5 <- (Y[,5]-mean(Y[,5]))/sd(Y[,5])
y6 <- (Y[,6]-mean(Y[,6]))/sd(Y[,6])
y7 <- (Y[,7]-mean(Y[,7]))/sd(Y[,7])
datastandar <- data.frame(x1=c(y1),x2=c(y2),x3=c(y3),x4=c(y4),x5=c(y5),x6=c(y6),x7=c(y7))
#Matriks X
X <- as.matrix(datastandar)
X

#----- MCD -----#
library(rrcov)
set.seed(123)
mcd = covMcd(datastandar, alpha=0.90)
mcd$center    
mcd$cov  

#----- Robust Distance -----#
rd = mahalanobis(datastandar, center = mcd$center,cov = mcd$cov);rd
#Deteksi Outlier
cutoff = qchisq(p = 0.95 , df = ncol(datastandar));cutoff
datastandar[rd > cutoff ,]

#----- Robust SVD -----#
#Nilai Eigen dan Vektor Eigen
Cov_MCD.e = eigen(mcd$cov)
Cov_MCD.e

#Matriks L
Eigen_values = Cov_MCD.e$values;Eigen_values
Eigen_sqrt = sqrt(Eigen_values);Eigen_sqrt
diag_L = Eigen_sqrt;diag_L
L = cbind(diag_L,diag_L,diag_L,diag_L,diag_L,diag_L,diag_L);L
L[1,2:7] = 0
L[2,3:7] = 0
L[3,4:7] = 0
L[4,5:7] = 0
L[5,6:7] = 0
L[5,7] = 0
L[2,1] = 0
L[3,1:2] = 0
L[4,1:3] = 0
L[5,1:4] = 0
L[6,1:5] = 0
L[7,1:6] = 0
L[6,7] = 0
L
colnames(L) = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7")
L

#Matriks A
A = Cov_MCD.e$vectors
A

#Matriks U
##Mencari matriks 1/L
K = 1/L
K[1,2:7] = 0
K[2,3:7] = 0
K[3,4:7] = 0
K[4,5:7] = 0
K[5,6:7] = 0
K[5,7] = 0
K[2,1] = 0
K[3,1:2] = 0
K[4,1:3] = 0
K[5,1:4] = 0
K[6,1:5] = 0
K[7,1:6] = 0
K[6,7] = 0
K

##Hasil Matriks U
U = X %*% A %*% K
U

#----- Identifikasi Persentase Keragaman Data -----#
#Komponen Utama
Ut = t(U);Ut
Tr = mcd$center;Tr
PCr = Ut%*%(X-Tr)
PCr

#Variasi Data
V = Eigen_values
V

#Persentase Keragaman PCr
PV1 = V[1]/sum(V)
PV2 = V[2]/sum(V)
PV3 = V[3]/sum(V)
PV4 = V[4]/sum(V)
PV5 = V[5]/sum(V)
PV6 = V[6]/sum(V)
PV7 = V[7]/sum(V)
data.frame(PV1,PV2,PV3,PV4,PV5,PV6,PV7)

#Persentase Keragaman Kumulatif
PC1 = PV1
PC2 = PC1+PV2
PC3 = PC2+PV3
PC4 = PC3+PV4
PC5 = PC4+PV5
PC6 = PC5+PV6
PC7 = PC6+PV7
data.frame(PC1,PC2,PC3,PC4,PC5,PC6,PC7)

#----- Membentuk Matriks G dan Matriks H -----#
#Matriks G
alpha = 0.5
G = U%*%(L^alpha)
round(G,3)

#Matriks H
Ht = L^(1-alpha)%*%t(A)
H = t(Ht)
round(H,3)

#----- Analisis Biplot -----#
#Matriks G3 dan Matriks H3
G3 = G[,c(1:3)]
round(G3,3)
H3 = H[,c(1:3)]
round(H3,3)

#Visualisasi Biplot 2D
##PC1 dan PC2
G12 <- G[,c(1,2)]
H12 <- H[,c(1,2)] 
biplot(G12,H12,cex=0.8,main="PC1 vs PC2",xlab = "PC 1",ylab="PC 2")
abline(h=0)
abline(v=0)

##PC1 dan PC3
G13 <- G[,c(1,3)]
H13 <- H[,c(1,3)]
biplot(G13,H13,cex=0.8,main="PC1 vs PC3",xlab = "PC 1",ylab="PC 3")
abline(h=0)
abline(v=0)

##PC2 dan PC3
G23 <- G[,c(2,3)]
H23 <- H[,c(2,3)]
biplot(G23,H23,cex=0.9,main="PC2 vs PC3",ylab = "PC 3",xlab="PC 2")
abline(h=0)
abline(v=0)

#Visualisasi Biplot 3D
##Kelompok Kabupaten/Kota
G3_C <- as.data.frame(G3)
G3_C[,4] <- NA
G3_C[,4] <- ifelse((G3_C[,1] > 0 & G3_C[,2] > 0 & G3_C[,3] > 0), 1,
                   ifelse((G3_C[,1] > 0 & G3_C[,2] < 0 & G3_C[,3] > 0), 2,
                          ifelse((G3_C[,1] > 0 & G3_C[,2] > 0 & G3_C[,3] < 0), 3,
                                 ifelse((G3_C[,1] > 0 & G3_C[,2] < 0 & G3_C[,3] < 0), 4,
                                        ifelse((G3_C[,1] < 0 & G3_C[,2] > 0 & G3_C[,3] > 0), 5,
                                               ifelse((G3_C[,1] < 0 & G3_C[,2] > 0 & G3_C[,3] < 0), 6,
                                                      ifelse((G3_C[,1] < 0 & G3_C[,2] < 0 & G3_C[,3] > 0), 7,
                                                             0)))))))
rownames(G3_C) = data[,1]
G3_C[order(G3_C$V4),]
##Visualisasi 3D
library(rgl)
K <- c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10",
       "K11","K12","K13","K14","K15","K16","K17","K18","K19","K20",
       "K21","K22","K23","K24","K25","K26","K27")
a <- c("UPLM1","UPLM2","UPLM3","UPLM4","UPLM5","UPLM6","UPLM7")
plot3d(G3)
points3d(G3+0.09,color=G3_C[,4]+1,size=7)
text3d(G3, texts=K,col=G3_C[,4]+1)
points3d(H3+0.09,color='red', size=6)
text3d(H3, texts=a,col="red")
coords <- NULL
for (i in 1:nrow(H3)) {
  coords <- rbind(coords,rbind(c(0,0,0),H3[i,1:3])) 
}
lines3d(coords,col="red",lwd=2)
grid3d(c("x", "y", "z"))
abclines3d(0,0,0,a=diag(3),col="black",lwd=1)

#----- Identifikasi Hasil PCA Biplot -----#
#Matriks Korelasi Kosinus antar Variabel
korelasi_cosinus <- function(X) {
  library(geometry)
  p <- nrow(X)
  r <- matrix(,nrow=p,ncol=p)
  for (i in 1:nrow(X)) {
    for (j in 1:nrow(X)) {
      r[i,j] <- dot(X[i,],X[j,]) /
        (sqrt(sum(X[i,]^2)) *
           sqrt(sum(X[j,]^2)))
    }
  }
  print(r)
}
korvar3 <- korelasi_cosinus(H3)
korvar_3 <- round(korvar3,3);korvar_3

#Nilai Variabel pada Suatu Objek
cosine_correlation <- function(G, H) {
  # Inisialisasi jumlah objek dan variabel
  n <- nrow(G)  # Jumlah objek (baris di matriks G)
  p <- nrow(H)  # Jumlah variabel (baris di matriks H)
  # Matriks korelasi kosinus, diinisialisasi dengan nilai 0
  rho <- matrix(0, n, p)
  # Loop untuk menghitung korelasi kosinus antara setiap objek dan variabel
  for (i in 1:n) {
    for (j in 1:p) {
      # Menghitung pembilang (dot product antara baris i G dan baris j H)
      numerator <- sum(G[i, ] * H[j, ])
      # Menghitung penyebut (akar dari jumlah kuadrat elemen pada baris i G dan baris j H)
      denominator <- sqrt(sum(G[i, ]^2) * sum(H[j, ]^2))
      # Menghitung korelasi kosinus
      rho[i, j] <- numerator / denominator
    }
  }
  return(rho)
}
obvar3<- cosine_correlation(G3, H3)
obvar_3 <- round(obvar3,3);obvar_3

#----- Analisis Klaster Hierarki -----#
#Menghitung Jarak Cosine
# Fungsi Cosinus Similarity
cosine_similarity <- function(x) {
  # Menghitung dot product antar semua baris
  similarity_matrix <- x %*% t(x)
  # Menghitung norma (panjang vektor) untuk setiap baris
  norms <- sqrt(rowSums(x^2))
  # Membagi dot product dengan hasil kali norma untuk setiap pasangan baris
  similarity_matrix <- similarity_matrix / outer(norms, norms)
  return(similarity_matrix)
}
# Menghitung Cosine Similarity Matriks G3
cosine_sim_G3 <- cosine_similarity(G3); cosine_sim_G3
# Cosine Distance
jarakcosine <- as.dist(1 - cosine_sim_G3)
round(jarakcosine,3)

#Metode Klaster Hierarki
##Single Linkage
single_clus <- hclust(d=jarakcosine, method="single")
plot(single_clus, main="Dendogram Single Linkage")
rect.hclust(single_clus,18)
data.frame(data$KAB.KOTA, cutree(single_clus, k=18))
cophen.single <- cophenetic(single_clus)
print(cor(cophen.single, jarakcosine))

##Complete Linkage
complete_clus <- hclust(d=jarakcosine, method="complete")
plot(complete_clus, main="Dendogram Complete Linkage")
rect.hclust(complete_clus,18)
data.frame(data$KAB.KOTA, cutree(complete_clus, k=18))
cophen.complete <- cophenetic(complete_clus)
print(cor(cophen.complete, jarakcosine))

##Average Linkage
avg_clus <- hclust(d=jarakcosine, method="average")
plot(avg_clus, main="Dendogram average Linkage")
rect.hclust(avg_clus,18)
data.frame(data$KAB.KOTA, cutree(avg_clus, k=18))
cophen.average <- cophenetic(avg_clus)
print(cor(cophen.average, jarakcosine))