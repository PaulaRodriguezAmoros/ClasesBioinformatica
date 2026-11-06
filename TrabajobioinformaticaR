# TRABAJO R - BIOINFORMÁTICA

# colores que uso en todo el trabajo
azul <- "steelblue"
rojo <- "tomato"

### 1. Cargar el archivo
datos <- read.table("archivo-entrada.txt", header = TRUE, sep = "\t", check.names = FALSE)
head(datos)
tail(datos)
dim(datos)   # 619 genes y 7 columnas

# columnas 2-4 = control, columnas 5-7 = condición
control <- datos[, 2:4]
condicion <- datos[, 5:7]

### 2. Media y desviación estándar de cada gen
datos$media_control <- apply(control, 1, mean)
datos$de_control <- apply(control, 1, sd)
datos$media_cond <- apply(condicion, 1, mean)
datos$de_cond <- apply(condicion, 1, sd)

head(datos)

### 3. Histogramas (en una sola imagen)
# uso log2(x+1) porque hay genes con valores muy altos y si no no se ve nada
par(mfrow = c(1, 2))
hist(log2(datos$media_control + 1), col = azul, main = "Media Control",
     xlab = "log2(media + 1)", ylab = "Frecuencia")
hist(log2(datos$media_cond + 1), col = rojo, main = "Media Condición",
     xlab = "log2(media + 1)", ylab = "Frecuencia")
par(mfrow = c(1, 1))

### 4. Boxplot
boxplot(log2(datos$media_control + 1), log2(datos$media_cond + 1),
        names = c("Control", "Condición"), col = c(azul, rojo),
        main = "Control vs Condición", ylab = "log2(media + 1)")
legend("topleft", legend = c("Control", "Condición"), fill = c(azul, rojo))

### 5. Normalidad y comparación de medias
shapiro.test(unlist(control))
shapiro.test(unlist(condicion))
# los dos p-values son mucho menores que 0.05, así que no son normales
# por eso uso el test de Wilcoxon (no paramétrico) en cada gen

pvalue <- c()
for (i in 1:nrow(datos)) {
  x <- as.numeric(control[i, ])
  y <- as.numeric(condicion[i, ])
  pvalue[i] <- suppressWarnings(wilcox.test(x, y)$p.value)
}
datos$pvalue <- pvalue

head(datos)

### 6. Genes diferencialmente expresados
diferenciales <- datos[datos$pvalue < 0.05, ]
nrow(diferenciales)
# no sale ningún gen. Con solo 3 réplicas por grupo, el test de Wilcoxon no puede
# dar un p-value menor de 0.1 (2/20 combinaciones posibles), así que nunca llega a 0.05
# por tanto hay 0 sobreexpresados y 0 reprimidos
min(datos$pvalue)

### 7. Heatmap
# como no hay genes diferenciales, hago el heatmap con los 50 genes con menor p-value
top <- datos[order(datos$pvalue), ][1:50, ]
matriz <- as.matrix(top[, 2:7])
rownames(matriz) <- top[, 1]
matriz <- log2(matriz + 1)

colores <- colorRampPalette(c("blue", "white", "red"))(100)
heatmap(matriz, scale = "row", col = colores,
        ColSideColors = c(rep(azul, 3), rep(rojo, 3)),
        main = "50 genes con menor p-value", cexRow = 0.6)
