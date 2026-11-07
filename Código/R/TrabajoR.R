####################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R       ##################
#####################                                       			##################
####################################################################################

# Claudia Peña Vega
# 4º GBME Grupo 2 TRAB 1
# NP: 154119


# Comprobación de que estamos en el directorio de trabajo correcto:
getwd()


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)
archivo<-read.table("archivo-entrada.txt",header=TRUE, sep="\t")  # Lee el archivo

head(archivo)     # Muestra las primeras 6 filas
tail(archivo)     # Muestra las últimas 6 filas
dim(archivo)      # Da las dimensiones de la tabla (619 filas y 7 columnas)

head(archivo,10)  # Muestra las primeras 10 filas
tail(archivo,3)   # Muestra las últimas 3 filas


##################################################################################################################################


### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

archivo$media_control <- apply(archivo[,2:4],1,mean)     # Media control
archivo$de_control <- apply(archivo[,2:4],1,sd)          # DE control
archivo$media_condicion <- apply(archivo[,5:7],1,mean)   # Media condicion
archivo$sd_condicion <- apply(archivo[,5:7],1,sd)        # DE condicion

# Comprobación:
head(archivo)  # Muestra las primeras 6 filas
dim(archivo)   # Muestra las dimensiones de la tabla: 619 filas y 11 columnas (se suman las 4 nuevas)


##################################################################################################################################


### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

rango <- range(c(archivo$media_control, archivo$media_condicion))  # Evalúa el valor máximo y mínimo para ajustar los ejes y que los histogramas sean comparables

par(mfrow=c(1,2))                             # Coloca los histogramas en la misma imagen en dos paneles uno al lado de otro (1 fila y 2 columnas)
hist(archivo$media_control,                   # Histograma de la media del control
     main = "Control", col = "thistle",       # Título y color
     xlab = "Media", ylab = "Frecuencia")     # Títulos de los ejes

hist(archivo$media_condicion,                 # Histograma de la media de la condición
     main = "Condición", col = "tomato3",     # Título y color
     xlab = "Media", ylab = "Frecuencia")     # Títulos de los ejes
par(mfrow=c(1,1))                             # Devuelve la configuración de la ventana de gráficos a 1 panel para que los siguientes salgan bien

# Para visualizar mejor podemos aplicar la transformación logarítmica mediante log2().
# Sumamos 1 porque hay valores = 0 y cualquier log de 0 es -infinito (no existe, perderíamos datos)
# Añadimos los valores como columnas nuevas en la tabla:
archivo$log_control   <- log2(archivo$media_control+1)    # Añadimos el log2 de la media del control +1
archivo$log_condicion <- log2(archivo$media_condicion+1)  # Añadimos el log2 de la media de la condición +1

# Volvemos a representar los histogramas con la transformación logarítmica:
rango <- range(c(archivo$media_control, archivo$media_condicion)) 

par(mfrow=c(1,2))                                    
hist(archivo$log_control,                          # Histograma de log2 media del control + 1
     main = "Control", col = "thistle",            # Título y color
     xlab = "log2 (Media)", ylab = "Frecuencia")   # Títulos de los ejes
 
hist(archivo$log_condicion,                        # Histograma de log2 media de la condición + 1
     main = "Condición", col = "tomato3",          # Título y color
     xlab = "log2 (Media)", ylab = "Frecuencia")   # Títulos de los ejes
par(mfrow=c(1,1))


##################################################################################################################################


### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

boxplot(archivo$log_control, archivo$log_condicion,    # Boxplot para control y condición       
        names = c("Control", "Condición"),             # Nombres eje X
        col = c("thistle", "tomato3"),                 # Colores
        main = "Boxplot", ylab = "log2(Media)",        # Título y nombre eje Y
        ylim = c(0, 20))                               # Límite eje Y para ajustar visualización y que no se tape con la leyenda
legend("topleft", legend = c("Control", "Condición"),  # Leyenda: localización (arriba izq.), títulos y
       fill = c("thistle", "tomato3"))                 # colores correspondientes


##################################################################################################################################


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

shapiro.test(archivo$media_control)     # Test de normalidad control
shapiro.test(archivo$media_condicion)   # Test de normalidad condición

# El test de normalidad es muy significativo en ambos casos; estamos fuera de la norma, ninguno sigue una distribución normal
# Necesitamos un test no paramétrico -> test de Wilcoxon-Mann-Whitney porque son datos independientes no pareados:

archivo$p_value <- apply(archivo[, 2:7], 1,                                  # Creamos la columna p value. A los datos de las columnas 2 a 7
                         function(f) wilcox.test(f[1:3], f[4:6])$p.value)    # aplicamos el test por grupos (control columnas 2,3,4 y condición columnas 5,6,7)

head(archivo)     # Comprobación de que se ha añadido correctamente la nueva columna


##################################################################################################################################


### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Si no te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

significativo_005 <- archivo[archivo$p_value <= 0.05,]  # Flitra la tabla quedándose con las filas con p value < 0,05
nrow(significativo_005)                                 # Contamos las filas que quedan (genes significativos)

# Con un punto de corte para el p value de 0,05 no tenemos ningún gen significativo. 

# Justificación:
# n = 3 en cada grupo suele ser el mínimo de réplicas requeridas para tener potencia estadística, 
# pero es muy reducido para poder encontrar diferencias significativas con el test empleado. 
# El test combina y ordena las réplicas para evaluar las posiciones (rangos)
# Con 3 controles y 3 condiciones solo existen 20 combinaciones posibles de orden.
# El p value es la probabilidad de que salga uno de los dos extremos posibles 
# separando los 3 controles por un lado y los 3 valores de la condición por otro:
# Ct, Ct, Ct, Cn, Cn, Cn o viceversa (Ct = control, Cn = condición)
# Podría salir significativo con un p value de 0,05 para contraste unilateral porque veríamos uno de los extremos (1/20) pero aquí
# necesitamos usar contraste bilateral (hipótesis de dos colas) para ver tanto los genes sobreexpresados como los reprimidos (por encima o por debajo)
# El p-value más pequeño que este test puede calcular matemáticamente incluyendo ambos extremos es 0.1 (2/20).

# Aumentando el tamaño muestral podría solucionarse, pero como los datos que tenemos son estos, 
# probamos a reducir la restricción del p value a 0,1 a ver si encontramos algo con lo que trabajar.

significativo_01 <- archivo[archivo$p_value <= 0.1, ]
nrow(significativo_01)

# Ahora tenemos 429 genes diferencialmente expresados para un p value menos restrictivo (0,1).
# Comprobamos cuántos están sobreexpresados y cuántos reprimidos:

archivo$log2FC <- rowMeans(archivo[, 5:7]) - rowMeans(archivo[, 2:4])  # Añadimos columna log2FC (Fold Change) = media tratamiento - media control

# Añadimos columna que clasifica según p-valor y cambio mínimo (umbral para el fold change, FC):
archivo$estado <- "no significativo"
archivo$estado[which(archivo$p_value <= 0.1 & archivo$log2FC >= 1)]  <- "sobreexpresado"
archivo$estado[which(archivo$p_value <= 0.1 & archivo$log2FC <= -1)] <- "reprimido"

table(archivo$estado)  # Crea una tabla con el nº de genes en cada categoría (no significativo, sobreexpresado o reprimido)

# Hay 190 genes no significativos, 360 reprimidos y 69 sobreexpresados


##################################################################################################################################


### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

# Primero generamos una matriz con los genes diferencialmente expresados:
filtrado <- as.matrix(significativo_01[, 2:7])
rownames(filtrado) <- significativo_01[, 1]

# Luego creamos los dendrogramas haciendo un cluster de las columnas y los genes:
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtrado))))) # Traspone la tabla, correlaciona columnas y transforma en distancia para clusterizar -> dendrograma
colv = as.dendrogram(hclust(as.dist(1-cor(filtrado))))    # Lo mismo pero sin trasponer

# Generamos el heatmap:
heatmap(filtrado, Rowv=rowv, Colv=colv,               # Colv y Rowv = dendogramas
        cexCol=0.7, labRow=FALSE,                     # Tamaño de letra del eje x, quita nombres de los genes (no caben)
        col=hcl.colors(50, palette = "Blue-Red 3"),   # Color: paleta de 50 tonos (rojo = sobreexpresado / azul = reprimido)
        margins = c(8, 2))                            # Ajusta márgenes


# Otro heatmap más elaborado con leyenda.
# Necesitamos instalar los paquetes gplots y RColorBrewer:
install.packages("gplots")		        
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

# Creamos el heatmap:
heatmap.2(filtrado, Rowv=rowv, Colv=colv, cexCol=0.7,
          col = rev(redblue(256)),                     # Invierte los colores para que sea rojo sobreexpresado (++) y azul reprimido (--)
          scale = "row",                               # Escala por filas calculando z score = (valor - media) / DE
          labRow=FALSE, margins = c(8, 2))             # Quita nombres de genes y ajusta márgenes


# Si lo queremos más limpio sin las líneas cian de la figura y la leyenda:
heatmap.2(filtrado, Rowv=rowv, Colv=colv, cexCol=0.7,
          col = rev(redblue(256)),                     # Invierte los colores para que sea rojo sobreexpresado (++) y azul reprimido (--)
          scale = "row",                               # Escala por filas calculando z score = (valor - media) / DE
          trace = "none", density.info = "none",       # Limpia la figura quitando las líneas cian del heatmap y la leyenda respectivamente
          labRow=FALSE, margins = c(8, 2))             # Quita nombres de genes y ajusta márgenes
