#TFM
#DAVID RIUS CARRETERO
# Calculo de qx interpoladas con GP 

################################### GP ACUM ##################################

# lx ; estandarizamos las variables de interes (lx)
z_aux = as.numeric(c(lpasem20ux2))
z <- (z_aux-mean(z_aux))/sd(z_aux)
y <- z
sd(y)
x <- c(seq(0,120, by =1))
# Creamos el conjunto de observaciones Input
obs <- data.frame(x = x,
                  y = z,
                  l_x = z_aux)

plot_valores_cohort <- ggplot(obs,aes(x=x,y=l_x)) +
  geom_point() +
  ggtitle("Curva Valores de Cohorte por edad anual: PASEM Unisex 2020, 2o orden")+
  scale_y_continuous(name="valores cohorte l(x)") +  xlab("edad actuarial, x")

plot_valores_cohort
# Creamos el vector con los valores de x interpolados en 12 periodos anuales
x_predict <- round(seq(0,120,by = 1/12),6)

#Definicion parametros para los tipo de kernel
l <- 1
sigma.n <- 0.15
a <- 1 #10 #0.1 #para observar la acentuación de la desviación cuadrática
per <- 1

############ Definición de todas las funciones de los Kernels ################
# Filtramos el kernel de interés a través de # 
# outer = producto exterior de dos matrices
############ GP exponencial ##################################################
SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
# cov <- function(X, Y) outer(X, Y, SE, l)
# COV <- cov(x_predict, x_predict)

############ GP MATTERN 3_2 ##################################################
MT_3_2 <- function(Xi,Xj, l) (1+ (sqrt(3)*abs(Xi - Xj)/ l))*exp(-sqrt(3)*abs(Xi - Xj)/ l)
# cov <- function(X, Y) outer(X, Y, MT_3_2, l)
# COV <- cov(x_predict, x_predict)


############ GP MATTERN 5_2 ##################################################
MT_5_2 <- function(Xi,Xj, l) (1+ ((sqrt(5)*abs(Xi - Xj)/ l)+(5*(Xi - Xj)^2)/ 3*l^2))*exp(-sqrt(5)*abs(Xi - Xj)/ l)
# cov <- function(X, Y) outer(X, Y, MT_5_2, l)
# COV <- cov(x_predict, x_predict)


############ GP rational quadratic ###########################################
QR <- function(Xi,Xj, l, a) (1 + (Xi - Xj)^2 / (2 * a * l^2))^(-a)
# cov <- function(X, Y) outer(X, Y, QR, l, a)
# COV <- cov(x_predict, x_predict)

############ GP periodic kernel################################################
# PK <- function(Xi,Xj, l, per) exp(-(2 * (sin(pi * abs(Xi - Xj) / per))^2) / l^2)
# cov <- function(X, Y) outer(X, Y, PK, l, per)
# COV <- cov(x_predict, x_predict)



############################## Prior #########################################

# Creamos con la función mvrnorm, del paquete MASS, nos permite crear una 
# muestra de datos con distribución normal multivariante con parametros 
# mu y sigma definidos, donde la media es cero y la varianza la matriz Kernel
set.seed(795795)
distr_prior <- mvrnorm(5, rep(0, length=length(x_predict)), COV)

# creamos un data.frame con las simulaciones realizadas
dat_aux <- data.frame(x=x_predict, t(distr_prior)) %>% # unificamos observaciones
  tidyr::pivot_longer(-x, names_to = "rep", values_to = "value") %>% # ordenamos
  mutate(rep = as.numeric(as.factor(rep))) # arreglamos formatos

# Visualización de las distribuciones multivariantes a prior definidas solo por 
# la matriz de covarianzas kernels
prior_plot <- ggplot(dat_aux,aes(x=x,y=value)) +
  geom_line(aes(group=rep), color =  rgb(0.7, 0.1, 0.4), alpha = 0.4) +
  scale_y_continuous(name="output, f(x) = Interpolación de Z(lx)") +  xlab("input, x")

# plot_prior_SE <- prior_plot + ggtitle("GP - Distr. A Priori-Kernel: SE")
# plot_prior_M_3_2<- prior_plot+ ggtitle("GP - Distr. A Priori-Kernel: M. 3/2")
# plot_prior_M_5_2<- prior_plot+ ggtitle("GP - Distr. A Priori-Kernel: M. 5/2")
# plot_prior_QR<- prior_plot+ ggtitle("GP - Distr. A Priori-Kernel: QR")

# extraemos la media y la matriz de covarianzas estimadas de nuestra distribuciones
# a prior
c <- cov(obs$x, obs$x)
det(c) 
cov_xx_inv <- solve(cov(obs$x, obs$x)) # inversa de la matriz covarianzas datos observados
# %*% = Multiplicación Matricial
# Las estimaciones siguen la formulacion reflejada en el paper
Ef <- cov(x_predict, obs$x) %*% cov_xx_inv %*% obs$y  #Estimaci?m media
Cf <- cov(x_predict, x_predict) - cov(x_predict, obs$x)  %*% cov_xx_inv %*% cov(obs$x, x_predict) #Estimaci?n matrices covarianza

############################ Posterior #######################################

########################### sin ruido ########################################
# 1er escenario sin integracion de ruido (media igual a cero) ruido blanco

# Creamos de nuevo nuestras distribuciones multivariantes normales pero con las
# medias y las matrices de covarianzas esperadas como parametros
set.seed(794794)
distr_post <- mvrnorm(50, Ef, Cf)

dat <- data.frame(x=x_predict, t(distr_post)) %>%
  tidyr::pivot_longer(-x, names_to = "rep", values_to = "value") %>% 
  mutate(rep = as.numeric(as.factor(rep)))

# Visualización al 95% de confianza de nuestro proceso gaussiano a posteriori
gp <- data.frame(x = x_predict, Ef = Ef, sigma = 1.96*sqrt(diag(Cf)) )

plot_gp <- ggplot(dat,aes(x=x,y=value)) + 
  geom_line(aes(group=rep), color =  rgb(0.7, 0.1, 0.4), alpha = 0.2) + #REPLICATES
  geom_ribbon(data = gp, aes(x, y = Ef, ymin = Ef - sigma, ymax = Ef + sigma),
              fill="grey", alpha = 0.4) +
  geom_line(dat = gp, aes(x=x,y=Ef), size=1) + # Visualización media
  geom_point(data=obs,aes(x=x,y=y)) +  # Integración Inputs
  scale_y_continuous(name="output, f(x) = Interpolación de Z(lx)") +  xlab("input, x")
  

# plot_gp_SE <- plot_gp + ggtitle("GP A posterori sin ruido-Kernel: SE")
# plot_gp_M_3_2<- plot_gp + ggtitle("GP A posterori sin ruido-Kernel: M. 3/2")
# plot_gp_M_5_2<- plot_gp + ggtitle("GP A posterori sin ruido-Kernel: M. 5/2")
# plot_gp_QR<- plot_gp + ggtitle("GP A posterori sin ruido-Kernel: QR")


yhat_norm <- c(as.numeric(gp$Ef)) # Extracción de la media esperada
yhat <- round(yhat_norm*sd(z_aux)+mean(z_aux),4) # Des-normalizamos los datos y (lx)
data_gp <- data.frame(x_predict, yhat)

# data_gp_SE <- data_gp
# data_gp_M_3_2 <- data_gp
# data_gp_M_5_2 <- data_gp
# data_gp_QR <- data_gp

########################### con ruido #########################################
# Creamos de nuevo nuestras distribuciones multivariantes normales pero con las
# medias y las matrices de covarianzas esperadas como parametros y con 
# integración de ruido (sigmma.n)

cov_xx_inv <- solve(cov(obs$x, obs$x) + sigma.n^2 * diag(1, length(obs$x)))
Ef <- cov(x_predict, obs$x) %*% cov_xx_inv %*% obs$y
Cf <- cov(x_predict, x_predict) - cov(x_predict, obs$x)  %*% cov_xx_inv %*% cov(obs$x, x_predict)

set.seed(793793)
distr_post_noise <- mvrnorm(50, Ef, Cf)

dat <- data.frame(x=x_predict, t(distr_post_noise))
dat <- reshape2::melt(dat, id="x")

gp_noise <- data.frame(x_predict = x_predict, Ef = Ef, sigma = 1.96*sqrt(diag(Cf)) ) #95% confianza

plot_gp_noise <- ggplot(dat,aes(x=x,y=value)) +
  geom_ribbon(data=gp_noise, aes(x=x_predict, y=Ef, ymin=Ef-sigma, ymax=Ef+sigma),
              fill="grey80") + # Variabilidad 
  geom_line(aes(group=variable), alpha=0.3, col=rgb(0.7, 0.1, 0.4, 0.1)) +
  geom_line(data=gp_noise,aes(x=x_predict,y=Ef), size=1) +
  geom_point(data=obs,aes(x=x,y=y)) + scale_y_continuous(lim=c(-2.5,2), 
              name="output, f(x) = Interpolación de Z(lx)") + xlab("input, x")

yhat_norm_noise <- c(as.numeric(gp_noise$Ef))
yhat_noise <- round(yhat_norm_noise*sd(z_aux)+mean(z_aux),4)
data_gp_noise <- data.frame(x_predict, yhat_noise)


# plot_gp_noise_SE <- plot_gp_noise+ggtitle("GP A posterori con ruido-Kernel: SE")
# plot_gp_noise__M_3_2 <- plot_gp_noise+ggtitle("GP A posterori con ruido-Kernel: M. 3/2")
# plot_gp_noise_M_5_2 <- plot_gp_noise+ggtitle("GP A posterori con ruido-Kernel: M. 5/2")
# plot_gp_noise_QR <- plot_gp_noise+ggtitle("GP A posterori con ruido-Kernel: QR")



# data_gp_noise_SE <- data_gp_noise
# data_gp_noise_M_3_2<- data_gp_noise
# data_gp_noise_M_5_2 <- data_gp_noise
# data_gp_noise_QR <- data_gp_noise

# grid.arrange(plot_gp_SE,plot_gp_noise_SE,
#           plot_gp_M_3_2,plot_gp_noise__M_3_2,
#           plot_gp_M_5_2,plot_gp_noise_M_5_2,
#           plot_gp_QR,plot_gp_noise_QR,
#           nrow=2,ncol=1)

# Una vez obtenidos los valores por interpolacion, procedemos a relizar un 
# ejemplo de su aplicación y posbiles beneficios

###################### calculo TAR seguro vida temporal 1 a?o ##################

# extraccion probabilidades fallecimiento (qx)
new_data_gp_aux <- data_gp %>% 
  mutate(qx = (lag(yhat,1)-yhat)/lag(yhat,1))

# reemplazamos el retardo por diferenciación
new_data_gp_aux$qx[1:1428] <-new_data_gp_aux$qx[2:1429]

# filtramos por las edades de interés (15 edad minima trabajar y 70 edad max cobertura)
new_data_gp <- new_data_gp_aux %>% 
    filter(x_predict>=15 & x_predict<=70) %>% 
    arrange(desc(yhat)) #ordenamos por valores qx

################### Bases Técnicas escenario inicial ###########################
ggi <- 0.10 # gastos internos
gge <- 0.10 # gastos externos
rate <- 0.0046 #BOE: tipo interes max para provisiones actualizacion 2022, anual


fecha_efecto <- as.Date("01-05-2022", format="%d-%m-%Y") # efecto calculo

# rescatamos el listado de asegurados random
lista_aseg <- lista_aseg[,0:3]
#Extraemos la edad actuarial tomando en cuenta a?o bisiesto (365.25)
lista_aseg$edad_act <- matrix(as.numeric(round((fecha_efecto-lista_aseg$birth_dates)/365.25,2)), ncol = 1)
colnames(lista_aseg)= c('id_aseg','birth_dates','Capital',"edad_actuarial")
aseg <- data.frame(lista_aseg, ggi, gge, rate) # Unificación BT's por asegurado

# Por otro lado, creamos data.frame con las lx originales (no interpoladas)
Tablas_lx <- data.frame(x = c(seq(0,120, by =1)),
                  y = lpasem20ux2)

Tablas_tasas <- Tablas_lx %>% 
  mutate(qx_real = (lag(y,1)-y)/lag(y,1)) %>%  # extraccion de las qx anuales
  filter(qx_real!=is.na(qx_real))
Tablas_tasas$qx_real[1:110] <-Tablas_tasas$qx_real[2:111]
colnames(Tablas_tasas) <- c("edad_real", "lx_real", "qx_real")

Tablas_int_aux <- new_data_gp %>% 
  mutate(edad_actuarial = round(x_predict,2)) %>%  # redondeamos edad actuarial a 2 decimales
  mutate(edad_real=rep(15:70, each =12, length.out = 661)) %>% # extraccion edad real anual
  data.frame()

# Unificación de las qx reales con las qx interpoladas
Tablas_int_aux_2 <- merge(Tablas_int_aux,Tablas_tasas, by = "edad_real")

# como la variación interpolada de qx será muy baja se ha propuesta el siguiente
# tipo de calculo, donde la qx interpolada tomara el siguiente valor
# sera el (1-1/12)*(qx_interpolada+qx real anual)
Tablas_int <- Tablas_int_aux_2 %>%
  mutate(qx_int = (qx + qx_real*(1-1/12))) %>%
      arrange((edad_actuarial))
      # *(1-1/12))
# filtramos por esos valores a mes vencido
Tabla_Aseg_aux <- fuzzy_left_join(aseg,Tablas_int, by = c("edad_actuarial"="edad_actuarial", 
                                                      "edad_actuarial"="edad_real"),
                              match_fun = list(`<=`, `>=`))
# filtramos por el valor mas cercano, dentro del mes vencido
# el matcheo perfecto no existe por decimales
Tabla_Aseg <- Tabla_Aseg_aux %>% 
  group_by(id_aseg) %>% 
  filter(abs(edad_actuarial.x-edad_actuarial.y) == min(abs(edad_actuarial.x-edad_actuarial.y))) %>% 
  data.frame()
 
# Calculo de primas metodologia interpolada
prima <- Tabla_Aseg %>% 
  mutate(tasa_pura = (qx_int*(1+rate)^(1/12)*1000)) %>% 
  mutate(tasa_tarifa = tasa_pura/(1-ggi-gge)) %>% 
  mutate(prima_neta = tasa_tarifa*Capital/1000)
prima <- data.frame(prima)
prima
sum(prima$prima_neta)

# prima_SE <- sum(prima$prima_neta)
# prima_M_3_2 <- sum(prima$prima_neta)
# prima_M_5_2 <- sum(prima$prima_neta)
# prima_QR <- sum(prima$prima_neta)

# Calculo de primas metodologia anual (real)

aseg$edad_act_real <- round(aseg$edad_actuarial,0)
Tablas_int$edad_act_real <- round(Tablas_int$edad_real,0)


Tabla_Aseg_real <- left_join(aseg,Tablas_int, by = c("edad_act_real"))
Tabla_Aseg_real<-distinct(Tabla_Aseg_real, qx_real, .keep_all = TRUE)

prima_real <- Tabla_Aseg_real %>% 
mutate(tasa_pura_real = (qx_real*(1+rate)^(1/12)*1000)) %>% 
  mutate(tasa_tarifa_real = tasa_pura_real/(1-ggi-gge)) %>% 
  mutate(prima_neta_real = Capital/1000*tasa_tarifa_real)
prima_real <- data.frame(prima_real)
prima_real
sum(prima_real$prima_neta_real)

# Grafica comparativa de probabilidades de muerte

colors <- c("qx'' interpolada GP" = "orange", "qx real anual" = "black")
general <-ggplot(Tablas_int, aes(x = edad_actuarial)) +
  geom_point(aes(y = qx_real, color = "qx real anual"), size = 0.75) +
  geom_line(aes(y = qx_int, color = "qx'' interpolada GP"), size = 1) +
    labs(x = "Edad_Actuarial", y = "(qx)", color = "Legend") +
  scale_color_manual(values = colors)+
  ggtitle("GP- qx interpolada contra qx anual")
#zoom muestro
zoom_int <-ggplot(Tablas_int[450:550,], aes(x = edad_actuarial)) +
  geom_point(aes(y = qx_int, color = "qx'' interpolada GP"), size = 2) +
  # geom_point(aes(y = qx_int[250:550], color = "qx interpolada GP"), size = 2.5) +
  geom_point(aes(y = qx_real, color = "qx real anual"), size = 1.5) +
  labs(x = "Edad_Actuarial", y = "(qx)", color = "Legend") +
  scale_color_manual(values = colors)+
  ggtitle("GP- qx interpolada contra qx anual (tramo ampliado)")
zoom_int
grid.arrange(general,zoom_int, ncol=2)

