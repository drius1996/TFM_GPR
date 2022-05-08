#TFM
#DAVID RIUS CARRETERO
# Calculo de qx interpoladas con GP 

################ grafica a priori y posteriori inferencia bayesiana ##########
Proporcion =seq (0,1, length =1000) # proporcion ratio
x_0 =15 #moneda cara
x_1 =10 # moneda cruz
a =2
b =5
prior = dbeta(Proporcion,a,b)
like=dbeta(Proporcion,x_0+1,x_1+1)
post=dbeta(Proporcion,a+x_0,b+x_1)
plot(Proporcion,post,type="l", ylab="Densidad", lty=2, lwd=3, col=4, main = "Ejemplo Inferencia Bayesiana")
 lines(Proporcion,like,lty=1, lwd=3, col=2)
 lines(Proporcion,prior,lty=3, lwd=3, col=6)
 legend(0.7,4,c("Distr.A priori", "Verosimilitud", "Distr.A posteriori"),
 lty=c(3,1,2) ,lwd=c(3,3,3), col=c(6,2,4))

### graficas datos random normalizados para la visualizacion de los kernels ####
x_predict <- seq(-1, 1, length.out = 201)  # x-coordenadas
l <- 1 # sd(x_predict) # si el kernel sigue una normal, l = sd (x)
sigma.n <- 1
a <- 1
per <- 1

  ############ GP exponencial 
 SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
 cov <- function(X, Y) outer(X, Y, SE, l)
 COV <- cov(x_predict, x_predict)
 
 ############ GP MATTERN 3_2 
 MT_3_2 <- function(Xi,Xj, l) (1+ (sqrt(3)*abs(Xi - Xj)/ l))*exp(-sqrt(3)*abs(Xi - Xj)/ l)
 # cov <- function(X, Y) outer(X, Y, MT_3_2, l)
 # COV <- cov(x_predict, x_predict)
 
 
 ############ GP MATTERN 5_2 
 MT_5_2 <- function(Xi,Xj, l) (1+ ((sqrt(5)*abs(Xi - Xj)/ l)+(5*(Xi - Xj)^2)/ 3*l^2))*exp(-sqrt(5)*abs(Xi - Xj)/ l)
 # cov <- function(X, Y) outer(X, Y, MT_5_2, l)
 # COV <- cov(x_predict, x_predict)
 
 
 ############ GP rational quadratic
 QR <- function(Xi,Xj, l, a) (1 + (Xi - Xj)^2 / (2 * a * l^2))^(-a)
 # cov <- function(X, Y) outer(X, Y, QR, l, a)
 # COV <- cov(x_predict, x_predict)
 
 ############ GP periodic kernel
 PK <- function(Xi,Xj, l, per) exp(-(2 * (sin(pi * abs(Xi - Xj) / per))^2) / l^2)
 # cov <- function(X, Y) outer(X, Y, PK, l, per)
 # COV <- cov(x_predict, x_predict)
 
set.seed(712712)
 values <- mvrnorm(5, rep(0, length=length(x_predict)), COV)
 
 # creamos un data.frame con las simulaciones realizadas
 dat <- data.frame(x=x_predict, t(values)) %>% # unificamos observaciones
   tidyr::pivot_longer(-x, names_to = "rep", values_to = "value") %>% # ordenamos
   mutate(rep = as.numeric(as.factor(rep))) # arreglamos formatos
 
 # Visualización de las distribuciones multivariantes a prior definidas solo por 
 # la matriz de covarianzas kernels
 priori <- ggplot(dat,aes(x=x,y=value)) +
   geom_line(aes(group=rep), color =  rgb(0.7, 0.1, 0.4), alpha = 1) +
   scale_y_continuous(name="Valores random normalizados") +  xlab("input, x")
l
a
per
priori
# SE_l_0.01 <- priori + ggtitle("GP - Distr. A Priori- Kernel: SE - length-scale (l)=0.01")
# SE_l_0.5   <- priori + ggtitle("GP - Distr. A Priori- Kernel: SE - length-scale (l)=0.5")
# SE_l_1 <- priori + ggtitle("GP - Distr. A Priori- Kernel: SE - length-scale (l)=1")
# 
# MT_3_2_l_0.01  <- priori + ggtitle("GP - Distr. A Priori- Kernel: MT_3/2 - length-scale (l)=0.01")
# MT_3_2_l_0.5 <- priori + ggtitle("GP - Distr. A Priori- Kernel: MT_3/2 - length-scale (l)=0.5")
# MT_3_2_l_1 <- priori + ggtitle("GP - Distr. A Priori- Kernel: MT_3/2 - length-scale (l)=1")
# 
# MT_5_2_l_0.01 <- priori + ggtitle("GP - Distr. A Priori- Kernel: MT_5/2 - length-scale (l)=0.01")
# MT_5_2_l_0.5 <- priori + ggtitle("GP - Distr. A Priori- Kernel: MT_5/2 - length-scale (l)=0.5")
# MT_5_2_l_1 <- priori + ggtitle("GP - Distr. A Priori- Kernel: MT_5/2 - length-scale (l)=1")
# 
# grid.arrange(SE_l_0.01, SE_l_0.5,SE_l_1,
#              MT_3_2_l_0.01, MT_3_2_l_0.5, MT_3_2_l_1,
#              MT_5_2_l_0.01, MT_5_2_l_0.5, MT_5_2_l_1,
#              nrow = 3, ncol = 3)
# # 
# SE_l_0.01; SE_l_0.5; SE_l_1
# MT_3_2_l_0.01; MT_3_2_l_0.5; MT_3_2_l_1
# MT_5_2_l_0.01; MT_5_2_l_0.5; MT_5_2_l_1


# QR_l_0.01_aplha_0.5 <- priori + ggtitle("GP - Distr. A Priori- Kernel: QR - length-scale (l)=0.01, alpha=0.5")
# QR_l_0.5_aplha_0.5 <- priori + ggtitle("GP - Distr. A Priori- Kernel: QR - length-scale (l)=0.5, alpha=0.5")
# QR_l_1_aplha_0.5 <- priori + ggtitle("GP - Distr. A Priori- Kernel: QR - length-scale (l)=1, alpha=0.5")
# 
# QR_l_0.01_aplha_1 <- priori + ggtitle("GP - Distr. A Priori- Kernel: QR - length-scale (l)=0.01, alpha=1")
# QR_l_0.5_aplha_1 <- priori + ggtitle("GP - Distr. A Priori- Kernel: QR - length-scale (l)=0.5, alpha=1")
# QR_l_1_aplha_1 <- priori + ggtitle("GP - Distr. A Priori- Kernel: QR - length-scale (l)=1, alpha=1")
# 
# PK_l_0.01_per_0.5 <- priori + ggtitle("GP - Distr. A Priori- Kernel: PK - length-scale (l)=0.01, periodic=0.5")
# PK_l_0.5_per_0.5 <- priori + ggtitle("GP - Distr. A Priori- Kernel: PK - length-scale (l)=0.5, periodic=0.5")
# PK_l_1_per_0.5 <- priori + ggtitle("GP - Distr. A Priori- Kernel: PK - length-scale (l)=1, periodic=0.5")
# 
# PK_l_0.01_per_1 <- priori + ggtitle("GP - Distr. A Priori- Kernel: PK - length-scale (l)=0.01, periodic=1")
# PK_l_0.5_per_1 <- priori + ggtitle("GP - Distr. A Priori- Kernel: PK - length-scale (l)=0.5, periodic=1")
# PK_l_1_per_1 <- priori + ggtitle("GP - Distr. A Priori- Kernel: PK - length-scale (l)=1, periodic=1")
# 

# QR_l_0.01_aplha_0.5;QR_l_0.5_aplha_0.5;QR_l_1_aplha_0.5
# QR_l_0.01_aplha_1;QR_l_0.5_aplha_1;QR_l_1_aplha_1
# PK_l_0.01_per_0.5;PK_l_0.5_per_0.5;PK_l_1_per_0.5
# PK_l_0.01_per_1;PK_l_0.5_per_1;PK_l_1_per_1
# grid.arrange(QR_l_0.01_aplha_0.5,QR_l_0.5_aplha_0.5,QR_l_1_aplha_0.5,
#              QR_l_0.01_aplha_1,QR_l_0.5_aplha_1,QR_l_1_aplha_1,
#              PK_l_0.01_per_0.5,PK_l_0.5_per_0.5,PK_l_1_per_0.5,
#              PK_l_0.01_per_1,PK_l_0.5_per_1,PK_l_1_per_1,
#              nrow = 4, ncol = 3)


################## visualización gráfica introduccion inputs ###################

x_train <- c(-0.81,-0.58,
             -0.2,0.33,0.51,0.9
             )
z <- c(3.21,1.8,
  0.2,-1.9,1.75,2.1
       )

y_train <- z

obs <- data.frame(x = x_train,
                  y = y_train)

x_predict <- seq(-1, 1, length.out = 101)  # x-coordenadas
l <- 0.1
sigma.n <- 0.2


############ Definición de todas las funciones de los Kernels ################
############ GP exponencial ##################################################
SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
cov <- function(X, Y) outer(X, Y, SE, l)
COV <- cov(x_predict, x_predict)

############################## Prior #########################################
set.seed(790790)
values <- mvrnorm(100, rep(0, length=length(x_predict)), COV)

# creamos un data.frame con las simulaciones realizadas
dat <- data.frame(x=x_predict, t(values)) %>% # unificamos observaciones
  tidyr::pivot_longer(-x, names_to = "rep", values_to = "value") %>% # ordenamos
  mutate(rep = as.numeric(as.factor(rep))) # arreglamos formatos

hist_media <- ggplot(dat, aes(x=value)) + 
  geom_density() + geom_density(color="darkblue", fill="lightblue")+
  ggtitle("Curva de Densidad (Distr. Prior)")+
  geom_vline(aes(xintercept=mean(Ef)),
             color="blue", linetype="dashed", size=1)
# Visualización de las distribuciones multivariantes a prior definidas solo por 
# la matriz de covarianzas kernels
priori <- ggplot(dat,aes(x=x,y=value)) +
  geom_line(aes(group=rep), color =  rgb(0.7, 0.1, 0.4), alpha = 0.4) +
  scale_y_continuous(name="valores de y en base al kernel SE") +  xlab("input, x")
plot_prior <- priori + ggtitle("GP - Distribuciones a Priori (Simulaciones)")
plot_prior
# extraemos la media y la matriz de covarianzas estimadas de nuestra distribuciones
# a prior
cov_xx_inv <- solve(cov(obs$x, obs$x)) # inversa de la matriz covarianzas datos observados
# Las estimaciones siguen la formulacion reflejada en el paper
Ef <- cov(x_predict, obs$x) %*% cov_xx_inv %*% obs$y  #Estimaci?m media
Cf <- cov(x_predict, x_predict) - cov(x_predict, obs$x)  %*% cov_xx_inv %*% cov(obs$x, x_predict) #Estimaci?n matrices covarianza

############################ Posterior #######################################

########################### sin ruido ########################################
# 1er escenario sin integracion de ruido (media igual a cero) ruido blanco

# Creamos de nuevo nuestras distribuciones multivariantes normales pero con las
# medias y las matrices de covarianzas esperadas como parametros
set.seed(791791)
values <- mvrnorm(100, Ef, Cf)

dat <- data.frame(x=x_predict, t(values)) %>%
  tidyr::pivot_longer(-x, names_to = "rep", values_to = "value") %>% 
  mutate(rep = as.numeric(as.factor(rep)))

# Visualización al 95% de confianza de nuestro proceso gaussiano a posteriori
gp <- data.frame(x = x_predict, Ef = Ef, sigma = 1.96*sqrt(diag(Cf)) )

post <- ggplot(dat,aes(x=x,y=value)) + 
  geom_line(aes(group=rep), color =  rgb(0.7, 0.1, 0.4), alpha = 0.2) + #REPLICATES
  geom_ribbon(data = gp, aes(x, y = Ef, ymin = Ef - sigma, ymax = Ef + sigma),
              fill="grey", alpha = 0.4) +
  geom_line(dat = gp, aes(x=x,y=Ef), size=1) + # Visualización media
  geom_point(data=obs,aes(x=x,y=y)) +  # Integración Inputs
  scale_y_continuous(name="output, f(x) = y'") +  xlab("input, x")

# plot_train_medias <- post + ggtitle("GP - Distribuciones a Posteriori sin Ruido - Punto 1")

plot_post <- post + ggtitle("GP - Distribuciones a Posteriori sin Ruido")
plot_post

########################### con ruido #########################################
# Creamos de nuevo nuestras distribuciones multivariantes normales pero con las
# medias y las matrices de covarianzas esperadas como parametros y con 
# integración de ruido (sigmma.n)

cov_xx_inv <- solve(cov(obs$x, obs$x) + sigma.n^2 * diag(1, length(obs$x)))
Ef <- cov(x_predict, obs$x) %*% cov_xx_inv %*% obs$y
Cf <- cov(x_predict, x_predict) - cov(x_predict, obs$x)  %*% cov_xx_inv %*% cov(obs$x, x_predict)

set.seed(792792)
values <- mvrnorm(100, Ef, Cf)

dat <- data.frame(x=x_predict, t(values))
dat <- reshape2::melt(dat, id="x")

gp_noise <- data.frame(x_predict = x_predict, Ef = Ef, sigma = 1.96*sqrt(diag(Cf)) ) #95% confianza

post_noise <- ggplot(dat,aes(x=x,y=value)) +
  geom_ribbon(data=gp_noise, aes(x=x_predict, y=Ef, ymin=Ef-sigma, ymax=Ef+sigma),
              fill="grey80") + # Variabilidad 
  geom_line(aes(group=variable), alpha=0.3, col=rgb(0.7, 0.1, 0.4, 0.1)) +
  geom_line(data=gp_noise,aes(x=x_predict,y=Ef), size=1) +
  geom_point(data=obs,aes(x=x,y=y)) + scale_y_continuous(lim=c(-3,4), 
                  name="output, f(x) = y''") + xlab("input, x")

# plot_train_medias_noise <- post_noise + ggtitle("GP - Distribuciones a Posteriori con Ruido - Punto 1")

plot_post_noise <- post_noise + ggtitle("GP - Distribuciones a Posteriori con Ruido")
plot_post_noise


grid.arrange(plot_prior, plot_train_medias, plot_train_medias_noise,
             hist_media,plot_post, plot_post_noise,
             nrow = 2, ncol = 3)
