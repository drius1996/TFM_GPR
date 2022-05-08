#TFM
#DAVID RIUS CARRETERO
# Calculo de Provisiones IBNR con CL y GP

# Claims
Claims 
######################### TRATAMIENTO INPUT ##################################

# Pasamos nuestros datos a una matriz triangular
cum.triangle <- with(Claims, {
  M <- matrix(nrow=n,
              ncol=n,
              dimnames=list(origin=levels(originf),
                            dev=1:n))
  
  M[cbind(originf, dev)] <- inc.paid
  
  M
})
cum.triangle

my_triangle <- as.triangle(cum.triangle) 
plot(my_triangle)
plot(my_triangle, lattice = TRUE)

# Calculamos los ultimos pagos realizados (diagonal inversa)
latest.paid <- cum.triangle[row(cum.triangle) == n - col(cum.triangle) + 1]
latest.paid

# Creamos nueva variable con cambio de nombre
Claims$cum.paid <- cum.triangle[with(Claims, cbind(originf, dev))]

################ metodologia chain-ladder determinista ########################

# En primer lugar extraemos los factores de desarrollo
f <- sapply((n-1):1, function(i) {
  sum(cum.triangle[1:i, n-i+1]) / sum(cum.triangle[1:i, n-i])
})
tail <- 1
f <- c(f, tail)
round(f, 6)

# Extraemos las predicciones con los factores de desarrollo sobre las
# observaciones reales
full.triangle <- cum.triangle
for(k in 1:(n-1)){
  full.triangle[(n-k+1):n, k+1] <- full.triangle[(n-k+1):n, k] * f[k]
}
round(full.triangle, 0)

# Calculamos de nuevo los ultimos pagos a realizar (ultima columna acumulado)
ultimate.paid <- full.triangle[,n]
round(ultimate.paid, 0)

# Calculamos las provisiones totales restando la ultima columna de pagos 
# pendientes con los valores de la diagonal inversa (pagos ya realizados)
round(ultimate.paid - rev(latest.paid), 0)

# Reserves (total)
Reserva_total_CL <- round(sum(ultimate.paid - rev(latest.paid)), 0)
Reserva_total_CL
# Alternativa de cáclulo que nos permite una mejor visualización gráfica del output
# con paquete CL --> aplica una regresi?n lineal simple (MCO)
model.cl <- chainladder(cum.triangle)
model.cl.p <- predict(model.cl)
round(model.cl.p, 0)

my_triangle <- as.triangle(full.triangle) 
plot(my_triangle, lattice = TRUE)

####### Calculo de Provisiones según GP - metodología en Lally y Hartman #######
# Obtención Ouput a través del lenguaje STAN (RStan)

####################### Tratamiento de Input descrito en paper ################

Claims <- Claims_aux
colnames(Claims) <- c('AY', 'dev', 'cum')
Claims$AY <- as.numeric(Claims$AY)
Claims <- Claims[order(Claims$dev),]
dat <- Claims[order(Claims$dev),]
dat$origin <- dat$AY-min(dat$AY)+1
dat <- dat %>% 
  mutate(cum_stand =(cum-mean(cum))/sd(cum)) #estandarizamos importe
dat <-   dat %>% # normalizamos los dos inputs
  mutate(x_2 = (origin-min(origin))/(max(origin)-min(origin))) %>% 
  mutate(x_1 = (dev-min(dev))/(max(dev)-min(dev)))

# mean(dat$cum)
# sd(dat$cum)
# mean(dat$cum_stand)
# sd(dat$cum_stand) # siguen una normal de media cero y sd 1

# creamos un vector con los años a predecir
nyears <- 10
newdat <- data.frame(
  origin=rep(1:10, each=nyears),
  AY=rep(sort(unique(dat$AY)), each=nyears),  
  dev=rep(seq(from=1, to=nyears, by=1), 1)
)
newdat <- merge(dat, newdat, all=TRUE) # unificacmos para predecir NA's
newdat <- newdat[order(newdat$dev),]

newdat_pred <- newdat %>% # normalizamos de nuevo en input existente y a predecir
  mutate(x_2 = (origin-min(origin))/(max(origin)-min(origin))) %>% 
  mutate(x_1 = (dev-min(dev))/(max(dev)-min(dev)))

newdat_pred_2<-  newdat_pred[is.na(newdat_pred$cum),]

input_data <- rbind(dat,newdat_pred_2) # agrupamos por vector ordenado dejando 
# en ultimo lugar los valores NA que queremos predecir


######## Definicion Inferencia Bayesiana por tipologia de Kernel ##############
#Las funciones definidas en Stan son las mismas que facilitan en el codigo
# Lally y Hartman en el anexo Code Stan
########################### mattern 3/2 ######################################
write("functions {
  matrix L_cov_matern32(vector x1, vector x2, real eta_sq,
                        real psi1, real psi2,
                        real a1, real b1, real a2, real b2,
                        real sigma_sq, real delta, int N) {
    // covariance matrix
    matrix[N, N] K;
    // warped inputs
    vector[N] wx1;
    vector[N] wx2;
    // warp input variables
    for(i in 1:N){
      wx1[i] = beta_cdf(x1[i], a1, b1);
    }
    for(i in 1:N){
      wx2[i] = beta_cdf(x2[i], a2, b2);
    }
    // construct covariance matrix
    for (i in 1:(N - 1)) {
      K[i,i] = eta_sq + sigma_sq + delta;
      for (j in (i + 1):N) {
        real dsq;
        dsq = psi1*(wx1[i]-wx1[j])^2 + psi2*(wx2[i]-wx2[j])^2;
        K[i, j] = eta_sq*((1 + sqrt(3)*sqrt(dsq))*exp(-sqrt(3)*sqrt(dsq)));
        K[j, i] = K[i, j];
      }
    }
    K[N,N] = eta_sq + sigma_sq + delta;
    return cholesky_decompose(K);
  }
 }
 data {
  int<lower=1> N; // sample size
  int<lower=1> N1; // training sample size
  int<lower=1> N2; // test sample size
  vector[N1] z1; // target (standardized losses)
  vector[N] x1; // development lag
  vector[N] x2; // accident year
 }
 transformed data {
  vector[N] mu; // mean vector for GP prior
  mu = rep_vector(0, N);
 }
 parameters{
  // bandwidth, signal and noise variance
  vector<lower=0>[2] psi;
  real<lower=0> eta_sq;
  real<lower=0> sigma_sq;
  // input warping parameters
  real<lower=0> a1;
  real<lower=0> b1;
  real<lower=0> a2;
  real<lower=0> b2;
  // test set predictions (lower triangle)
  vector[N2] zmissing;
 }
 transformed parameters {
  // target + predictions
  vector[N] z;
  // computations
  for (n in 1:N1) z[n] = z1[n];
  for (n in 1:N2) z[N1 + n] = zmissing[n];
 }
 model {
  // Cholesky decomposed covariance matrix
  matrix[N,N] L_K;
  L_K = L_cov_matern32(x1, x2, eta_sq,
                       psi[1], psi[2],
                       a1, b1, a2, b2,
                       sigma_sq, 0.01,N);
  // priors on warping functions
  a1 ~ lognormal(0, 0.5);
  b1 ~ lognormal(0, 0.5);
  a2 ~ lognormal(0, 0.5);
  b2 ~ lognormal(0, 0.5);
  // priors on bandwidth, signal and noise variance
  psi ~ gamma(4,4);
  sigma_sq ~ student_t(4,0,1);
  eta_sq ~ student_t(4,0,1);
  // GP
  z ~ multi_normal_cholesky(mu, L_K);
}", "MT3_2_model.stan")


MT3_2model <- "MT3_2_model.stan" # Grabamos las funciones en nuestro environment


########################### mattern 5/2 ######################################


write("functions {
    matrix L_cov_matern52(vector x1, vector x2, real eta_sq,
                          real psi1, real psi2,
                          real a1, real b1, real a2, real b2,
                          real sigma_sq, real delta, int N) {
      // covariance matrix
      matrix[N, N] K;
      // warped inputs
      vector[N] wx1;
      vector[N] wx2;
      // warp input variables
      for(i in 1:N){
        wx1[i] = beta_cdf(x1[i], a1, b1);
      }
      for(i in 1:N){
        wx2[i] = beta_cdf(x2[i], a2, b2);
      }
      // construct covariance matrix
      for (i in 1:(N - 1)) {
        K[i,i] = eta_sq + sigma_sq + delta;
        for (j in (i + 1):N) {
          real dsq;
          dsq = psi1*(wx1[i]-wx1[j])^2 + psi2*(wx2[i]-wx2[j])^2;
          K[i, j] = eta_sq*((1 + sqrt(5)*sqrt(dsq) + (1.666667)*dsq)*exp(-sqrt(5)*sqrt(dsq)));
          K[j, i] = K[i, j];
        }
      }
      K[N,N] = eta_sq + sigma_sq + delta;
      return cholesky_decompose(K);
    }
  }
 data {
  int<lower=1> N; // sample size
  int<lower=1> N1; // training sample size
  int<lower=1> N2; // test sample size
  vector[N1] z1; // target (standardized losses)
  vector[N] x1; // development lag
  vector[N] x2; // accident year
 }
 transformed data {
  vector[N] mu; // mean vector for GP prior
  mu = rep_vector(0, N);
 }
 parameters{
  // bandwidth, signal and noise variance
  vector<lower=0>[2] psi;
  real<lower=0> eta_sq;
  real<lower=0> sigma_sq;
  // input warping parameters
  real<lower=0> a1;
  real<lower=0> b1;
  real<lower=0> a2;
  real<lower=0> b2;
  // test set predictions (lower triangle)
  vector[N2] zmissing;
 }
 transformed parameters {
  // target + predictions
  vector[N] z;
  // computations
  for (n in 1:N1) z[n] = z1[n];
  for (n in 1:N2) z[N1 + n] = zmissing[n];
 }
 model {
  // Cholesky decomposed covariance matrix
  matrix[N,N] L_K;
  L_K = L_cov_matern52(x1, x2, eta_sq,
                       psi[1], psi[2],
                       a1, b1, a2, b2,
                       sigma_sq, 0.01,N);
  // priors on warping functions
  a1 ~ lognormal(0, 0.5);
  b1 ~ lognormal(0, 0.5);
  a2 ~ lognormal(0, 0.5);
  b2 ~ lognormal(0, 0.5);
  // priors on bandwidth, signal and noise variance
  psi ~ gamma(4,4);
  sigma_sq ~ student_t(4,0,1);
  eta_sq ~ student_t(4,0,1);
  // GP
  z ~ multi_normal_cholesky(mu, L_K);
}", "MT5_2_model.stan")


MT5_2model <- "MT5_2_model.stan" # Grabamos las funciones en nuestro environment

########################### SE KERNEL ######################################
write("functions {
    matrix L_cov_sqexp(vector x1, vector x2, real eta_sq,
                          real psi1, real psi2,
                          real a1, real b1, real a2, real b2,
                          real sigma_sq, real delta, int N) {
      // covariance matrix
      matrix[N, N] K;
      // warped inputs
      vector[N] wx1;
      vector[N] wx2;
      // warp input variables
      for(i in 1:N){
        wx1[i] = beta_cdf(x1[i], a1, b1);
      }
      for(i in 1:N){
        wx2[i] = beta_cdf(x2[i], a2, b2);
      }
      // construct covariance matrix
      for (i in 1:(N - 1)) {
        K[i,i] = eta_sq + sigma_sq + delta;
        for (j in (i + 1):N) {
          real dsq;
          dsq = psi1*(wx1[i]-wx1[j])^2 + psi2*(wx2[i]-wx2[j])^2;
          K[i, j] = eta_sq*exp(-dsq);
          K[j, i] = K[i, j];
        }
      }
      K[N,N] = eta_sq + sigma_sq + delta;
      return cholesky_decompose(K);
    }
  }
 data {
  int<lower=1> N; // sample size
  int<lower=1> N1; // training sample size
  int<lower=1> N2; // test sample size
  vector[N1] z1; // target (standardized losses)
  vector[N] x1; // development lag
  vector[N] x2; // accident year
 }
 transformed data {
  vector[N] mu; // mean vector for GP prior
  mu = rep_vector(0, N);
 }
 parameters{
  // bandwidth, signal and noise variance
  vector<lower=0>[2] psi;
  real<lower=0> eta_sq;
  real<lower=0> sigma_sq;
  // input warping parameters
  real<lower=0> a1;
  real<lower=0> b1;
  real<lower=0> a2;
  real<lower=0> b2;
  // test set predictions (lower triangle)
  vector[N2] zmissing;
 }
 transformed parameters {
  // target + predictions
  vector[N] z;
  // computations
  for (n in 1:N1) z[n] = z1[n];
  for (n in 1:N2) z[N1 + n] = zmissing[n];
 }
 model {
  // Cholesky decomposed covariance matrix
  matrix[N,N] L_K;
  L_K = L_cov_sqexp(x1, x2, eta_sq,
                       psi[1], psi[2],
                       a1, b1, a2, b2,
                       sigma_sq, 0.01,N);
  // priors on warping functions
  a1 ~ lognormal(0, 0.5);
  b1 ~ lognormal(0, 0.5);
  a2 ~ lognormal(0, 0.5);
  b2 ~ lognormal(0, 0.5);
  // priors on bandwidth, signal and noise variance
  psi ~ gamma(4,4);
  sigma_sq ~ student_t(4,0,1);
  eta_sq ~ student_t(4,0,1);
  // GP
  z ~ multi_normal_cholesky(mu, L_K);
}" , "SE_model.stan")


SE_model <- "SE_model.stan" # Grabamos las funciones en nuestro environment

# Para poder importar el codigo Stan en R, primero lo definimos como caracteres,
# y con la funcion del paquete stan: stan, la lectura del caracter permite leerlo
# como la compilatoria de la funcion

############################# definicion del input ############################

# la introducción de datos en las funciones creadas en stan deben estar definidas
# en formato lista

acum_data_stan <- list(N = 100, #sample size
                       N1=length(dat$x_2), #Numero FILAS MUESTRA
                       N2=length(newdat_pred_2$x_1), #Numero FILAS A PREDECIR
                       z1=dat$cum_stand, # Observaciones reales de 
                                        # los pagos acumulados ya estandarizados
                       x1=input_data$x_1, # Input con los años de desarrollo
                       x2=input_data$x_2) # Input con los años de origen

############################# fit del modelo ##################################
# Seleccionamos el modelo que queremos calcular, ponenido # en los que no

select_model <- (MT3_2model
                 #MT5_2model
                 #SE_model
)

options(mc.cores = parallel::detectCores()) # deteccion de caracteres
# funcion stan para lectura de codigo stan, nos permite de finir el numero de iteraciones
# warmup = primeras 1000 iteraciones que no se usaran pero se utilizan a modo de 
# calentamiento del modelo
# chains = cadenas a definir de iter= iteraciones ; 4 cadenas de Markov (predefinido)
# como hemos definido 1000 warmup, las iteraciones reales que se realizan son 1000
# iter = 2000 - 1000 warmp = 1000, con lo cual, nos quedan 1000 iteraciones para cada cadena
# que son las asignadas para la muestra de la distribuciones posterior bayesiana
set.seed(797797)
fit1 <- stan(file = select_model, data = acum_data_stan, warmup = 1000, iter = 2000, chains = 4, cores = 4, thin = 1)
monitor(fit1)
fit1
# con la visualización de fit, el output que obtenemos es el siguiente:
# mean es la media esperada, se_mean el error de la media, sd=sd, definidos
# intervalos de confianza para la media 
# también nos devuelve el valor esperado de los hiperparámetros estimados

# extraccion de todas las estimaciones realizadas predichas
posterior <- rstan::extract(fit1)
str(posterior)

# grafico de interés de los valores de todos los parametros optimos
plot(fit1, pars =c("eta_sq", "sigma_sq", "psi[1]", "psi[2]", "a1", "b1", "a2", "b2") )
# tabla con los valores de los parametros optimos con el detalle
print(fit1, c("eta_sq", "sigma_sq", "psi[1]", "psi[2]", "a1", "b1", "a2", "b2", "lp__"),
      probs=c(0.5, 0.75, 0.975))

# también nos permite extraer los valores de los parametros para cada iteracion
ult <- as.data.frame(rstan::extract(fit1, paste0("psi[", 1:2, "]")))
# estadisticos de estos parametros
summary(apply(ult, 1, sum))

# grafico con los valores optimos que toman los parametros a cada iteracion
# plot_SE_PARAM <- rstan::traceplot(fit1, c("eta_sq", "sigma_sq", "psi[1]", "psi[2]", "a1", "b1", "a2", "b2","lp__"), inc_warmup = TRUE)
# plot_M_3_2_PARAM <- rstan::traceplot(fit1, c("eta_sq", "sigma_sq", "psi[1]", "psi[2]", "a1", "b1", "a2", "b2","lp__"), inc_warmup = TRUE)
# plot_M_5_2_PARAM <- rstan::traceplot(fit1, c("eta_sq", "sigma_sq", "psi[1]", "psi[2]", "a1", "b1", "a2", "b2","lp__"), inc_warmup = TRUE)


#################### visualizaci?n output #####################################

# combinacion de todos los chains definidos (4) = cadenas de Markov
sampler_params <- get_sampler_params(fit1, inc_warmup = TRUE)
# estadisticos de la muestra
summary(do.call(rbind, sampler_params), digits = 2)
# cada chain separado
lapply(sampler_params, summary, digits = 2)

# pairs(fit1, pars = c("eta_sq", "sigma_sq", "psi[1]", "psi[2]", "a1", "b1", "a2", "b2","lp__"), las = 1)

# En el gr?fico anterior, la distribuci?n marginal de cada par?metro seleccionado se incluye como un histograma a lo largo 
# de la diagonal. De forma predeterminada, los sorteos con accept_stat__ por debajo de la mediana (tasa de aceptaci?n de 
# propuestas de MCMC) se trazan debajo de la diagonal y aquellos con accept_stat__ por encima de la mediana se trazan sobre 
# la diagonal (esto se puede cambiar usando el argumento de condici?n). Cada cuadrado fuera de la diagonal representa una 
# distribuci?n bivariada de los sorteos para la intersecci?n de la variable de fila y la variable de columna. Idealmente,
# la intersecci?n por debajo de la diagonal y la intersecci?n por encima de la diagonal de las mismas dos variables deber?an
# tener distribuciones que sean im?genes especulares entre s?. Cualquier punto amarillo indicar?a transiciones en las que 
# se alcanz? la profundidad m?xima del ?rbol__, y los puntos rojos indican una transici?n divergente.

# Extraccion Ouput para nuestras Provisiones Estimadas
Y_mean <- rstan::extract(fit1, "z") # extraccion media (coeficientes modelo)
Y_mean_cred <- apply(Y_mean$z, 2, quantile, c(0.025, 0.975)) # extraccion IC media
Y_mean_mean <- apply(Y_mean$z, 2, mean) #calculamos la media total de cada obs
Y_pred <- rstan::extract(fit1, "zmissing") # extraccion solo de la media predicha
Y_pred_cred <- apply(Y_pred$zmissing, 2, quantile, c(0.025, 0.975))  # extr IC media pred
Y_pred_mean <- apply(Y_pred$zmissing, 2, mean) #calculamos la media total de cada obs

#Unificamos en una sola base de datos los valores predichos 
newdat_pred_3 <- cbind(newdat_pred_2, Y_pred_mean, Y_pred_cred[1,], Y_pred_cred[2,])
colnames(newdat_pred_3)<- c("AY", "dev", "origin", "cum", "cum_stand", "x_2","x_1","y", "Y_pred_cred5", "Y_pred_cred95")
# Añadimos nuestras predicciones a la totalidad de nuestra muestra input
dat$y <- dat$cum_stand
dat$Y_pred_cred5 <- dat$cum_stand
dat$Y_pred_cred95 <- dat$cum_stand
dat2 <- rbind(dat, newdat_pred_3)
dat2 <- dat2 %>% 
  arrange(dev)

####################### grafica sin estandarizar los datos ####################

key <- list(rep=FALSE, 
  lines=list(col=c("#00526D", "purple"), type=c("p","l"), pch=1),
  text=list(lab=c("Observation","Mean Estimate")),
  rectangles = list(col=adjustcolor("yellow", alpha.f=0.5), border="grey"),
  text=list(lab="95% Prediction credible interval"))

xyplot( Y_pred_cred5 + Y_pred_cred95 + y + cum_stand ~ dev | factor(AY), 
        data=dat2, as.table=TRUE,
        xlab="dev", ylab="Y_mean", 
        main="Curva de Crecimiento con mu esperada e IC",
        scales=list(alternating=1), layout=c(5,2), key=key,
        panel=function(x, y){
          n <- length(x)
          k <- n/2
          upper <- y[(k/2+1):k]
          lower <- y[1:(k/2)]
          x <- x[1:(k/2)]
          panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                        col = adjustcolor("yellow", alpha.f = 0.5), 
                        border = "grey")
          panel.lines(x, y[(k+1):(k+n/4)], col="purple")
          panel.points(x, y[(n*3/4+1):n], lwd=2, col="#00526D")
        })


###################### grafica con datos estandarizados #################

# Pasamos de nuevo las variables normalizadas al formato origianl
sigma <-c(rep(sd(dat$cum),100))
media <-c(rep(mean(dat$cum),100))

dat2$Y_pred_mean_aux <- dat2$y*sigma+media
dat2$Y_pred_cred5_aux <- dat2$Y_pred_cred5*sigma+media
dat2$Y_pred_cred95_aux <- dat2$Y_pred_cred95*sigma+media

key <- list(
  rep=FALSE, 
  lines=list(col=c("#00526D", "purple"), type=c("p","l"), pch=1),
  text=list(lab=c("Observation","Mean Estimate")),
  rectangles = list(col=adjustcolor("yellow", alpha.f=0.5), border="grey"),
  text=list(lab="95% Prediction credible interval"))
plot_output <- xyplot( Y_pred_cred5_aux + Y_pred_cred95_aux + Y_pred_mean_aux + cum ~ dev | factor(AY), 
        data=dat2, as.table=TRUE,
        xlab="dev", ylab="Y_mean", 
        main="Curva de Crecimiento con media esperada e IC; GP - Matern 5/2",
        scales=list(alternating=1), layout=c(5,2), key=key,
        panel=function(x, y){
          n <- length(x)
          k <- n/2
          upper <- y[(k/2+1):k]
          lower <- y[1:(k/2)]
          x <- x[1:(k/2)]
          panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                        col = adjustcolor("yellow", alpha.f = 0.5), 
                        border = "grey")
          panel.lines(x, y[(k+1):(k+n/4)], col="purple")
          panel.points(x, y[(n*3/4+1):n], lwd=2, col="#00526D")
        })

# plot_output_SE <- plot_output
# plot_output_M_3_2 <- plot_output
# plot_output_M_5_2 <- plot_output

# Una vez definido el output, volvemos al formato triangulo para por comparar
# con el observado y el metodo CL

####################### Triangulo GP ########################################

#pasar los datos de nuevo a triangulo
triangle_aux <- dat2 %>%
  select(AY,dev,Y_pred_mean_aux) # seleccionamos solo las variables de interés

inc.triangle <- as.triangle(triangle_aux, origin="AY", dev="dev", "Y_pred_mean_aux")
inc.triangle

# calculo diagonal (pagado hasta dia de hoy)
latest.paid <- inc.triangle[row(inc.triangle) == n - col(inc.triangle) + 1]
latest.paid
# Calculo de la ultima columna, pagos esperados
ultimate.paid <- inc.triangle[,n]
round(ultimate.paid, 0)
# Calculo de Provisiones por año consiguiente
round(ultimate.paid - rev(latest.paid), 0)

# Total de provisiones esperadas por el metodo GP
Reserva_total_GP <- round(sum(ultimate.paid - rev(latest.paid)), 0)

# Reserva_total_GP_SE <-Reserva_total_GP
# Reserva_total_GP_3_2 <- Reserva_total_GP
# Reserva_total_GP_5_2 <- Reserva_total_GP
################### RMSE OUTPUTS #########################
T_OBS <- 100

# RMSE_GP_SE <- ((Reserva_total_GP_SE-Reserva_total_observado)^2)/T_OBS;RMSE_GP_SE
# RMSE_GP_3_2 <- ((Reserva_total_GP_3_2-Reserva_total_observado)^2)/T_OBS;RMSE_GP_3_2
# RMSE_GP_5_2 <- ((Reserva_total_GP_5_2-Reserva_total_observado)^2)/T_OBS;RMSE_GP_5_2
RMSE_CL <- ((Reserva_total_CL-Reserva_total_observado)^2)/T_OBS;RMSE_CL
