# David Rius Carretero
# Regresion de Procesos Gaussianos: 2 aplicaciones actuariales

# Librerias a utilizar

library(MASS) # Simulaciones distribución normal multivariante
library(ChainLadder) #Visualizacion grafica y calculo CL
library(dplyr) # Manipulacion de datos
library(rstan) #inferencia bayesiana
library(bayesplot) #visualizacion grafica 
library(lattice) #visualizacion grafica
library(plot3Drgl)#visualizacion grafica
library(plot3D)#visualizacion grafica
library(fuzzyjoin) #join aproximado
library(plgp) #funcion distance
library(gridExtra)#visualizacion grafica

# remotes::install_github("Ractuary/casdata")
library(casdata) #base datos provisiones

########################## BBDD PROVISIONES ###################################

########### paqueteria con triangulos enteros ########################
# proviene de RActuary )paqueterias y bbdd para proyectos actuariales) #

# comauto
# View(comauto)
# unique(comauto$GRNAME)
# data("comauto")
bbdd <- comauto %>% 
  filter(GRCODE  == 388) %>% # y/o ==353
  select(AccidentYear, DevelopmentLag, IncurLoss_C, BulkLoss_C) %>% 
  mutate(inc.paid = IncurLoss_C-BulkLoss_C) %>% 
  data.frame()
colnames(bbdd) <- c("originf","dev","IncurLoss_C","BulkLoss_C",'inc.paid')

# bbdd <- comauto %>% 
#   filter(GRCODE  == 'Canal Ins Co Grp') %>% 
#   select(AccidentYear, DevelopmentLag, CumPaidLoss_C) %>% 
#   data.frame()
# colnames(bbdd) <- c("originf","dev", 'inc.paid')

bbdd$originf <- as.factor(bbdd$originf)
bbdd$dev <- as.numeric(bbdd$dev)
bbdd$CumPaidLoss <- as.numeric(bbdd$inc.paid)


# Pasamos nuestros datos a una matriz triangular

# De la base de datos CLAIMS extraeremos las predicciones
n <- 10

data <- data.frame(originf = factor(rep(1988:1997, n:1)),
                   dev = sequence(n:1))
data$dev <- as.numeric(data$dev)

Claims <- left_join(data, bbdd, by= c("originf","dev"))
# Claims$inc.paid <- ave(Claims$inc.paid, Claims$originf, FUN=cumsum)
Claims_aux <- data.frame(as.numeric(Claims$originf),Claims$dev,Claims$inc.paid)
colnames(Claims_aux) <- c("originf","dev", 'inc.paid')

########################## REAL OBSERVADO #####################################
# bbdd$inc.paid <- ave(bbdd$inc.paid, bbdd$originf, FUN=cumsum)
cum.triangle <- with(bbdd, {
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
plot(my_triangle, lattice = TRUE, main="Total Provisiones reales acumuladas observadas: Cod.id: 388")

# latest paiyments
latest.paid <- cum.triangle[row(cum.triangle) == n - col(cum.triangle) + 1]
latest.paid
# Ultimate payments
ultimate.paid <- cum.triangle[,n]
round(ultimate.paid, 0)

# Reserves (by origin period)
round(ultimate.paid - rev(latest.paid), 0)

# Reserves (total)
Reserva_total_observado <- round(sum(ultimate.paid - rev(latest.paid)), 0)
Reserva_total_observado


########################## BBDD TAR VIDA-RIESGO PASEM 2020 UNI 2o orden ########

#Garantía afectada: Fallecimiento por cualquier causa

# Tablas PASEM2020_General_2ndo.orden para los seguros de vida-riesgo, excluido el seguro de decesos, 
# que se recogen en el anexo 1.2. El a?o central base de estas tablas es el a?o natural 2019.

lpasem20ux2 <- c(1000000.0000,998294.9800,998172.9796,998071.7897,997989.2776,997922.0762,997865.6898,997816.1495, 
                 997770.1353, 997725.9314,997681.9204,997635.5265,997583.8501,997522.4432,997446.8016,997352.3745,
                 997235.1733,997093.8154,996929.4620,996768.3741,996606.4382,996442.2636,996275.0989,996103.9219,
                 995928.9431,995750.4300,995569.2640,995387.3690,995206.4395,995027.9277,994852.8305,994679.9432,
                 994507.6367,994332.8940,994150.9921,993955.6374, 993740.0856,993499.8624,993231.9956,992935.1680,
                 992608.3522,992249.0964,991853.6040,991396.1530,990861.2715,990234.6901,989503.4373,988652.7000,
                 987663.1791,986518.5574,985199.7037,983689.9765,982028.0109,980205.6311, 978216.0219,976052.8141,
                 973709.0689,971176.2670,968444.0306,965500.0795,962330.2539,958918.7854,955248.2256, 951298.6781,
                 947047.4077,942468.2475,937531.0709,932201.1720,926437.7495,920192.3061,913407.1144,906013.2149,
                 897928.8256,889058.7631,879296.4481,868527.2726,856631.7741,843485.7437,828956.4785,812899.3402,
                 795155.4903,775554.0347,753899.6505,729998.4232,703672.3834,674776.1679,643216.4272,608972.6978,
                 572116.4980,532829.3492, 491418.2697,448322.2945,404108.8986,359453.6736,315109.4793,271866.1309,
                 230506.7036,191768.0144,156296.9004, 124603.8172, 97020.8579, 73674.4719, 55012.3539, 38907.0536,
                 25233.3969, 14272.7591,6550.7101,2236.8656,505.9667,63.9154, 3.5218, 0.0000, 0.0000, 0.0000, 
                 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000)

####################### listado de asegurados random ##########################

# Creamos el listado de asegurados random
id_aseg <- c("aseg_1","aseg_2","aseg_3","aseg_4","aseg_5")
birth_dates <- c("07-09-1981","01-04-1959","30-01-1971","15-07-1992","14-12-1969")
Capital <- c(25000,100000,50000,35000,150000)
# Prepamos formatos adecuados
lista_aseg <- data.frame(cbind(id_aseg, birth_dates, Capital))
lista_aseg$birth_dates <- as.Date(lista_aseg$birth_dates,format="%d-%m-%Y")
lista_aseg$Capital<- as.numeric(lista_aseg$Capital)

