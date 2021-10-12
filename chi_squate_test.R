
library("MASS")


directory = dirname(rstudioapi::getActiveDocumentContext()$path)
dataset1="/results/twobox.csv"
par = read.csv(paste(directory,dataset1,sep = ""))


n_par = 1000


#Atomo de Argon
escalaR <- 5e6
escalaT <- 1e9

#Se multiplica por un factor 1x10^7 por la escala de la simulaci�n

r <- 2e-10*escalaR
mu <- 1.660539040e-27
atomic_mass <- 39.948  # masa atomica del Argon
m <- mu*atomic_mass
constR <- 8.314472
k_botlz <- 1.380649e-23
Temp <- 273.15

tam_box <- 1

#Velocidad de la simulaciÃ³n
totalv <- sqrt((par$Vx^2)+(par$Vy^2))
totalv <- tail( totalv, n_par)*(escalaT/escalaR)


teorical <- function(x){
  sigm <- sqrt( k_botlz*Temp / m  )
  y <-  ( x / (sigm^2) )* exp( - (  x^2 / (2*sigm^2) ) ) 
  return(y)
}

pteorical <- function(x){
  sigm <- sqrt( k_botlz*Temp / m  )
  y <-  1 -  exp( - (  x^2 / (2*sigm^2) ) ) 
  return(y)
}
  

number_bins = 100
min_x = min(totalv)
max_x = max(totalv)


# truehist(totalv*(escalaT/escalaR), nbins=number_bins)
# curve( teorical(x), from = min_x, to=max_x, add=TRUE )

probs = hist(totalv, plot=FALSE, breaks=number_bins)


ejex = probs$mids
expected =  teorical(ejex) 
expected = expected / sum(expected)

observed = probs$density
observed = observed / sum(observed)


plot(ejex, cumsum(expected), type="l", col="blue", main="Cumulative Distribution")
lines(ejex, cumsum(observed), col="red")
legend(0, 1, legend=c("expected", "observed") , col=c("blue", "red") , lty=1:2, cex=0.8)


d <- density(totalv, from = min_x, to = max_x, adjust = 2, kernel = "gaussian")

truehist(totalv, nbins=number_bins, main="Density Distribution", col="white")
lines(ejex, teorical(ejex), type="l", col="blue", add=TRUE)
lines(d$x, d$y, type="l", col="red", add=TRUE)

legend(170, 0.015, legend=c("expected", "observed") , col=c("blue", "red") , lty=1:2, cex=0.8)

ks.test(totalv, pteorical)


