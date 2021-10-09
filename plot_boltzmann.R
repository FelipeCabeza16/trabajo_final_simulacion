#install.packages('transformr')

library(ggplot2)
library(gganimate)
library(transformr)
library("MASS")

#cambiar la ruta absoluta por la de su sistema
# y las dos simulaciones están en los datasets (particles_n_n.csv) o (particles_2_box.csv)
par = read.csv("D:\\Proyects\\R\\simulacion\\trabajo_final_simulacion\\particles_2_box.csv")
n_par = 1000
escalaR <- 1e7
escalaT <- 1e9
m <- mu*atomic_mass
k_botlz <- 1.380649e-23
Temp <- 273.15


totalv <- sqrt((par$Vx^2)+(par$Vy^2))
totalv <- tail( totalv, n_par)

x_graf <- seq(min(totalv*(escalaT/escalaR)),max(totalv*(escalaT/escalaR)),by=0.01)


sigm <- sqrt( k_botlz*Temp / m  )

y_graf <-  ( x_graf / (sigm^2) )* exp( - (  x_graf^2 / (2*sigm^2) ) ) 

datosjaja <- data.frame(cbind(x_graf,y_graf))
lines(x_graf,y_graf,type="l",col="red")

#primero graficamos el histograma 

X_v <- matrix(x_graf, ncol = 1)
Y_v <- matrix(y_graf, ncol = 1)

anim_plot1 <- ggplot(par, aes((sqrt((par$Vx^2)+(par$Vy^2)))*(escalaT/escalaR))) + geom_histogram(col = "black",fill = "blue", aes(y=..density..)) + geom_line(aes(x=x_graf,y=y_graf),data=datosjaja) + view_follow(fixed_x = TRUE) + transition_time(par$time) 
anim_plot1

      