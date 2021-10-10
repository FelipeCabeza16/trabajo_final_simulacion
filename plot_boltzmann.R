<<<<<<< HEAD
#install.packages('transformr')

=======
>>>>>>> bb5e9a3d59d0cab4bd5921c96c14b0528d69b910
library(ggplot2)
library(gganimate)
library(magick)


#cambiar la ruta absoluta por la de su sistema, chmod 777
# y las dos simulaciones están en los datasets (particles_n_n.csv) o (particles_2_box.csv)
<<<<<<< HEAD
par = read.csv("D:\\Proyects\\R\\simulacion\\trabajo_final_simulacion\\particles_2_box.csv")
=======
par = read.csv(file = "/home/felipe/Documentos/Programming/R/Simulación/Trabajo_Final/Trabajo_final/dataset/particles_n_n.csv")

>>>>>>> bb5e9a3d59d0cab4bd5921c96c14b0528d69b910
n_par = 1000
escalaR <- 1e7
escalaT <- 1e9
mu <- 1.660e-24
atomic_mass <- 1.00784
m <- mu*atomic_mass
k_botlz <- 1.380649e-23
Temp <- 273.15

tam_box <- 1
r <- 120e-12*escalaR
#Velocidad de la simulación
totalv <- sqrt((par$Vx^2)+(par$Vy^2))
totalv <- tail( totalv, n_par)

<<<<<<< HEAD
x_graf <- seq(min(totalv*(escalaT/escalaR)),max(totalv*(escalaT/escalaR)),by=0.01)

=======
>>>>>>> bb5e9a3d59d0cab4bd5921c96c14b0528d69b910

#Velocidad analítica
x_graf <- seq(min(totalv*(escalaT/escalaR)),max(totalv*(escalaT/escalaR)),by=0.01)
sigm <- sqrt( k_botlz*Temp / m  )
y_graf <-  ( x_graf / (sigm^2) )* exp( - (  x_graf^2 / (2*sigm^2) ) ) 

<<<<<<< HEAD
datosjaja <- data.frame(cbind(x_graf,y_graf))
lines(x_graf,y_graf,type="l",col="red")
=======
data <- data.frame(cbind(x_graf,y_graf))


graph <- ggplot(par, aes((sqrt((par$Vx^2)+(par$Vy^2)))*(escalaT/escalaR))) + 
              geom_histogram(col = "black",fill = "purple", aes(y=..density..)) + 
              view_follow(fixed_x = TRUE, fixed_y = TRUE) + transition_time(par$time)  + 
              geom_line(aes(x=x_graf,y=y_graf),data=data) + 
              ggtitle("Simulación vs solución analítica", ) + 
              labs(x = "Velocidad [m/s]", y = "f(v)", col="Solución analítica") +
              #scale_color_manual(labels = c("f(v)"), values = c("red"))
              theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
                    legend.title = element_blank(), legend.text = element_text(size = 14),
                    plot.title = element_text(size = 12, face = "bold"))



simu <- ggplot(par, aes(x=unlist(par["Sx"]),y=unlist(par["Sy"])))+ 
            geom_point(aes(colour = factor(particle)),size=r,alpha = 0.7, show.legend = F) + theme_bw()+
            scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
            labs(title = 'Time: {frame_time}', x = 'Posición [m]', y = 'Posición [m]') +
            transition_time(unlist(par["time"])) + 
            geom_line(data=data.frame("x"=c(0,0),"y"=c(0,tam_box)),aes_string(x="x",y="y"),linetype="dashed")+
            geom_line(data=data.frame("x"=c(0,tam_box),"y"=c(tam_box,tam_box)),aes_string(x="x",y="y"),linetype="dashed")+
            geom_line(data=data.frame("x"=c(0,tam_box),"y"=c(0,0)),aes_string(x="x",y="y"),linetype="dashed")+
            geom_line(data=data.frame("x"=c(tam_box,tam_box),"y"=c(0,tam_box)),aes_string(x="x",y="y"),linetype="dashed") +
            ggtitle("Colisión de partículas") +
            theme(plot.title = element_text(size = 12, face = "bold"))
>>>>>>> bb5e9a3d59d0cab4bd5921c96c14b0528d69b910


a_gif <- animate(graph, width = 240, height = 240, fps = 20, duration = 15,
                 end_pause = 60, res = 100 ,renderer = gifski_renderer())
b_gif <- animate(simu, width = 240, height = 240, fps = 20, duration = 15,
                 end_pause = 60, res = 100 ,renderer = gifski_renderer())

<<<<<<< HEAD
anim_plot1 <- ggplot(par, aes((sqrt((par$Vx^2)+(par$Vy^2)))*(escalaT/escalaR))) + geom_histogram(col = "black",fill = "blue", aes(y=..density..)) + geom_line(aes(x=x_graf,y=y_graf),data=datosjaja) + view_follow(fixed_x = TRUE) + transition_time(par$time) 
anim_plot1

      
=======
a_mgif <- image_read(a_gif)
b_mgif <- image_read(b_gif)

new_gif <- image_append(c(a_mgif[1], b_mgif[1]))
for(i in 2:100){
  combined <- image_append(c(a_mgif[i], b_mgif[i]))
  new_gif <- c(new_gif, combined)
}

new_gif
anim_save("simulation.gif", animation = new_gif)
>>>>>>> bb5e9a3d59d0cab4bd5921c96c14b0528d69b910
