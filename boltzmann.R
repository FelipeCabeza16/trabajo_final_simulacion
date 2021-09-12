#Parámetros:
#           time -> Vector de intervalos de tiempo
#           sx -> Vector de posiciones en x
#           sy -> Vector de posiciones en y
#           vx -> Vector de velocidades en x
#           vy -> Vector de velocidades en y
#Return:
#           particles -> dataframe de partículas con los datos ingresados como parámetro,
#                        con un identificador único para cada partícula
init_particles = function(time,sx,sy,vx,vy){
  particles = data.frame(time, c(seq(1,n_particles,by=1)), sx, sy,vx,vy)
  names(particles) <- c('time', 'particle','Sx', 'Sy','Vx','Vy')
  return(particles)
}

#Parámetros: 
#           vx -> Vector de velocidades en x
#           vy -> Vector de velocidades en y
#Return:
#           speeds <- dataframe de velocidades
init_speeds = function(vx,vy){
  speeds = data.frame(c(seq(1,n_particles,by=1)),sx, sy)
  names(speeds) <- c('particle','Vx', 'Vy')
  return(speeds)
}

# Charge libraries:
library(ggplot2)
library(gganimate)
library(gifski)
library(png)

#inicializar los valores 

n_particles = 100

#Atomo de HidrÃ³geno 
#Se multiplica por un factor 1x10^7 por la escala de la simulación
r = 120e-12*1e7
m = 1e-24

#Tamaño de la caja n*n
tam_box = 1

#Temperatura constante
T = 300
FPS = 30
dt = 1/FPS


#Generar de manera aleatoría las posiciones iniciales
sx = runif(n_particles,0,tam_box)
sy = runif(n_particles,0,tam_box)


#Velocidad mÃ¡xima
v_max = 1.5

#Generar de manera aleatoría las velocidades iniciales
thetha = runif(n_particles) * 2 * pi
vx = v_max * cos(thetha)
vy = v_max * sin(thetha)


#Tiempo máximo en seg
time_max = 5
#Diferencial de tiempo 
dt = 0.1


#Parámetros:
#           dt -> Double: Diferencial de tiempo
#           time_max -> Tiempo máximo
#           particles -> dataframe con información de partículas
#Return:
#           particles -> dataframe de partículas con los datos después de realizar
#                        interacciones, ordenado por partícula
forward = function(dt, time_max, particles){
  time = 0
  pos = array()
    while(time <= time_max){
      time = time + dt  
      sx <- sx + vx*dt
      sy <- sy + vy*dt
      for (i in 1:length(sx)){
        for(j in 1:length(sx)){
          # Si dos partículas diferentes se encuentran en el radio una de otra
          if ((abs(sx[j]-sx[i]) <= 4*r) && (abs(sy[j]-sy[i]) <= 40*r) && j!=i){
            print("Choque")
            
            # Se almacenan las posiciones en variables auxiliares
            auxPosY = sy[i]
            auxPosX = sx[i]  
            
            # Se almacenan las velocidades en variables auxiliares
            auxX = vx[i]
            auxY = vy[i]
            
            # Se asume un choque perfectamente elástico
            # por tanto se intecambian las velocidades
            vx[i] = vx[j]
            vy[i] = vy[j]
            
            vx[j] = auxX
            vy[j] = auxY
            
            # Se regresa mueven las partículas de forma que no queden superpuestas
            # Se realiza la operación vx[i]/abs(vx[i]) con el fin de calcular la dirección del movimiento
            sx[i] = sx[j]  + 2*r*(vx[i]/abs(vx[i]))
            sx[j] = auxPosX + 2*r*(vx[j]/abs(vx[j]))
            
            sy[i] = sy[j]  + 2*r*(vy[i]/abs(vy[i]))
            sy[j] = auxPosY + 2*r*(vy[j]/abs(vy[j]))
          }
        }
        
        # Se buscan partículas que se encuentren en uno de los bordes
        # En caso de encontrarse en alguno de ellos, se invierte la dirección de su movimiento
        # (cambiando el signo a la velocidad), y se pone en una posición en la que no se
        # superponga
        if(sx[i] >= tam_box-(2*r)){
          sx[i] = tam_box - (2*r) - abs(tam_box - sx[i]) 
          vx[i] = - vx[i]
        }else if(sx[i] <= 0+(2*r)){
          sx[i] = 0 + (2*r) + abs(sx[i]) 
          vx[i] = - vx[i]
        }
        
        if(sy[i] >= tam_box-(2*r)){
          sy[i] = tam_box-(2*r) - abs(tam_box - sy[i])
          vy[i] = - vy[i]
        }else if(sy[i] <= 0+(2*r)){
          sy[i] = 0 + (2*r) + abs(sy[i]) 
          vy[i] = - vy[i]
        }
      }
      
      #Se agregan nuevas filas al dataframe de partículas
      particles <- rbind(particles,init_particles(rep(time, n_particles), sx, sy,vx,vy))
    }
  
  #Se ordena la partícula en base a un identificador
  particles = particles[order(particles$particle), ]
  return (particles)
} 

particles = init_particles(rep(0, n_particles), sx,sy,vx,vy)
speeds = init_speeds(vx,vy)
particles <- forward(dt,time_max,particles)

# Se realiza la gráfica con los datos del dataframe, utilizando cada intervalo de tiempo
anim <- ggplot(particles, aes(x=unlist(particles["Sx"]),y=unlist(particles["Sy"])))+ 
  geom_point(aes(colour = factor(particle),size = r/100),alpha = 0.7, show.legend = F) + theme_bw()+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  labs(title = 'Time: {frame_time}', x = 'X', y = 'Y')+
  transition_time(unlist(particles["time"])) 
# Se generan n imágenes correspondientes a los movimientos

# A partir de estas imágenes se genera la animación
animate(anim, height = 500, width = 600, fps = 30, duration = 10,
        end_pause = 60, res = 100 ,renderer = gifski_renderer())

# Se exporta la animación cómo un gif
anim_save("graph.gif")
