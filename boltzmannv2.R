#Par�metros:
#           time -> Vector de intervalos de tiempo
#           sx -> Vector de posiciones en x
#           sy -> Vector de posiciones en y
#           vx -> Vector de velocidades en x
#           vy -> Vector de velocidades en y
#Return:
#           particles -> dataframe de part�culas con los datos ingresados como par�metro,
#                        con un identificador �nico para cada part�cula
init_particles = function(time,sx,sy,vx,vy){
  particles = data.frame(time, c(seq(1,n_particles,by=1)), sx, sy,vx,vy)
  names(particles) <- c('time', 'particle','Sx', 'Sy','Vx','Vy')
  return(particles)
}

#Par�metros: 
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

#VARIOS ESCENARIOS

#inicializar los valores 

n_particles <- 1000

#Atomo de Hidrógeno
escalaR <- 1e7
escalaT <- 1e9

#Se multiplica por un factor 1x10^7 por la escala de la simulaci�n

r <- 120e-12*escalaR

mu <- 1.660e-24
atomic_mass <- 1.00784
m <- mu*atomic_mass
constR <- 8.314472
k_botlz <- 1.380649e-23

#Tama�o de la caja n*n
tam_box <- 1

#Temperatura constante
Temp <- 273.15
FPS <- 30
dt <- 1/FPS


#Generar de manera aleator�a las posiciones iniciales
sx <- runif(n_particles,0,tam_box)
sy <- runif(n_particles,0,tam_box)


#Velocidad máxima
v_max <-  sqrt(  2*k_botlz*Temp / m  )*(escalaR/escalaT)

#Generar de manera aleator�a las velocidades iniciales
thetha <- runif(n_particles) * 2  * pi
v_real <- abs( rnorm(n_particles, mean=v_max) )

vx <- cos(thetha) * v_real 
vy <- sin(thetha) * v_real


#Tiempo m�ximo en seg
time_max = 10
#Diferencial de tiempo 
dt = 0.01


#Par�metros:
#           dt -> Double: Diferencial de tiempo
#           time_max -> Tiempo m�ximo
#           particles -> dataframe con informaci�n de part�culas
#Return:
#           particles -> dataframe de part�culas con los datos despu�s de realizar
#                        interacciones, ordenado por part�cula
box_n_n = function(dt, time_max, particles){
  time = 0
    while(time <= time_max){
      time = time + dt  
      sx <- sx + vx*dt
      sy <- sy + vy*dt
      for (i in 1:length(sx)){
        for(j in 1:length(sx)){
          # Si dos part�culas diferentes se encuentran en el radio una de otra
          if ((abs(sx[j]-sx[i]) <= 2*r) && (abs(sy[j]-sy[i]) <= 2*r) && j!=i){
            print("Choque")
            
            rel_posx = sx[i] - sx[j]
            rel_posy = sy[i] - sy[j]
            
            rel_velx = vx[i] - vx[j]
            rel_vely = vy[i] - vy[j]
            
            v_rel = rel_posx*rel_velx + rel_posy*rel_vely
            r_rel = (rel_posx^2)+(rel_posy^2)
            # Se almacenan las posiciones en variables auxiliares
            auxPosY = sy[i]
            auxPosX = sx[i]  
            
            # Se almacenan las velocidades en variables auxiliares
            auxX = vx[i]
            auxY = vy[i]
            
            # Se asume un choque perfectamente el�stico
            # por tanto se intecambian las velocidades
            vcmx = (vx[i] + vx[j])/2
            vcmy = (vy[i] + vy[j])/2
            
            vx[i] = vx[i] - (v_rel/r_rel)*(sx[i]-sx[j])
            vy[i] = vy[i] - (v_rel/r_rel)*(sy[i]-sy[j])
            #vx[i] = vx[j]
            #vy[i] = vy[j]
            
            vx[j] = vx[j] - (-v_rel/-r_rel)*(sx[j]-sx[i])
            vy[j] = vy[j] - (-v_rel/-r_rel)*(sy[j]-sy[i])
            
            #vx[j] = auxX
            #vy[j] = auxY
            
            # Se regresa mueven las part�culas de forma que no queden superpuestas
            # Se realiza la operaci�n vx[i]/abs(vx[i]) con el fin de calcular la direcci�n del movimiento
            sx[i] = sx[j]  + 2*r*(vx[i]/abs(vx[i]))
            sx[j] = auxPosX + 2*r*(vx[j]/abs(vx[j]))
            
            sy[i] = sy[j]  + 2*r*(vy[i]/abs(vy[i]))
            sy[j] = auxPosY + 2*r*(vy[j]/abs(vy[j]))
          }
        }
        
        # Se buscan part�culas que se encuentren en uno de los bordes
        # En caso de encontrarse en alguno de ellos, se invierte la direcci�n de su movimiento
        # (cambiando el signo a la velocidad), y se pone en una posici�n en la que no se
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
      
      #Se agregan nuevas filas al dataframe de part�culas
      particles <- rbind(particles,init_particles(rep(time, n_particles), sx, sy,vx,vy))
    }
  
  #Se ordena la part�cula en base a un identificador
  #particles = particles[order(particles$particle), ]
  return (particles)
} 


twobox = function(dt, time_max, particles){
  time = 0
  dimpasarelay <- c(1/3,2/3)
  dimpasarelax <- c(0,1/5)
  dim2dacajay <- c(1/6,5/6)
  dim2dacajax <- c(dimpasarelax[2],4/6)
  pos = c(rep(0, n_particles))
  while(time <= time_max){
    time = time + dt  
    sx <- sx + vx*dt
    sy <- sy + vy*dt
    for (i in 1:length(sx)){
      for(j in 1:length(sx)){
        # Si dos part�culas diferentes se encuentran en el radio una de otra
        if ((abs(sx[j]-sx[i]) <= 2*r) && (abs(sy[j]-sy[i]) <= 2*r) && j!=i){
          print("Choque")
          
          rel_posx = sx[i] - sx[j]
          rel_posy = sy[i] - sy[j]
          
          rel_velx = vx[i] - vx[j]
          rel_vely = vy[i] - vy[j]
          
          v_rel = rel_posx*rel_velx + rel_posy*rel_vely
          r_rel = (rel_posx^2)+(rel_posy^2)
          # Se almacenan las posiciones en variables auxiliares
          auxPosY = sy[i]
          auxPosX = sx[i]  
          
          # Se almacenan las velocidades en variables auxiliares
          auxX = vx[i]
          auxY = vy[i]
          
          # Se asume un choque perfectamente el�stico
          # por tanto se intecambian las velocidades
          vcmx = (vx[i] + vx[j])/2
          vcmy = (vy[i] + vy[j])/2
          
          vx[i] = vx[i] - (v_rel/r_rel)*(sx[i]-sx[j])
          vy[i] = vy[i] - (v_rel/r_rel)*(sy[i]-sy[j])
          #vx[i] = vx[j]
          #vy[i] = vy[j]
          
          vx[j] = vx[j] - (-v_rel/-r_rel)*(sx[j]-sx[i])
          vy[j] = vy[j] - (-v_rel/-r_rel)*(sy[j]-sy[i])
          
          #vx[j] = auxX
          #vy[j] = auxY
          
          # Se regresa mueven las part�culas de forma que no queden superpuestas
          # Se realiza la operaci�n vx[i]/abs(vx[i]) con el fin de calcular la direcci�n del movimiento
          sx[i] = sx[j]  + 2*r*(vx[i]/abs(vx[i]))
          sx[j] = auxPosX + 2*r*(vx[j]/abs(vx[j]))
          
          sy[i] = sy[j]  + 2*r*(vy[i]/abs(vy[i]))
          sy[j] = auxPosY + 2*r*(vy[j]/abs(vy[j]))
        }
      }
      
      # Se buscan part�culas que se encuentren en uno de los bordes
      # En caso de encontrarse en alguno de ellos, se invierte la direcci�n de su movimiento
      # (cambiando el signo a la velocidad), y se pone en una posici�n en la que no se
      # superponga
      if(pos[i]==0){
        if(sx[i] >= tam_box-(2*r)){
          if(sy[i]>= dimpasarelay[2]*tam_box || sy[i]<= dimpasarelay[1]*tam_box){
            sx[i] = tam_box - (2*r) - abs(tam_box - sx[i]) 
            vx[i] = - vx[i]
          }else{
            pos[i]=1
          }
          
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
      if (pos[i]==1){
        if(sx[i] < tam_box && vx[i]<0){
            pos[i]=0
        }else if(sx[i] >= (1+dimpasarelax[2])*tam_box && vx[i]>0){
          pos[i]=2
        }
        
        if(sy[i] >= dimpasarelay[2]*tam_box-(2*r)){
          sy[i] = dimpasarelay[2]*tam_box-(2*r) - abs(dimpasarelay[2]*tam_box - sy[i])
          vy[i] = - vy[i]
        }else if(sy[i] <= dimpasarelay[1]*tam_box+(2*r)){
          sy[i] = dimpasarelay[1]*tam_box + (2*r) + abs(sy[i]) 
          vy[i] = - vy[i]
        }
      }
      
      if (pos[i]==2){
        if(sx[i] <= (1+dimpasarelax[2])*tam_box+(2*r)){
          if(sy[i]>= dimpasarelay[2]*tam_box || sy[i]<= dimpasarelay[1]*tam_box){
            sx[i] = (1+dimpasarelax[2])*tam_box + (2*r) + abs(sx[i]-(1+dimpasarelax[2])*tam_box) 
            vx[i] = - vx[i]
          }else{
            pos[i]=1
          }
          
        }else if(sx[i] >= ((1+dimpasarelax[2])+(dim2dacajax[2]))*tam_box-(2*r)){
          sx[i] = (((1+dimpasarelax[2])+(dim2dacajax[2]))*tam_box) - (2*r) - abs(((1+dimpasarelax[2])+(dim2dacajax[2]))*tam_box-sx[i]) 
          vx[i] = - vx[i]
        }
        
        if(sy[i] >= (dim2dacajay[2])*tam_box-(2*r)){
          sy[i] = (dim2dacajay[2])*tam_box-(2*r) - abs(tam_box - sy[i])
          vy[i] = - vy[i]
        }else if(sy[i] <= (dim2dacajay[1])*tam_box+(2*r)){
          sy[i] = (dim2dacajay[1])*tam_box + (2*r) + abs((dim2dacajay[1])*tam_box-sy[i]) 
          vy[i] = - vy[i]
        }
      }
    }
    
    #Se agregan nuevas filas al dataframe de part�culas
    particles <- rbind(particles,init_particles(rep(time, n_particles), sx, sy,vx,vy))
  }
  
  #Se ordena la part�cula en base a un identificador
  #particles = particles[order(particles$particle), ]
  return (particles)
}


particles = init_particles(rep(0, n_particles), sx,sy,vx,vy)
speeds = init_speeds(vx,vy)
particles <- box_n_n(dt,time_max,particles)
#particles <- twobox(dt,time_max,particles)

# Se realiza la gr�fica con los datos del dataframe, utilizando cada intervalo de tiempo
anim <- ggplot(particles, aes(x=unlist(particles["Sx"]),y=unlist(particles["Sy"])))+ 
  geom_point(aes(colour = factor(particle)),size=r,alpha = 0.7, show.legend = F) + theme_bw()+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  labs(title = 'Time: {frame_time}', x = 'X', y = 'Y') +
  transition_time(unlist(particles["time"])) + 
  geom_line(data=data.frame("x"=c(0,0),"y"=c(0,tam_box)),aes_string(x="x",y="y"),linetype="dashed")+
  geom_line(data=data.frame("x"=c(0,tam_box),"y"=c(tam_box,tam_box)),aes_string(x="x",y="y"),linetype="dashed")+
  geom_line(data=data.frame("x"=c(0,tam_box),"y"=c(0,0)),aes_string(x="x",y="y"),linetype="dashed")+
  geom_line(data=data.frame("x"=c(tam_box,tam_box),"y"=c(0,tam_box)),aes_string(x="x",y="y"),linetype="dashed")
  #geom_line(data=data.frame("x"=c(tam_box,tam_box),"y"=c((2/3)*tam_box,tam_box)),aes_string(x="x",y="y"),linetype="dashed")+
  #geom_line(data=data.frame("x"=c(tam_box,tam_box),"y"=c(0,(1/3)*tam_box)),aes_string(x="x",y="y"),linetype="dashed")+
  #geom_line(data=data.frame("x"=c(tam_box,tam_box*1.2),"y"=c((1/3)*tam_box),(1/3)*tam_box),aes_string(x="x",y="y"),linetype="dashed")+
  #geom_line(data=data.frame("x"=c(tam_box,tam_box*1.2),"y"=c((2/3)*tam_box),(2/3)*tam_box),aes_string(x="x",y="y"),linetype="dashed")+
  #geom_line(data=data.frame("x"=c(tam_box*1.2,tam_box*1.2),"y"=c(((1/6)*tam_box),(1/3)*tam_box)),aes_string(x="x",y="y"),linetype="dashed")+
  #geom_line(data=data.frame("x"=c(tam_box*1.2,tam_box*1.2),"y"=c(((2/3)*tam_box),(5/6)*tam_box)),aes_string(x="x",y="y"),linetype="dashed")+
  #geom_line(data=data.frame("x"=c(tam_box*1.2,tam_box*(1.2+(4/6))),"y"=c((1/6)*tam_box),(1/6)*tam_box),aes_string(x="x",y="y"),linetype="dashed")+
  #geom_line(data=data.frame("x"=c(tam_box*1.2,tam_box*(1.2+(4/6))),"y"=c((5/6)*tam_box),(5/6)*tam_box),aes_string(x="x",y="y"),linetype="dashed")
# Se generan n im�genes correspondientes a los movimientos

# A partir de estas im�genes se genera la animaci�n
animate(anim, height = 500, width = 600, fps = 30, duration = 20,
        end_pause = 60, res = 100 ,renderer = gifski_renderer())

# Se exporta la animaci�n c�mo un gif
anim_save("graph.gif")

totalvel <- sqrt((particles$Vx^2)+(particles$Vy^2))
totalvel <- tail( totalvel, n_particles)


library("MASS")

#k_botlz <- 8.617333262e-5

truehist(totalvel*(escalaT/escalaR),nbins = 100)
x_grafica <- seq(min(totalvel),max(totalvel),by=0.1)

sigma <- sqrt( k_botlz*T / m  )

y_grafica <-  ( x_grafica / (sigma^2) )* exp( - (  x_grafica^2 / (2*sigma^2) ) ) 

#y_grafica <- sqrt((M/(2*pi*k_botlz)))*(exp((-M*(x_grafica^2))/(2*k_botlz*T)))

lines(x_grafica,y_grafica,type="l",col="red")
print(totalvel)
