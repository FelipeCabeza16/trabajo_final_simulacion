

# Charge libraries:
library(ggplot2)
library(gganimate)
#inicializar los valores 

n_particles = 1000
#Atomo de Hidrógeno 
r = 120e-12
#Temeratura constante
T = 300
FPS = 30
dt = 1/FPS

sx = runif(n_particles)
sy = runif(n_particles)

#Velocidad máxima
v_max = 1.5

thetha = runif(n_particles) * 2 * pi
vx = v_max * cos(thetha)
vy = v_max * sin(thetha)

time_max = 5
dt = 0.1


#retorna el dataset con los datos iniciales 
init_particles = function(){
  particles = data.frame(sx, sy, vx, vy)
  names(particles) <- c('Sx', 'Sy','Vx', 'Vy')
  return(particles)
}

forward = function(dt, time_max, particles){

time = 0

pos = array()
  while(time <= time_max){
    time = time + dt  
    sx = sx + vx*dt
    sy = sy + vy*dt
  }
} 

particles = init_particles()  
#forward(dt, time_max)