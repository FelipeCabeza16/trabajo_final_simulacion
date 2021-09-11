

# Charge libraries:
library(ggplot2)
library(gganimate)
#inicializar los valores 

particles = 1000
#Atomo de Hidrógeno 
r = 120e-12
#Temeratura constante
T = 300
FPS = 30
dt = 1/FPS

sx = runif(particles)
sy = runif(particles)

s = array(c(sx, sy), dim = c(particles, 2))

#Velocidad máxima
v_max = 1.5

thetha = runif(particles) * 2 * pi
vx = v_max * cos(thetha)
vy = v_max * sin(thetha)

v = array(c(vx, vy), dim = c(particles, 2))

time_max = 5
dt = 0.1

forward = function(dt, time_max){

time = 0

pos = array()
  while(time <= time_max){
    time = time + dt  
    s = s + v*dt
  }
} 
  
forward(dt, time_max)