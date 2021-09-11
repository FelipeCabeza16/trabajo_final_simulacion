

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
init_particles = function(time,sx,sy){
  particles = data.frame(time,c(seq(1,n_particles,by=1)),sx, sy, vx, vy)
  names(particles) <- c('time','particle','Sx', 'Sy','Vx', 'Vy')
  return(particles)
}

forward = function(dt, time_max, particles){
  time = 0
  n_steps <-0
  pos = array()
    while(time <= time_max){
      time = time + dt  
      sx <- sx + vx*dt
      sy <- sy + vy*dt
      n_steps <- n_steps + 1
      particles <- rbind(particles,init_particles(rep(n_steps,n_particles),sx,sy))
    }
  return (particles)
} 

particles = init_particles(c(integer(n_particles)),sx,sy)
particles <- forward(dt,time_max,particles)
ggplot(particles, aes(x=unlist(particles["Sx"]),y=unlist(particles["Sy"])))+ 
  geom_point() +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  # gganimate specific bits:
  labs(title = 'Year: {frame_time}', x = 'GDP per capita', y = 'life expectancy') +
  transition_time(particles['time']) +
  ease_aes('linear')
#forward(dt, time_max)