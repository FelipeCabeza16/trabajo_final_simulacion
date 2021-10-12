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
  
  
  # Charge libraries:
  library(ggplot2)
  library(gganimate)
  library(gifski)
  library(png)
  
  #VARIOS ESCENARIOS
  
  #inicializar los valores 
  
  n_particles <- 1000
  
  
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
  
  #Tama�o de la caja n*n
  tam_box <- 1
  
  
  FPS <- 30
  dt <- 1/FPS
  # dt = 0.01
  
  #Generar de manera aleator�a las posiciones iniciales
  sx <- runif(n_particles,0,tam_box)
  sy <- runif(n_particles,0,tam_box)
  
  
  #Velocidad máxima
  v_max <-  sqrt( 2*k_botlz*Temp / m  )*(escalaR/escalaT)*1.7
  
  #Generar de manera aleator�a las velocidades iniciales
  thetha <- runif(n_particles) * 2  * pi
  
  # v_real <- abs( rnorm(n_particles, mean=v_max) )
  
  v_real <- runif(n_particles )*v_max
  
  vx <- cos(thetha) * v_real 
  vy <- sin(thetha) * v_real
  
  
  #Tiempo m�ximo en seg
  time_max = 20
  #Diferencial de tiempo 

  
  
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
            
            dist = sqrt( (sx[j] - sx[i])**2 + (sy[j] - sy[i])**2 )
            
            if ( dist <= 2*r && i < j){
              print("Choque")
              rel_posx = sx[i] - sx[j]
              rel_posy = sy[i] - sy[j]
              
              rel_velx = vx[i] - vx[j] 
              rel_vely = vy[i] - vy[j]
              
              
              r_rel = (rel_posx**2)+(rel_posy**2)
              v_rel = rel_posx*rel_velx + rel_posy*rel_vely
              
              # Se asume un choque perfectamente elastico
              # por tanto se intercambian las velocidades
              vcmx = (vx[i] + vx[j])/2
              vcmy = (vy[i] + vy[j])/2
              
              change_x = 2*rel_posx*v_rel / r_rel - rel_velx
              change_y = 2*rel_posy*v_rel / r_rel - rel_vely
              
              
              vx[i] = vcmx - change_x/2
              vy[i] = vcmy - change_y/2
              
              vx[j] = vcmx + change_x/2
              vy[j] = vcmy + change_y/2
  
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
          
          dist = sqrt( (sx[j] - sx[i])**2 + (sy[j] - sy[i])**2 )
          
          if ( dist <= 2*r && i < j){
            print("Choque")
            rel_posx = sx[i] - sx[j]
            rel_posy = sy[i] - sy[j]
            
            rel_velx = vx[i] - vx[j] 
            rel_vely = vy[i] - vy[j]
            
            
            r_rel = (rel_posx**2)+(rel_posy**2)
            v_rel = rel_posx*rel_velx + rel_posy*rel_vely
            
            # Se asume un choque perfectamente elastico
            # por tanto se intercambian las velocidades
            vcmx = (vx[i] + vx[j])/2
            vcmy = (vy[i] + vy[j])/2
            
            change_x = 2*rel_posx*v_rel / r_rel - rel_velx
            change_y = 2*rel_posy*v_rel / r_rel - rel_vely
            
            
            vx[i] = vcmx - change_x/2
            vy[i] = vcmy - change_y/2
            
            vx[j] = vcmx + change_x/2
            vy[j] = vcmy + change_y/2
            
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
  #crear el dataset
  
  speeds = init_speeds(vx,vy)
  # particles <- box_n_n(dt,time_max,particles)
  particles <- twobox(dt,time_max,particles)
  write.csv(particles, "C:\\Users\\BRAYANMONROY\\Documents\\2021-I\\simulacion\\trabajo_final_simulacion\\results\\twobox.csv")
  
  
  
  
  
  
  
  
