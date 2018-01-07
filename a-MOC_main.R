##############################################
#           Jorge A. Parraga-Alava
#           jorge.parraga@usach.cl
############################################## 
setwd("C:/Users/jbele/Dropbox/Papers/2018/comparacion_biologia/a-MOC/Algoritmo")

source("a-MOC_librerias.R")
source("a-MOC_parametros.R")
source("a-MOC_distancias.R")
source("a-MOC_funciones.R")


funcion_poblacion_inicial=dget("a-MOC_poblacion_inicial.R") 
funcion_poblacion_final=dget("a-MOC_obtiene_pareto_soluciones.R") 

options(warn=-1) #poner 0 para volver activar advertencias 


cat("----------------------------------------","\n", 
    "Dataset:",nameBD, "with",numgenes, "x", ncol(matriz_exp), "\n",
    "Distance Exp:",nombre_dist_exp, "\n",
    "Distance Bio:",nombre_dist_bio, "\n",
    "Objectives:",nombre_objetivos, "\n",
    "Num. Generations:",generation, "\n",
    "Population Size:",popSize, "\n",
    "Random Population:", poblacion_inicial_aleatorio, "\n",
    "Crossover Rate:",ratCruz, "\n",
    "Mutation Rate:",ratMuta, "\n",
    "----------------------------------------","\n","\n",
    "---------------------------------------- ", "\n",
    "--------------- * Start * -------------- ","\n",
    "---------------------------------------- ","\n","\n")


    lista_funciones=c(lsf.str())
    lista_paquetes=c(.packages())  


    
  
    
    #archivo_log = paste("log_", max_iteration, "_alfa_", serie_alfa,"_", nameBD, ".txt",sep = "")
    #cat(c(""), file =archivo_log, append=FALSE) 
    
    #writeLines(c(""), paste("log_", max_iteration, "_alfa_", serie_alfa,"_", nameBD, ".txt",sep = ""))
    
    
    for  (iteration in min_iteration:max_iteration) { 
      
          cat("\n"," Iteration ",iteration,"of",max_iteration,"\n")
       
          #------------------------------------#
          #------------------------------------#
                  cl <- makeCluster(nucleos)
                  registerDoParallel(cl)
          #------------------------------------#
          #------------------------------------#
      
           #for  (alfa in serie_alfa) { 
           foreach (alfa = serie_alfa, .export=lista_funciones, .packages=lista_paquetes, .errorhandling="stop") %dopar% {
            
           #     
              
                #sink(paste("log_", max_iteration, "_alfa_", serie_alfa,"_", nameBD, ".txt",sep = ""), append=TRUE)
             
                cat("\n","  alfa=",alfa)#, file =archivo_log, append=TRUE)
                
                ##########################################################################################
                ############################## MATRIZ DISTANCIA PONDERADA ################################
                
                matriz_dist_ponderada<-generaMatrizDistanciaPonderada(alfa, distancia_exp, distancia_bio)
                
                ##########################################################################################            
                ##########################################################################################
                
                
                      #------------------------------------#
                      #------------------------------------#
                      cl_k <- makeCluster(nucleos)
                      registerDoParallel(cl_k)
                      #------------------------------------#
                      #------------------------------------#
                      
                  foreach (ik = serie_k, .export=lista_funciones, .packages=lista_paquetes, .errorhandling="stop") %dopar% {
                      #for  (ik in serie_k){
                      
                        
                                TiempoInicio <- proc.time()
                        
                        
                                num_k=varNo=ik
                                
                                #cat("     k=",num_k,"\n")#, file =archivo_log, append=TRUE)
                                
                                ##########################
                                #Inicializa poblacion P
                                ##########################
                                populationP = funcion_poblacion_inicial(nameBD, nameAlgorithm, objDim, alfa, ik, numgenes, matriz_dist_ponderada, 
                                                                      nombre_dist_exp, nombre_dist_bio, poblacion_inicial_aleatorio, 
                                                                      generation, popSize, ratCruz, ratMuta, tourSize)
                            
                                
                                promedio_poblacion_historia=vector()
                                respuesta=FALSE
                                g=1
                    
                                while (g<=generation) {
                                  
             
                                          
                                      ######################## Poblacion P ########################
                                      verficaSingletonsP<-repararSingletons(populationP, distancia_exp, distancia_bio, num_k, alfa) 
                                      tablaGruposP<-verficaSingletonsP[[1]] #1 corresponde al primer argumento devuelvo por repararSingletons y corresponde a los grupos formados
                                      populationP<- as.matrix(verficaSingletonsP[[2]]) #1 corresponde al segundo argumento ... es la poblacion sin singletons
                                      
                                      populationP<-calculaJerarquiasDensidad(popSize, populationP, tablaGruposP, objDim)
                                      
                                      #Imprime objetivos y soluciones de clustering de la poblacion inicial (semilla)
                                      if (iteration==1 && g==1){
                                          write.table(populationP[, -((varNo+objDim+1):(varNo+objDim+2))], file = paste("../Resultados/",nameBD,"/1_",nameBD, "_",nameAlgorithm,"_k_",ik,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_0_poblacion_semilla_clustering_objetivos.csv",sep = ""),sep = " ",col.names = TRUE,row.names = FALSE)
                                          #Guardar soluciones de poblacion P
                                          clustering_solution<-as.data.frame(tablaGruposP); colnames(clustering_solution)<-1:length(tablaGruposP)
                                          write.table(clustering_solution, file = paste("../Resultados/",nameBD,"/1_",nameBD, "_",nameAlgorithm,"_k_",ik,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_0_poblacion_semilla_clustering_solutions.csv",sep = ""),sep = " ",col.names = TRUE)
                                      }
                                      
                                      ######################## Seleccion, Cruzamiento, Mutacion  ########################                  
                                      #Seleccion
                                      matingPool <- tournamentSelection(populationP,popSize,tourSize)
                                      populationQ <- t(sapply(1:popSize, function(u) array(rep(0,num_k))))
                                      #Cruzamiento y Mutacion
                                      cruza <-cruzamiento_k_puntos(popSize,num_k,matingPool, populationQ, ratCruz)
                                      populationQ<-controlaFactibilidad(cruza,populationQ)
                                      muta <-mutacion_controller_random(popSize,num_k,populationQ, ratMuta)
                                      populationQ<-controlaFactibilidad(muta,populationQ)
                                      
                                      ######################## Poblacion Q  ########################
                                      verficaSingletonsQ<-repararSingletons(populationQ, distancia_exp, distancia_bio, num_k, alfa) 
                                      tablaGruposQ<-verficaSingletonsQ[[1]]
                                      populationQ<- as.matrix(verficaSingletonsQ[[2]])
                                      populationQ<-calculaJerarquiasDensidad(popSize, populationQ, tablaGruposQ, objDim)
                                      
                                      ######################## Poblacion R   ########################
                                      #como P y Q han sido corregidas para no tener singletons, ya no debo verificar eso en poblacion R, ya que es la union de ambas
                                      populationR <- rbind(populationP,populationQ)
                                      rownames(populationR)<-1:nrow(populationR) 
                                      populationR<-populationR[, -((ik+1):(varNo+objDim+2))] #deja solo cromosomas, quita objetivos y las dos ultimas columnas o sea Jerarq y Densidad
                                      #recalcula Jerarquia Pareto y densidad en poblacion R
                                      tablaGruposR<-generaGrupos(nrow(populationR), populationR, distancia_exp, distancia_bio, alfa)
                                      populationR<-calculaJerarquiasDensidad(popSize*2, populationR, tablaGruposR, objDim)
                            
                                       
                                      verficaSingletonsR<-repararSingletons(populationR, distancia_exp, distancia_bio, num_k, alfa) 
                                      tablaGruposR<-verficaSingletonsR[[1]]
                                      populationR<- as.matrix(verficaSingletonsR[[2]])
                                      populationR<-populationR[,1:num_k]
                                      rownames(populationR)<-1:nrow(populationR)
                                      
                                      populationR<-calculaJerarquiasDensidad(popSize, populationR, tablaGruposR, objDim)
                                      
                                 
                                      populationR<-as.data.frame(populationR)
                                  
 
                                      ####################### CRITERIO DE TÉRMINO #############################
                                      #populationPareto= subset(populationR, rnkIndex=='1') 
                                      #cat("k:", num_k," g: ", g, " obj1: min. ", min(populationPareto$obj1), "avg. ", mean(populationPareto$obj1), "\n")
                                      #cat("k:", num_k," g: ", g, " obj2: min. ", min(populationPareto$obj2), "avg. ", mean(populationPareto$obj2), "\n \n")
                                    
                                      #promedio_poblacion_actual=apply(populationPareto[,(varNo+1):(varNo+objDim)],2, mean)
                                      #promedio_poblacion_historia=rbind(promedio_poblacion_historia, promedio_poblacion_actual)
                                      #rownames(promedio_poblacion_historia)=1:nrow(promedio_poblacion_historia)
                                      #promedio_poblacion_repetidos=plyr::count(promedio_poblacion_historia)$freq
                                      #se termina la ejecucion cuando hay umbral veces repetido el promedio
                                      #terminos=sapply(promedio_poblacion_repetidos, function (x) if (x>=umbral) respuesta=TRUE else respuesta=FALSE )
                                      #if(sum(terminos)>0) respuesta=TRUE else respuesta=FALSE
                                      
                                      #if(respuesta==TRUE){
                                       # cat("criterio termino en generacion : ", g, " num. sol.", nrow(populationR), "\n")
                                        #g=generation
                                      #}
                                      #########################################################################
                                  
                                        
                                      ######################################################################################
                                      #imprime poblacion , objetivo y ordenada por frente en la primera y ultima generacion
                                      if (g==1 || g==generation) {
                                        
                                          #Se guardan valores positivos. Ya que trato el problema como minimizacion a pesar de ser de maximizacion
                                          populationR$obj1<-populationR$obj1*-1
                                          populationR$obj2<-populationR$obj2*-1
                                          
                                          
                                          d<-as.data.frame(populationR[,])
                                          D = d[order(d$rnkIndex, -d$densidad), ]
                                          #las que son Pareto
                                          paretos<-subset(D, rnkIndex=='1') 
                                          
                                          nondom<-paretos[!duplicated(paretos[,(1):(varNo+objDim)]),]
                                          
                                          if(g==1) write.table(nondom[1:nrow(nondom),], file = paste("../Resultados/",nameBD,"/",iteration, "_",nameBD, "_",nameAlgorithm,"_k_",num_k,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_0_generation_",g,".csv",sep = ""),sep = " ",col.names = TRUE,row.names = FALSE)
                                          if(g==generation) write.table(nondom[1:nrow(nondom),], file = paste("../Resultados/",nameBD,"/",iteration, "_",nameBD, "_",nameAlgorithm,"_k_",num_k,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_1_generation_",g,".csv",sep = ""),sep = " ",col.names = TRUE,row.names = FALSE)
                                      }
                                      #####################################################################################
                                  
                                      
                                      
                                      ######################## Poblacion P+1  ########################                        
                                      populationP<- as.matrix(populationR[1:popSize,1:num_k])
                                  
                                      
                                      g=g+1

                                      
                                }
                    
                                
                                genera_soluciones_finales=funcion_poblacion_final(nameBD, nameAlgorithm, objDim, alfa, num_k, nombre_dist_exp, 
                                                                                  nombre_dist_bio, generation, popSize, ratCruz, ratMuta, 
                                                                                  populationR, varNo, iteration, distancia_exp, distancia_bio)
                        
                                
                                
                                TiempoFin<- proc.time()- TiempoInicio   # Detiene el cronometro
                                cat("\n k=", num_k, (TiempoFin))#, file =archivo_log, append=TRUE)
                                
                                
                                
                      }
            
                      #------------------------------------#
                      #------------------------------------#                    
                       stopCluster(cl_k) 
                      #------------------------------------#
                      #------------------------------------#  
            
          }
              #Procesa plot e hypervolumen y guarda en carpeta finales
              #source("Dual-MOC_calcula_plot_soluciones_pareto_vs_mono.R")
              #source("Dual-MOC_calcula_hypervolume.R")
           #------------------------------------#
           #------------------------------------#                    
                   stopCluster(cl) 
           #------------------------------------#
           #------------------------------------#  
              
      
    }
    
 
    
    
    
  

    
  cat(" ---------------------------------","\n")
  cat("--------------End----------------","\n")
  cat("---------------------------------","\n")



print(TiempoFin)
 

