
poblacion.inicial=function( nameBD, nameAlgorithm, objDim, alfa, num_k, numgenes, matriz_dist_ponderada,
                            nombre_dist_exp, nombre_dist_bio, poblacion_inicial_aleatorio, generation, popSize, ratCruz, ratMuta, tourSize){

  
  source("a-MOC_funciones.R")
  
  
  ik=num_k
  
  
  if (poblacion_inicial_aleatorio==TRUE) {
    
        #######################################
        # genera poblacion inicial aleatoria
        #######################################
        set.seed(1987)
        populationP = t(sapply(1:popSize, function(u)  sample(1:numgenes, num_k ,replace=F)  ))
        rm(.Random.seed, envir=globalenv()) 
        write.table(populationP, file = paste("Resultados/",nameBD,"/1_",nameBD, "_",nameAlgorithm,"_k_",ik,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_0_poblacion_semilla_cromosomas.csv",sep = ""),sep = " ",col.names = TRUE,row.names = FALSE)
        
        #hacerTodo<-TRUE
    
  }else{
  
        #############################################################
        # genera poblacion inicial usando soluciones mono objetivo
        ############################################################
          #1, para que independiente de la iteracion seimpre se usen la misma pobliacion inicial
          #si no existe el archivo con los cromosomas de la poblacion los creo.
          if(!file.exists(paste("../Resultados/",nameBD,"/1_",nameBD, "_",nameAlgorithm,"_k_",ik,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_0_poblacion_semilla_cromosomas.csv",sep = ""))){
               
     
                        ######################################################
                        # funcion para encontrar elemento medoide en grupos  #
                          clust.medoid = function(i, distmat, clusters) {
                            ind = (clusters == i)
                            names(which.min(rowSums( distmat[ind, ind] )))
                          }
                        ################################################### 
            
           
                          #lee la ruta de los archivos con soluciones (particiones) mono-objetivo
                          files_semillaAlg1<- unname(paste("../../Dual-MOC_mono/Resultados/",nameBD,"/",nameBD, "_k_",ik,"_",nombre_dist_exp,"_",nombre_dist_bio,"_alfa_",alfa,"_clustering_solutions_kmeans.csv",sep = ""))
                          files_semillaAlg2<- unname(paste("../../Dual-MOC_mono/Resultados/",nameBD,"/",nameBD, "_k_",ik,"_",nombre_dist_exp,"_",nombre_dist_bio,"_alfa_",alfa,"_clustering_solutions_PAM.csv",sep = ""))      
                          files_semillaAlg3<- unname(paste("../../Dual-MOC_mono/Resultados/",nameBD,"/",nameBD, "_k_",ik,"_",nombre_dist_exp,"_",nombre_dist_bio,"_alfa_",alfa,"_clustering_solutions_UPGMA.csv",sep = ""))      
                          
                          
                          #inicializo matriz de poblacion inicial con 0
                          populationP_semilla <- t(sapply(1:popSize, function(u)  rep(0, num_k)))
                          
                          #controlo la columna actual de las particiones de k-means
                          col_sol<-1
                          
                          for (j in 1:popSize) {
                            
                                    matriz= matriz_dist_ponderada
                                    
                                    if((j==(popSize-1))){
                                            #TOMO LA UNICA SOLUCION UPGMA
                                            solucion<-read.table(files_semillaAlg3,header = T) 
                                            nombres<-rownames(solucion)
                                            solucion<-as.vector(solucion[,1])
                                            names(solucion)<-nombres
                                    }else if(j==popSize){
                                            #TOMO LA UNICA SOLUCION PAM
                                            solucion<-read.table(files_semillaAlg2,header = T)
                                            nombres<-rownames(solucion)
                                            solucion<-as.vector(solucion[,1])
                                            names(solucion)<-nombres
                                    }else{
                                            #TOMO LA SOLUCIONES K-means
                                            solucion<-read.table(files_semillaAlg1,header = T)
                                            nombres<-rownames(solucion)
                                            solucion<-as.vector(solucion[,col_sol])
                                            names(solucion)<-nombres
                                            col_sol<-col_sol+1 
                                    }
                                    
                                    mydist = as.matrix(dist(matriz))
                                    clusters_sol<-solucion
                                    
                                    nombre_medoide<-sapply(unique(clusters_sol), clust.medoid, mydist, quitarSingletons(clusters_sol))
                                    indice_medoide<-which(rownames(matriz) %in% nombre_medoide)
                                    
                                    #Si hay grupos sin genes, genero un nuevo medoide para ese cluster.
                                    #O sea es como reemplazar el individuo actual por uno parecido donde se mantienen los demas medoides.
                                    
                                    if (length(indice_medoide)<num_k) { 
                                        for(ii in (length(indice_medoide)+1):num_k) {
                                          indice_medoide[ii]<-sample(1:numgenes,1, replace=F)
                                        }
                                    }
                                    
                                    populationP_semilla[j, ]<-indice_medoide
                            
                          }
                          
            
                    populationP <- populationP_semilla
                    write.table(populationP, file = paste("../Resultados/",nameBD,"/1_",nameBD, "_",nameAlgorithm,"_k_",ik,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_0_poblacion_semilla_cromosomas.csv",sep = ""),sep = " ",col.names = TRUE,row.names = FALSE)
      
          }else{
                  
                  #si existe el archivo lo leo   
                  population_semilla_leida<-read.table(paste("../Resultados/",nameBD,"/1_",nameBD, "_",nameAlgorithm,"_k_",ik,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_0_poblacion_semilla_cromosomas.csv",sep = ""),header = T)
                  populationP <- as.matrix(population_semilla_leida)
          }
        
          #hacerTodo<-FALSE
          
  
  }
  
  
  return(populationP)
  
  
  
}
