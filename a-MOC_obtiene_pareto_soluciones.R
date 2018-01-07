
genera.paretos.finales=function(nameBD, nameAlgorithm, objDim, alfa, num_k, 
                                nombre_dist_exp, nombre_dist_bio, generation, 
                                popSize, ratCruz, ratMuta, populationR, varNo, 
                                iteration, distancia_exp, distancia_bio ){
  
  source("a-MOC_funciones.R")
  
  
  ############ obtener Pareto de la ultima generacion
  d<-as.data.frame(populationR[1:popSize,])
  D = d[order(d$rnkIndex, -d$densidad), ]
  #las que son Pareto
  paretos<-subset(D, rnkIndex=='1') 
  
  nondom<-paretos[!duplicated(paretos[,(1):(varNo+objDim)]),]
  
  
  
  write.table(nondom, 
              file = paste("../Resultados/",nameBD,"/",iteration, "_",nameBD, "_",nameAlgorithm,"_k_",num_k,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_1_poblacion_pareto_clustering_objetivos.csv",sep = ""),sep = " ", row.names = FALSE)
 
  #Forma grupos de cada solucion Pareto
  populationPARETO<-as.matrix(nondom[,1:num_k])
  tablaGruposPARETO<-generaGrupos(nrow(populationPARETO), populationPARETO, distancia_exp, distancia_bio, alfa)

  #Guardar soluciones pareto
  pareto_solutions<-as.data.frame(tablaGruposPARETO)
  colnames(pareto_solutions)<-1:length(tablaGruposPARETO)

 
  write.table(pareto_solutions, 
              file = paste("../Resultados/",nameBD,"/",iteration, "_",nameBD, "_",nameAlgorithm,"_k_",num_k,"_",nombre_dist_exp,"_",nombre_dist_bio,"_P_",popSize,"_TG_",generation,"_C_",ratCruz,"_M_",ratMuta,"_alfa_",alfa,"_1_poblacion_pareto_clustering_solutions.csv",sep = ""),sep = ",", col.names = TRUE)




}