
generaMatrizDistanciaPonderada<-function(par_alfa, par_matriz_exp,par_matriz_bio){
    
    matriz_distancia_ponderada <- ( par_alfa * par_matriz_exp ) + ( (1 - par_alfa) * par_matriz_bio )
    return(matriz_distancia_ponderada)
  
}

generaGrupos<-function(par_popsize, par_medoides, par_matriz_exp, par_matriz_bio, alfa){

  
  grupos<-list(1:par_popsize)
  formacion<-rep(NA,nrow(par_matriz_exp))

  
  for (p in 1:par_popsize) {
 
          for (i in 1:nrow(par_matriz_exp)) {
            grupo<-0
            gen<-i
         
            valor<-(alfa * par_matriz_exp [gen, par_medoides[p,]]) + ( (1-alfa) * par_matriz_bio[gen, par_medoides[p,]])
            grupo<-unname(which.min(valor))
            #grupo<-unname(which.min(par_matrizdistancia[gen,par_medoides[p,]]))
            formacion[i]<-grupo
          }
    
  grupos[[p]] <- formacion
  names(grupos[[p]])<-rownames(par_matriz_exp)

  }
  
  return(grupos)
  
}
 



calculaIndiceSilhouette<-function(par_grupos, par_matriz_distancia){
	#silueta, cercano a 1 es mejor
  fobj<-rep(0,length(par_grupos))
  
      for (x in 1:length(par_grupos)) {
        fobj[[x]]<-(summary(cluster::silhouette(x=unname(unlist(par_grupos[[x]])), dmatrix=par_matriz_distancia))$avg.width)
      } 
  
  return(unlist(fobj))
             
}



calculaIndiceBHI<-function(par_popsize,  par_grupos, par_anotacion, par_categoria){

 
   #se multiplica -1 porque las librerias de NSGA-II estan programadas para trabajar con minizacion  
   fobj<-rep(0,par_popsize)
   
		   for (x in 1:length(par_grupos)) {
			  fobj[x]<-BHI(unlist(par_grupos[[x]]), annotation=par_anotacion, category=par_categoria)*-1  #Lo convierto a problema de minimizacion
		   }
   
   return(unlist(fobj))
 
 
}

calculaIndiceCI <- function(par_popsize, par_grupos, par_matriz_exp) {

  #se multiplica -1 porque las librerias de NSGA-II estan programadas para trabajar con minizacion  
  fobj<-rep(0,par_popsize)
  
  for(x in 1:length(par_grupos)){
    
            CI<-sapply(unique(unlist(par_grupos[[x]])), CoexpressionIndicator, unlist(par_grupos[[x]]), par_matriz_exp)
            CI<-unlist(Filter(Negate(is.null), CI))
            CIsol<-mean(CI)
            CIsol<-CIsol*-1 #Lo convierto a problema de minimizacion
            fobj[x]=CIsol
  }   
  
  #return(CIsol)
  return(unlist(fobj))
  
}


CoexpressionIndicator<- function(i, clusters, par_matriz_expresion) { 
 
  
  ind <- (clusters == i)
  lista_genes_cluster_actual<-names(clusters[mapply(function(X) {if(X==i) return(TRUE) else return(FALSE) }, clusters)==TRUE])
  
  cardinalidad_cluster_actual<-sum(clusters == i)
  
  if (length(lista_genes_cluster_actual)>1) {
    
    lista_genes_cluster_actual<-which(ind)
    
    exp_cluster_actual<-as.matrix(par_matriz_expresion[lista_genes_cluster_actual,])
    
    
    prom_expr_gen<-unname(rowMeans(exp_cluster_actual, na.rm = FALSE, dims = 1))
    sd_expr_gen<-unname(rowSds(exp_cluster_actual, na.rm = FALSE, dims = 1))
    
    par_genes_clu_act<-t(combn(length(lista_genes_cluster_actual), 2)) #genera combinacion sin repeticion
    
    coexp<- (
      (exp_cluster_actual[par_genes_clu_act[,1],] - prom_expr_gen[par_genes_clu_act[,1]])/(sd_expr_gen[par_genes_clu_act[,1]]) *
        (exp_cluster_actual[par_genes_clu_act[,2],] - prom_expr_gen[par_genes_clu_act[,2]])/(sd_expr_gen[par_genes_clu_act[,2]])
    )
    
   
    #promedio_condiciones<-sum(abs(coexp))/ncol(exp_cluster_actual)
    A=1-(coexp)
    promedio_condiciones<-sum((2-A))/ncol(exp_cluster_actual)
    
    Fcard<-2/(cardinalidad_cluster_actual*(cardinalidad_cluster_actual-1))
    CI_i<-Fcard*promedio_condiciones
  }
  
  
}



########################## Distancias  ##########################
#Todas las distancias usadas para procesar silueta son cuadradas
#para hacerlas comparables con la literatura ya que alli usan
#la distancia euclideana cuadrada.

calculaDistanciasResultados<-function(par_matriz_exp, par_nombre_distancia){
  

  
  switch(par_nombre_distancia, 
         pearson={
           sq_distancia <- amap::Dist(par_matriz_exp, method="correlation") ^ 2 #correlation o pearson, es lo mismo. 1-p
         },
         euclidean={
           sq_distancia <- amap::Dist(par_matriz_exp, method="euclidean") ^ 2
         },
         spearman={
           sq_distancia <- amap::Dist(par_matriz_exp, method="sperman") ^ 2
         },
         kendall={
           sq_distancia <- amap::Dist(par_matriz_exp, method="kendall") ^ 2
         }
         
  )
  
  sq_distancia=as.matrix(sq_distancia)
  
  return(sq_distancia)
  
}

cruzamiento_k_puntos<-function(par_popsize, par_num_k, par_pob_Origen, par_pob_Destino, par_cruzRate){
  
  p=0.50
  
  for (i in seq(1, par_popsize, by = 2)) {
    
    parejas<-sample(1:par_popsize,2, replace=F)
    
    cruce=runif(1,0,1)
    
    
    if (cruce < par_cruzRate) {
      
      for (j in 1:par_num_k) {
        
        aleatorio=runif(1,0,1)
        if (aleatorio<=p) {
          par_pob_Destino[i,j]=as.matrix(par_pob_Origen[parejas[1],j])
          par_pob_Destino[i+1,j]=as.matrix(par_pob_Origen[parejas[2],j])
        }else{
          par_pob_Destino[i,j]=as.matrix(par_pob_Origen[parejas[2],j])
          par_pob_Destino[i+1,j]=as.matrix(par_pob_Origen[parejas[1],j])
        }  
      }
      
      
    }else{
      
      par_pob_Destino[i,1:par_num_k]<-as.matrix(par_pob_Origen[parejas[1],1:par_num_k])
      par_pob_Destino[i+1,1:par_num_k]<-as.matrix(par_pob_Origen[parejas[2],1:par_num_k])
      
    }  
    
    
  }
  
  return(par_pob_Destino)
  
}



mutacion_controller_random<-function(par_popsize, par_num_k, par_poblacion, par_cruzRate){
  
  for (p in 1:par_popsize) {
     
    muta=runif(1,0,1)
    
        if (muta < par_cruzRate) {
            posicion_mutar = sample(1:num_k,1, replace=F)
            medoide=sample(1:numgenes,1, replace=F)
            #length(medoide)
            #length(posicion_mutar)
            par_poblacion[p,posicion_mutar]=medoide
        }
     
  }
  
  return(par_poblacion)
  
}

controlaFactibilidad<-function(par_poblacion, par_destino){
  
	#recorre cada individuo, si encuentra un cromosoma repetido(medoide) lo reemplaza aleatoriamente por otro.
	  for (p in 1:nrow(par_poblacion)) {
		
			while(length(which(duplicated(par_poblacion[p,])))>0){
			  id_rep<-which(duplicated(par_poblacion[p,]))
			  num_rep<-length(id_rep)
			  
				  for (i in 1:num_rep) {
					par_poblacion[p,id_rep[i]]=sample(1:numgenes,1, replace=F) 
				  }
			  
			}

	  }
  
  par_destino<-par_poblacion

  return(par_destino)
  
}

quitarSingletons<-function(grupos_formados) {
  #se usa para el calculo de BHI, si hay genes singletons, entonces se quitan de la lista de analisis
  
  i=1
  cluster_singleton<-0
  
  for (variable in 1:max(grupos_formados)) {
    
    if(length(which(grupos_formados == variable))<3){
      cluster_singleton[i]<-variable 
      i<-i+1
    }
    
  }
  
  grupos_sin_singletons<-subset(grupos_formados, !(grupos_formados %in% cluster_singleton))
  
  
  return(grupos_sin_singletons)
}



repararSingletons<-function(par_poblacion_reparar, par_distancia_exp, par_distancia_bio, num_k, alfa) {
 

  par_poblacion_reparar=par_poblacion_reparar[,1:num_k]
  
  
  tablaGruposReparar<-generaGrupos(nrow(par_poblacion_reparar), par_poblacion_reparar, par_distancia_exp, par_distancia_bio, alfa)
  
  #si encuentra solucion con cluster singletons, quita el medoide que produce eso y lo reemplazo por uno aleatorio.
  #el reemplazo se verifica para evitar que sea un medoide ya incluido y ademas que no produzca otros singletons
 
   genes_total <- c(1:numgenes)
  
  for (individuo in 1:nrow(par_poblacion_reparar)) {
    
      #puede suceder que alterar por ejemplo el ultimo medoide del cromosoma hago que el medoide 1 del mismo
      #cromosoma ahora produzca singletons. Por eso cada vez que hay singletons vuelve a verificar desde el primer medoide del cromosoma
      columna <- 1
      
      while(columna <= ncol(par_poblacion_reparar)){
        
              todoOK=FALSE
                
              while(todoOK==FALSE){
                
                    if(length(which(unlist(tablaGruposReparar[[individuo]]) == columna))<3){
  
                        #cat("individuo ", individuo, " medoide ", columna, " produce singletons", "\n")
                        #obtiene los medoides del individuo donde se encontraron singletons
                        genes_cromosoma <- par_poblacion_reparar[individuo,1:ncol(par_poblacion_reparar)]
                        #identifica cuales de todo el conjunto de genes no estan en el cromosoma,
                        #los cuales serviran para el reemplazo del medoide que es singleton
                        posibles_reemplazos<-setdiff(genes_total,intersect(genes_total, genes_cromosoma))
                        #de los genes reciente identificados escojo uno aleatoriamente
                        reemplazo<-sample(posibles_reemplazos,1, replace=F) 
                        
                        #reemplazo tiene el gene que sera ahora usado como medoide
                        par_poblacion_reparar[individuo,columna]=reemplazo
                        #re forma los grupos con el nuevo medoide.
                        tablaGruposReparar[[individuo]]<-generaGrupos(1, as.matrix(t(par_poblacion_reparar[individuo,])), par_distancia_exp, par_distancia_bio, alfa)
                        todoOK=FALSE
                        columna=1
                        
                    }else{
                      todoOK=TRUE
                      columna=columna+1
                    }
                      
              }
              
              
      }
      
      #el individuo i debe estar mejorado cuando haya salido de todos los while, o sea
      #el par_poblacion_reparar[individuo, ] y tablatablaGruposReparar[[individuo]] debe estar con dos genes en cada grupo
 }
 
   #grupos_sin_singletons<-subset(grupos_formados, !(grupos_formados %in% cluster_singleton))

   

  return(list("integer" = tablaGruposReparar, "integer"=par_poblacion_reparar))
  
}




calculaJerarquiasDensidad<-function(par_popsize, par_population, par_grupos, par_objDim){
 
 
  
  #Calcula Objetivos
  obj1<-calculaIndiceCI(par_popsize, par_grupos, matriz_exp)
  obj2<-calculaIndiceBHI(par_popsize,par_grupos, matriz_clases_func,"all")
  
  par_population <- cbind(par_population,obj1,obj2)
  
  #calcula ordenamiento no-dominado
  ranking <- fastNonDominatedSorting(par_population[,(varNo+1):(varNo+par_objDim)])
  #calcula ranking no-dominado
  rnkIndex <- calculaJerarquiaDeFrentes(par_popsize,ranking)
  par_population <- cbind(par_population,rnkIndex)
  #calcula densidad poblacional. Maximos y minimos
  objRange <- calculaRangoObj(par_population,varNo,par_objDim)
  #calcula crowding distance
  cd <- crowdingDist4frnt(par_population,ranking,objRange)
  
  par_population <- cbind(par_population,densidad=c(apply(cd,1,sum)))
  par_population<-as.data.frame(par_population)
  par_population<-par_population[order(par_population$rnkIndex, -(par_population$densidad)), ]
  par_population<-as.matrix(par_population)
  
  #devuelve la poblacion ordenada segun Jerarquia y Densidad. 
  
  return(par_population)
  
  #par_poblacion es un objeto matrix
  
}

calculaJerarquiaDeFrentes<-function(par_popsize, par_ranking){
  
  rankIndex <-integer(par_popsize)
  i <- 1
  while (i <= length(par_ranking)) {
    rankIndex[par_ranking[[i]]] <- i
    i <- i + 1
  }
  
  return(rankIndex)
}

calculaRangoObj<-function(par_poblacion, par_var, par_obj){

  rangoObj<- apply(par_poblacion[,(par_var+1):(par_var+par_obj)], 2, max) - apply(par_poblacion[,(par_var+1):(par_var+par_obj)], 2, min)
  
  return(rangoObj)
  
}



