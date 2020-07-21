
  pathBD<-"/Datasets/arabidopsis/"
  nameBD<-"arabidopsis" #nombre de la BD, debe ser el mismo nombre de la carpeta
  organismoNAME<-"arabidopsis"  #organismo a usar por la distancia semantica para obtener relaciones entre genes
  organismoBD<-"org.At.tair.db" #

  #arabidopsis      arabidopsis   org.At.tair.db
  #cell_cycle       yeast         org.Sc.sgd.db
  #sporulation      yeast         org.Sc.sgd.db
  #serum            human         org.Hs.eg.db

  nameAlgorithm<-"a-MOC"
 
  nucleos<-3
  #Parametros distancias y k
  objDim <-2
  min_iteration<-1
  max_iteration<-1
  serie_alfa<-seq(0.10,0.0,-0.05)
  serie_k<-c(5:6)
  ###############################
  nombre_dist_exp<-"pearson"
  nombre_dist_bio<-"wang"
  nombre_objetivos<-c("CI","BHI")
  #Parametros geneticos
  poblacion_inicial_aleatorio<-TRUE #SI ES FALSO, ENTONCES USO VERSIONES KMEANS, PAM, UPGMA COMO POBLACION INICIAL QUE DEBEN SER EJECUTADAS PREVIAMENTE
  generation<-10
  popSize <- 6 #Debe ser numero par
  ratCruz<-0.80
  ratMuta<-0.01
  tourSize <-2 #ceiling(popSize*0.20)
  #umbral <-10 # criterio terminio, num. veces repeticiones permitidas de promedios de CI o BHI




