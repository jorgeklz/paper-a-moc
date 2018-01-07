
  pathBD<-"../../Datasets/sporulation/"
  nameBD<-"sporulation" #nombre de la BD, debe ser el mismo nombre de la carpeta
  organismoNAME<-"yeast"  #organismo a usar por la distancia semantica para obtener relaciones entre genes
  organismoBD<-"org.Sc.sgd.db" #

  #arabidopsis      arabidopsis   org.At.tair.db
  #cell_cycle       yeast         org.Sc.sgd.db
  #sporulation      yeast         org.Sc.sgd.db
  #serum            human         org.Hs.eg.db
  #drosophila       fly     org.Dm.eg.db
  #alzheimer        human   org.Hs.eg.db
  #seminoma         human   org.Hs.eg.db

  nameAlgorithm<-"a-MOC"
 
  nucleos<-3
  #Parametros distancias y k
  objDim <-2
  min_iteration<-3
  max_iteration<-10
  serie_alfa<-seq(1,0,-0.05)
  serie_k<-c(4:6)
  ###############################
  nombre_dist_exp<-"pearson"
  nombre_dist_bio<-"resnik"
  nombre_objetivos<-c("CI","BHI")
  #Parametros geneticos
  poblacion_inicial_aleatorio<-TRUE #SI ES FALSO, ENTONCES USO VERSIONES KMEANS, PAM, UPGMA COMO POBLACION INICIAL
  generation<-100
  popSize <- 50
  ratCruz<-0.80
  ratMuta<-0.01
  tourSize <-2 #ceiling(popSize*0.20)
  #umbral <-10 # criterio terminio, num. veces repeticiones permitidas de promedios de CI o BHI




