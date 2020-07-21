    ##########################################
    # leer datos de expresion y clases funcionales
    ##########################################
    
    matriz_exp<-as.matrix(read.table(paste(getwd(),pathBD,"BD.csv",sep=""), header=T, sep=","))
    #matriz_exp<-normalize(matriz_exp, byrow=TRUE)
    matriz_clases_func<-as.matrix(read.table(paste(getwd(),pathBD,"BIO.csv",sep=""), header=T, sep=";"))
    numgenes<-nrow(matriz_exp)
    
    ##########################################
    # calcula matriz distancia por Expresi?n
    ##########################################
    if (!file.exists(paste(getwd(),pathBD,nameBD, "_",nameAlgorithm,"_matriz_distancia_expresion_",nombre_dist_exp,".csv",sep = ""))){
      
      switch(nombre_dist_exp, 
             pearson={
               distancia_exp <- amap::Dist(matriz_exp, method="correlation")
             },
             spearman={
               distancia_exp <- amap::Dist(matriz_exp, method="sperman")
             },
             euclidean={
               distancia_exp <- amap::Dist(matriz_exp, method="euclidean")
             },                
             kendall={
               distancia_exp <- amap::Dist(matriz_exp, method="kendall")
             }#,{
             #distancia_exp <- amap::Dist(matriz_exp, method="abspearson")#(1-abs(cor(t(matriz_exp), method="pearson"))
             #}
      )
      
      write.table(as.matrix(distancia_exp), file = paste(getwd(),pathBD,nameBD, "_",nameAlgorithm,"_matriz_distancia_expresion_",nombre_dist_exp,".csv",sep = ""),sep = ",",quote =FALSE)
    }else{
      #si ya existe el archivo matriz de distancia expresion lo lee
      distancia_exp<-stats::as.dist(read.table(paste(getwd(),pathBD,nameBD, "_",nameAlgorithm,"_matriz_distancia_expresion_",nombre_dist_exp,".csv",sep=""), header=T, sep=","))
    }
    
    
    
    
    
    ##########################################
    # calcula matriz distancia por Biologia
    ##########################################
    
    if (!file.exists(paste(getwd(),pathBD,nameBD, "_",nameAlgorithm,"_matriz_distancia_biologia_",nombre_dist_bio,".csv",sep = ""))){
      
      #carga lista de genes a mapear en caso de medulloblastoma
      lista_genes<-row.names(matriz_exp)
      
      # WANG requiere que los genes sean ENTREZID y mi bd los tiene como symbol, por eso los convierto
      if (organismoBD=="org.Hs.eg.db" && nameBD=="medulloblastoma" || nameBD=="serum") {
        mapearID<-AnnotationDbi::select(org.Hs.eg.db, keys=lista_genes,keytype="SYMBOL", columns=c("ENTREZID"))
        lista_genes<-mapearID[,2]
        lista_genes<-lista_genes[!is.na(lista_genes)]
      }
      
      
      switch(nombre_dist_bio, 
             wang={
               distancia <- mgeneSim(genes=lista_genes, ont = "BP", organism = organismoNAME, measure = "Wang",verbose = "TRUE")
             },
             lin={
               distancia<- mgeneSim(genes=lista_genes, ont = "BP", organism = organismoNAME, measure = "Lin",verbose = "TRUE")
             },
             jiang={
               distancia<- mgeneSim(genes=lista_genes, ont = "BP", organism = organismoNAME, measure = "Jiang",verbose = "TRUE")
             },                
             resnik={
               distancia<- mgeneSim(genes=lista_genes, ont = "BP", organism = organismoNAME, measure = "Resnik",verbose = "TRUE")
             }
      )
      
      
      distancia_bio<-stats::as.dist(apply(distancia,1,function(x) (1-(x)))) #1 mas cerca genes, convierto a 0 mas cerca genes
      distancia_bio<-as.matrix(distancia_bio)
      
      #vuelvo a asignar los nombres de genes originales por los que fueron mapeados  
      if (organismoBD=="org.Hs.eg.db" && nameBD=="medulloblastoma" || nameBD=="serum") {
        remapearID<-na.omit(mapearID)
        remapearID<-remapearID[which(!duplicated(remapearID[,2])),]
        remapearID<-remapearID[which(remapearID[,2] %in% rownames(distancia)),1]
        rownames(distancia_bio)<- remapearID
        colnames(distancia_bio)<- remapearID
      }
      
      #write.table(distancia_bio, file = paste(pathBD,nameBD, "_mono_matriz_distancia_biologia_",nombre_dist_bio,".csv",sep = ""),sep = ",",quote =FALSE)
      
      distancia_bio<-stats::as.dist(distancia_bio)
      
      
      ############################################################################
      # agregar genes que no estan reportados en la matriz de distancia biologica
      ############################################################################
      matriz_bio <- cmdscale(distancia_bio, k=ncol(matriz_exp)) # escala la matriz de distancia a puntos en ncol(expresion) coordenadas
      numero_genes_anotados<-nrow(as.matrix(distancia_bio))
      lista_genes_anotados<-row.names(as.matrix(distancia_bio))
      
      #si numero genes anotados en matriz distancia biologica es menor que los usados en la consulta
      if (nrow(matriz_exp)>numero_genes_anotados) {
        
        # esta linea solo interesa obtener la expresion, aunque despues eso lo convierto a 1 que es la peor distancia por biologia
        matriz_exp_rezagados<-matriz_exp[!(rownames(matriz_exp) %in% lista_genes_anotados), ] 
        
        original.points <- as.matrix(matriz_bio)                #genes con distancia biologica
        d0 <- as.matrix(distancia_bio)                          #matriz distancia genes con dist biologica
        extra.points <- matriz_exp_rezagados                    #genes sin distancia biologica
        
        inner.dist   <- rdist(extra.points)                     #simula distancia entre genes sin distancia biologica
        inner.dist[inner.dist > 0] <- 1 
        
        outer.dist   <- rdist(extra.points, original.points)    #simula distancia entre genes con y sin distancia biologica
        outer.dist[outer.dist > 0] <- 1                        
        #asigno 1 en todos los casos ya que eso indica que no hay relacion biologica entre los genes
        #con distancia biologica y los que no la tienen.
        
        d1 <- rbind(cbind(d0, t(outer.dist)), cbind(outer.dist, inner.dist))
        distancia_bio<-as.matrix(distancia_bio)
        row.names(d1)<-c(row.names(matriz_bio),row.names(matriz_exp_rezagados))
        colnames(d1)<-c(row.names(matriz_bio),row.names(matriz_exp_rezagados))
        
        #asigno la informacion unida como distancia biologica
        distancia_bio<-as.dist(d1)
        
      }
      
      ############################################################################
      # asigno datos finales a ambas matrices: expresion y biologia
      ############################################################################             
      
      
      distancia_exp<-as.matrix(distancia_exp)
      distancia_bio<-as.matrix(distancia_bio)
      
      #ordenar matriz de distancia biologica para que queden en el mismo orden los genes de ella y los de expresion.
      #es decir, que el gen1 en exp coincida con el gen1 en bio.
      distancia_bio<-distancia_bio[rownames(distancia_exp),rownames(distancia_exp),drop=FALSE]
      write.table(distancia_bio, file = paste(getwd(),pathBD,nameBD, "_",nameAlgorithm,"_matriz_distancia_biologia_",nombre_dist_bio,".csv",sep = ""),sep = ",",quote =FALSE)
      
      
      

    }else{
      #si ya existe el archivo matriz de distancia biologica lo lee
      distancia_bio<-as.dist(read.table(paste(getwd(),pathBD,nameBD, "_",nameAlgorithm,"_matriz_distancia_biologia_",nombre_dist_bio,".csv",sep=""), header=T, sep=","))
    }
    
    
    
    ############################################################################
    # asigno datos finales a ambas matrices: expresion y biologia
    ############################################################################             
    
    
    distancia_exp<-as.matrix(distancia_exp)
    distancia_bio<-as.matrix(distancia_bio)
    
    #ordenar matriz de distancia biologica para que queden en el mismo orden los genes de ella y los de expresion.
    #es decir, que el gen1 en exp coincida con el gen1 en bio.
    distancia_bio<-distancia_bio[rownames(distancia_exp),rownames(distancia_exp),drop=FALSE]
