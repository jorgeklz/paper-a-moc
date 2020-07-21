# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

#if (!require("Biobase"))      { BiocManager::install("Biobase")       ;     library("Biobase")        }   #carga funciones para manejo genico y anotaciones.
#if (!require("annotate"))     { BiocManager::install("annotate")      ; library("annotate")       }   #carga AnnotationDbi, stats4, Iranges, S4vectors, Xml
if(!require("GO.db"))         { BiocManager::install("GO.db")         ;     library("GO.db")          }   #carga DBI
if(!require("clValid"))       { BiocManager::install("clValid")       ;     library("clValid")        }   #implementa funcion BHI
if(!require("org.Sc.sgd.db")) { BiocManager::install("org.Sc.sgd.db") ;     library("org.Sc.sgd.db")  }   #anotaciones de Saccharomyces cerevisiae o sea yeast cell cycle 
#if(!require("org.Dm.eg.db"))  { BiocManager::install("org.Dm.eg.db")  ;     library("org.Dm.eg.db")   }   #anotaciones de drosophila o sea mosca
if(!require("org.Hs.eg.db"))  { BiocManager::install("org.Hs.eg.db")  ;     library("org.Hs.eg.db")   }   #anotaciones de homosapiens o sea humano
if(!require("org.At.tair.db")){ BiocManager::install("org.At.tair.db");     library("org.At.tair.db") }   #anotaciones de homosapiens arabidopsis thaliana
#if(!require("Rcpp"))      { install.packages('Rcpp', dependencies = TRUE)     ;     library("Rcpp")       }   #implementa similaridad semantica segun GO (Resnik, Schlicker, Jiang, Lin and Wang)
if(!require("GOSemSim"))      { BiocManager::install("GOSemSim")      ;     library("GOSemSim")       }   #implementa similaridad semantica segun GO (Resnik, Schlicker, Jiang, Lin and Wang)


if(!require("mco"))           { BiocManager::install("mco")           ;     library("mco")            } 
if(!require("emoa"))          { install.packages("emoa")  ;     library("emoa")           }   #calculo de hypervolume
if(!require("nsga2R"))        { BiocManager::install("nsga2R")        ;     library("nsga2R")         }  

if(!require("cluster"))       { BiocManager::install("cluster")       ;     library("cluster")        }   # implementa jerarquicos (diana, mona, agnes) y particionales (pam, clara, fanny). !!! no permite pearson
if(!require("amap"))          { BiocManager::install("amap")          ;     library("amap")           }   # implementa k-means y hclust y calcula pearson, euclidea, kendall y otras
if(!require("e1071"))         { BiocManager::install("e1071")         ;     library("e1071")          }   # implementa funciones varias, entre ellas fuzzy c-means. !!! no permite pearson
#if(!require("fpc"))           { biocLite("fpc")           ;     library("fpc")            }   # implementa DBSCAN, kmeans, pam, clara y metodos de estimaci?n de k.

#if(!require("pheatmap"))      { biocLite("pheatmap")      ;     library("pheatmap")       }   # implementa heatmaps buenos dise?os

#remove.packages(c("ggplot2")
#install.packages('Rcpp', dependencies = TRUE)
#install.packages('ggplot2', dependencies = TRUE)


if(!require("ggplot2"))       { install.packages("ggplot2", dependencies = TRUE)       ;     library("ggplot2")        }   # implementa plots buenos dise?os
if(!require("ggthemes"))      { install.packages("ggthemes");   library("ggthemes")     }   # ciertas configuracions de color ggplot

if(!require("reshape2"))      { BiocManager::install("reshape2")      ;     library("reshape2")       }   # implementa funciones para trasnformar formato de datos 
if(!require("miscTools"))     { BiocManager::install("miscTools")     ;     library("miscTools")      }   # implementa sumas, promedios por filas y columnas.

if(!require("matrixStats"))   { BiocManager::install("matrixStats")   ;     library("matrixStats")    }   # para calcular desv estandar por filas
if(!require("fields"))        { BiocManager::install("fields")        ;     library("fields")         }   
#if(!require("dendextend"))    { BiocManager::install("dendextend")    ;     library("dendextend")     }
if(!require("som"))           { BiocManager::install("som")           ;     library("som")            }   # implementa funcion para normalizar datos

if(!require("Rmisc"))         { BiocManager::install("Rmisc")         ;     library("Rmisc")          }
if(!require("clusterCrit"))   { BiocManager::install("clusterCrit")   ;     library("clusterCrit")    }   # impementa 50 indices de clustering, todos los que requiero
# getCriteriaNames(TRUE)
if(!require("fclust"))        { BiocManager::install("fclust")        ;     library("fclust")         }   
if(!require("vegan"))         { install.packages("vegan") ;     library("vegan")          }   # pa plot de soluciones de clustering
if(!require("ggfortify"))     { install.packages("ggfortify") ; library("ggfortify")  }   

if(!require("foreach"))       { install.packages("foreach")   ; library("foreach")    }   # implementa foreach para usar con paralelismo
if(!require("doParallel"))    { install.packages("doParallel"); library("doParallel") }   # implementa paralelismo

if(!require("plyr"))          { install.packages("plyr");       library("plyr") }  #para contar filas repetidas data.frame, lo uso par verificar convergencia
 

 
#if(!require("pheatmap")) biocLite("pheatmap")         # implementa heatmaps buenos dise?os
 
 