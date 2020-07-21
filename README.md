# a-moc
It is the implementation of a-MOC (Multi-Objective Gene Clustering Algorithm) which is used to cluster gene based on the combination of similarities of expression profile (EP) and biological functional profile (BP) between a pair of gene, using an alpha parameter.

This is a R project, to replicate experiments of paper {INFLUENCE OF THE GO-BASED SEMANTIC SIMILARITY MEASURES IN MULTI-OBJECTIVE GENE CLUSTERING ALGORITHM PERFORMANCE} you should run a-MOC_main.R 

If you want changing any parameter or configuration such as dataset, GO-based measure, number of cluster and others, you should modify a-MOC_parametros.R

The repository includes the dataset, expression and biological similarities matrix, and functional classes. So there is no need to recalculate this. At the end of the execution, the algorithm will save the results in the Results folder.
