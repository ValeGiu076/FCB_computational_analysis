# FCB_computational_analysis
This is a method to perform immunophenotyping on high throughput Fluorescent Cell Barcoding data, and for quantification and minimization of inter-operator variability.
Debris are first removed using rectangleGate function in flowCore, and then flowClust applied for clustering. Lymphocyte and monocyte clusters are identified based on SSC.A values, because monocytes are bigger and more complex compared to lymphocytes. Then, unsupervised deconvolution is performed using FCB dye channels by flowClust. 
On each barcoded lymphocyte population, an unsupervised clustering is carried out for identification of CD3+ and CD20+ cells. The population with the highest Maximum value of PerCP.Cy5.5.A is assigned as CD3+ cells, and CD4+ and CD8+ cells are further gated. 
Similarly, after deconvolution of monocytes by flowClust, CD14+ cells are identified using linear parameter and CD14 expression by rectangleGate in flowCore by fixing the cut-off at the 1st quartile of the PE.A parameter. 
