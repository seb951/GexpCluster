
library(cluster)
library(DT)
library(wesanderson)



wes_colors = c(wes_palettes$GrandBudapest1,wes_palettes$GrandBudapest2,wes_palettes$Zissou1,wes_palettes$Rushmore)
###


#' render  twss
#'
#' render  twss
#' @return render
#' @examples
#' eg = renderstwss();
#' @export
renderstwss <- function(multi_k =  clustering_statistics(data[[2]]),
                        colors = wes_colors)
{
  plot(multi_k[,1],multi_k[,3], type='b', col=colors, pch = 19, lwd = 3,main="Total within sum of square \n gene expression dataset", xlab='Number of clusters (k)', ylab='Total Within Sum of Square', frame=FALSE)
}


#' render  silhouette
#'
#' render  silhouette
#' @return render
#' @examples
#' eg = rendersilhouette();
#' @export
rendersilhouette <- function(multi_k =  clustering_statistics(data[[2]]),
                             colors = wes_colors)
  {
  plot(multi_k[,1],multi_k[,2], type='b', col=colors, pch = 19, lwd = 3,main="Silhouette score \n gene expression dataset", xlab='Number of clusters (k)', ylab='Average Silhouette Scores', frame=FALSE)
}



#' render detailed silhouette
#'
#' render detailed silhouette
#' @return render
#' @examples
#' eg = render_detailed_silhouette();
#' @export
render_detailed_silhouette = function(input_data = data[[2]],k = 2,colors = wes_colors,type=c("kmeans","dendrogram")){
  if(type == "kmeans"){
    input_pca_data.pca = prcomp(input_data[,1:ncol(input_data)], center = TRUE,scale. = TRUE)
    km <- kmeans(input_pca_data.pca$x[,1:10], centers = k, nstart=25)
    ss <- silhouette(km$cluster, dist(input_pca_data.pca$x))
    rownames(ss) = names(km$cluster)
  }

  if(type == "dendrogram"){
    cl <- hclust(dist(input_data))
    ss <- silhouette(cutree(cl, k=k) ,dist(input_data), title=title(main = 'Good'))
    rownames(ss) = cl$labels
  }
  plot(ss,col = colors[1:k],max.strlen=20,nmax.lab = 200,cex.names = 0.2,main = paste0("Detailed Silhouette plot for k = ",k))
}



#' render pca
#'
#' render pca
#' @return render
#' @examples
#' eg = renderpca();
#' @export
renderpca = function(input_data=data[[2]],pcX=1,pcY=2,k = 2,colors = wes_colors,type="kmeans")
  {
  #kmeans
  input_pca_data.pca = prcomp(input_data[,1:ncol(input_data)], center = TRUE,scale. = TRUE)
  if(type == "kmeans"){
    km <- kmeans(input_pca_data.pca$x[,1:10], centers = k, nstart=25)
  }
  if(type == "dendrogram"){
    cl <- hclust(dist(input_data))
    km <- list(cluster = cutree(cl, k=k))
  }


  colors= colors[1:length(unique(km$cluster))]

  colors_val <- lapply(data.frame(km$cluster), function(x) as.character(factor(km$cluster, labels=colors)))
  colors_val = as.character(colors_val$km.cluster)
  eigen_vector = input_pca_data.pca$sdev^2
  pve= round(eigen_vector/sum(eigen_vector)*100,2)
  plot(input_pca_data.pca$x[,pcX],input_pca_data.pca$x[,pcY],xlab = paste0("PC",pcX," (",pve[pcX]," %)"),ylab =paste0("PC",pcY," (",pve[pcY]," %)"),col=colors_val,pch=19,lwd=3,main = paste0("Principal Component Analysis \n (grouping according to ",type,")"))
  #biplot(data,choices=c(pcX,pcY))
  }




#' render summary data
#'
#' render summary data
#' @return render
#' @examples
#' eg = render_summary_data();
#' @export
render_summary_data = function(clinical = clinical,variable = colnames(clinical),colors = wes_colors) {
  if(variable == c("survival"))
  {
    hist(clinical[,variable],breaks =10,col = colors,xlab = 'days of survival',ylab = "Number of patients",main = "Survival rates since diagnosis")
  }

  if(variable == c("age"))
  {
    hist(clinical[,variable],breaks =10,col = colors,xlab = 'Age (years)',ylab = "Number of patients",main = "Age at diagnosis")
  }

  if(variable == c("stage"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(stage = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ stage, data = summary,main ="stage",col = colors,ylab = "Number of patients")
  }

  if(variable == c("histology"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(histology = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ histology, data = summary,main ="histology",col = colors,ylab = "Number of patients")

  }

  if(variable == c("gender"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(gender = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ gender, data = summary,main ="gender",col = colors,ylab = "Number of patients")

  }

  if(variable == c("kmeans"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(kmeans = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ kmeans, data = summary,main ="Clusters kmeans",col = colors,ylab = "Number of patients")

  }

  if(variable == c("hierarchical"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(hierarchical = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ hierarchical, data = summary,main ="Clusters hierarchical",col = colors,ylab = "Number of patients")

  }
}

#' Hierarchical clustering plot
#'
#' Convert a dataframe of Gexpr to a heatmap
#' @param max_genes max nb of genes to show
#' @param gexp gexp df
#' @return the heatmap
#' @examples
#' heatmap = hc();
#' @export
hc = function(gexp = data[[2]],max_genes=200){
  gexp_max = t(gexp)[1:max_genes,]
  stats::heatmap(gexp_max,main=list("Heatmap (genes X samples)",cex = 1.1))
  }

#' switch clustering plots
#'
#' switch clustering plots
#' @return plots
#' @examples
#' eg = switch_clustering_plots();
#' @export
switch_clustering_plots = function(multi_k = multi_k, data = data,clustering_metrics = '',k = 2) {
  switch(clustering_metrics,
       silhouette = rendersilhouette(multi_k = multi_k),
       detailed = render_detailed_silhouette(input_data = data[[2]],k = k,type = 'kmeans'),
       twss = renderstwss(multi_k = multi_k),
       silhouette_hc = rendersilhouette(multi_k = multi_k),
       detailed_hc = render_detailed_silhouette(input_data = data[[2]],k = k,type = 'dendrogram'),
       twss_hc = renderstwss(multi_k = multi_k)
  )

}

