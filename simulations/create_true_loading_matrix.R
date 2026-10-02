# Create factor loading matrix
lambda = matrix(0, 12, 3)
lambda[1:6, 1] = 1
lambda[7:10, 2] = 1
lambda[11:12, 3] = 1

library(pheatmap)
#plotting function
plot_pheatmap = function(m, name='heatmap', breaks=NULL, save_name=""){
  if(!save_name==""){
    try(dev.off())
    pdf(save_name)
  }
  if(is.null(breaks)){
    m = m + matrix(rnorm(prod(dim(m)),0, 0.00001), dim(m)[1])
    pheatmap(m, cluster_rows = FALSE, cluster_cols=FALSE, main=name)}
  else{
    pheatmap(m, breaks=breaks, cluster_rows = FALSE, cluster_cols=FALSE, main=name)}
  if(!save_name==""){
    try(dev.off())
  }
}

plot_pheatmap(lambda, name = 'True Loading Matrix for simulations 1-2.', save_name='true_factor_loading_matrix.pdf')
