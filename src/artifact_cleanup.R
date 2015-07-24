### Translated by Timothy Tickle and Brian Haas
### Note, code is translated from matlab code written by Nir Yosef
### Special thanks to Nir Yosef


# Constants
i_var = 90
### Clean up covers 90% of the variance in the covariate matrix
i_number_bins = 10
### Number of bins for global-scaling normalization
i_number_of_pc = 10
### Number of PC considered


compute_pc_k = function(
### Calculate the number of PCs to select
df_covariate_matrix,
### (complexity, 3-5' bias etc.), rows are cells, columns are different covariates
){
  # perform principle components analysis
  results_princomp = princomp( df_covariate_matrix )
  mtrx_loadings = results_princomp$loadings
  vctr_variance = results_princomp$sdev^2

  # Get the projection of the data
  mtrx_proj = df_covariate_matrix * mtrx_loadings

  # Get K where cumulative sum is above 90 %
  d_cumsum = ( cumsum( vctr_variance ) / sum( vctr_variance ) ) * 100
  i_K = which( d_cumsum > 90 )[1]

  # Reduce projection down to the selected components
  mtrx_proj = mtrx_proj[,1:i_K]

  # Output to user
  print( paste( "Covering ",d_cumsun[ i_k ], " percent of variance (K=",i_K,")" ) )

  return( list( K=i_K, proj=mtrx_proj ) )
}


current_status_pca = function(
df_data,
mtrx_covariates
){
  # Calculate standardized values
  df_zscore_data = scale( df_data )

  # perform principle components analysis
  results_princomp = princomp( df_zscore_data )
  mtrx_loadings = results_princomp$loadings

  # project data to the predefined number of components
  mtrx_projected = df_zscore_data * mtrx_loadings[,1:3]

  # Get the minimum correlation
  vectr_cor = cor( mtrx_projected , mtrx_covariates )
  return( min( min( vctr_cor ) )
}


artifact_cleanup = function(
### Perform Global scaling normalization or PCA shaving for normalizing RNA-seq data
df_data,
### Data data frame, columns are cells, rows are genes
df_covariate_matrix,
### (complexity, 3-5' bias etc.), rows are cells, columns are different covariates
i_normalization_option = 1,
### Normalization option 1 (PCA shaving) or 2 (Global-scaling normalization)
f_use_mean = FALSE
### Toggles between using mean or median for measures of centrality; used in global-scaling normalization
){
  # First select the number of K beased on the covariates
  results_pca_k = compute_pc_k( df_covariate_matrix ) 

  # Select cleanup method
  if( i_normalization_option == 1 )
  {
    return( pca_cleanup( df_data, results_pca_k$proj, 0.00001 ) )
  } else if ( i_normalization_option == 2 ) {
    return( global_scaling( ) )
  } 
}


global_scaling = function(
df_data,
mtrx_covariates,
i_K,
i_number_bins
){
  # Select function for determining centrality
  func_eval = mean
  if( ! f_use_mean )
  {
    func_eval = median
  }

  # Data copies
  df_best = df_data
  df_raw = df_data

  # Get the min p
  d_best_p = current_status_pca( df_best, mtrx_covariates[,1:i_K] )
  print( paste( "0: ", d_best_p ) )

  for( i_bin in i_number_bins )
  {
    # Make the sequence of bins and chenge the endigs to inf to catch any
    # Unexpected values
    vctr_quantiles = quantile( mtrx_covariates, seq( 1/i_bin, 1,1/i_bin ) )
    vctr_quantiles[0] = -Inf
    vctr_quantiles[ length( vctr_quantiles ) ] = Inf
    i_quantile_length = length( vctr_quantiles )

    # PC in selected PCs 
    for( i_index_k in 1:i_K )
    {
      vctr_centers = repmat( func_eval( df_data, 2 ), 1, nrow( df_data ) )
      for( i_index_quantiles in 1: ( i_quantile_length-1 ) )
      {
        vctr_selected = mtrx_covariates[,i_index_k] >= vctr_quantiles[i_index_i,i_index_k]
        vctr_selected = vctr_selected and ( mtrx_covariates[,i_index_k] < vctr_quantiles[i_index+1, i_index_k]  )
        i_count_selected = length( which( vctr_selected ) )
        df_data[,vctr_selected] = df_data[,vctr_selected] - repmat( func_eval( df_data[,vctr_selected],1, i_count_selected ) + vctr_centers[,1:i_count_selected] ) 
      }
    }
    
  } 
  i_min = current_status_pca( df_data, mtrx_covariates[,1:i_K] )
  if ( i_min > i_best )
  {
    d_best_p = i_min
    d_best_p = df_data
  }
  print( paste( i_bin,": ", i_min ) )
  return( df_best )
} 


pca_cleanup = function(
### Clean up by PCA
df_data,
### Data data frame, columns are cells, rows are genes
df_covariate_matrix,
### (complexity, 3-5' bias etc.), rows are cells, columns are different covariates
d_cutoff = 0.001,
### P-value cutoff of correlation with covariate vectors below which the PC will be removed.
i_number_PC = 10
### Number of PCs considered
){
  # Calculate mean
  vctr_mean = mean( df_data )
  # Calculate standard deviation
  vctr_std = sd( df_data )

  # Calculate standardized values
  df_zscore_data = t( scale( df_data ) )

  # perform principle components analysis
  results_princomp = princomp( df_zscore_data )
  mtrx_loadings = results_princomp$loadings
  mtrx_scores = results_princomp$scores
  vctr_eigen_values = eigenvals( results_princomp )

  # project data to the predefined number of components
  mtrx_projected = df_zscore_data * mtrx_loadings[,1:i_number_of_pc]

  # Get the p-values of the correlations with projection and covariates
  vctr_pvalues = c()
  for( i_index in 1:i_number_PC)
  {
    vctr_pvalues = c( vctr_pvalues, cor.test( mtrx_projected[ ,i_index ], df_covariate_matrix )$p.value )
  }
  vctr_sig_pvalues = vctr_pvalues < d_cutoff
  i_sig_pvalues_count = length( which( vctr_sig_pvalues ) )

  # Remove what can be explianed in the relationship between the loadings and the zscored data.
  df_data_return = df_score_data - t( ( matrix_loadings[, vctr_sig_pvalues ]*solve( matrix_loadings[ , vctr_sig_pvalues ], t( df_zscore_data ) ) ) ) 
  print( paste( "Removed", i_sig_pvalues, "PC with", 100*( 1-explained_var( matrix_loadings, df_data_return ), "percent variance." ) )

  df_data_return = t( df_data_return )
  df_data_return = ( df_data_return * repmat( vctr_std, nrow( df_data_return ), 1  ) ) + repmat( vctr_mean, nrow( df_data ), 1 )
  print( paste( "After transformation:", 100*explained_var( df_data , df_data_return ), " remaining." ) )
  return list( data=df_data_return, sig_count=i_sig_pvalues )
}

repmat = function(X,m,n){
### R equivalent of repmat found at http://haky-functions.blogspot.com/2006/11/repmat-function-matlab.html
  mx = dim(X)[1]
  nx = dim(X)[2]
  return( matrix( t( matrix(X,mx,nx*n ) ), mx*m, nx*n, byrow=TRUE ) )
}
