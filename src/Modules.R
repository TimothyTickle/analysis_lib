# Constant methods
# Note if these change you have to change the switches
# at the bottom that select modules based on keywords
STR_ANOVA = "anova"
STR_SAMPLE_DENDROGRAM = "dendrogram"
STR_PCA = "pca"
STR_NMDS = "nmds"


# Key words for modules listed below
vctr_batch_methods = c( "none", "combat", "edgeR", "regression", "PCA_skimming", "ruvseq" )
vctr_describe_feature_groupings_functionally = c()
vctr_discriminate_genes_methods = c( "none", STR_ANOVA )
vctr_filters_methods = c( "none", "percentile", "sd", "downsample", "top", "occurence", STR_PCA, "sparsity" )
vctr_impute_methods = c( "none", "complete_case", "mean", "median", "multiple_imputation" )
vctr_normalization_methods = c( "none", "cpm", "upper_quartile" )
vctr_ordination_methods = c( "none", STR_PCA, STR_NMDS, "scde" )
vctr_select_sample_groups_methods = c( "none", STR_SAMPLE_DENDROGRAM, "mclust" )
vctr_select_feature_group_methods = c( "none", STR_SAMPLE_DENDROGRAM )
vctr_select_transform = c( "none", "asin_sqrt", "center", "log", "tss", "zscale" )


##################################
# Description of data
##################################


func_describe_metadata <- function(
### Describes the matrix of metadata
df_metadata,
### Output project directory
str_output_dir
){
  loginfo( "Describe metadata" )

  # Manage output directory for PDF
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }
  pdf( file.path( str_output_dir, "describe_metadata.pdf"), useDingbats=FALSE )

  # Counts of metadata
  vctr_missing = apply( df_metadata, 2, function(x) length( union(which(is.na(x)), which(is.nan(x)))))
  vstr_labels = paste( names(df_metadata), "(", nrow( df_metadata ), ")", sep = " ")
  barplot( vctr_missing, names.arg = vstr_labels, ylab = "Count of missing values", xlab = "Metadata", main = "Missing Values per Metadata" )
  dev.off()

  # Metadata associations
  if( sum( vctr_missing ) > 0 )
  {
    i_ncol = ncol( df_metadata )
    i_nrow = nrow( df_metadata )
    df_metadata_discrete = as.data.frame( matrix(rep(0, i_ncol * i_nrow), nrow = i_nrow, ncol = i_ncol ) )
    names( df_metadata_discrete ) = names( df_metadata )
    row.names( df_metadata_discrete ) = row.names( df_metadata )
    df_metadata_discrete[ is.na( df_metadata ) ] = 1
    dist_cor = 1 - func_get_custom_distance_matrix( t(df_metadata_discrete), "jaccard" )
    df_metadata_discrete =  data.frame( t(combn( names( df_metadata_discrete ), 2 ) ), as.numeric( dist_cor ))
    func_write_text( "Similarity of Missing Data\n\n", file.path( str_output_dir, "describe_metadata.txt" ))
    func_write_frame( df_metadata_discrete, file.path( str_output_dir, "describe_metadata.txt" ), TRUE)
  }
}


func_describe_count_data_matrix <- function(
### Make basic plots to describe the data
obj_analysis,
### Holds the data associated with the current analysis
str_output_dir,
### Directory to output data and visualizations
str_cor_metric = "euclidean",
### Correlation metric to correlate samples
vctr_color_groups,
### Coloring for sample metadata
...
){
  loginfo( "Describe Count Data. Please note assumptions are made here that the data is count." )

  # Manage output directory for PDF
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }
  pdf( file.path( str_output_dir, "describe_count_data.pdf"), useDingbats=FALSE )

  # Basic counts
  plot.new()
  func_describe_data_frame_values( obj_analysis )

  # Heat map raw data distribution
  # Add metadata as columns
  func_heatmap( obj_analysis, str_linkage = "average", vctr_factor_groups = NULL, plt_colors = func_monochromatic_palette( ),... ) 

  # Heat map correlation matrix
  # TODO diag = NA for copy matrix
  func_cor_matrix( obj_analysis, str_linkage = "average", f_correlate_samples = TRUE, f_log = TRUE, ... )
  func_cor_matrix( obj_analysis, str_linkage = "average", f_correlate_samples = FALSE, f_log = TRUE, ... )

  # Sample distributions
  boxplot( get.logged.matrix( obj_analysis ), main = "Sample Distributions (log2)" )

  # SD vs Mean
  func_feature_sd_vs_mean( obj_analysis )

  # Distribution by sparsity
  func_feature_distributions_by_sparsity( obj_analysis )

  # Plot by metadata
  # TODO Use metadata properly
  if( ! is.null( metadata( obj_analysis ) ) )
  {
    par( mar = c( 5.1, 4.1, 4.1, 8.1 ), xpd = TRUE )
    for( str_metadata in names( metadata( obj_analysis ) ) )
    {
      # Make factor data
      vctr_groups = metadata( obj_analysis )[[ str_metadata ]]
      vctr_grouping_colors = rep( NA, length( vctr_groups ) )
      if( is.numeric( vctr_groups ) || is.integer( vctr_groups ) )
      {
        vctr_groups = cut( vctr_groups, min( 10, length( vctr_groups ) ) )
      } else {
        vctr_groups = factor( vctr_groups )
      }

      # Change to colors
      lret_grouping_colors = func_factor_to_metadata_color( vctr_groups )

      # Add color to the plots
      boxplot( get.logged.matrix( obj_analysis ), main="Sample Distribution (log2)", col = lret_grouping_colors$vctr_grouping_colors, las = 2 )
      legend( "topright", inset = c( -0.2, 0 ), title = str_metadata, legend = lret_grouping_colors$vctr_levels, col = lret_grouping_colors$vctr_colors, pch = 16 )
      func_feature_sd_vs_mean( obj_analysis, lret_grouping_colors$vctr_colors )

      # Heat map raw data distribution
      # Add metadata as columns
      func_heatmap( obj_analysis, str_linkage = "average", plt_colors = func_monochromatic_palette( ), vctr_grouping = vctr_groups, f_log = TRUE, ... ) 
      # Heat map correlation matrix
      func_cor_matrix( obj_analysis, str_linkage = "average", f_correlate_samples = TRUE, vctr_grouping = vctr_groups, f_log = TRUE, ... )
      func_cor_matrix( obj_analysis, str_linkage = "average", f_correlate_samples = FALSE, f_log = TRUE, ... )
    }
  }
  dev.off()
}


func_describe_data_frame_values = function(
### Prints basic counts and stats about a data frame / matrix
obj_analysis
### Holds the data associated with the current analysis
){
  # Calculate basic descriptions
  i_nrows = nrow( data( obj_analysis ) )
  i_cols = ncol( data( obj_analysis ) )
  vctr_row_counts = apply( data( obj_analysis ), 1, sum )
  vctr_col_counts = apply( data( obj_analysis ), 2, sum )
  i_zero_rows = length( which( vctr_row_counts == 0 ) )
  i_zero_cols = length( which( vctr_col_counts == 0 ) )

  # Report values
  i_line_increment = 0
  mtext( "Describing count data:", side=3, adj = 0 ,line = i_line_increment )
  i_line_increment = i_line_increment - 2
  str_line = paste( "Number of Rows:", i_nrows )
  mtext( str_line, side = 3, adj = 0, line = i_line_increment )
  str_line = paste( "Number of Rows with no measurments:", i_zero_rows )
  i_line_increment = i_line_increment - 1
  mtext( str_line, side = 3, adj = 0, line = i_line_increment )
  str_line = paste( "Number of columns:", i_cols )
  i_line_increment = i_line_increment - 1
  mtext( str_line, side = 3, adj = 0, line = i_line_increment )
  str_line = paste( "Number of columns with no measurements:", i_zero_cols )
  i_line_increment = i_line_increment - 1
  mtext( str_line, side = 3, adj = 0, line = i_line_increment )
}


######################################
# Data transformation
######################################

func_zscale <- function(
### Mean centers and standardizes by sd to make norm values
### Transforms within columns
df_data
### Data to be transformed
){
  return( data.frame( scale( df_data, center = TRUE, scale = TRUE ) ) )
}


func_center <- function(
### Mean centers data
### Transforms within columns
df_data
### Data to be transformed
){
  return( data.frame( scale( df_data, center = TRUE, scale = FALSE ) ) )
}


func_asin_sqrt <- function(
### Asin sqrt transform for proportional data exhibiting heteroscedasticity
### Expects propotional data
### Transforms within columns
df_data
### Data to be transformed
){
  return( asin( sqrt( df_data )  ) )
}


######################################
# Filters
######################################

func_filter_downsample <- function(
### Performs downsampling on a full data frame
### Downsamples all samples to the lowest sample depth among them
obj_analysis,
### Holds the data associated with the current analysis
str_output_dir = NA,
### Output directory to write results to
...
){
  loginfo( "Downsample counts" )

  df_cur_data = data( obj_analysis )

  # Downsample to the lowest sample depth
  i_min_depth_sample = min( colSums( df_cur_data ) )

  # If plotting is needed, make pdf
  f_plotting = !is.na( str_output_dir )
  if( f_plotting )
  {
    pdf( file.path( str_output_dir, "downsampled_populations.pdf" ), useDingbats = FALSE )
  }

  for( i_index in 1:ncol( df_cur_data ) )
  {
    # Plot before distribution
    if( f_plotting )
    {
      vctr_cur_data = df_cur_data[[ i_index ]]
      barplot( vctr_cur_data, main=paste( "Column", i_index, "(Before and After)"), col = STR_COLOR_BEFORE )
    }
    # Downsample each column (sample)
    df_cur_data[[ i_index ]] = func_downsample( df_cur_data[[ i_index ]], i_min_depth_sample )
    # Plot after distribution
    if( f_plotting )
    {
      barplot( df_cur_data[[ i_index ]], main="", col = STR_COLOR_AFTER, add=TRUE )
      legend( "topright", c("Before","After"), bg = "white", fill = c( STR_COLOR_BEFORE, STR_COLOR_AFTER ))
    }
  }

  if( f_plotting )
  {
    dev.off()
  }

  data( obj_analysis, islogged = is.logged( obj_analysis ) ) = df_cur_data
  return( obj_analysis )
}


func_downsample <- function(
### Down samples counts while attempting to keep the distribution shape
### This will not up sample a sample size; in that scenario it will return the original distribution.
### This works on one sample distribution
vctr_i_sample_counts,
### Original feature ( distrbution ) to be downsampled
i_target_population_sum = 100,
### Target total counts in distribution after sampling
...
){
  loginfo( paste( "Downsample to ", i_target_population_sum, " total counts." ) )

  # Get the sum of the population vector
  i_population_sum = sum( vctr_i_sample_counts )
  i_population_length = length( vctr_i_sample_counts )

  # If more counts are requested than exist just return the original distribution
  if( i_target_population_sum > i_population_sum )
  {
    return( vctr_i_sample_counts )
  }

  # Check and make sure the length of the sample vector is more than
  # one measurement or the sample method will act up
  # If there is only one value in the distribution return it or the amount of sampling requeted
  if( i_population_length == 1 )
  {
    return( min( i_target_population, vctr_i_sample_counts ) )
  }

  # Sample each instance and count instances (positionally)
  vctr_downsampling = rep( 0, i_population_length)
  vi_index_counts = c()
  
  for(i_index in 1:length( vctr_i_sample_counts ))
  {
    vi_index_counts = c(vi_index_counts,rep( i_index, floor(vctr_i_sample_counts[ i_index ] ) ) )
  }
  vi_counts = table( sample( vi_index_counts, size = i_target_population_sum, replace = FALSE ) )
  vctr_downsampling[ as.integer( names( vi_counts ) ) ] = vi_counts
  return( vctr_downsampling )
}


func_filter_extreme_by_pca <- function(
### Select by PCA
obj_analysis,
### Holds the data associated with the current analysis 
i_PC_Components = NULL,
### Number of PC components
i_count = NULL,
### The number of top ranked features to select
str_output_dir = NA
### Output directory to document to
){
  loginfo( "Feature by PCA (selecting extreme features)" )
  # Row center and log
  mtrx_scale = t( scale( t( get.logged.matrix( obj_analysis ) ), center=TRUE, scale=TRUE ) )
  # Remove constant rows
  mtrx_scale = mtrx_scale[ !is.na( mtrx_scale[, colnames( mtrx_scale )[ 1 ] ] ), ]
  # Perfrom PCA
  results_pca = princomp( mtrx_scale, cor = TRUE )
  # For each component get the top genes
  mtrx_scores = results_pca$scores
  # How many genes exist
  i_gene_counts = nrow( mtrx_scores )
  # Get percent variance
  vctr_i_variance = results_pca$sdev^2
  vctr_i_percent_variance = ( vctr_i_variance ) / sum( vctr_i_variance )
  i_eigen_count = length( results_pca$sdev )

  # If check output directory path
  if( is.na( str_output_dir ) )
  {
    str_output_dir = "."
  }

  # If no count is given us the upper quartile expressed
  if( is.null( i_count ) )
  {
    i_count = floor( i_gene_counts * 0.25 )
  }

  # Get the stick break distribution
  if( is.null( i_PC_Components ) )
  {
    result_stick = func_select_pca_components_stick_breaking( vctr_i_percent_variance = vctr_i_percent_variance, str_output_dir = str_output_dir )
    i_PC_Components = result_stick$count
  }

  # Holds the selection per component
  vctr_return_feature_list = c()

  # For each component select the top genes
  for( i_index_component in 1:i_PC_Components )
  {
    vctr_genes = row.names( mtrx_scores[ order( abs( mtrx_scores[,i_index_component] ) ), ] )[ 1 : i_count ]
    if( ! is.na( str_output_dir ) )
    {
      func_write_frame( vctr_genes, file.path( str_output_dir, paste( "Extreme",i_count, "genes","by","PCA_component",i_index_component,".txt", sep="_" ) ) )
    }
    vctr_return_feature_list = unique( c( vctr_return_feature_list, vctr_genes ) )
  }

  data( obj_analysis, islogged = is.logged( obj_analysis ) ) = data( obj_analysis )[vctr_return_feature_list,]
  return( obj_analysis )
  ### Keep only the top extreme features in loadings in PCA
  ### Return the analysis object with only selected feature in the data
}


func_filter_by_SD <- function(
### Feature filter. Filter by standard deviation
obj_analysis,
### Holds the data associated with the current analysis
str_output_dir = NA,
### Directory to save plots and such
...
){
  loginfo( "Filter by SD" )  

  # Filter rows by variance
  d_min_sd = 0.1

  # Filter
  vctr_sd = apply( data( obj_analysis ), 1, function( x ){ sd( x, na.rm = TRUE ) } )
  vctr_sd_no_zero = apply( data( obj_analysis ), 1, function( x ) { sd( x[ x > 0 ], na.rm = TRUE) } )
  vctr_f_keep = vctr_sd >= d_min_sd

  if( !is.na( str_output_dir ) )
  {
    vctr_colors = rep( STR_COLOR_BEFORE, length( vctr_sd ) )
    vctr_colors[ vctr_f_keep ] = STR_COLOR_STANDARD

    pdf( file.path( str_output_dir, "sd_filter.pdf" ), useDingbats=FALSE )
    plot( sort( vctr_sd ), main = "SD ( with zeros )", col = vctr_colors[ order( vctr_sd ) ] )
    plot( sort( vctr_sd_no_zero ), main = "SD ( with out zeros )", col = vctr_colors[ order( vctr_sd ) ] )
    dev.off()
  }

  data( obj_analysis, islogged = is.logged( obj_analysis ) ) = data( obj_analysis )[ vctr_f_keep ,]
  return( obj_analysis )
}


func_filter_by_first_rows <- function(
### Feature filter. Filter by giving the first N rows sorted by expression (top)
obj_analysis,
### Holds the data associated with the current analysis
d_top_rows = 50,
### Top number of rows ( in order ) to keep
str_output_dir = NA,
### Directory to save plots
...
){
  loginfo( "Filter by top (n first top rows )" )

  # Head data frame
  vctr_ordered_by_expression = order( apply( data( obj_analysis ), 1, sum ), decreasing = TRUE )[ 1:d_top_rows ]
  df_ret = data( obj_analysis )[ vctr_ordered_by_expression , ]

  plt_colors = func_monochromatic_palette()

  # Optionally plot
  if( !is.na( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "top_filter.pdf" ) )
    # Heatmap before
    func_heatmap( obj_analysis, str_linkage = "average", plt_colors = plt_colors, str_title = "Data before top filter", f_log = TRUE, ... )
  }

  # Save abridged data frame
  data( obj_analysis, islogged = is.logged( obj_analysis ) ) = df_ret

  # Optionally plot
  if( !is.na( str_output_dir ) )
  {
    # Heatmap after
    func_heatmap( obj_analysis, str_linkage = "average", plt_colors = plt_colors, str_title = "Data after topN filter", f_log = TRUE, ... )
    dev.off() 
  }
  return( obj_analysis )
}


func_filter_by_min_value_occurence <- function(
### Feature filter. Filter by min value in min number samples
obj_analysis,
### Holds the data associated with the current analysis
d_min_value = 10,
### Minimum value a sample must have
d_min_occurence = NULL,
### Minimum times the value must occur in a feature to be kept
str_output_dir = NA,
### Directory to save plots
...
){
  loginfo( "Filter by min value occurence" )

  # If the min occurence is not given, make it 10% of the current data frame size (min 1)
  if( is.null( d_min_occurence ) )
  {
    d_min_occurence = max(1, ceiling( ncol( data( obj_analysis ) ) * .1 ))
  }
  # True of false that a feature passes the filter
  vctr_f_keep = apply( data( obj_analysis ), 1, function( x ){ return( length( which( x >= d_min_value ) ) >= d_min_occurence ) } )

  # Optionally plot
  if( !is.na( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "value_occurence_filter.pdf" ) )
    plt_colors = func_monochromatic_palette()

    # Calculate values
    vctr_min_value_occurence_keep = apply( data( obj_analysis )[ vctr_f_keep,], 1, function( x ){ return( length( which( x >= d_min_value ) ) ) })
    vctr_min_value_keep = apply( data( obj_analysis )[ vctr_f_keep,], 1, function( x ){ return( min( x[x >= d_min_value] ) ) })
    vctr_min_value_occurence = apply( data( obj_analysis )[ ! vctr_f_keep,], 1, function( x ){ return( length( which( x >= d_min_value ) ) ) })
    vctr_min_value = apply( data( obj_analysis )[ ! vctr_f_keep,], 1, function( x ){ return( min( x[x >= d_min_value] ) ) })
    i_length_keep = length( vctr_min_value_occurence_keep ) 
    i_length_no_keep = length( vctr_min_value_occurence ) 

    # Plot distributions
    if( ( i_length_keep == 0 ) && ( i_length_no_keep > 0 ) )
    {
       hist( vctr_min_value, col = STR_COLOR_BEFORE,
          main = "Filtering by min value occurence ( All Filtered ) ", xlab = "Min value", ylab = "Occurence across samples" )
    } else if( ( i_length_keep > 0 ) && ( i_length_no_keep == 0 ) )
    {
       hist( vctr_min_value_keep, col = STR_COLOR_STANDARD,
          main = "Filtering by min value occurence ( None Filtered ) ", xlab = "Min value", ylab = "Occurence across samples" )
    } else {
      ret_hist1 = hist( vctr_min_value_keep, plot = FALSE )
      ret_hist2 = hist( vctr_min_value, plot = FALSE )
      d_xlim = c( min( vctr_min_value_keep, vctr_min_value ), max( vctr_min_value_keep, vctr_min_value ) )
      d_ylim = c( 0, max( max( ret_hist1$counts ), max( ret_hist2$counts ) ) )
      vctr_d_breaks = hist( c( vctr_min_value_keep, vctr_min_value ), plot = FALSE )$breaks
      hist( vctr_min_value_keep, col = STR_COLOR_STANDARD, breaks = vctr_d_breaks, xlim = d_xlim, ylim = d_ylim,
          main = "Filtering by min value occurence", xlab = "Min value", ylab = "Occurence across samples" )
      hist( vctr_min_value, col = STR_COLOR_BEFORE, add = TRUE )
    }
    legend( "topright", title = "Filtering", legend = c("Filtered","Kept"), fill = c( STR_COLOR_BEFORE, STR_COLOR_STANDARD ) ) 

    # Plot the values of min value vs min occurence, indicating who was kept
    if( i_length_keep > 0 )
    {
      plot( vctr_min_value_keep, vctr_min_value_occurence_keep, col = STR_COLOR_STANDARD, main = "Filtering by min value occurence", xlab = "Min value", ylab = "Occurence across samples", pch = 16 )
      points( vctr_min_value, vctr_min_value_occurence, col = STR_COLOR_BEFORE, pch = 16 )
    } else {
      plot( vctr_min_value, vctr_min_value_occurence, col = STR_COLOR_BEFORE, main = "Filtering by min value occurence", xlab = "Min value", ylab = "Occurence across samples", pch = 16 )
    }
    legend( "topright", title = "Filtering", legend = c("Filtered","Kept"), fill = c( STR_COLOR_BEFORE, STR_COLOR_STANDARD ) ) 

    # Heatmap before
    func_heatmap( obj_analysis, str_linkage = "average", plt_colors = plt_colors, str_title = "Data before occurence filter", f_logged = TRUE, ... )
  }

  # Return kept data frame
  data( obj_analysis, islogged = is.logged( obj_analysis ) ) = data( obj_analysis )[ vctr_f_keep, ]
  if( !is.na( str_output_dir ) )
  {
    if( i_length_keep > 0 )
    {
      # Heatmap after
      func_heatmap( obj_analysis, str_linkage = "average", plt_colors = plt_colors, str_title = "Data after occurence filter", f_logged = TRUE, ... )
    }
    dev.off() 
  }

  return( obj_analysis )
}


func_filter_by_percentile_percentage <- function(
### Feature filter. Filter rows that do not have a percentage of values in a top percentile of each group
obj_analysis,
### Holds the data associated with the current analysis
d_min_percentile = 0.5,
### Minimum percentile of a sample a feature must be in
d_min_occurence = NULL,
### Minimum times the feature must be in the given percentile to be kept
str_output_dir = NA,
### Directory to save plots
...
){
  loginfo( "Filter by percentile" )

  # Set min occurence if not given
  if( is.null( d_min_occurence ) )
  {
    d_min_occurence = ncol( data( obj_analysis ) ) * 0.1
  }

  # Get the min value for needed percentile for each sample
  vctr_min_values_at_percentile = apply( data( obj_analysis ), 2, function(x){ return( quantile( x, d_min_percentile ) ) } )

  # Find features that have the min occurence above the given percentile value.
  vctr_f_keep = apply( data( obj_analysis ), 1, func_min_occurence_at_min_value, vctr_min_values = vctr_min_values_at_percentile ) >= d_min_occurence

  # Optionally plot
  if( !is.na( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "value_percentile_filter.pdf" ) )
    plt_colors = func_monochromatic_palette()

    # Calculate values
    vctr_at_percentile = apply( data( obj_analysis ), 1, func_min_occurence_at_min_value, vctr_min_values = vctr_min_values_at_percentile )
    vctr_keep = vctr_at_percentile[ vctr_f_keep ]
    vctr_no_keep = vctr_at_percentile[ ! vctr_f_keep ]
    i_length_keep = length( vctr_keep )
    i_length_no_keep = length( vctr_no_keep )

    # Plot distributions
    if( ( i_length_keep == 0 ) && ( i_length_no_keep > 0 ) )
    {
       hist( vctr_at_percentile, col = STR_COLOR_BEFORE,
          main = "Filtering by percentile ( All Filtered ) ", xlab = "Occurence in top percentile across samples", ylab = "Density" )
    } else if( ( i_length_keep > 0 ) && ( i_length_no_keep == 0 ) )
    {
       hist( vctr_at_percentile, col = STR_COLOR_STANDARD,
          main = "Filtering by percentile ( None Filtered ) ", xlab = "Occurence in top percentile across samples", ylab = "Density" )
    } else {
      ret_hist1 = hist( vctr_keep, plot = FALSE )
      ret_hist2 = hist( vctr_no_keep, plot = FALSE )
      d_xlim = c( min( vctr_no_keep, vctr_keep ), max( vctr_keep, vctr_no_keep ) )
      d_ylim = c( 0, max( max( ret_hist1$counts ), max( ret_hist2$counts ) ) )
      vctr_d_breaks = hist( c( vctr_keep, vctr_no_keep ), plot = FALSE )$breaks
      hist( vctr_keep, col = STR_COLOR_STANDARD, breaks = vctr_d_breaks, xlim = d_xlim, ylim = d_ylim,
          main = "Filtering by percentile", xlab = "Occurence in top percentile across samples", ylab = "Density" )
      hist( vctr_no_keep, col = STR_COLOR_BEFORE, add = TRUE )
    }
    legend( "topright", title = "Filtering", legend = c("Filtered","Kept"), fill = c( STR_COLOR_BEFORE, STR_COLOR_STANDARD ) )

    # Heatmap before
    func_heatmap( obj_analysis, str_linkage = "average", plt_colors = plt_colors, str_title = "Data before percentile filter", f_log = TRUE, ... )
  }

  # Return kept data frame
  data( obj_analysis, islogged = is.logged( obj_analysis ) ) = data( obj_analysis )[ vctr_f_keep, ]

  if( !is.na( str_output_dir ) )
  {
    if( i_length_keep > 0 )
    {
      # Heatmap after
      func_heatmap( obj_analysis, str_linkage = "average", plt_colors = plt_colors, str_title = "Data after percentile filter", f_log = TRUE, ... )
    }
    dev.off()
  }

  return( obj_analysis )
}


func_min_occurence_at_min_value <- function(
### Returns if a vector has a minimum occurence of values at or equal to the given min value
### This is a helper function for func_filter_by_percentile_percentage
vctr_values_to_evaluate,
### Values to check if the are greater or equal to the min value
vctr_min_values,
### Minimum values (positionally) a value must be to be counted
d_min_occurence,
### Minimum counts a feature must have above the value
...
){
  return( length( which( vctr_values_to_evaluate >= vctr_min_values ) ) )
}


func_filter_by_sparsity <- function(
### Sample Filter. Filter by percent zeros ( sparsity )
obj_analysis,
### Holds the data associated with the current analysis
d_percentile = 0.9,
### Samples must have this max percent of zeros
str_output_dir = NA,
### Directory to save plots and such
...
){
  loginfo( "Filter by Sparsity" )

  i_nrow = nrow( data( obj_analysis ) )
  i_ncol = ncol( data( obj_analysis ) )

  # Filter rows
  vctr_percent_zero_rows = apply( data( obj_analysis ), 1, function(x) length( which(x == 0 ) ) ) / i_ncol
  vctr_row_keep = which( vctr_percent_zero_rows < d_percentile )
  vctr_row_col = rep( STR_COLOR_BEFORE, i_nrow )
  vctr_row_col[ vctr_row_keep ] = STR_COLOR_STANDARD

  # Reduce data frame rows
  data( obj_analysis, islogged = is.logged( obj_analysis ) ) = data( obj_analysis )[ vctr_row_keep, ]

  # Filter columns
  vctr_percent_zero_cols = apply( data( obj_analysis ), 2, function(x) length( which( x==0 ) ) ) / i_nrow
  vctr_col_keep = which( vctr_percent_zero_cols < d_percentile )
  vctr_col_col = rep( STR_COLOR_BEFORE, i_ncol )
  vctr_col_col[ vctr_col_keep ] = STR_COLOR_STANDARD

  # Filter out samples
  subset.samples( obj_analysis, vctr_col_keep )

  # Plot
  if( !is.na( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "sparsity_filter.pdf" ), useDingbats=FALSE )
    # Plot filtered second so that it is not hidden
    vctr_cur = vctr_percent_zero_rows[vctr_row_keep]
    plot( rnorm( length( vctr_cur ), 1, .1 ), vctr_cur, axes = FALSE, col = STR_COLOR_STANDARD, ylab = "", xlab = "", xlim = c(0,3), pch = 16, cex = .85, ylim = c( min( vctr_percent_zero_rows ) * 0.9, max( vctr_percent_zero_rows ) * 1.1  ) )
    vctr_cur = vctr_percent_zero_rows[ -1 * vctr_row_keep ]
    points( rnorm( length( vctr_cur ), 1, .1 ), vctr_cur, col = STR_COLOR_BEFORE, ylab = "", xlab = "", pch = 16, cex = .85 )

    vctr_cur = vctr_percent_zero_cols[ vctr_col_keep ]
    points( rnorm( length( vctr_cur ), 2, .1 ), vctr_cur, col = vctr_col_col, ylab = "", xlab = "", pch = 16, cex = .85 )
    vctr_cur = vctr_percent_zero_cols[ -1 * vctr_col_keep ]
    points( rnorm( length( vctr_cur ), 2, .1 ), vctr_cur, col = STR_COLOR_BEFORE, ylab = "", xlab = "", pch = 16, cex = .85 )
    boxplot( vctr_percent_zero_rows, vctr_percent_zero_cols, names = c("Row","Columns"), main = "Percent Zeros (columns after row filtering)", ylab = "Percent Zero", add = TRUE )
    legend( "topright", title = "Filtering", legend = c("Filtered","Kept"), col = c( STR_COLOR_BEFORE, STR_COLOR_STANDARD ), pch = 16 ) 
    dev.off()
  }
  return( obj_analysis )
}


###########################################
# Normalization methods
###########################################


func_cpm <- function(
### Normalize counts to counts per million
### Normalizes within columns
df_data,
### Count data to be transformed
...
){
  return( tss( df_data ) * 1e6 )
}


func_tss <- function(
### Total Sum Scaled
### Normalize columns of the data set buy dividing each observation by the total of the sample
df_data,
### Data to be transformed
...
){
  return( sweep( df_data, 2, colSums( df_data), "/" ) )
}


func_upper_quartile <- function(
### Normalize values using the upper quartile value
### http://vinaykmittal.blogspot.com/2013/10/fpkmrpkm-normalization-caveat-and-upper.html
df_data,
### Count data to be transformed
f_scale = FALSE,
...
){
  loginfo( "Upper quartile normalization" )

  df_normalized = df_data

  # Remove transcripts / genes/ rows that have 0 expression in all samples
  vctr_zero_rows = which( apply( df_data, 1, sum ) == 0 )
  if( length( vctr_zero_rows ) )
  {
    df_normalized = df_normalized[ -1 * vctr_zero_rows, ]
  }

  # For each column / sample find the 75th percentile value (upper quartile)
  vctr_upper_quartile = apply( df_normalized, 2, quantile, probs = c(0.75) ) 

  # Divide all the expression values by the upper quartile
  df_normalized = sweep( df_normalized, MARGIN = 2, vctr_upper_quartile, "/" ) 

  # Optionally scale all the values by the mean of the 75th percentile quartiles
  # To increase very low values created by the normalization
  if( f_scale )
  {
    return( df_normalized * mean( vctr_upper_quartile ) )
  }

  return( df_normalized )
}


############################################
# Quality Control
############################################

func_QC <- function(
### Perform quality control on matrix of data
obj_analysis,
### Holds the data associated with the current analysis
func_handle_outliers = func_do_nothing,
### Function to handle outliers in metadata
func_imputation = func_do_nothing,
### Function for imputing data
func_batch_correct = func_do_nothing,
### The method for batch correction
str_output_dir = NULL
###
){
  loginfo( "func_QC" )

  df_cleaned_data = data( obj_analysis )
  df_cleaned_metadata = metadata( obj_analysis )

  # Handle outliers
  # TODO
  
  # Handle imputation
  if( !is.null( metadata( obj_analysis ) ) && ( ! is.null( func_imputation ) ) ) 
  {
    list_ret = func_imputation( data( obj_analysis ), metadata( obj_analysis ) )
    df_cleaned_data = list_ret$data
    df_cleaned_metadata = list_ret$metadata
  }

  # Handle batch correction
  # TODO
  data( obj_analysis, islogged = is.logged( obj_analysis ) ) = df_cleaned_data
  metadata( obj_analysis ) = df_cleaned_metadata

  return( obj_analysis )
}


#--- Metadata imputation


func_complete_case <- function(
### Handle missing data by reducing the metadata and data observations
### to those that do not have missing data
df_data,
### Data to reduce based on metadata (assumed to be the same observations as metadata)
df_metadata
### Data to reduce to observation with no missing values
){
  loginfo( "Handle missing values by reducing metadata and data to observations without missing values" )
  df_complete_metadata = na.omit( df_metadata )
  df_complete_data = df_data[ ,row.names( df_complete_metadata ) ]
  return( list( data = df_complete_data, metadata = df_complete_metadata ) )
}


func_impute_mean <- function(
### Impute missing values with a mean for numerics and integers
### binary and factor data make a sperate NA category
df_data,
### Data measurements (will also be imputed)
df_metadata
### Metadata measurements to impute missing data
){
  loginfo( "Impute with mean" )

  # Reset the metadata
  df_metadata = func_impute_with_center_helper( df_metadata, mean )
  # Reset the data
  df_data = func_impute_with_center_helper( df_data, mean )
  # Return
  return( list( data=df_data, metadata=df_metadata ) )
}


func_impute_median <- function(
### Impute missing values with a median for numerics and integers
### binary and factor data make a sperate NA category
df_data,
### Data measurements (will also be imputed)
df_metadata
### Metadata measurements to impute missing data
){
  loginfo( "Impute with median" )

  # Reset the metadata
  df_metadata = func_impute_with_center_helper( df_metadata, median )
  # Reset the data
  df_data = func_impute_with_center_helper( df_data, median )
  # Return
  return( list( data=df_data, metadata=df_metadata ) )
}


func_impute_with_center_helper <- function(
df_matrix,
### Values to impute
func_center
### Function to use to make the value to imput with
){
  # Reset the data
  for( str_feature in names( df_matrix ) )
  {
    vctr_values = df_matrix[[str_feature]]
    vctr_omit = union( which( is.na( vctr_values ) ),which( is.nan( vctr_values ) ) )
    if( is.numeric( vctr_values ) )
    {
      df_matrix[ vctr_omit, str_feature ] = func_center( vctr_values, na.rm = TRUE )
    } else if ( is.integer( vctr_values ) )
    {
      df_matrix[ vctr_omit, str_feature ] = round( func_center( vctr_values, na.rm = TRUE ) )
    } else if ( is.factor( vctr_values ) || is.ordinal( vctr_values ) )
    {
      levels( vctr_values ) = c( levels(vctr_values), "NA" )
      vctr_values[ vctr_omit ] = "NA"
      df_matrix[, str_feature] = vctr_values
    }
  }
  return( df_matrix )
}


# TODO
# func_impute_multiple <- function(
#


#---- Batch affect correction


func_PCA_skimming <- function(
### Removes batch affect by removing PCs correclated with the metadata
df_frame,
### Data to correct
df_metadata
### Metadata to control for
){
  loginfo( "func_PCA_skimming" )
  df_cleaned = df_frame
  return( df_cleaned )
}


func_combat <- function(
### Removes batch affect by using the combat method
df_frame,
### Data to correct
df_metadata
### Metadata to control for
){
  loginfo( "func_combat" )

  df_cleaned = df_frame
  return( df_cleaned )
}


func_control_for_batching <- function(
### Removes batch affects by regession
df_frame,
### Data to correct
df_metadata
### Metadata to control for
){
  loginfo( "func_control_for_batching" )

  df_cleaned = df_frame
  return( df_cleaned )
}


func_edgeR_batch_control <- function(
### Removes batch affects with edgr
df_frame,
### Data to correct
df_metadata
### Metadata to control for
){
  loginfo( "func_edgeR_batch_control" )
  library( edgeR )

  print( removeBatchEffect( x = df_frame, covariates = df_metadata ) )
}

func_ruvseq_batch_control <- function(
### Removes batch affect with ruvseq as per their tutorial
### http://www.bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf
df_frame,
### Counts
df_metadata,
### Metadata to control for
...
){
  loginfo( "func_ruvseq_batch_control" )
  library( RUVseq )

  # Filtering out non-expressed genes should have already occurred.
  # The tutorial requires 5 reads in atleat 2 samples for each gene.
  # Needs to be count data
  # 4 methods for defining the control genes depending if you have:
  # 1. Spike-ins or known housekeeping genes unaffected by the batch affects
  # 2. Invariant genes in the data set ("in-silico empirical" genes)
  # 3. Estimating the factors of unwanted variation using replicate/negative control samples
  # 4. Computing residuals from a GLM file, without RUVSeq normalization but potentially using standard normalizations

  # Get gene and spike names
  ls_genes = rownames()
  ls_spikes = rownames()

  # Make SeqExpression object
  seq_ex_cur = 

  # 
  
  
}


####################################################
# Generic selection of groups of samples or features
####################################################

func_do_dendrogram_selection <- function(
### Selection of groups given a dendrogram
mtrx_distance,
### Distance matrix to select from
str_key
### Key to indicate methodology to use
){
  return( switch( str_key,
    gmd = func_do_dendrogram_selection_gmd( mtrx_distance ),
    wgcna = func_do_dendrogram_selection_wgcna( mtrx_distance )))
}


func_do_dendrogram_selection_gmd <- function(
### Select selection of groups given a dendrogram
### This uses the generalized minimum distance
mtrx_distance
### Distance matrix to select from
){
  # Get HCL
  rtrn_hclust = hclust( mtrx_distance )

  # Evaluate clusters
  rtrn_elbow = elbow.batch( css.hclust( mtrx_distance, rtrn_hclust ) )

  # Select groups
  vctr_factor_groups = as.factor( cutree( rtrn_hclust, k=rtrn_elbow$k ) )
  names( vctr_factor_groups ) = names( mtrx_distance )

  # return named factor vector
  return( vctr_factor_groups )
}


func_do_dendrogram_selection_wgcna <- function(
### Select selection of groups given a dendrogram
### using Weighted gene correlation matrices
mtrx_distance
### Distance matrix to select from
){
 pass
 #TODO
}


###########################################
# Sample Selection
###########################################

func_do_dendrogram_sample_selection <- function(
### Select sample groups using dendrogram cutting
obj_analysis,
### Holds the data associated with the current analysis
str_output_dir,
### Directory to output data and visualizations
...
){
  loginfo( "Sample selection by automated dendrogram cutting" )

  # Create distance matrix
  mtrx_distance = sample.dist( obj_analysis )

  # Get groupings
  sample.groups( obj_analysis )[[ STR_SAMPLE_DENDROGRAM ]] = func_do_dendrogram_selection_gmd( mtrx_distance )

  # If plotting, plot the HCL with the groupings
  if( !is.null( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "sample_selection_dendrogram.pdf") , useDingbats=FALSE )
    func_heatmap( obj_analysis, vctr_grouping = sample.groups( obj_analysis )[[ STR_SAMPLE_DENDROGRAM ]] )
    dev.off()
  }
  
  # Return analysis object
  return( obj_analysis )
}


func_do_mclust <- function(
### Select sample groups using mixture gaussians
obj_analysis,
### Holds the data associated with the current analysis
str_output_dir,
### Directory to output data and visualizations
...
){
  loginfo( "Sample selection with mclust" )

  print( ordinations( obj_analysis ) )

  # Select groups from each ordination
  for( str_ordinations in names( ordinations( obj_analysis ) ) )
  {
    # If plotting open a pdf
    if( !is.null( str_output_dir ) )
    {
      pdf( file.path( str_output_dir, paste( "sample_selection_mclust_using_",str_ordinations,".pdf",sep="")), 
           useDingbats=FALSE )
    }

    # Call mclust on two dimensions
    rtrn_mclust = Mclust( ordinations( obj_analysis )[[ str_ordinations ]][c(1,2)] )

    # Get groups
    sample.groups( obj_analysis )[[ paste( str_ordinations, "mclust", sep="_" ) ]] = rtrn_mclust$classification

    # Optionally plot
    if( !is.null( str_output_dir ) )
    {
      plot( rtrn_mclust )
    }

    # Close pdf if opened
    if( !is.null( str_output_dir ) )
    {
      dev.off()
    }
  }

  # Return analysis object
  return( obj_analysis )
}


###########################################
# Discriminate features
###########################################


func_discriminate_by_anova <- function(
### Select feature groups using ANOVA
obj_analysis,
### Holds the data associated with the current analysis
str_output_dir
### Output directory
){
  loginfo( "Feature selection by ANOVA" )

  for( vctr_grouping in sample.groups( obj_analysis ) )
  {
    print("****" )
    print( vctr_grouping )
    print( "***" )
    i_group_count = length( unique( vctr_grouping ) )
    list_results = list()
    vctr_pvalues = c()
    vctr_names = c()
    vctr_contrast = c()
    vctr_grouping = as.factor( vctr_grouping )

    if( i_group_count < 2 )
    {
      print( paste( "Can not perform anova or t-test, number of groups found are ", i_group_count, "." ) )
    } else {

      ### Perform ANOVA
      for( row_name in row.names( df_frame ) )
      {
        results_anova = aov( as.vector( unlist(df_frame[ row_name, ])) ~ vctr_grouping )
        vctr_pvalues = c( vctr_pvalues, summary(results_anova)[[1]][["Pr(>F)"]][1] )
        vctr_names = c( vctr_names, row_name )
        result_pt = TukeyHSD( results_anova )
        vctr_index_sig = which( result_pt$'as.factor(vctr_grouping)'[,4] <= 0.05 )
        if( length( vctr_index_sig ) == 1 )
        {
          vctr_contrast = c( vctr_contrast, result_pt$'as.factor(vctr_grouping)'[,4][ vctr_index_sig ] )
        } else {
          vctr_contrast = c( vctr_contrast, NA )
        }
      }

      ### Correct for multiple testing
      ### Filter for discrimination in one group
      vctr_qvalues = p.adjust( vctr_pvalues, "BH" )
      i_index_keep_pvalue = which( vctr_qvalues <= i_q_value_cutoff )
      i_index_keep_contrast = which( !is.na( vctr_index_sig ) )
      i_index_keep = i_index_keep_pvalue && i_index_keep_constrast
      df_corrected = cbind( p_value = vctr_pvalues[ i_index_keep ], q_value = vctr_qvalues[ i_index_keep ], contrast = names( vctr_contrast[ i_index_keep ] ), constrast_adj_pvalue = vctr_contrast[ i_index_keep ] )
      row.names( df_corrected ) = vctr_names[ i_index_keep ]

      ### Add results to discriminant information
      discriminant( list_results )[[ paste( names(vctr_grouping),STR_ANOVA,sep="_")]] = df_corrected
    }
  }

  return( list_results )
}


#############################
# Feature groupings
#############################

func_select_features_by_dendrogram <- function(
### Select feature groups using dendrogram cutting 
obj_analysis,
### Holds the data associated with the current analysis
str_output_dir,
### Directory to output data and visualizations
...
){
  loginfo( "Feature selection by automated dendrogram cutting" )

  # Create distance matrix
  mtrx_distance = feature.dist( obj_analysis )

  # Get groupings
  feature.groups( obj_analysis )[[ STR_SAMPLE_DENDROGRAM ]] = func_do_dendrogram_selection_gmd( mtrx_distance )

  # If plotting, plot the HCL with the groupings
  if( !is.null( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "sample_selection_dendrogram.pdf") , useDingbats=FALSE )
# TODO    func_heatmap( obj_analysis, vctr_grouping = feature.groups( obj_analysis )[[ STR_SAMPLE_DENDROGRAM ]] )
    dev.off()
  }

  # Return analysis object
  return( obj_analysis )
}


func_select_features_by_pca <- function(
### Select by PCA
df_frame,
### Data to perform selection from
i_PC_Components = NULL,
### Number of PC components
i_count = NULL,
### The number of top ranked features to select
str_output_dir = NULL
### Output directory to document to
){
  loginfo( "Feature selection by PCA" )
  # Row center and log
  mtrx_scale = t( scale( t( as.matrix( func_log( df_frame ) ) ), center=TRUE, scale=TRUE ) )
  # Remove constant rows
  mtrx_scale = mtrx_scale[ !is.na( mtrx_scale[, colnames( mtrx_scale )[ 1 ] ] ), ]
  # Perfrom PCA
  results_pca = princomp( mtrx_scale, cor = TRUE )
  # For each component get the top genes
  mtrx_scores = results_pca$scores
  # How many genes exist
  i_gene_counts = nrow( mtrx_scores )
  # Get percent variance
  vctr_i_variance = results_pca$sdev^2
  vctr_i_percent_variance = ( vctr_i_variance ) / sum( vctr_i_variance )
  i_eigen_count = length( results_pca$sdev )

  # If check output directory path
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }

  # If no count is given us the upper quartile expressed
  if( is.null( i_count ) )
  {
    i_count = floor( i_gene_counts * 0.25 )
  }

  # Get the stick break distribution
  if( is.null( i_PC_Components ) )
  {
    result_stick = func_select_pca_components_stick_breaking( vctr_i_percent_variance = vctr_i_percent_variance, str_output_dir = str_output_dir )
    i_PC_Components = result_stick$count
  }

  # Holds the selection per component
  vctr_return_feature_list = c()

  # For each component select the top genes
  for( i_index_component in 1:i_PC_Components )
  {
    vctr_genes = row.names( mtrx_scores[ order( abs( mtrx_scores[,i_index_component] ) ), ] )[ 1 : i_count ]
    if( ! is.null( str_output_dir ) )
    {
      print( str_output_dir )
      func_write_frame( vctr_genes, file.path( str_output_dir, paste( "Extreme",i_count, "genes","by","PCA_component",i_index_component,".txt", sep="_" ) ) )
    }
    vctr_return_feature_list = unique( c( vctr_return_feature_list, vctr_genes ) )
  }
  return( vctr_return_feature_list )
  ### Returns a list of selected genes
}


###########################################
# Functional analysis
###########################################

func_david <- function(
### Perform functional analysis by interfacting with the david website
vctr_genes,
### Genes to evaluate
str_output_dir
### Output directory
){
  library("RDAVIDWebService")

  david_service = DAVIDWebService$new(email=str_email)
  david_list = addList( david_service, vctr_genes, idType=str_annotation_type, listName= paste( format(Sys.time(), "%b_%d_%Y_%X"), "automated_analysis",sep="_"), listType="Gene")
  setAnnotationCategories( david_service, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
  termCluster = getClusterReport( david_service, type="Term" )
  getClusterReportFile( david_service, type="Term", fileName=path.join( str_output_dir, "termClusterReport1.tab" ) )
  
}


func_gsea_preranked <- function(
### Performs GSEA with preranked gene lists

){
# From the limma package
wilcoxGST()
}


func_gsea_response <- function(
### Performs GSEA with response

){
GSEABase
}


func_stringdb <- function(
### Interface with String DB to make gene networks
vctr_genes,
### Genes to make a network with
str_output_dir
### Output directory
){

}


###########################################
# Generic visualizations
###########################################


#func_matrix_plot <- function(
#### This function plots raw dataframe
#df_frame,
#### Data frame oto plot
#title = ""
#### Title to figure
#){
#  heatmap.3b( df_frame, dendrogram = "both", scale = "center", col = func_monochromatic_palette( ), 
#              trace = "none", tracecol = C_STR_DENSITY_COLOR, density.info = "histogram", main = title )
#}


###########################################
# Module Infrastructure
###########################################

func_do_nothing <- function(
### This function does nothing
df_frame,
### Data frame to act on
df_metadata = null
### Metadata data frame
){
  return( df_frame )
  ### Returns the identical data frame
}


func_get_batch_module <- function(
### This function take a keyword and returns the function
### that would perform batch analysis associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    combat = func_combat,
    edgeR = func_edgeR_batch_control,
    none, func_do_nothing,
    pca_skimming = func_PCA_skimming,
    regression = func_control_by_regression )
  # Return function
  return( func_ret )
  ### Returns batch function to perform the needed analysis 
}


func_get_discriminant_module <- function(
### This function take a keyword and returns the function
### that would perform that discriminate feature analysis associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    anova = func_discriminate_by_anova,
    none = func_do_nothing )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed discriminant feature analysis 
}


func_get_filter_module <- function(
### This function take a keyword and returns the function
### that would perform the filtering associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    downsample = func_filter_downsample,
    top = func_filter_by_first_rows,
    none = func_do_nothing,
    occurence = func_filter_by_min_value_occurence,
    pca = func_filter_extreme_by_pca,
    percentile = func_filter_by_percentile_percentage,
    sd = func_filter_by_SD,
    sparsity = func_filter_by_sparsity )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed filtering 
}


func_get_impute_module <- function(
### This function take a keyword and returns the function
### that would perform the imputation associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    complete_case = func_complete_case,
    mean = func_impute_mean,
    median = func_impute_median,
    multiple_imputation = func_impute_multiple,
    none = func_do_nothing )
  # Return function
  return( func_ret )
  ### Returns function to perform the imputation needed
}


func_get_normalization_module <- function(
### This function taks a keyword and returns the function
### that woudl perform the imputation associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    cpm = func_cpm,
    none = func_do_nothing,
    upper = func_upper_quartile )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed transformation
}


func_get_ordinate_module <- function(
### This function take a keyword and returns the function
### that would perform the ordination associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    none = func_do_nothing,
    pca = func_PCA,
    nmds = func_PCoA )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed ordination
}


func_get_select_sample_groups_module <- function(
### This function take a keyword and returns the function
### that would perform the sample group selection associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    dendrogram = func_do_dendrogram_sample_selection,
    mclust = func_do_mclust,
    none = func_do_nothing )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed sample group selection 
}


func_get_select_feature_groups_module <- function(
### This function take a keyword and returns the function
### that would perform the feature group associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    none = func_do_nothing,
    dendrogram = func_select_features_by_dendrogram )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed feature group selection 
}


func_get_transform_module <- function(
### This function take a keyword and returns the function
### that would perform the data transform associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    asin_sqrt = func_asin_sqrt,
    center = func_center,
    log = func_log,
    none = func_do_nothing,
    tss = func_tss,
    zscale = func_zscale )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed data transform
}
