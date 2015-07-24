#############################
# Color
#############################

# Constants
C_STR_DENSITY_COLOR = "cyan"
### Color for lengend density lines
STR_COLOR_AFTER = "#00660075" # green
### Used for after a method (used as a secondary to STR_COLOR_BEFORE)
STR_COLOR_BEFORE = "#DD640075" # orange
### USed for the before color in a before / after picture.
### Also used as the selection color for plots.
STR_COLOR_STANDARD = "#00000075" # transparent black ( grey )
### Color to be used for normal plotting

# Colors for quantiles
STR_MIN_HEX = "#ff0000"
STR_Q1_HEX = "#EE82EE"
STR_MEDIAN_HEX = "#00ffff"
STR_Q3_HEX = "#0000ff"
STR_MAX_HEX = "#00ff00"


func_factor_to_metadata_color <- function(
vctr_factors
### Factors to turn to metadata colors
){
  vctr_grouping_colors = rep("NA", length( vctr_factors ) )
  names( vctr_grouping_colors ) = as.character( vctr_factors )
  vctr_levels = levels( vctr_factors )
  vctr_colors = func_metadata_palette( length( vctr_levels ) )
  for( i_color_index in 1:length( vctr_levels ) )
  {
    vctr_grouping_colors[ which( vctr_factors == vctr_levels[ i_color_index ] ) ] = vctr_colors[ i_color_index ]
  }
  return( list( vctr_grouping_colors = vctr_grouping_colors, vctr_levels = vctr_levels, vctr_colors = vctr_colors ) )
}


################################
# Feature Plots
################################

func_feature_distributions_by_sparsity <- function(
### Show distributions of data given different levels of sparsity
obj_analysis
){
  loginfo( "Visualize distributions by sparsity." )

  # Make sure these are logged or logg the data matrix +1 before density
  df_frame = data( obj_analysis )
  if( ! is.logged( obj_analysis ) )
  {
    df_frame = get.logged.matrix( obj_analysis )
  }

  # The number of features per quantile to plot
  i_number_features = min( 10, round( nrow( df_frame )/5 ) )

  # Get the sorted order for the genes
  feature.sum.order <- order( apply( df_frame, 1, sum) )
  feature.summary <- summary( 1:length( feature.sum.order ))

  # Get the 10 most sparse genes
  index_min <- feature.sum.order[ 1:i_number_features ]

  # 1 quartile sparsity genes
  index_q1 <- feature.sum.order[ floor(feature.summary[[2]]):(floor(feature.summary[[2]])+i_number_features) ]
  # Median sparsity genes
  index_median <- feature.sum.order[ floor(feature.summary[[3]]):(floor(feature.summary[[3]])+i_number_features) ]
  # 3rd quartile sparsity genes
  index_q3 <- feature.sum.order[ floor(feature.summary[[5]]):(floor(feature.summary[[5]])+i_number_features) ]
  # Get the least sparse genes
  index_max <- feature.sum.order[ (length(feature.sum.order)-i_number_features):length(feature.sum.order) ]

  # Get double the max density for the max group
  i_height = 0
  for( i_max_plot in index_max ){
    i_height = max( i_height, density(as.matrix( df_frame[ i_max_plot,]))$y )
  }
  i_height = i_height * 2

  # Plot
  plot( x=0,y=0,type="p", xlim=c(0,max(df_frame)), ylim=c(0,i_height), main="Feature Count Distributions by Sparsity", ylab="Density of feature counts", xlab="Feature count value (Log)" )
  for( i_q3_plot in index_q3 ){ lines( density( as.matrix( df_frame[ i_q3_plot, ] )), col = paste(STR_Q3_HEX,"75",sep=""))}
  for( i_median_plot in index_median ){ lines( density( as.matrix( df_frame[ i_median_plot, ])), col = paste(STR_MEDIAN_HEX,"75",sep=""))}
  for( i_min_plot in index_min ){ lines( density( as.matrix( df_frame[ i_min_plot, ])), col = paste(STR_MIN_HEX,"75",sep=""))}
  for( i_q1_plot in index_q1 ){ lines( density( as.matrix( df_frame[ i_q1_plot, ])), col = paste(STR_Q1_HEX,"75",sep=""))}
  for( i_max_plot in index_max ){ lines( density( as.matrix( df_frame[ i_max_plot, ])), col = paste(STR_MAX_HEX,"75",sep=""))}
  legend( "topright", c("Min","Q1","Median","Q3","MAX"), fill=c( STR_MIN_HEX, STR_Q1_HEX, STR_MEDIAN_HEX, STR_Q3_HEX, STR_MAX_HEX ), title="Sparsity group" )
}


func_feature_sd_vs_mean <- function(
###
obj_analysis,
vctr_colors = NULL
){
  loginfo( "Visualize data with Mean vs SD." )

  # Make sure these are logged or logg the data matrix +1 before density
  df_frame = data( obj_analysis )

  if( ! is.logged( obj_analysis ) )
  {
    df_frame = get.logged.matrix( obj_analysis )
  }

  # Make colors default unless given
  if( is.null( vctr_colors ) )
  {
    vctr_colors = rep( STR_COLOR_STANDARD, nrow( df_frame ) )
  }

  vctr_sd = apply( df_frame, 1, sd )
  vctr_mean = apply( df_frame, 1, mean )
  plot( vctr_mean, vctr_sd, main="Feature SD vs Mean", xlab="Mean (log)", ylab="SD (log)", col = vctr_colors )

  # TODO add coeffient of variation line
}


func_heatmap <- function(
### Create a heatmap of data frame
obj_analysis,
### The collected analysis run
str_linkage = "average",
### Any valid linkage method for hclust
vctr_grouping = NULL,
### Optional Groups the samples
str_output_dir = NULL,
### The directory to output figure and such
plt_colors = NULL,
### Color pallette, if not given ploychromatic default will be used
str_title = "Data Heatmap",
### Title of plot
f_log = FALSE,
### Indicates the data should be logged if TRUE
...
){

  ## TODO Tranform (center and such) should happen before this

  df_data = NULL
  if( f_log )
  {
    df_data = data( obj_analysis )
  } else {
    df_data = get.logged.matrix( obj_analysis )
  }

  loginfo( "Visualize data with heatmap " )
  print( dim( df_data ) )
  # Get a consistent color scheme for the plotting
  if( is.null( plt_colors ) )
  {
    plt_colors = func_polychromatic_palette( )
  }

  # Get distance matrix
  dist_row = make.feature.distance.matrix( obj_analysis, f_log = f_log )
  print( dim( dist_row ) )
  dist_col = make.sample.distance.matrix( obj_analysis, f_log = f_log )
  print( dim( dist_col ) )
  # Visualize heatmap and possibly record in pdf file
  if( ! is.null( str_output_dir ) )
  {
    pdf( paste( str_output_dir, "data_heatmap.pdf" ) )
  }

  if( !is.null( vctr_grouping ) )
  {
    # Set the groupings to colors
    lret_coloring = func_factor_to_metadata_color( vctr_grouping )
    vctr_grouping_colors = matrix( lret_coloring$vctr_grouping_colors, nrow = 1 )
    names( vctr_grouping_colors ) = names( lret_coloring$vctr_grouping_colors )
    heatmap.3b( x= df_data, dendrogram='both', 
              Rowv=as.dendrogram( hclust( dist_row ), method = str_linkage ), 
              Colv=as.dendrogram( hclust( dist_col ), method = str_linkage ), 
              scale='none', symm=TRUE, key=TRUE, col = plt_colors,
              density.info='histogram', denscol = C_STR_DENSITY_COLOR, trace='none', 
              symkey=FALSE, margins=c(10,10), cexCol=1, ColSideColors = vctr_grouping_colors,
              cexRow=1, cex.main=0.75, main = str_title, ...)
  } else {
    heatmap.3b( x=df_data, dendrogram='both',
              Rowv=as.dendrogram( hclust( dist_row ), method = str_linkage ),
              Colv=as.dendrogram( hclust( dist_col ), method = str_linkage ),
              scale='none', symm=TRUE, key=TRUE, col = plt_colors,
              density.info='histogram', denscol = C_STR_DENSITY_COLOR, trace='none',
              symkey=FALSE, margins=c(10,10), cexCol=1,
              cexRow=1, cex.main=0.75, main = str_title, ...)
  }

  if( ! is.null( str_output_dir ) )
  {
    dev.off()
  }
}


func_cor_matrix <- function(
### Show the correlation of the date in the data frame (samples or features)
obj_analysis,
### The analysis object containing the state of the analysis run
vctr_grouping = NULL,
### Grouping of samples or features by color
str_linkage = "average",
### Any valid linkage method for hclust
str_output_dir = NULL,
### The directory to output figure and such
f_correlate_samples = TRUE,
### Indicates if samples (TRUE; columns) or features (FALSE; rows) are correlated
str_title = NULL,
### Title for plot
f_log = FALSE,
### Log data before calculating
...
){
  loginfo( "Visualize correlation" )
  if( f_correlate_samples )
  {
    loginfo( "Sample Correlation" )
  } else {
    loginfo( "Feature Correlation" )
  }

  # Get a consistent color scheme for the plotting
  plt_colors = func_polychromatic_palette( )

  # Create correlation matrix
  mtrx_cor = NA
  if( f_correlate_samples )
  {
    mtrx_cor = make.sample.distance.matrix( obj_analysis, f_log = f_log )
  } else {
    mtrx_cor = make.feature.distance.matrix( obj_analysis, f_log = f_log )
  }

  # Remove diag
  mtrx_cor = as.matrix( mtrx_cor )
  diag( mtrx_cor ) = NA

  # Distance between items
  dgrm_samples = as.dendrogram( hclust( as.dist( 1-mtrx_cor ), method = str_linkage ) )

  # Visualize heatmap and possibly record in pdf file
  if( ! is.null( str_output_dir ) )
  {
    str_file_name = "sample_correlation_matrix.pdf"
    if( ! f_correlate_samples )
    {
      str_file_name = "feature_corrleation_matrix.pdf"
    }
    pdf( paste( str_output_dir, str_file_name ) )
  }

  if( is.null( str_title ) )
  {
    str_title = "Sample Correlation Matrix"
    if( ! f_correlate_samples )
    {
      str_title = "Feature Correlation Matrix"
    }
  }

  if( !is.null( vctr_grouping ) )
  {
      lret_grouping_colors = func_factor_to_metadata_color( vctr_grouping )
      vctr_grouping_colors = matrix( lret_grouping_colors$vctr_grouping_colors, nrow = 1 )
      names( vctr_grouping_colors ) = names( lret_grouping_colors$vctr_grouping_colors )
      heatmap.3b( mtrx_cor, dendrogram = "both", Rowv = dgrm_samples, Colv = dgrm_samples, col = plt_colors, scale = "none", revC = TRUE, key = TRUE, density.info = "histogram", trace = "none", margin.for.labRow = 10, margin.for.labCol = 10, cexCol = 1, cexRow = 1, cex.main = 0.75, main = str_title, sep.color = "black", sep.lwd = 1, denscol = C_STR_DENSITY_COLOR, ColSideColors = vctr_grouping_colors, RowSideColors = t( vctr_grouping_colors ), ... ) 
  } else {
      heatmap.3b( mtrx_cor, dendrogram = "both", Rowv = dgrm_samples, Colv = dgrm_samples, col = plt_colors, scale = "none", revC = TRUE, key = TRUE, density.info = "histogram", trace = "none", margin.for.labRow = 10, margin.for.labCol = 10, cexCol = 1, cexRow = 1, cex.main = 0.75, main = str_title, sep.color = "black", sep.lwd = 1, denscol = C_STR_DENSITY_COLOR, ... ) 
  }

  if( ! is.null( str_output_dir ) )
  {
    dev.off()
  }
}


func_weighted_bar_plot <- function(
  ### Plots a bar plot but the color is weighted based on a given weighting
  vctr_values,
  ### Values to plot
  vctr_weights,
  ### Weighted values, should be between 0 and 1.
  ### Will be normalized so that the largest value is 1
  vctr_colors,
  ### Hex colors to guide the bar plot coloring
  ...
  ){
  vctr_weights = round( vctr_weights / max( vctr_weights ) * 99 )
  # Normalize weights
  
  vctr_updated_colors = c()
  for( i_index in 1:length(vctr_colors) )
  {
    vctr_updated_colors = c( vctr_updated_colors, func_change_hex_color_transparency( vctr_colors[ i_index ], vctr_weights[ i_index ]  ) )
  }
  # Update the colors
  
  barplot( vctr_values, col = vctr_updated_colors, ... )
}


func_change_hex_color_transparency <- function(
  ### Change the transparency of a hex color
  str_hex_color,
  ### The hex color to change
  int_transparency
  ### Transparency
){
  return( paste( substr( str_hex_color, 0, 7 ), sprintf( "%02d", int_transparency ), sep = "" ) )
}


################################
# Sample plots
################################

func_PCA = function(
### Perform ordination with PCA
obj_analysis,
### Holds the data for the analysis run
str_output_dir = NULL,
### Output directory
vctr_str_color_groups = NULL,
...
){
  loginfo( "Ordination by PCA" )

  # Row center and log
  mtrx_scale = t( scale( t( get.logged.matrix( obj_analysis ) ), center=TRUE, scale=TRUE ) )
  # Remove constant rows
  mtrx_scale = mtrx_scale[ !is.na( mtrx_scale[, colnames( mtrx_scale )[ 1 ] ] ), ]
  # Perfrom PCA
  results_pca = prcomp( mtrx_scale, retx = TRUE )
  # Get percent variance
  vctr_i_variance = results_pca$sdev^2
  vctr_i_percent_variance = ( vctr_i_variance ) / sum( vctr_i_variance )
  i_eigen_count = length( results_pca$sdev )
  # Get the stick break distribution
  result_stick = func_select_pca_components_stick_breaking(vctr_i_percent_variance = vctr_i_percent_variance, str_output_dir = str_output_dir )
  vctr_break_distribution = result_stick$distribution
  i_max_pc = result_stick$count

  # Store the loadings of interest
  ordinations( obj_analysis )[[STR_PCA]] = as.data.frame( results_pca$rotation[,1:i_max_pc] )

  # Write loadings
  func_write_frame( results_pca$rotation, file.path( str_output_dir, "pca_ordination_loadings.tsv" ) )
  # write scores
  func_write_frame( results_pca$x, file.path( str_output_dir, "pca_ordination_scores.tsv" ) )

  # Plot in cur directory if output directory not indicated
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }

  # create colorings
  vctr_str_colors = rep( STR_COLOR_STANDARD, ncol( data( obj_analysis ) ) )

  # Plot scree vs stick break
  # Plot selected components
  pdf( file.path( str_output_dir, "pca.pdf" ), useDingbats=FALSE )
  for( i_index in 1:i_max_pc )
  {
    plot( results_pca$rotation[, i_index], results_pca$rotation[, i_index + 1 ], 
        xlab = paste("PC", i_index, "(", round( vctr_i_percent_variance[ i_index ]*100, digits = 2 ), "% variance )" ),
        ylab = paste( "PC", i_index + 1, "(", round(vctr_i_percent_variance[ i_index + 1 ]*100, digits = 2 ), "% variance )" ),
        col=vctr_str_colors, main = "Ordination by PCA", pch = 16 )
  }

  # Plot scree plot
  vctr_str_comp_colors = rep( "grey",i_eigen_count )
  vctr_str_comp_colors[ 1: i_max_pc ] = "red"
  plot( results_pca, pch = "16", main = "Scree plot", col = vctr_str_comp_colors, xlab = "Component")

  # Paint with metadata
  for( str_name in names( metadata( obj_analysis ) ) )
  {
    vctr_metadata = metadata( obj_analysis )[[ str_name ]]

    vctr_meta_colors = NULL
    vctr_meta_levels = NULL
    vctr_str_colors = NULL

    if( is.numeric( vctr_metadata ) )
    {
      vctr_str_colors = func_polychromatic_palette( length( vctr_metadata ) )
      vctr_meta_levels = c( paste( max(vctr_metadata), "(Max)", sep = " "), paste( min(vctr_metadata), "(Min)", sep = " ") )
      vctr_meta_colors = c( vctr_str_colors[ length( vctr_str_colors ) ], vctr_str_colors[ 1 ] )
    } else {
      if( ! is.factor( vctr_metadata ) )
      {
        vctr_metadata = as.factor( vctr_metadata )
      }
      ret_color_plt = func_factor_to_metadata_color( metadata( obj_analysis )[[ str_name ]] )
      vctr_meta_colors = ret_color_plt$vctr_colors
      vctr_meta_levels = ret_color_plt$vctr_levels
      vctr_str_colors = ret_color_plt$vctr_grouping_colors
    }

    for( i_index in 1:i_max_pc )
    {
      plot( results_pca$rotation[, i_index], results_pca$rotation[, i_index + 1 ], 
        xlab = paste("PC", i_index, "(", round( vctr_i_percent_variance[ i_index ]*100, digits = 2 ), "% variance )" ),
        ylab = paste( "PC", i_index + 1, "(", round(vctr_i_percent_variance[ i_index + 1 ]*100, digits = 2 ), "% variance )" ),
        col=vctr_str_colors, main = paste( "Ordination by PCA (showing metadata ",str_name,")", sep ="" ), pch = 16)
      legend( "topright", vctr_meta_levels, fill = vctr_meta_colors, title = str_name )
    }
  }

  dev.off()
  return( obj_analysis )
}


func_PCoA <- function(
obj_analysis,
### Holds the data for the analysis run
i_dimensions = 2,
### The number of dimensions to estimate
mtrx_dist = NULL,
### The distance matrix to use (must be given or str_metric must be given)
str_output_dir = NULL,
### Output directory
...
){

  # Distance metric
  str_metric = sample.corr.metric( obj_analysis )

  ###
  # Note:
  # Much thanks to the following post which inspired much of this code. Well done!
  # http://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
  ##

  if( is.null( str_metric ) )
  {
    loginfo( paste( "NMDS (PCoA) Ordination using a custom distance/dissimilarity matrix." ) )
  } else {
    loginfo( paste( "NMDS (PCoA) Ordination using", str_metric, "distance/dissimilarity", sep = " " ) )
  }

  # Plot ordination
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }

  # Start pdf
  pdf( file.path( str_output_dir, paste( "pcoa_", str_metric,".pdf", sep = "" ) ), useDingbats=FALSE )

  # Check that we can make or have a distance matrix
  if( is.null( mtrx_dist ) )
  {
    if( is.null( str_metric ) )
    {
      logerror( "NMDS (PCoA) Ordination. Please give the either a distance/dissimilarity metric or matrix. Ordination was not performed." )
      return( NULL )
    }

    # Generate a distance matrix
    mtrx_dist = make.sample.distance.matrix( obj_analysis, f_log = is.logged( obj_analysis ) )
  }

  # Perform NMDS
  ret_NMDS = metaMDS( mtrx_dist, k=i_dimensions, autotransfor = FALSE, trymax = 100 )

  # Plot ordination but only samples
  # Plot sample ordination with sample names
  for( str_plot_type in c( "points", "text" ) )
  {
    ordiplot( ret_NMDS, choices = c(1,2), display = "sites", type = str_plot_type, xlab = "NMDS 1", ylab = "NMDS 2",
            main = paste( "NMDS Ordination (",str_metric,") Stress=", ret_NMDS$stress, sep = " " ))
  }

  # Plot stress
  stressplot( ret_NMDS )

  # Store the ordination
  ordinations( obj_analysis )[[STR_NMDS]] = as.data.frame( ret_NMDS$points )

  # If no metadata is given return at this point
  if( is.null( metadata( obj_analysis ) ) )
  {
    dev.off()
    return( obj_analysis )
  }

  # Plot sample ordination painted by metadata
  for( str_name in names( metadata( obj_analysis ) ) )
  {
    vctr_metadata = metadata( obj_analysis )[[ str_name ]]
    if( is.numeric( vctr_metadata ) )
    {
      # If continuous, plot with gradient
      ordisurf( ret_NMDS, vctr_metadata, choices = c(1,2), col = STR_COLOR_AFTER, xlab = "NMDS1", ylab = "NMDS2",
                main = paste( "NMDS Ordination(",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )
      plot( ret_NMDS$points[,1], ret_NMDS$points[,2], col = STR_COLOR_STANDARD, add = TRUE )
      legend( "topright", c( paste( max(vctr_metadata), "(Max)", sep = " " ), paste( min(vctr_metadata), "(Min)", sep = " " )),
              title = str_name ) 

      # Plot samples as names
      ordisurf( ret_NMDS, vctr_metadata, choices = c(1,2), col = STR_COLOR_AFTER, xlab = "NMDS1", ylab = "NMDS2",
                main = paste( "NMDS Ordination(",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )
      orditorp( ret_NMDS, choices = c(1,2), display = "sites", col = STR_COLOR_STANDARD )
      legend( "topright", c( paste( max(vctr_metadata), "(Max)", sep = " " ), paste( min(vctr_metadata), "(Min)", sep = " " )),
              title = str_name ) 

    } else {
      # If not numeric make factor and plot
      if( !is.factor( vctr_metadata ) )
      {
        vctr_metadata = as.factor( vctr_metadata )
      }

      # Make color
      ret_meta_color = func_factor_to_metadata_color( vctr_metadata )
      vctr_metadata_levels = ret_meta_color$vctr_levels
      vctr_plt_colors = ret_meta_color$vctr_colors
      vctr_colors = ret_meta_color$vctr_grouping_colors

      # If discrete data, plot with ellipse
      plot( ret_NMDS$points[,1], ret_NMDS$points[,2], col = vctr_colors, add = TRUE,
            main = paste( "NMDS Ordination (",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )
      legend( "topright", vctr_metadata_levels, fill = vctr_plt_colors, title = str_name )

      # If discrete data, plot with ellipse
      plot( ret_NMDS$points[,1], ret_NMDS$points[,2], col = vctr_colors, add = TRUE,
            main = paste( "NMDS Ordination (",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )
      legend( "topright", vctr_metadata_levels, fill = vctr_plt_colors, title = str_name )

      # Plot each group to change color for groups
      for( i_index in 1:length( vctr_colors ) )
      {
        str_cur_level = vctr_metadata_levels[ i_index ]
        str_cur_color = vctr_plt_colors[ i_index ]
        ordiellipse( ret_NMDS, vctr_metadata, choices = c(1,2), display = "sites", kind = "se", conf = 0.95 ,
                   draw = "line", w = NULL, col = str_cur_color, label = FALSE, show.groups = str_cur_level )
        ordispider( ret_NMDS, vctr_metadata, choices = c(1,2), display = "sites", col = str_cur_color, label = FALSE, show.groups = str_cur_level )
      }
      legend( "topright", vctr_metadata_levels, fill = vctr_plt_colors, title = str_name )

      # Plot samples as names
      ordiplot( ret_NMDS, choices = c(1,2), type = "none", xlab = "NMDS1", ylab = "NMDS2",
                main = paste( "NMDS Ordination (",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )

      # Plot each group to change color for groups
      for( i_index in 1:length( vctr_colors ) )
      {
        str_cur_level = vctr_metadata_levels[ i_index ]
        str_cur_color = vctr_plt_colors[ i_index ]
        ordiellipse( ret_NMDS, vctr_metadata, choices = c(1,2), display = "sites", kind = "se", conf = 0.95 ,
                   draw = "line", w = NULL, col = str_cur_color, label = FALSE, show.groups = str_cur_level )
        ordispider( ret_NMDS, vctr_metadata, choices = c(1,2), display = "sites", col = str_cur_color, label = FALSE, show.groups = str_cur_level )
      }
      orditorp( ret_NMDS, choices = c(1,2), display = "sites", col = vctr_colors )
      legend( "topright", vctr_metadata_levels, fill = vctr_plt_colors, title = str_name )
    }
  }

  dev.off()
  return( obj_analysis )
}
