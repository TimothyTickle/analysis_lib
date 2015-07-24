#!/usr/bin/env Rscript


# This is a generic analysis engine for performing analysis on matrices.
# The modules this will access will focus on sparse high-dimensional data for single-cell rnaseq analysis
# Modules will be able to generically be pulled in and out of the analysis pipeline and modified though commandline.


# Install libraries if not installed
# Load libraries either way
library( optparse )

# A way to automatically install and load packages
# This reduces the headache of a user running the script
# Needs a little troubleshooting and not priority
#vctr_libraries = c( "logging" )
#for( str_library in vctr_libraries )
#{
#  if( !require( str_library, character.only=TRUE ) )
#  {
#    install.packages( pkgs=str_library, repos="http://cran.us.r-project.org" )
#  }
#  library( str_library )
#}



# Describe data options
STR_BEFORE = "before"
STR_BOTH = "both"
STR_AFTER = "after"
STR_NEVER = "never"


run <- function(
### This runs the analysis script based on the command line arguments given by the user
pArgs
### Contains all the command line given by the user.
 ){

  # Parse the arguments
  lsArgs <- parse_args( pArgs, positional_arguments = TRUE )

  # Read in libraries and code
  library( caret )
  library( GMD )
  library( gplots )
  library( logging )
  library( mclust )
  library( RColorBrewer )
  library( tools )
  library( vegan )

  # Source scripts
  source( file.path( "src", "analysis_run.R" ) )
  source( file.path( "src", "heatmap.3b.R" ) )
  source( file.path( "src", "Modules.R" ) )
  source( file.path( "src", "Plots.R" ) )
  source( file.path( "src", "Utilities.R" ) )

  # Data file name
  str_data_file = lsArgs$args[ 1 ]
  # Metadata file name
  str_metadata_file = lsArgs$options$str_metadata_file

  # Manage the output directory
  # Create a name if it is not given and manage creating it if needed.
  if( is.null( lsArgs$options$str_output_dir ) )
  {
    lsArgs$options$str_output_dir = file.path( getwd(), file_path_sans_ext( basename( str_data_file ) ) )
  }
  if( ! file.exists( lsArgs$options$str_output_dir ) )
  {
    dir.create( lsArgs$options$str_output_dir )
  }

  # Will eventually become a proper object to hold results of the run
  obj_analysis = analysis.run( lsArgs$options$str_output_dir ) 

  # Set correlations if given
  sample.corr.metric( obj_analysis ) = lsArgs$options$str_sample_corr
  feature.corr.metric( obj_analysis ) = lsArgs$options$str_sample_corr
  if( !is.null( lsArgs$options$str_feature_corr ) )
  {
    feature.corr.metric( obj_analysis ) = lsArgs$options$str_feature_corr
  }

  # Check the arguments
  # Should be added at some point but not priority yet.

  # Normalize string arguments
  lsArgs$options$str_describe_when = tolower( lsArgs$options$str_describe_when )

  # Create the logger
  logr_single_cell = getLogger( "singe_cell" )
  addHandler( writeToConsole, logr_single_cell )
  setLevel( lsArgs$options$str_verbosity, logr_single_cell )

  # Document arguments
  func_document_arguments( lsArgs, file.path( output.dir( obj_analysis ), "run_arguments.txt" ) )

  # Get the function to perform the analysis step
  func_batch_correct = func_get_batch_module( lsArgs$options$str_batch_method )
  func_discriminate = func_get_discriminant_module( lsArgs$options$str_discriminate_method )
  func_feature_selection = func_get_select_feature_groups_module( lsArgs$options$str_feature_selection )
  func_impute_metadata = func_get_impute_module( lsArgs$options$str_impute_metadata )
  func_normalization = func_get_normalization_module( lsArgs$options$str_normalization )
  func_ordination = func_get_ordinate_module( lsArgs$options$str_ordination_method )
  func_sample_group_selection = func_get_select_sample_groups_module( lsArgs$options$str_select_sample_groups_method )
  func_feature_group_selection = func_get_select_feature_groups_module( lsArgs$options$str_select_feature_groups_method )

  # Read in data table
  loginfo( "Reading data.", logr_single_cell )
  data( obj_analysis, islogged = FALSE ) = func_read_frame( str_data_file )
  str_data_file = func_update_name_to_dir( str_data_file, output.dir( obj_analysis ) )

  # Keeps track of steps performed as analysis occurs
  # Numbering is put in the name to help guide user in looking through results
  i_analysis_step = 1

  # Read in metadata table
  df_metadata = NULL
  if( !is.null( str_metadata_file ) )
  {
    loginfo( "Reading metadata.", logr_single_cell )
    metadata( obj_analysis ) = func_read_frame( str_metadata_file )
    str_metadata_file = func_update_name_to_dir( str_metadata_file, output.dir( obj_analysis ) )
  }

  # Reduce data to a targeted list if given
  if( !is.null( lsArgs$options$str_targeted_features_file ) )
  {
    vctr_features = func_read_list_file( lsArgs$options$str_targeted_features_file )
    vctr_features = vctr_features[ vctr_features %in% row.names( obj_analysis$data ) ]
    data( obj_analysis, islogged = is.logged( obj_analysis ) ) = data( obj_analysis )[vctr_features,]
  }

  if( lsArgs$options$str_describe_when %in% c( STR_BOTH, STR_BEFORE ) )
  {
    # Describe the matrix before filtering
    loginfo( "Describing RAW matrix/s." )
    ret_path = func_make_analysis_directory( output.dir( obj_analysis ), "describe_matrix", i_analysis_step )
    i_analysis_step = ret_path$next_step
    func_describe_count_data_matrix( obj_analysis=obj_analysis, str_output_dir=ret_path$path, str_cor_metric="euclidean" )
    if( !is.null( str_metadata_file ) )
    {
      func_describe_metadata( df_metadata = metadata( obj_analysis ), str_output_dir = ret_path$path )
    }
  }

  # Perform filters
  # If describing is not to happen util after filtering do no give filtering an output dir to plot into
  # If one does the filters will plot the data in similar ways to the describing data step which may not be possible
  # depending on the size of the plots
  loginfo( "Start filtering.", logr_single_cell )
  if( ! is.null( lsArgs$options$str_filter_method ) )
  {
    for( str_filter in unlist( strsplit( lsArgs$options$str_filter_method, ";" ) ) )
    {
      func_filter = func_get_filter_module( str_filter )
      if( is.null( func_filter ) )
      {
        logerror( paste( "The following filter is not known: ", str_filter, " analysis stopped."), logr_single_cell )
        stop( "Abnormal termination, analysis is incomplete. Please DO NO RELY on the results of this run. Filter unknown." )
      }
      ret_path = func_make_analysis_directory( output.dir( obj_analysis ), paste( str_filter, "filter", sep="_"), i_analysis_step )
      i_analysis_step = ret_path$next_step
      loginfo( paste( "Data dimension before filtering:", nrow( data( obj_analysis ) ), "x", ncol( data( obj_analysis ) ) ) )
      obj_analysis = func_filter( obj_analysis, str_output_dir = ifelse ( lsArgs$options$str_describe_when == STR_AFTER, NA, ret_path$path ) )
      loginfo( paste( "Data dimension AFTER filtering:", nrow( data( obj_analysis ) ), "x", ncol( data( obj_analysis ) ) ) )
    }
    str_data_file = func_add_keyword_to_file_name( str_data_file, "filtered" )
    func_write_frame( data( obj_analysis ), str_data_file )
  }

  if( lsArgs$options$str_describe_when %in% c( STR_BOTH, STR_AFTER ) )
  {
    # Describe the matrix
    loginfo( "Describing QCed matrix/s." )
    ret_path = func_make_analysis_directory( output.dir( obj_analysis ), "describe_qced_matrix", i_analysis_step )
    i_analysis_step = ret_path$next_step
    func_describe_count_data_matrix( obj_analysis=obj_analysis, str_output_dir=ret_path$path )
    # Describe metadata
    if( !is.null( str_metadata_file ) )
    {
      func_describe_metadata( df_metadata = metadata( obj_analysis ), str_output_dir = ret_path$path )
    }
  }

  # Perform checks and
  # Perform quality control
  loginfo( "Start quality control.", logr_single_cell )
  if( lsArgs$options$f_qc )
  {
    ret_path = func_make_analysis_directory( output.dir( obj_analysis ), "qc", i_analysis_step )
    i_analysis_step = ret_path$next_step
    obj_analysis = func_QC( obj_analysis=obj_analysis, func_handle_outliers=func_handle_outliers, 
                        func_imputation = func_impute_metadata, func_batch_correct=func_batch_correct, ret_path$path )
    str_data_file = func_add_keyword_to_file_name( str_data_file, "data_qc" )
    str_metadata_file = func_add_keyword_to_file_name( str_metadata_file, "metadata_qc" )
    func_write_frame( data( obj_analysis ), str_data_file )
    func_write_frame( metadata( obj_analysis ), str_metadata_file )
  }

  # Visualize using ordinations
  if( ! is.null( lsArgs$options$str_ordination_method ) )
  {
    loginfo( "Start ordinations.", logr_single_cell )
    ret_path = func_make_analysis_directory( output.dir( obj_analysis ), "ordination", i_analysis_step )

    i_analysis_step = ret_path$next_step

    for( str_ordination in unlist( strsplit( lsArgs$options$str_ordination_method, ";" ) ) )
    {
      func_ordination = func_get_ordinate_module( str_ordination )
      if( is.null( func_ordination ))
      {
        logerror( paste( "Did not understand the request for the ordination", str_ordination, ", stopped analysis.", sep = " " ) )
      }
      obj_analysis = func_ordination( obj_analysis=obj_analysis, str_output_dir = ret_path$path )

      # Write ordination to file
      func_write_frame( ordinations( obj_analysis )[[str_ordination]],
                        func_add_keyword_to_file_name( str_data_file, paste( str_ordination,"ordination","dims",sep="_" )))
    }
  }

  # Find sample groups
  loginfo( "Start sample group selection.", logr_single_cell )
  ret_path = func_make_analysis_directory( output.dir( obj_analysis ), "select_sample_groups", i_analysis_step )
  i_analysis_step = ret_path$next_step
  obj_analysis = func_sample_group_selection( obj_analysis, ret_path$path )
  func_document_groupings_selection( sample.groups( obj_analysis ), ret_path$path )

  # Find discriminant features (genes) associated with sample groups
  loginfo( "Start discriminate genes.", logr_single_cell )
  ret_path = func_make_analysis_directory( output.dir( obj_analysis ), "discriminate_features", i_analysis_step )
  i_analysis_step = ret_path$next_step
  obj_analysis = func_discriminate( obj_analysis, ret_path$path )
  # TODO document discriminate
iiiiii
  # group discriminant features (make gene groups)
  loginfo( "Start gene group selection.", logr_single_cell )
  ret_path = func_make_analysis_directory( output.dir( obj_analysis ), "select_feature_groups", i_analysis_step )
  i_analysis_step = ret_path$next_step
  func_feature_group_selection( obj_analysis, ret_path$path )
  func_document_groupings_selection( feature.groups( obj_analysis ), ret_path$path )

zzzzzzzzzzzzzz
  # Find functional descriptions of feature groups
  loginfo( "Start functional analysis", logr_single_cell )
  ret_path = func_make_analysis_directroy( output.dir( obj_analysis ), "functional_analysis", i_analysis_step )
  i_analysis_step = ret_path$next_step
#  func_functionally_describe_feature_groups( obj_analysis, ret_path$path )

  loginfo( "Successfully completed.", logr_single_cell )
}

# Parse commandline
pArgs <- OptionParser( usage = "%prog [options] <input_matrix.txt>" )
pArgs <- add_option( pArgs, c("-a","--describe"), type="character", action="store", dest="str_describe_when", metavar="Describe_when", default=STR_AFTER, help= paste( "Indicates when to desribe the data. Either before filtering, after filtering, or both. [Use the values ", STR_NEVER, STR_BEFORE, ", ", STR_AFTER, ", or ", STR_BOTH, "]", sep="" ) )
pArgs <- add_option( pArgs, c("-b","--batch_control"), type="character", action="store", dest="str_batch_method", metavar="Batch_method", default=NULL, help="Method for controlling for batch affects." )
pArgs <- add_option( pArgs, c("-c","--correlation" ), type="character", action="store", dest="str_sample_corr", metavar="correlation", default="spearman", help="Correlation used for samples, use for features as well if the feature correlation is not indicated. These correlation metrics will be used as need for clustering and ordination." )
pArgs <- add_option( pArgs, c("-d","--discriminate"), type="character", action="store", dest="str_discriminate_method", metavar="Discriminant_method", default="anova", help="Method to use for discrimination of genes for sample groups. [Default %default]")
pArgs <- add_option( pArgs, c("-e","--feature_correlation" ), type="character", action="store", dest="str_feature_corr", metavar="feature_correlation", default=NULL, help="Correlation used only for features. These correlation metrics will be used as need for clustering and ordination." )
pArgs <- add_option( pArgs, c("-f","--filter"), type="character", action="store", dest="str_filter_method", metavar="Filter", default="percentile;occurence;sparsity;sd;pca", help="Preprocessing filter to be ran on data. [Default %default]")
pArgs <- add_option( pArgs, c("-i","--impute"), type="character", action="store", dest="str_impute_metadata", metavar="impute_metadata", default = NULL, help="The method to use to handle NA and missing information in metadata." )
pArgs <- add_option( pArgs, c("-l","--limit"), type="character", action="store", dest="str_targeted_features_file", metavar="Limit_to_list", default=NULL, help="Limits the data to the features in this file. Should be a path to a file." )
pArgs <- add_option( pArgs, c("-m","--metadata"), type="character", action="store", dest="str_metadata_file", metavar="Metadata", default=NULL, help="Metadata file. [Default %default]" )
pArgs <- add_option( pArgs, c("-n","--normalization"), type="character", action="store", dest="str_normalization", metavar="Normalization", default=NULL, help="Normalization. Default %default]" )
pArgs <- add_option( pArgs, c("-r","--ordination"), type="character", action="store", dest="str_ordination_method", metavar="Ordination", default="pca;nmds", help="Ordinationmethod. [Default %default]")
pArgs <- add_option( pArgs, c("-s","--sample_group"), type="character", action="store", dest="str_select_sample_groups_method", metavar="Sample_groups", default="mclust", help="Select samples groups. [Default %default]")
pArgs <- add_option( pArgs, c("-o","--out_dir"), type="character", action="store", dest="str_output_directory", metavar="Output_Dir", default=NULL, help="Directory to place output. [Default path is made from the input file name and placed in the current working directory]")
pArgs <- add_option( pArgs, c("-q","--qc"), type="logical", action="store_false", dest="f_qc", metavar="QC_off", default=TRUE, help="Turn off QC. [Default %default]")
pArgs <- add_option( pArgs, c("-t","--feature"), type="character", action="store", dest="str_select_feature_groups_method", metavar="Feature_groups", default="dendrogram", help="Method to select features to use." )
pArgs <- add_option( pArgs, c("-v","--verbosity"), type="character", action="store", dest="str_verbosity", metavar="Verbosity", default="DEBUG", help="Sets the logging level, controlling more or less logging about the internal processes. [Default %default]")

run( pArgs )
