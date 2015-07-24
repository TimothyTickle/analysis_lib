###########################
# File and File I/O
###########################


func_read_frame <- function(
str_path_to_frame
){
  return( read.csv( str_path_to_frame, row.names = 1, header = TRUE, sep = "\t" ) )
}


func_read_list_file <- function(
### Read in a file which is a list
str_path_to_file
){
  return( as.vector( as.matrix( read.table( str_path_to_file ) ) ) )
}


func_add_keyword_to_file_name <- function(
### Tacks on a keyword to the end of a file name before the extension
str_path,
### Path to file
str_key
### Keyword to add
){
  if( is.null( str_path ) )
  {
    return( str_path )
  }
  return( paste( file_path_sans_ext( str_path ), "_", str_key, ".", file_ext( str_path ), sep = "" ) )
}


func_update_name_to_dir <- function(
str_file_name,
str_directory
){
  if( is.null( str_file_name ) || is.null( str_directory ) )
  {
    return( str_file_name )
  }
  return( file.path( str_directory, basename( str_file_name ) ) )
}


func_write_frame <- function(
df_frame,
str_path_to_save_file,
f_append = FALSE
){
  if( ( !is.null( df_frame ) ) && ( !is.null( str_path_to_save_file ) ) )
  {
    write.table( df_frame, str_path_to_save_file, quote = FALSE, sep = "\t", col.names = NA, append = f_append )
  }
}


func_write_text <- function(
### Writes text to a file
str_text,
### The text to write
str_file = NULL,
### The file to write to
f_append = FALSE
){
  if( !is.null( str_file ) )
  {
    write( x = str_text, file = str_file, ncolumns = 1, append = f_append )
  }
}


func_make_analysis_directory <- function(
### Creates a directory to place analysis in
str_parent_path,
### Path to place the dirctory in
str_analysis_keyword,
### Keyword indicating the type of analysis placed in directory
i_analysis_step
### Which step in the pipeline the analysis is
){
  str_dir_path = file.path( str_parent_path, paste( i_analysis_step, str_analysis_keyword, sep="_" ) )
  if( ! file.exists( str_dir_path ) )
  {
    dir.create( str_dir_path )
  }
  i_analysis_step = i_analysis_step + 1
  return( list( path = str_dir_path, next_step = i_analysis_step ) )
}


##############################
# Standardized color schemes
##############################


func_monochromatic_palette <- function(
### Return the standard monochromatic color scheme
i_number_colors = 100
){
  return( colorpanel( i_number_colors, "white", "violet", "purple" ) )
}


func_polychromatic_palette <- function(
### Return the standard polychromatic color scheme
i_number_colors = 100
){
  return( colorpanel( i_number_colors, "purple", "black", "yellow" ) )
}


func_metadata_palette <- function(
### Return the standard color scheme for metadata
i_number_colors = 100
){
  return( rainbow( i_number_colors ) )
}


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


###############################
# Self-documentation
###############################


func_document_arguments <- function(
### Documents the run arguments to a file
args_run,
### The commandline arguments for the run
str_file_path
### The file to write the argument selections
){
  func_write_text( "Run arguments:", str_file_path )
  func_write_text( paste( "Input data matrix: ", args_run$args[1] ), str_file_path, f_append = TRUE )
  for( str_name in names( args_run$options ) )
  {
    func_write_text( paste( str_name, args_run$options[[str_name]], sep = ": " ), str_file_path, f_append = TRUE )
  }
}


func_document_groupings_selection <- function(
### List of groupings to document, list should contain named vectors
list_groupings,
### Holds the data associated with the current analysis
str_output_dir
### Directory to output data and visualizations
){
  for( str_selection in names( list_groupings ) )
  {
    func_write_frame( list_groupings[[ str_selection ]], 
                      file.path( str_output_dir, paste( "groupings_",str_selection,".txt", sep="" ) ) )
  }
}


###############################
# Statistical Methods
###############################


func_select_pca_components_stick_breaking <- function(
### Select important pca components by comparing to a stick breaking distribution
vctr_i_percent_variance,
str_output_dir = NULL
){
  # Get the stick break distribution
  i_max_pc = 1
  f_record_break = TRUE
  vctr_break_distribution = c()
  i_variance_count = length( vctr_i_percent_variance )
  for( i_k_break in 1:i_variance_count )
  {
    vctr_break_distribution = c(vctr_break_distribution, (1/i_variance_count) * sum( 1 / ( i_k_break : i_variance_count ) ) )
    if( f_record_break && ( vctr_break_distribution[ i_k_break ] < vctr_i_percent_variance[ i_k_break ] ) )
    {
      i_max_pc = i_k_break
    } else {
      f_record_break = FALSE
    }
  }
  i_max_pc = max( i_max_pc, 2 )

  # If an output directroy is given plot the stick breaking distribution.
  if( ! is.null( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "Stick_breaking_distribution.pdf"), useDingbats = FALSE )
    vctr_str_comp_colors = rep( "grey",i_variance_count )
    vctr_str_comp_colors[ 1: i_max_pc ] = "red"
    vctr_positions = barplot( vctr_i_percent_variance, ylim = c(0, max( vctr_break_distribution )*1.1), xlim = c(0, i_variance_count + 3), main = "Percent variance per component", col = vctr_str_comp_colors, xlab = "Component", ylab = "Percent variance" )
    points( vctr_positions , vctr_break_distribution, pch = 21, bg = vctr_str_comp_colors, col ="black", add=TRUE )
    lines( vctr_positions, vctr_break_distribution )
    dev.off()
  }
  return( list( count = i_max_pc, distribution = vctr_break_distribution ) )
}

