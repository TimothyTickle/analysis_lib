################
# This is an object to keep the state of an analysis run.
#
# This object allows one to understand if the data is currently
# transformed (like logged) or if it already has certain information
# computed so that it does not have to happen again (like distance matrices).
################
# Class specification follows recommendations from the following post
# http://www.johnmyleswhite.com/notebook/2009/12/14/object-oriented-programming-in-r-the-setter-methods/
################

# Constructor
analysis.run <- function( 
### Constructor method for analysis runs
str_out_dir = NULL
### Output directory of the analysis run
### If not given, the current working directory will be used
){
  # Create object
  object = list( data = NULL, metadata = NULL, ordinations = list(), sample.groups = list(),
                 feature.groups = list(), discriminant.features = list(), is.logged = FALSE,
                 output.dir = NULL, sample.corr.metric = NULL, feature.corr.metric = NULL,
                 sample.dist = NULL, feature.dist = NULL )
  class( object ) = "analysis.run"

  # Set the output directory
  if( is.na( str_out_dir ) || is.null( str_out_dir ) )
  {
    str_out_dir = getwd()
  }
  object$output.dir = str_out_dir

  # Return the created object instance
  return( object )
}


## Correlation metrics
# Get feature correlation
feature.corr.metric <- function( x ) UseMethod( 'feature.corr.metric', x )
feature.corr.metric.analysis.run <- function( x ){
  return( x$feature.corr.metric )
}

# Set feature correlation
'feature.corr.metric<-' <- function( x, value ) UseMethod( 'feature.corr.metric<-', x )
'feature.corr.metric<-.analysis.run' <- function( x, value )
{
  x$feature.corr.metric = value
  x$feature.dist = NULL
  return( x )
}


# Get sample correlation
sample.corr.metric <- function( x ) UseMethod( 'sample.corr.metric', x )
sample.corr.metric.analysis.run <- function( x ){
  return( x$sample.corr.metric )
}

# Set sample correlation
'sample.corr.metric<-' <- function( x, value ) UseMethod( 'sample.corr.metric<-', x )
'sample.corr.metric<-.analysis.run' <- function( x, value ) {
  x$sample.corr.metric = value
  x$sample.dist = NULL
  return( x )
}


### Data
### Log data
### Transform the given data frame to log2
### Expects values of ) or greater
### 1 is added to the value before transformation
### Transforms within columns
log <- function( x ) UseMethod( "log", x )
log.analysis.run <- function( x ){
  if( !is.logged( x ) )
  {
    data( x, islogged = TRUE ) = func_log( data( x ) )
    x$is.logged = TRUE
    x$sample.dist = NULL
    x$feature.dist = NULL
    ordinations( x ) = list()
    sample.groups( x ) = list()
    feature.groups( x ) = list()
    discriminant.features( x ) = list()
  }
  return( data( x ) )
}

# Is logged
# Get is logged
is.logged <- function( x ) UseMethod( "is.logged", x )
is.logged.analysis.run <- function( x ){
  return( x$is.logged )
}

# Set is logged
# Should not have to set if the data i logged except when the data is set

# Get data
data <- function( x ) UseMethod( "data", x )
data.analysis.run <- function( x ){
  return( x$data )
}

# Set data
# When setting data the distance matrices should be reset
# Make sure to get if the data is logged or not 
'data<-' <- function( x, value, islogged ) UseMethod( 'data<-', x )
'data<-.analysis.run' <- function( x, value, islogged ){
  # Update data
  x$data = value
  x$is.logged = islogged

  # Invalidate distance matrices
  x$sample.dist = NULL
  x$feature.dist = NULL
  ordinations( x ) = list()
  sample.groups( x ) = list()
  feature.groups( x ) = list()
  discriminant.features( x ) = list()

  return( x )
}

# Return a matrix of the data that is logged but does not update
# the analysis run
get.logged.matrix <- function( x ) UseMethod( 'get.logged.matrix', x )
get.logged.matrix.analysis.run <- function( x ){
  df_return = data( x )
  if( ! is.logged( x ) )
  {
    df_return = log2( df_return + 1 )
  }
  return( df_return )
}

# Return a matrix of data that is NOT logged , this does not update the analysis run object
get.unlogged.matrix <- function( x ) UseMethod( 'get.unlogged.matrix', x )
get.unlogged.matrix.analysis.run <- function( x ){
  df_return = data( x )
  if( is.logged( x ) )
  {
    df_return = ( 2^df_return ) - 1
  }
  return( df_return ) 
}


## Sample distance matrix
# Get sample distance matrix
sample.dist <- function( x ) UseMethod( 'sample.dist', x )
sample.dist.analysis.run <- function( x ){
  # Check to see if we have a stored matrix
  # If not, make one and store
  if( is.null( x$sample.dist ) )
  {
    x$sample.dist = make.sample.distance.matrix( x, f_log = is.logged( x ) )
  }
  return( x$sample.dist )
}

# Set sample distance matrix
# No need to set the sample distance matrix
# Should just get and one will be created as needed. 

# Create a sample distance matrix
make.sample.distance.matrix <- function(
obj_analysis,
f_log
){
  # Try pulling from or storing in cache first if appropriate
  if( f_log && is.logged( obj_analysis ) )
  {
    if( is.null( obj_analysis$sample.dist ) )
    {
      obj_analysis$sample.dist = make.custom.distance.matrix( get.logged.matrix( obj_analysis ), sample.corr.metric( obj_analysis ) )
    }
    return( sample.dist( obj_analysis ) )
  }

  if( !f_log && !is.logged( obj_analysis ) )
  {
    if( is.null( obj_analysis$sample.dist ) )
    {
      obj_analysis$sample.dist = make.custom.distance.matrix( get.unlogged.matrix( obj_analysis ), sample.corr.metric( obj_analysis ) )
    }
    return( sample.dist( obj_analysis ) )
  }

  # If not storing
  if( f_log )
  {
    return( make.custom.distance.matrix( get.logged.matrix( obj_analysis ), sample.corr.metric( obj_analysis ) ) )
  } else {
    return( make.custom.distance.matrix( get.unlogged.matrix( obj_anlaysis ), sample.corr.metric( obj_analysis ) ) )
  }
}


## Feature distance matrix
# Get feature distance matrix
feature.dist <- function( x ) UseMethod( "feature.dist", x )
feature.dist.analysis.run <- function( x ){
  # Check to see if we have a stored matrix
  # If not, make one and store
  if( is.null( x$feature.dist ) )
  {
    x$feature.dist = make.feature.distance.matrix( x, f_log = is.logged( x ) )
  }
  return( x$feature.dist )
}

## Set feature distance matrix
# Should not set the distance matrix, just get it and it will be created as needed.

# Create a feature distance matrix
make.feature.distance.matrix <- function(
obj_analysis,
f_log = FALSE
){
  # Try pulling from or storing in cache first if appropriate
  if( f_log && is.logged( obj_analysis ) )
  {
    if( is.null( obj_analysis$feature.dist ) )
    {
      obj_analysis$feature.dist = make.custom.distance.matrix( t( get.logged.matrix( obj_analysis ) ), feature.corr.metric( obj_analysis ) )
    }
    return( feature.dist( obj_analysis ) )
  }

  if( !f_log && !is.logged( obj_analysis ) )
  {
    if( is.null( obj_analysis$feature.dist ) )
    {
      obj_analysis$feature.dist = make.custom.distance.matrix( t( get.unlogged.matrix( obj_analysis ) ), feature.corr.metric( obj_analysis ) )
    }
    return( feature.dist( obj_analysis ) )
  }

  if( f_log )
  {
    return( make.custom.distance.matrix( t( get.logged.matrix( obj_analysis ) ), feature.corr.metric( obj_analysis ) ) )
  } else {
    return( make.custom.distance.matrix( t( get.unlogged.matrix( obj_analysis ) ), feature.corr.metric( obj_analysis ) ) )
  }
}


## Discriminant features
# Get discriminant features
discriminant.features <- function( x ) UseMethod( 'discriminant.features', x )
discriminant.features.analysis.run <- function( x, value ){
  return( x$discriminant.features )
}

# Set discriminant features
'discriminant.features<-' <- function( x, value ) UseMethod( 'discriminant.features<-', x )
'discriminant.features<-.analysis.run' <- function( x, value ){
  x$discriminant.features = value
  return( x )
}


## Feature.groups
# get feature.groups
feature.groups <- function( x ) UseMethod( 'feature.groups', x )
feature.groups.analysis.run <- function( x, value ){
  return( x$feature.groups )
}

# set feature.groups
'feature.groups<-' <- function( x, value ) UseMethod( 'feature.groups<-', x )
'feature.groups<-.analysis.run' <- function( x, value ){
  x$feature.groups = value
  return( x )
}


## Metadata
# Get metadata
metadata <- function( x ) UseMethod( 'metadata', x )
metadata.analysis.run <- function( x, value ){
  return( x$metadata )
}

# Set metadata
'metadata<-' <- function( x, value ) UseMethod( 'metadata<-', x )
'metadata<-.analysis.run' <- function( x, value ){
  x$metadata = value
  return( x )
}


## Ordination
# Get ordinations
ordinations <- function( x ) UseMethod( 'ordinations', x )
ordinations.analysis.run <- function( x, value ){
  return( x$ordinations )
}

# Set ordinations
'ordinations<-' <- function( x, value ) UseMethod( 'ordinations<-', x )
'ordinations<-.analysis.run' <- function( x, value ){
  x$ordinations = value
  return( x )
}


## Sample.groups
# Get sample.groups
sample.groups <- function( x ) UseMethod( 'sample.groups', x )
sample.groups.analysis.run <- function( x ){
  return( x$sample.groups )
}

# Set sample.groups
'sample.groups<-' <- function( x, value ) UseMethod( 'sample.groups<-', x )
'sample.groups<-.analysis.run' <- function( x, value ){
  x$sample.groups = value
  discriminant.features( x ) = list()
  return( x )
}


## Output directory
# Get output dir
output.dir <- function( x ) UseMethod( 'output.dir', x )
output.dir.analysis.run <- function( x ){
  return( x$output.dir )
}

# Set output dir
# This is set in the constructor and then should not change
'output.dir<-' <- function( x, value ) UseMethod( 'output.dir<-', x )
'output.dir<-.analysis.run' <- function( x, value ){
  print( "Not changing output directory." )
  return( x )
}


##########################
# Utility functions
##########################

subset.samples <- function( 
### Should be used to subset samples to make sure
### both the data and the metadata is updated together
obj_analysis,
### Analysis run object to update 
vctr_indices
### Sample indices to keep
){
  data( obj_analysis, islogged = is.logged( obj_analysis ) ) = data( obj_analysis )[,vctr_indices]
  if( !is.null( metadata( obj_analysis ) ) )
  {
    metadata( obj_analysis ) = metadata( obj_analysis )[vctr_indices,]
  }
}

make.custom.distance.matrix <- function(
### Will create a symmetric distance matrix given the str_metric
### Will return null if the distance metric is not handled.
df_data,
### Data to act on, it is assumed the columns are
### subjects / samples and the rows are observations / features
str_metric
### Distance metric to return
){
  print(str_metric)
  # Normalize metric keyword
  if( is.null( str_metric ) )
  {
    print( "str_metric is null")
    return( NULL )
  }
  str_key = tolower( str_metric )

  # Default value for return is NULL incase metric is unknown
  func_ret = NULL

  # Get function
  func_ret = switch( str_metric,
    altGower = vegdist( as.matrix( df_data ), method = str_metric ),
    binomial = vegdist( as.matrix( df_data ), method = str_metric ),
    bray = vegdist( as.matrix( df_data ), method = str_metric ),
    canberra = vegdist( as.matrix( df_data ), method = str_metric ),
    cao = vegdist( as.matrix( df_data ), method = str_metric ),
    chao = vegdist( as.matrix( df_data ), method = str_metric ),
    euclidean = vegdist( as.matrix( df_data ), method = str_metric ),
    gower = vegdist( as.matrix( df_data ), method = str_metric ),
    horn = vegdist( as.matrix( df_data ), method = str_metric ),
    jaccard = vegdist( as.matrix( df_data ), method = str_metric ),
    kulczynski = vegdist( as.matrix( df_data ), method = str_metric ),
    kendall = cor( df_data, use = "all.obs", method = "kendall" ),
    manhattan = vegdist( as.matrix( df_data ), method = str_metric ),
    morisita = vegdist( as.matrix( df_data ), method = str_metric ),
    mountford = vegdist( as.matrix( df_data ), method = str_metric ),
    none = NULL,
    pearson = cor( df_data, use = "all.obs", method = "pearson" ),
    raup = vegdist( as.matrix( df_data ), method = str_metric ),
#TODO    reciprocal = ,
    spearman = cor( df_data, use = "all.obs", method = "spearman" ))
  # Return function
  print("End Dim")
  print( dim( func_ret ) )
  print("endend dim")
  return( as.dist( func_ret ) )
  ### Returns the function to create the needed
}

