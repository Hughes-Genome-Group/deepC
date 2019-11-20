# FUNCTIONS FOR HiC DATA ANALYSIS AND VISUALISATION
# Ron Schweßinger
# 
library(tidyverse)
library(RColorBrewer)
library(grid, gridExtra)
library(dplyr)

#define theme for plotting
# science_theme <- theme(
#   panel.grid.major = element_line(size = 0.5, color = "grey"),
#   axis.line = element_line(size = 0.7, color = "black"),
#   text = element_text(size = 14, family="Arial")
# )

science_theme <- theme(
  panel.grid.major = element_line(size = 0.5, color = "grey"),
  panel.grid.minor = element_blank(),
  # text = element_text(size = 14, family="Arial"),
  axis.line = element_line(color="black", size = 0.7),
  axis.line.x = element_line(color="black", size = 0.7),
  axis.line.y = element_line(color="black", size = 0.7)
)

# define themes for metling plots
melting_theme_hic <- theme(axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks=element_blank(),
                           axis.line = element_blank(),
                           panel.border = element_blank(),
                           axis.line.x = element_blank(),
                           axis.line.y = element_blank())

melting_theme_signal <- theme(
  axis.text.x=element_blank(),
  axis.title.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.line.x=element_blank(),
  panel.border=element_blank()
)

melting_theme_signal_bottom <- theme(
  # plot.margin = unit(c(0,1,1,1), "lines")
)


# FUNCTIONS ====================================================================

ConvertBinToPos <- function(x, bin.size, start.pos, start.id){
  # Format bin ids to genomic postion
  x <- ((x - start.id) * bin.size) + start.pos
  
  return(x)
}

FormatBp <- function(bp, type="M"){
# format a single p positon to kilo (k) or Mega (M) given type=k/M
  
  if(type == "k"){
    bp <- round(bp/1000, digits=2) 
    bp <- paste0(bp, "k")
  }
  
  if(type == "M"){
    bp <- round(bp/1000000, digits=2)
    bp <- paste0(bp, "M")
  }
  
  return(bp)
}

GetBinsStart <- function(hic){
  # Function to retrieve the start of the first bin of a hic object
  x <- hic$start.pos - hic$bin.size/2
  return(x)
}

GetBinsEnd <- function(hic){
  # Function to retrieve the end of the last bin of a hic object
  x <- hic$coords$center[length(hic$coords$center)] + hic$bin.size/2
  return(x)
}

GetDifferenceMatrix <- function(a, b){
  # Helper function to unify and caluclate the difference of two matrices (must be normalised towards each other)
  # Arguments:
  #   a: matrix a
  #   b: matrix b
  # Returns: Matrix with delta = a - b column
  
  a <- a$matrix.df
  b <- b$matrix.df
  
  # combine bins to unique id
  a$id <- apply(a, 1, function(x){x <- paste0(x[1], ":", x[2]) })
  b$id <- apply(b, 1, function(x){x <- paste0(x[1], ":", x[2]) })
  # sort and combine unified
  comb.id <- unique(sort(c(a$id, b$id)))
  # initialize and empty df with 0 entries for every uniue bin pair id
  d <- data.frame(id=comb.id)
  d$y <- sapply(d$id, function(x) { x <- as.numeric(unlist(strsplit(as.character(x), ":"))[1]) })
  d$x <- sapply(d$id, function(x) { x <- as.numeric(unlist(strsplit(as.character(x), ":"))[2]) })
  d <- d[order(d$y, d$x),]
  
  # get interactions present in a and present b assigned to the respective column
  d$a <- rep(0, nrow(d))
  d$b <- rep(0, nrow(d))
  d$a[which(d$id %in% a$id)] <- a$value
  d$b[which(d$id %in% b$id)] <- b$value
  d$a <- as.numeric(d$a)
  d$b <- as.numeric(d$b)
  
  # make differential column and get rid of a and b and id column
  d$value <- d$a - d$b
  d <- d[,-c(1,4,5)]
  
  # return 
  return(d)
}

ImportHicproMatrix <- function(matrix, coords="empty", chr="empty", bin.size=0, nrow=-1){
  # Wrapper function for importing (pruned HiCPro) matrix and bin coords
  # Import a HiCPro matrix and the according bed coordinates assining the ids
  # Add a centric, single coord per bin region
  # Calculate bin size and extract first id and genomic position
  # Store as list including a ggplot2 friendly data frame
  #
  # Args:
  #   matrix: 3 column matrix as result frm HiCPro or pruned version (id.x id.y, interaction)
  #   coords: 4 column bed file (chr, start, end, id) 
  #   if left "empty" will assume coords can be read from the bins of the matrix and extract it from there
  #   chr: needs to specified to add a chr to a matrix without bed format
  #   nrows: number of rows to read in
  #
  # Returns:
  #   list: .$matrix.df .$bin.size .$start.pos .$start.id
  
  #read matrix
  matrix.in <- read.table(matrix, header=FALSE, colClasses = rep("numeric", 3), nrow=nrow)
  #make data frame
  # matrix.df <- as.data.frame(matrix.in)
  matrix.df <- as.tibble(matrix.in)
  colnames(matrix.df) <- c("y","x", "value")
  
  # if left "empty" will assume coords can be read from the bins of the matrix and extract it from there
  if(coords == "empty"){
    if(chr == "empty"){
      # check that a chromosome to name is supplied
      stop("Must supply a chromosome name when no bed coord file is present")
    }
    if(bin.size == 0){
      # check that a chromosome to name is supplied
      stop("Must supply a positive, matching bin.size if no coord file present")
    }
    # create pseudo coord file
    bins <- unique(sort(matrix.df$y))
    # # get bin size
    # bin.size <- abs(bins[1] - bins[2])
    
    # make ids
    bin.ids <- bins/bin.size
    # convert to id coord dataframe
    coords.in <- as.tibble(data.frame(
      chr = rep(chr, length(bins)),
      start = bins,
      end = bins + bin.size,
      id = bin.ids
    ))
    # reduce matrix to ids
    matrix.df$y <- matrix.df$y/bin.size
    matrix.df$x <- matrix.df$x/bin.size
    
  }else{
    #read and format coords
    coords.in <- as.tibble(read.table(coords, header=FALSE))
    colnames(coords.in) <- c("chr", "start", "end", "id")
    bin.size <- coords.in$end[1] - coords.in$start[1] # get bin size
  }

  # add centrier column
  coords.in$center <- coords.in$start + ((coords.in$end - coords.in$start)/2) #add a centric column

  #get genomic position and id from start bin
  start.id <- min(matrix.df$x)
  start.pos <- coords.in[coords.in$id == start.id, ]$center
  
  # add x.id and y.id column to identify across chromosomes
  matrix.df$y.id <- matrix.df$y
  matrix.df$x.id <- matrix.df$x
  
  #convert bin ids to center of bin region coords
  matrix.df$x <- ConvertBinToPos(matrix.df$x,
                                     bin.size=bin.size,
                                     start.pos=start.pos,
                                     start.id=start.id)

  matrix.df$y <- ConvertBinToPos(matrix.df$y,
                                     bin.size=bin.size,
                                     start.pos=start.pos,
                                     start.id=start.id)
  
  
  # sort after$y
  # matrix.df <- matrix.df[order(matrix.df$y, matrix.df$x), ]
  
  return(list(matrix.df=matrix.df, coords=coords.in, bin.size=bin.size, start.pos=start.pos, start.id=start.id))
  
}

IntersectBedDataframe <- function(df, chr, start, end){
  # Wrapper to subset/intersect a given BED format df with a region of interest
  #
  # Input:
  #   df: bedlike 4 or more column df
  #   chr: chromsome of interest
  #   start: start coord of interest
  #   end: end coord of interest
  
  df <- df[df[,1] == chr,] #chr

  df <- df[
    (df[, 2] >= start & df[, 2] < end) |
    (df[, 3] > start & df[, 3] <= end) |
    (df[, 2] <= start & df[, 3] >= end),] # coords
  
  return(df)
  
}

MakeAnnotationDataFrame <- function(df, name="annotation", row = "1"){
  # Adjust a bed like data frame for annotation ggplot2 like plotting
  #
  # Input:
  #   df: bed like data frame (min 3 columns(chr, start, end, (name)))
  #   name: name for the annotation to appear
  #   row: row in which annotation should appear 1 lowest ++ (default=1)
  #
  # Returns: 
  # Adjusted bed like data frame with annotation row, up and down border and annotation name column
  
  df$row <- rep(row, nrow(df))
  #   df$down <- rep(anno.row-0.5, nrow(df))
  #   df$up <- rep(anno.row+0.5, nrow(df))
  df$annotation <- rep(name, nrow(df))
  
  return(df)
  
}

MakeBedGraphDataFrameForSignal <- function(bg){
  
  bg2 <- bg[, c(1,3,4)]
  bg <- bg[, c(1,2,4)]
  colnames(bg2) <- c("chr", "pos", "value")
  colnames(bg) <- c("chr", "pos", "value")
  
  #correct bed coords
  bg$pos <- bg$pos + 1
  
  bg <- rbind(bg, bg2)  # join
  
  bg <- bg[order(bg$chr, bg$pos, decreasing=F),]  # sort
  
  return(bg)
  
}

MakeMirrorMatrix <- function(df){
  # Mirror a half full HiC(Pro) ggplot2 df toto get, plot a full matrix df to plot
  # taking a single df as input
  
  mirrored <- df[ , c(2,1,3)]
  colnames(mirrored) <- c("y", "x", "value")
  mirrored <- subset(mirrored, x != y )
  df <- rbind(df, mirrored)
  
  return(df)
  
}

MakeTriangleMatrix <- function(df){
  # Take a single HiC(Pro) (non mirrored) data frame as input 
  # convert it to triangular like polygon df over genomic locus
  
  # estimate bin size bin
  bin.size <- diff( sort( unique(df$y) ) [c(1,2)])
  
  # 1) initialise 4 matrices by copying from half square df & add polygon pos and bin id column
  tri <- df
  tri$mean <- apply(tri, 1, function(x){ mean(c(x[1], x[2])) })  # add mean column
  tri$cy <- tri$cx <- rep(0, nrow(tri))  #init empty coord corrected columns
  
  tri1 <- tri2 <- tri3 <- tri4 <- tri #quadruple for polygon positions
  
  tri1$pos <- rep(1, nrow(tri1))  #indicate polygon positions
  tri2$pos <- rep(2, nrow(tri2))
  tri3$pos <- rep(3, nrow(tri3))
  tri4$pos <- rep(4, nrow(tri4))
  
  tri1$id <- tri2$id <- tri3$id <- tri4$id <- c(1:nrow(tri1)) #add polygon group id
  
  # 2) add adjusted coordinate columns
  tri1$cy <- (tri1$x - tri1$y)/bin.size
  tri2$cy <- (tri2$x - tri2$y)/bin.size - 1
  tri3$cy <- (tri3$x - tri3$y)/bin.size
  tri4$cy <- (tri4$x - tri4$y)/bin.size + 1
  
  tri1$cx <- tri1$mean - bin.size/2
  tri2$cx <- tri2$mean
  tri3$cx <- tri3$mean + bin.size/2
  tri4$cx <- tri4$mean
  
  # 3) rbind and replace y and x coord columns
  tri <- rbind(tri1, tri2, tri3, tri4)
  tri$x <- tri$cx
  tri$y <- tri$cy
  tri <- tri[,c("y", "x", "value", "pos", "id")]
  tri <- tri[order(tri$id, tri$pos, decreasing = FALSE ), ] # sort
  
  tri$y[tri$y < 0] <- 0 #  for negative y
  
  return(tri)
  
}

PruneHicproMatrix <- function(hics, chrs, starts = 0, ends = 0){
  # Function to prune a matrix to zoom into a region of interest
  # Take a list object as imported with ImportHicproMatrix
  # and prune it to include only interactions between the genomix coordinates chosen
  
  if(starts == 0 & ends == 0)
    # only chromosome
    hics$coords <- hics$coords %>%
      filter(chr == chrs)
  else{
    # get bins in ROI
    hics$coords <- hics$coords %>%
      filter(chr == chrs & center >= starts & center <= ends)
    
  }
  
  # get interactions in ROI
  hics$matrix.df <- hics$matrix.df[((hics$matrix.df$y.id %in% hics$coords$id) & (hics$matrix.df$x.id %in% hics$coords$id)),]
  
  # adjust start.id and start.pos (for plotting)
  hics$start.id <- min(hics$coords$id)
  hics$start.pos <- min(hics$matrix.df$y)
  
  return(hics)
  
}

SetValueRange <- function(vector, min=-Inf, max=Inf){
  # Function to cap the values for a max and min value range
  # Take a vector and a min and max numeric argument as input. 
  # Return the vector with every value below min = min and above max = max
  
  vector <- sapply(vector, function(x){
    if(x < min){ x <- min;}
    if(x > max){ x <- max;}
    return(x)
  })
  return(vector)
  
}

# Plotting ==============================================================

PlotAnnotation <- function(
  anno.df,
  xmin,
  xmax,
  break.number=4,
  format="s",
  pal="Set1"
){
  # Create an ggplot2 annotation like plotting under a triangular HiC plot and for multiplot grobbing
  #
  # Input:
  #   anno.df: annotation version of a bed like data frame (result from MakeAnnotationDataFrame)
  #   xmin: minimum x value to plot for limit
  #   xmax: maximum x value to plot for limit
  #   break.number: number of breaks to print on axis default=5
  #   format: options to format the genomic postion (s=single, k=in kilo, M=in Mega) default="s"
  #   pal: RColorBrewer palette to use default="YlOrRd"
  
  # colors <- brewer.pal(9, pal)
  
  # make breaks and lables
  step <- ((xmax - xmin) - ((xmax - xmin) %% (break.number-1)))/(break.number-1)
  breaks <- seq(from=xmin, to=xmax, by=step)
  labels <- sapply(breaks, function(x){ FormatBp(x, type=format)}) 
  
  # make plot
  a <- ggplot(anno.df, aes(fill=annotation)) +
    # geom_hline(yintercept = c(1:max(row))) +
    geom_rect(aes(xmin=start, xmax=end, ymin=row-0.45, ymax=row+0.45)) +
    scale_fill_brewer(palette = pal) +
    scale_x_continuous(breaks=breaks, labels=labels) +
    # geom_text(data=anno.df[anno.df$name %in% c(""), ], aes(label=name, x=as.numeric(start), y= 0.2), size=3.5, col="black", fontface=2) +
    coord_cartesian(xlim=c(xmin, xmax), ylim=c(0, max(anno.df$row) + 1)) +
    science_theme + 
    theme(
      legend.position="bottom",
      axis.text.y=element_blank(), 
      axis.line.y=element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank()
      # plot.margin = unit(c(0,1,0,1), "lines")
    )
  
  return(a)
  
}

PlotSquareMatrix <- function(
  hic,
  break.number = 5,
  square=FALSE,
  format = "s",
  pal = "YlOrRd"
){
  # Wrapper function to plot a matrix given a 3 column dataframe as input
  # Print a Square Matrix plot given a three column (y, x, value) and grafik options
  #
  # Args:
  #   hic object list as above
  #   break.number: number of breaks to print on axis default=5
  #   square: FALSE/TRUE [default FALSE] indicate if to mirror a trianglualr matrix to get a full matrix plot
  #   format: options to format the gneomic postion (s=single, k=in kilo, M=in Mega) default="s"
  #   pal: RColorBrewer palette to use default="YlOrRd"
  #
  # Returns:
  #   ggplot2 plot
  # select color scheme colours
  colors <- brewer.pal(9, pal) 
  
  # get matrix
  matrix.df <- hic$matrix.df
  
  # match and get chromosomal position from coords id
  matrix.df$y <- hic$coords$center[match(matrix.df$y.id, hic$coords$id)]
  matrix.df$x <- hic$coords$center[match(matrix.df$x.id, hic$coords$id)]
  matrix.df <- matrix.df[,c(1:3)]
  
  #mirror data frame if selected
  if(square){
    matrix.df <- MakeMirrorMatrix(matrix.df)
  }
  
  # pre-set breaks
  step <- ((max(matrix.df$y) - min(matrix.df$y)) - ((max(matrix.df$y) - min(matrix.df$y)) %% (break.number-1)))/(break.number-1)
  breaks <- seq(from=min(matrix.df$y), to=max(matrix.df$y), by=step)
  labels <- sapply(breaks, function(x){ FormatBp(x, type=format)}) 
  
  p <- ggplot(data = matrix.df, aes(x=x, y=y, fill=value)) + 
    geom_tile(width=hic$bin.size, height=hic$bin.size) + 
    scale_x_continuous(breaks=breaks, label=labels) +
    scale_y_reverse(breaks=breaks, label=labels) +
    scale_fill_gradientn(colours=colors) + 
    science_theme + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  
  return(p)
  
}


PlotTriangleMatrix <- function(
  hic,
  break.number = 5,
  format = "s",
  pal = "YlOrRd"
){
  # Print a Triangular Matrix over genomic locus plot given a three column (y, x, value) and grafic options
  #
  # Args:
  #   matrix.df: 3 column df c(y, x, value)
  #   break.number: number of breaks to print on axis default=5
  #   format: options to format the gneomic postion (s=single, k=in kilo, M=in Mega) default="s"
  #   pal: RColorBrewer palette to use default="YlOrRd"
  #
  # Returns:
  #   ggplot2 plot
  
  # select color scheme colours
  colors <- brewer.pal(9, pal) 
  
  # get matrix
  matrix.df <- hic$matrix.df

  # match and get chromosomal position from coords id
  matrix.df$y <- hic$coords$center[match(matrix.df$y.id, hic$coords$id)]
  matrix.df$x <- hic$coords$center[match(matrix.df$x.id, hic$coords$id)]
  matrix.df <- matrix.df[,c(1:3)]

  # convert matrix to triangle plot polygon matrix format if not in right format yet
  if(!"pos" %in% colnames(matrix.df)){
    tri <- MakeTriangleMatrix(matrix.df)
  }else{
    tri <- matrix.df
  }
    
  # pre-set breaks (#TODO --> fix to round values rather then performing an actual rounding!!!)
  step <- ((max(matrix.df$x) - min(matrix.df$x)) - ((max(matrix.df$x) - min(matrix.df$x)) %% (break.number-1)))/(break.number-1)
  breaks <- seq(from=min(matrix.df$x), to=max(matrix.df$x), by=step)
  labels <- sapply(breaks, function(x){ FormatBp(x, type=format)}) 
  
  # Plot 
  t <- ggplot(data = tri, aes(x=x, y=y, fill=value, group=id)) + 
      geom_polygon() +
      scale_x_continuous(breaks=breaks, label=labels) +
      coord_cartesian(xlim=c(min(matrix.df$x), max(matrix.df$x)), ylim=c(0, max(tri$y))) +
      scale_fill_gradientn(colours=colors) +
      science_theme + 
      theme(
        legend.position = c(0.925, 0.75),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # strip.ba = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=45, vjust=0.5),
        strip.background = element_blank()
      )
  
  return(t)
    
}

PlotSignalTrack <- function(
  bg.df,
  label="",
  xlim,
  ylim="",
  break.number = 5,
  format = "s",
  color = "black"
){
  
  #convert dataframe for plotting
  bg.df <- MakeBedGraphDataFrameForSignal(bg.df)
  
  # add 2 times zero values flanking for AUC plotting
  bg.df <- rbind(bg.df[1, ], bg.df, bg.df[nrow(bg.df), ])
  bg.df[1, ]$pos <- bg.df[1, ]$pos - 1
  bg.df[1, ]$value <- 0
  bg.df[nrow(bg.df), ]$pos <- bg.df[nrow(bg.df), ]$pos + 1
  bg.df[nrow(bg.df), ]$value <- 0
  
  # pre-set breaks
  step <- ((max(bg.df$pos) - min(bg.df$pos)) - ((max(bg.df$pos) - min(bg.df$pos)) %% (break.number-1)))/(break.number-1)
  breaks <- seq(from=min(bg.df$pos), to=max(bg.df$pos), by=step)
  labels <- sapply(breaks, function(x){ FormatBp(x, type=format)}) 
  
  
  bg.p <- ggplot(bg.df, aes(x=pos, y=value)) + 
    # geom_path(colour=color, size=0.2) +
    geom_polygon(fill=color) + 
    labs(y=label) +
    expand_limits(x=xlim) +
    coord_cartesian(ylim=c(1, max(bg.df$value)), xlim=xlim) +
    scale_x_continuous(breaks=breaks, label=labels) +
    science_theme + 
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), panel.border=element_blank())
  
  # add custom ylim if desired
  if(ylim != ""){ bg.p <- bg.p + coord_cartesian(ylim=ylim, xlim=xlim) }

  return(bg.p)
  
}


# HELPER ===========================================
binInteractionValues <- function(hic, bins=9){
  # Take hic interactions and bin them into [bins] distinc classes equally split
  # will reserve one bin for true 0s
  # Arguments:
  #   bins: default=9, number of bins.quantiles to split the valus into
  bins <- bins -1
  hic$matrix.df$value <- ntile(hic$matrix.df$value[hic$matrix.df$value != 0], bins)
  return(hic)
}

leftHandNotate <- function(matrix){
  # helper Function to change matrix notation to left msot base of the bin notation
  matrix$x <- matrix$x - 500
  matrix$y <- matrix$y - 500
  return(matrix)
}

trimHicRange <- function(hic.obj, range=1000000){
  # Function to trim a hic object to only maintain interactions with a maximum interaction range
  hic.obj$matrix.df <- hic.obj$matrix.df[abs(hic.obj$matrix.df$x - hic.obj$matrix.df$y) <= range,]
  return(hic.obj)
}

trimHicCoords<- function(hic.obj, start=1000000, end=2000000){
  # Function to trim a hic object to only maintain interactions within start and end coordinate
  hic.obj$matrix.df <- hic.obj$matrix.df[hic.obj$matrix.df$x >= start,]
  hic.obj$matrix.df <- hic.obj$matrix.df[hic.obj$matrix.df$x <= end,]
  hic.obj$matrix.df <- hic.obj$matrix.df[hic.obj$matrix.df$y >= start,]
  hic.obj$matrix.df <- hic.obj$matrix.df[hic.obj$matrix.df$y <= end,]
  # hic.obj$matrix.df <- hic.obj$matrix.df[hic.obj$matrix.df$start >= start,]
  # hic.obj$matrix.df <- hic.obj$matrix.df[hic.obj$matrix.df$start <= end,]
  # hic.obj$matrix.df <- hic.obj$matrix.df[hic.obj$matrix.df$end >= start,]
  # hic.obj$matrix.df <- hic.obj$matrix.df[hic.obj$matrix.df$end <= end,]
  
  return(hic.obj)
}


# LEGACY ==========================================
PlotSquareMatrixOld <- function(
  matrix.df,
  break.number = 5,
  square=FALSE,
  format = "s",
  pal = "YlOrRd"
){
  # Wrapper function to plot a matrix given a 3 column dataframe as input
  # Print a Square Matrix plot given a three column (y, x, value) and grafik options
  #
  # Args:
  #   matrix.df: 3 column df c(y, x, value)
  #   break.number: number of breaks to print on axis default=5
  #   square: FALSE/TRUE [default FALSE] indicate if to mirror a trianglualr matrix to get a full matrix plot
  #   format: options to format the gneomic postion (s=single, k=in kilo, M=in Mega) default="s"
  #   pal: RColorBrewer palette to use default="YlOrRd"
  #
  # Returns:
  #   ggplot2 plot
  
  # select color scheme colours
  colors <- brewer.pal(9, pal) 
  
  #mirror data frame if selected
  if(square){
    matrix.df <- MakeMirrorMatrix(matrix.df)
  }
  
  # pre-set breaks
  step <- ((max(matrix.df$y) - min(matrix.df$y)) - ((max(matrix.df$y) - min(matrix.df$y)) %% (break.number-1)))/(break.number-1)
  breaks <- seq(from=min(matrix.df$y), to=max(matrix.df$y), by=step)
  labels <- sapply(breaks, function(x){ FormatBp(x, type=format)}) 
  
  p <- ggplot(data = matrix.df, aes(x=x, y=y, fill=value)) + 
    geom_raster() + 
    scale_x_continuous(breaks=breaks, label=labels) +
    scale_y_reverse(breaks=breaks, label=labels) +
    scale_fill_gradientn(colours=colors) + 
    science_theme + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  
  return(p)
  
}

PruneHicproMatrixOld <- function(hic, chr, start, end){
  # Function to prune a matrix to zoom into a region of interest
  # Take a list object as imported with ImportHicproMatrix
  # and prune it to include only interactions between the genomix coordinates chosen
  
  #get bins in ROI
  hic$coords <- hic$coords[((hic$coords$chr == chr) & (hic$coords$center >= start) & (hic$coords$center <= end)),]
  
  #get interactions in ROI
  hic$matrix.df <- hic$matrix.df[((hic$matrix.df$y %in% hic$coords$center) & (hic$matrix.df $x %in% hic$coords$center)),]
  
  #adjust start.id and start.pos (for plotting)
  hic$start.id <- min(hic$coords$id)
  hic$start.pos <- min(hic$matrix.df$y)
  
  return(hic)
  
}

PlotTriangleMatrixOld <- function(
  matrix.df,
  break.number = 5,
  format = "s",
  pal = "YlOrRd"
){
  # Print a Triangular Matrix over genomic locus plot given a three column (y, x, value) and grafic options
  #
  # Args:
  #   matrix.df: 3 column df c(y, x, value)
  #   break.number: number of breaks to print on axis default=5
  #   format: options to format the gneomic postion (s=single, k=in kilo, M=in Mega) default="s"
  #   pal: RColorBrewer palette to use default="YlOrRd"
  #
  # Returns:
  #   ggplot2 plot
  
  # select color scheme colours
  colors <- brewer.pal(9, pal) 
  
  # get matrix
  matrix.df <- hic$matrix.df
  
  # match and get chromosomal position from coords id
  matrix.df$y <- hic$coords$center[match(matrix.df$y, hic$coords$id)]
  matrix.df$x <- hic$coords$center[match(matrix.df$x, hic$coords$id)]
  
  # convert matrix to triangle plot polygon matrix format if not in right format yet
  if(!"pos" %in% colnames(matrix.df)){
    tri <- MakeTriangleMatrix(matrix.df)
  }else{
    tri <- matrix.df
  }
  
  # pre-set breaks (#TODO --> fix to round values rather then performing an actual rounding!!!)
  step <- ((max(matrix.df$x) - min(matrix.df$x)) - ((max(matrix.df$x) - min(matrix.df$x)) %% (break.number-1)))/(break.number-1)
  breaks <- seq(from=min(matrix.df$x), to=max(matrix.df$x), by=step)
  labels <- sapply(breaks, function(x){ FormatBp(x, type=format)}) 
  
  # Plot 
  t <- ggplot(data = tri, aes(x=x, y=y, fill=value, group=id)) + 
    geom_polygon() +
    scale_x_continuous(breaks=breaks, label=labels) +
    coord_cartesian(xlim=c(min(matrix.df$x), max(matrix.df$x)), ylim=c(0, max(tri$y))) +
    scale_fill_gradientn(colours=colors) +
    science_theme + 
    theme(
      legend.position = c(0.925, 0.75),
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # strip.ba = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle=45, vjust=0.5),
      strip.background = element_blank()
    )
  
  return(t)
  
}
