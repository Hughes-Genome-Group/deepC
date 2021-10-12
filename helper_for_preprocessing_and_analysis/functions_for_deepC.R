# FUNCTIONS FOR DEEPC DATA ANALYSIS AND VISUALISATION
# Ron Schwessinger
# Requires tidyverse and RColorBrewer

require(tidyverse)
require(RColorBrewer)

# THEMES ---------------------------------------------
# themes for plotting
upper_overlay_theme <- theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                             axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             legend.margin = margin(t= 1, l = 1, b = 1, r = 1, unit = "pt"))

lower_overlay_theme <- upper_theme <- theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                                            axis.title.x = element_blank(),
                                            legend.margin = margin(t= 1, l = 1, b = 1, r = 1, unit = "pt"))

# Processing Helper =======================================
#' deduplicate
#'
#' Remove duplicated rows from a tibble/data frame. "Unique" on tibbles works better and faster
#' @param df data frame / tibble
#' @return data frame with duplicate rows removed.
#' @examples
#' \dontrun{
#' deduplicate(df)
#' }
#' @export
deduplicate <- function(df){
  df <- df[!duplicated(df),]
  return(df)
}

# HELPER FUNCTION read shape variant deploy file

#' readDeepcVariantFile
#'
#' Read in a deepC variant file from deploy net output. Assumes a deepC header in the input file.
#' @param file Path to variant output file. Tab separated ("chr", "start", "end", predictions ...)
#' @param prediction.bins Number of bins/values in the output vector, depends on bin size and context window. Default = 101
#' @param bin.size Bin size in bp. default = 10000
#' @param gather TRUE/FALSE if to gather the data in long format for tidyverse/ggplot default = TRUE.
#' @param tag.col TRUE/FALSE if a tagging or id column is present ("chr", "start", "end", "tag", predictions ...) default = FALSE
#' @param zigzag TRUE/FALSE if to correc thte positions for zigzag pole encoding (shifting every second bin position) default = TRUE
#' @return Dataframe/tibble of predicted interactions/
#' @examples
#' \dontrun{
#' readDeepcVariantFile(file, prediction.bins = 201, bin.size = 5000, gather = TRUE, tag.col = FALSE, zigzag = TRUE)
#' }
#' @export
readDeepcVariantFile <- function(file, prediction.bins = 101, bin.size = 10000, gather=TRUE, tag.col = FALSE, zigzag = TRUE){

  # 1) Read Header
  con <- file(file,"r")
  header <- readLines(con,n=3)
  header.1 <- header[1]
  header.2 <- header[2]
  header.3 <- header[3]
  close(con)
  # 2) extract relative mapping positions
  header.1 <- strsplit(header.1, "\\s+")[[1]]
  header.2 <- strsplit(header.2, "\\s+")[[1]]
  header.3 <- strsplit(header.3, "\\s+")[[1]]
  variant.chrom <- header.1[4]
  variant.start <- as.numeric(header.1[5])
  variant.end <- as.numeric(header.1[6])
  relative.chrom <- header.2[7]
  relative.start <- as.numeric(header.2[8])
  relative.end <- as.numeric(header.2[9])
  bp.adjust <- as.numeric(header.3[7])
  # 3) load variant dataframe
  df <- as_tibble(read.table(file, header = F, skip = 3))
  if(tag.col == TRUE){
    names(df) <- c("chr", "start", "end", "tag", c(1:prediction.bins))
    df <- df %>% select(-tag)
  }else{
    names(df) <- c("chr", "start", "end", c(1:prediction.bins))
  }

  # 4) gather data frame over bins if specified
  if(gather == TRUE){
    df <- df %>%
      mutate(pos = start + (end - start)/2) %>%
      select(chr, pos, c(4:(prediction.bins+3))) %>%
      gather(bin, value, -chr, -pos) %>%
      mutate(bin = as.numeric(bin)) %>%
      filter(chr == relative.chrom)

    # correct position for zigzag binning
    if( zigzag == TRUE){
      df <- df %>% mutate(pos = if_else(bin %% 2 != 0, pos + bin.size/2, pos))
    }

  }

  newlist <- list(df = df,
                  variant.chrom = variant.chrom,
                  variant.start = variant.start,
                  variant.end = variant.end,
                  relative.chrom = relative.chrom,
                  relative.start = relative.start,
                  relative.end = relative.end,
                  bp.adjust = bp.adjust)
  return(newlist)
}

#' readDeepcVariantFileNoHeader
#'
#' Read in a deepC variant file from deploy net output. Assumes a no header present in the input file.
#' @param file Path to variant output file. Tab separated ("chr", "start", "end", predictions ...)
#' @param prediction.bins Number of bins/values in the output vector, depends on bin size and context window. Default = 101
#' @param bin.size Bin size in bp. default = 10000
#' @param gather TRUE/FALSE if to gather the data in long format for tidyverse/ggplot default = TRUE.
#' @param tag.col TRUE/FALSE if a tagging or id column is present ("chr", "start", "end", "tag", predictions ...) default = FALSE
#' @param zigzag TRUE/FALSE if to correc thte positions for zigzag pole encoding (shifting every second bin position) default = TRUE
#' @return Dataframe/tibble of predicted interactions/
#' @examples
#' \dontrun{
#' readDeepcVariantFileNoHeader(file, prediction.bins = 201, bin.size = 5000, gather = TRUE, tag.col = FALSE, zigzag = TRUE)
#' }
#' @export
readDeepcVariantFileNoHeader <- function(file, prediction.bins = 101, bin.size = 10000, gather=TRUE, tag.col = FALSE, zigzag = TRUE){

  # 1) load variant dataframe
  df <- as_tibble(read.table(file, header = F))
  if(tag.col == TRUE){
    names(df) <- c("chr", "start", "end", "tag", c(1:prediction.bins))
    df <- df %>% select(-tag)
  }else{
    names(df) <- c("chr", "start", "end", c(1:prediction.bins))
  }

  # 4) gather data frame over bins if specified
  if(gather == TRUE){
    df <- df %>%
      mutate(pos = start + (end - start)/2) %>%
      select(chr, pos, c(4:(prediction.bins+3))) %>%
      gather(bin, value, -chr, -pos) %>%
      mutate(bin = as.numeric(bin))

    # correct position for zigzag binning
    if( zigzag == TRUE){
      df <- df %>% mutate(pos = if_else(bin %% 2 != 0, pos + bin.size/2, pos))
    }

  }

  return(df)
}


#' readHicDataBins
#'
#' Read in a hic data file with interactions in bins in a comma separated string in the 4th column..use when start and end columns are genomic positions/coordinates.
#' @param file Path to variant output file. Tab separated ("chr", "start", "end", "qbins"). Where qbins are hic interactions in comma separated string.
#' @param prediction.bins Number of bins/values in the output vector, depends on bin size and context window. Default = 101
#' @param bin.size Bin size in bp. default = 10000
#' @param gather TRUE/FALSE if to gather the data in long format for tidyverse/ggplot default = TRUE.#' @param zigzag TRUE/FALSE if to correc thte positions for zigzag pole encoding (shifting every second bin position) default = TRUE
#' @return Dataframe/tibble of quantized hic interactions.
#' @examples
#' \dontrun{
#' readHicDataBins(file, prediction.bins = 50, gather = TRUE)
#' }
#' @export
readHicDataBins <- function(file, prediction.bins = 101, gather=TRUE){

  sdf <- as_tibble(read.table(file))
  sdf$V4 <- as.character(sdf$V4)
  names(sdf) <- c("chr", "start", "end", "qbins")

  df <- sdf[,c(1:3)]
  names(df) <- c("chr", "start", "end")
  df <- df %>%
    mutate(pos = (end - start)/2 + start) %>%
    select(c(chr, pos))

  # split and fill classes into numeric matrix
  mat <- matrix(data=0, nrow = dim(sdf)[1], ncol = prediction.bins)
  b <- as.character(sdf$qbins)
  for(i in c(1:length(b))){
    s <- as.numeric(strsplit(b[i], ",")[[1]])
    mat[i,] <- s
  }

  # combine and gather for plotting
  sdf <- as_tibble(cbind(df, mat))

  #gather if specified
  if(gather == TRUE){
    sdf <- sdf %>%
      gather(bin, value, -chr, -pos) %>%
      mutate(bin = as.numeric(bin))
  }

  return(sdf)
}


#' readDeepcInputHicBed
#'
#' Read in a hic data file with interactions in bins in a comma separated string in the 4th column. Use when coordinates are bin ids rather then genomic poistions/
#' @param file HiC/DeepC file tab separated ("chr", "start", "end", "qbins"). Where qbins are hic interactions in comma separated string.
#' @param prediction.bins Number of bins/values in the output vector, depends on bin size and context window. Default = 101
#' @param bin.size Bin size in bp. default = 10000
#' @param gather TRUE/FALSE if to gather the data in long format for tidyverse/ggplot default = TRUE.#' @param zigzag TRUE/FALSE if to correc thte positions for zigzag pole encoding (shifting every second bin position) default = TRUE
#' @return Dataframe/tibble of quantized interactions.
#' @examples
#' \dontrun{
#' readDeepcInputHicBed(df, prediction.bins = 201, bin.size = 5000, gather = TRUE, zigzag = TRUE)
#' }
#' @export
readDeepcInputHicBed <- function(file, prediction.bins = 101, bin.size = 10000, gather=TRUE, zigzag = TRUE){

  # 1) load variant dataframe
  sdf <- as_tibble(read.table(file, header = F))
  names(sdf) <- c("chr", "start", "end", "qbins")
  sdf$V4 <- as.character(sdf$qbins)

  #  subset to data frame and matrix add pos columnadd position
  df <- sdf[,c(1:3)]
  names(df) <- c("chr", "start", "end")
  if(gather == TRUE){
  df <- df %>%
    mutate(pos = (end - start)/2 + start) %>%
    select(c(chr, pos))
  }
  # split and fill classes into numeric matrix
  mat <- matrix(data=0, nrow = dim(sdf)[1], ncol = prediction.bins)
  b <- as.character(sdf$qbins)
  for(i in c(1:length(b))){
    s <- as.numeric(strsplit(b[i], ",")[[1]])
    mat[i,] <- s
  }

  # combine and gather for plotting
  sdf <- as_tibble(cbind(df, mat))
  if(gather == TRUE){
    sdf <- sdf %>%
      gather(bin, value, -chr, -pos) %>%
      mutate(bin = as.numeric(bin))
    
    # correct position for zigzag binning
    if( zigzag == TRUE){
      sdf <- sdf %>% mutate(pos = if_else(bin %% 2 == 0, pos - bin.size/2, pos))
    }
    
  }



  return(sdf)
}



#' plotDiffDeepC
#'
#' Plot difference between two hic/deepc data frames (df1 - df2). e.g reference - variant
#' @param df1 Interaction data frame (reference)
#' @param df2 Interaction data frame (variant)
#' @param bin.size Bin size in bp. default = 10000
#' @param threshold Set every interaction below "threshold" to zero
#' @return ggplot2 plot of (predicted) interaction difference
#' @examples
#' \dontrun{
#' plotDiffDeepC(df1, df2, bin.size = 5000, threshold = 2)
#' }
#' @export
plotDiffDeepC <- function(df1, df2, bin.size = 10000, threshold = 0){

  sub.df <- df1 %>% filter(pos %in% df2$pos)
  sub.df$var.value <- pull(df2 %>% filter(pos %in% df1$pos) %>% select(value))
  sub.df <- triangularize(sub.df, bin.size, extra = c("value", "var.value"))
  diffplot <- sub.df %>%
    mutate(diff = var.value - value) %>%
    mutate(diff = if_else(abs(diff) >= threshold, diff, 0)) %>%
    ggplot(aes(x = pos, y = (bin*bin.size)/1000,, fill = diff, group = polyid)) + geom_polygon() +
    scale_fill_gradientn(colours = brewer.pal(9, 'Spectral')) +
    labs(y = "Genomic Distance [kb]")

  return(diffplot)

}

#' getDiffDeepC
#'
#' Calcualte difference between two hic/deepc data frames (df1 - df2). e.g reference - variant
#' @param df1 Interaction data frame (reference)
#' @param df2 Interaction data frame (variant)
#' @return data frame of (predicted) interaction difference
#' @examples
#' \dontrun{
#' getDiffDeepC(df1, df2)
#' }
#' @export
getDiffDeepC <- function(df1, df2){
  sub.df <- df1 %>% filter(pos %in% df2$pos)
  sub.df$var.value <- pull(df2 %>% filter(pos %in% sub.df$pos) %>% select(value))
  sub.df <- sub.df %>% mutate(diff = var.value - value)
  return(sub.df)
}

#' getDistanceWindow
#'
#' Helper fuction to pick a distance window with the same distance from center position from a hic/deepC data frame
#' @param d hic data frame / tibble
#' @param center center position of region of interest
#' @param distance distance from center to extract
#' @return Hic/deepc data frame from region  of interest
#' @examples
#' \dontrun{
#' getDistanceWindow(d, center = 1000000, distance = 500000)
#' }
#' @export
getDistanceWindow <- function(d, center, distance){
  #
  a <- d %>%
    filter(abs(pos - center) <= distance & bin <= (distance - (abs(pos - center)))/(bin.size/2))
  return(a)
}


#' calcInsulationScoreChange
#'
#' Calculate changes in the interaction score that pass the diagonals starting from
#' the variant bin ---> insulation score changes. The "ins" snsulation score change.
#' Also calculate changes in the score that do not pass the diagonal but
#' are pot. affected by the centre of a domain. The "cont" containment score change.
#' @param diff.df differential data frame as from getDiffDeepC.
#' @param bin.size Bin size in bp. default = 10000
#' @return List of length two with the .$ins change and the .$cont change
#' @examples
#' \dontrun{
#' calcInsulationScoreChange(diff.df, bin.size = 5000)
#' }
#' @export
calcInsulationScoreChange <- function(diff.df, bin.size = 10000){
  center <- custom_floor(median(diff$pos), bin.size)
  diff.df <- diff.df %>%
    mutate(left.passing = if_else(pos <= center & (bin + 1) > (center - pos)/(bin.size/2), 1, 0)) %>%
    mutate(right.passing = if_else(pos > center & (bin + 1) > (pos - center)/(bin.size/2), 1, 0)) %>%
    mutate(passing = case_when(left.passing == 1 ~ 1, right.passing == 1 ~ 1, TRUE ~ 0)) %>%
    mutate(contained = 1-passing) %>%
    select(-left.passing, -right.passing)

  insulation.score.change <- sum(diff.df$value * diff.df$passing)
  containment.score.change <- sum(diff.df$value * (diff.df$contained))
  # normalize
  insulation.score.change <- insulation.score.change / sum(diff.df$passing)
  containment.score.change <- containment.score.change / sum(diff.df$contained)

  newlist <- list('ins' = insulation.score.change, 'cont' = containment.score.change)

  return(newlist)
}


#' triangularize
#'
#' Triangularize a hic/deepc dataframe for polygon style plotting/
#' Take a df with distance bins (deepCregr ZigZag)
#' convert it to triangular like polygon df over genomic locus
#' Specify which extra cloumn names to extract.
#' @param df Hic/deepc data frame with interactions/predictions in zig zag encoding via bins and genomic position.
#' @param bin Bin size in bp. default = 10000
#' @param extra Character vector specifying which columns to extract interactions etc. from. Default c("value"). Optional length
#' @return Data frame specifying polygon coordinates for HiC ploting.
#' @examples
#' \dontrun{
#' triangularize(df, bin.size = 5000, extra = c("value"))
#' }
#' @export
triangularize <- function(df, bin = 10000, extra = c("value")){
  # Take a df with distance bins (deepCregr ZigZag/Diamond)
  # convert it to triangular like polygon df over genomic locus
  # extra cloumn names to extract

  half_bin <- bin/2

  # 1) initialise 4 matrices by copying from half square df & add polygon pos and bin id column
  tri <- df
  tri$cy <- tri$cx <- rep(0, nrow(tri))  #init empty coord corrected columns

  tri1 <- tri2 <- tri3 <- tri4 <- tri #quadruple for polygon positions

  # tri1$polypos <- rep(1, nrow(tri1))  #indicate polygon positions
  # tri2$polypos <- rep(2, nrow(tri2))
  # tri3$polypos <- rep(3, nrow(tri3))
  # tri4$polypos <- rep(4, nrow(tri4))
  tri1$polyid <- tri2$polyid <- tri3$polyid <- tri4$polyid <- c(1:nrow(tri1)) # add polygon group id

  # 2) add adjusted coordinate columns
  tri1$cy <- tri$bin - 1
  tri2$cy <- tri$bin
  tri3$cy <- tri$bin + 1
  tri4$cy <- tri$bin

  tri1$cx <- tri1$pos
  tri2$cx <- tri1$pos - half_bin
  tri3$cx <- tri1$pos
  tri4$cx <- tri1$pos + half_bin

  # 3) rbind and replace y and x coord columns
  tri <- rbind(tri1, tri2, tri3, tri4)
  tri$pos <- tri$cx
  tri$bin <- tri$cy
  # need to make that more generalizable ....!
  if("start" %in% names(tri)){
    tri <- tri[,c("chr", "start", "end", "pos", "bin", "polyid", extra)]
  }else{
    tri <- tri[,c("chr", "pos", "bin", "polyid", extra)]
  }

  # tri <- tri[order(tri$polyid, tri$key, decreasing = FALSE ), ] # sort

  tri$bin[tri$bin < 0] <- 0 #  for negative y

  return(tri)

}

# Virtual 4C --------------------------------------------------------------
#' custom_round
#'
#' Round to a custom base (nearest base)
#' round(x/base)*base
#' @param x Input number to round
#' @param base Base to round to (base=1 for ordinary rounding)
#' @param extra Character vector specifying which columns to extract interactions etc. from. Default c("value"). Optional length
#' @return Number x rounded to nearest base.
#' @examples
#' custom_round(1503, 500)
#' @export
custom_round <- function(x, base) {
  round(x/base)*base
}


#' custom_floor
#'
#' Floor round to a custom base (nearest base)
#' @param x Input number to round
#' @param base Base to round to (base=1 for ordinary floor rounding)
#' @param extra Character vector specifying which columns to extract interactions etc. from. Default c("value"). Optional length
#' @return Number x rounded to nearest base floor.
#' @examples
#' custom_floor(1503, 500)
#' @export
custom_floor <- function(x, base) {
  r <- round(x/base)*base
  if(r > x){
    r <- r - base
  }
  return(r)
}


#' custom_ceil
#'
#' Ceiling round to a custom base (nearest base)
#' @param x Input number to round
#' @param base Base to round to (base=1 for ordinary ceiling rounding)
#' @param extra Character vector specifying which columns to extract interactions etc. from. Default c("value"). Optional length
#' @return Number x rounded to nearest base ceiling.
#' @examples
#' custom_ceil(1503, 500)
#' @export
custom_ceil <- function(x, base) {
  r <- round(x/base)*base
  if(r < x){
    r <- r + base
  }
  return(r)
}


#' vec_custom_round
#'
#' Vectorised custom_round
#' Round to a custom base (nearest base)
#' round(x/base)*base
#' @param x Input number to round
#' @param base Base to round to (base=1 for ordinary rounding)
#' @param extra Character vector specifying which columns to extract interactions etc. from. Default c("value"). Optional length
#' @return Vector x rounded to nearest base.
#' @examples
#' vec_custom_round(c(2505, 1503, 499), 500)
#' @export
vec_custom_round <- Vectorize(custom_round)


#' vec_custom_floor
#'
#' Vectorised custom_floor
#' Floor round to a custom base (nearest base)
#' @param x Input number to round
#' @param base Base to round to (base=1 for ordinary floor rounding)
#' @param extra Character vector specifying which columns to extract interactions etc. from. Default c("value"). Optional length
#' @return Vector x rounded to nearest base floor
#' @examples
#' vec_custom_floor(c(2505, 1503, 499), 500)
#' @export
vec_custom_floor <- Vectorize(custom_floor)


#' vec_custom_ceil
#'
#' Vectorised custom_ceil
#' CeilingFloor round to a custom base (nearest base)
#' @param x Input number to round
#' @param base Base to round to (base=1 for ordinary floor rounding)
#' @param extra Character vector specifying which columns to extract interactions etc. from. Default c("value"). Optional length
#' @return Vector x rounded to nearest base ceiling
#' @examples
#' vec_custom_ceil(c(2505, 1503, 499), 500)
#' @export
vec_custom_ceil <- Vectorize(custom_ceil)


#' triangularize
#'
#' Virtual 4C current implementation.
#' Map to floor and ceiling of nearby bins
#' For downstream interactions map between pos and - bin.size
#' For upstream interactions map between pos and + bin.size
#' @param idf input dataframe in deepC/hic format (bin, pos)
#' @param chr chromsome
#' @param pos genomic positon  of virtual 4C "viewpoint"
#' @param bin.size Bin size in bp default = 10000
#' @param window.size Window in which to calculate the virtual 4C. Distance from viewpoint. Set pred.start and pred.window to 0 if using this mode (default)
#' @param pred.start Alternatively. leave window.size = 0 and set pred.start as start of region of interested.
#' @param pred.start Alternatively. leave window.size = 0 and set pred.end as end of region of interested.
#' @return Data frame specifying polygon coordinates for HiC ploting.
#' @examples
#' \dontrun{
#' virtual4C(idf, 'chr16', 1000000, bin.size = 5000, window.size = 2000000)
#' }
#' @export
virtual4C <- function(idf, chr, pos, bin.size, window.size = 0, pred.start = 0, pred.end = 0){

  temp.window.size <- window.size - bin.size

  downstream.pred.pos <- custom_ceil(pos, bin.size) + bin.size/2
  upstream.pred.pos <- custom_floor(pos, bin.size) + bin.size/2

  if(pred.start == 0){

    down.df <- tibble(chr = chr,
                      pos = seq(downstream.pred.pos, downstream.pred.pos + (temp.window.size/2) - bin.size, bin.size),
                      value = 0)
    up.df <- tibble(chr = chr,
                    pos = rev(seq(upstream.pred.pos - (temp.window.size/2), upstream.pred.pos, bin.size)),
                    value = 0)

  }else if(window.size == 0){
    pred.start <- custom_ceil(pred.start, bin.size) + bin.size/2
    pred.end <- custom_floor(pred.end, bin.size) + bin.size/2

    down.df <- tibble(chr = chr,
                      pos = seq(downstream.pred.pos, pred.end, bin.size),
                      value = 0)
    up.df <- tibble(chr = chr,
                    pos = rev(seq(pred.start, upstream.pred.pos, bin.size)),
                    value = 0)
  }else{
    print("Specify window around position of interest or a pred.start and pred.end point!")
    return(NA)
  }


  # fill data frames
  for(i in c(1:nrow(down.df))){

    j <- i + 1
    temp.pos <- down.df$pos[i]
    # adjust query position to find value along diagonal
    query.pos <- temp.pos
    # add half bin if mapping to even bin distance
    corr.factor <- (j - (j %% 2))/2 - 1
    query.pos <- query.pos - corr.factor * bin.size

    if(j %% 2 == 0){
      query.pos <- query.pos - bin.size/2
    }

    to.put <- idf %>% filter(pos == query.pos & bin == j)
    if(nrow(to.put) > 0){
      down.df[i,"value"] <- pull(to.put %>% select(value))
    }else{
      down.df[i,"value"] <- 0
    }
  }

  for(i in c(1:nrow(up.df))){

    if(temp.pos >= 0){
      temp.pos <- up.df$pos[i]
      # add half bin if mapping to even bin distance
      query.pos <- temp.pos
      corr.factor <- (i - (i %% 2))/2
      query.pos <- query.pos + corr.factor * bin.size

      if(i %% 2 == 0){
        query.pos <- query.pos + bin.size/2
      }

      to.put <- idf %>% filter(pos == query.pos & bin == i)
      if(nrow(to.put) > 0){
        up.df[i,"value"] <- pull(to.put %>% select(value))
      }else{
        up.df[i,"value"] <- 0
      }
    }
  }


  # combine and remove negative positions (or above chrom_sizes)
  v4c.df <- rbind(down.df %>% mutate(direction = "ds"), up.df %>% mutate(direction = "us"))
  v4c.df <- v4c.df %>%
    filter(pos >= 0) %>%
    arrange(pos)

  return(v4c.df)
}


# HiC Processing ---------------------------------------------------

#' trimHicRange
#'
#' (duplicate from HiC helper scripts)
#' Function to trim a hic object to only maintain interactions with a maximum interaction range
#' @param hic.obj hic object (list) as imported.
#' @param range maximum linear distance range to maintain interactions
#' @return hic object with interactions trimmed for linear disance
#' @examples
#' \dontrun{
#' trimHicRange(hic, range=1000000)
#' }
#' @export
trimHicRange <- function(hic.obj, range=1000000){
  hic.obj$matrix.df <- hic.obj$matrix.df[abs(hic.obj$matrix.df$x - hic.obj$matrix.df$y) <= range,]
  return(hic.obj)
}


#' getBinnedChrom
#'
#' helper function to get a binned genome using bedtools 
#' Create a binned genome data frame on the specified chromosome from start to end positions.
#' @param chr chromosome
#' @param start start position
#' @param end end postion
#' @param window size of the window / bins
#' @param step The step size between windows/bins. The increment.
#' @param timestamp automatically produces a timestamp to uniquely identify temporary files
#' defaut = paste0(gsub("\\s+", "_", Sys.time()), "_", runif(n=1))
#' @return data frame with genomic bins
#' @examples
#' \dontrun{
#' getBinnedChrom(chr='chr1', start=0, end=10000, window=1000, step=500)
#' }
#' @export
getBinnedChrom <- function(chr='chr1', start=0, end=10000, window=1000, step=500, timestamp = paste0(gsub("\\s+", "_", Sys.time()), "_", runif(n=1))){
  # HELPER FUCNTION: get binned genome 
  
  # calc effective end point
  effective.end <- start + window * (floor((end - start)/window) + 1)
  #start with end points
  bg <- tibble(end = seq(start.pos + window, effective.end, step))
  bg <- bg %>% mutate(start = end - window, chr = chr) %>% select(chr, start, end)
  
  return(bg)
}



#' getBinSubMatrixfunction
#'
#' Get A pruned submatrix from hic object given a region of interest as well as an interaction data frame.
#'
#' @param chr chromosome
#' @param start start position
#' @param end end postion
#' @param step The step size between windows/bins. The increment.
#' @param bin.size Bin size in bp default = 10000
#' defaut = paste0(gsub("\\s+", "_", Sys.time()), "_", runif(n=1))
#' @return List with two elements "hicp" the pruned Hic object and "idf" the interaction data frame.
#' @examples
#' \dontrun{
#' getBinSubMatrixfunction(chr='chr1', start=0, end=10000, bin.size = 5000)
#' }
#' @export
getBinSubMatrixfunction <- function(hic.obj, chr='chr1', start=0, end=20000, bin.size=10000){

  chromosome <- chr
  coord.start <- start
  coord.end <- end
  coord.center <- (coord.end - coord.start)/2 + coord.start

  # get pruned matrix
  hic.pruned <- PruneHicproMatrix(hic.obj, chr=chromosome, start=coord.start, end=coord.end)$matrix
  hic.pruned <- leftHandNotate(hic.pruned)
  hic.pruned <- hic.pruned[order(hic.pruned$y, hic.pruned$x),]

  # create empty dataframe to assemble the interaction values
  idf <- data.frame(
    x=seq(coord.start, coord.center-bin.size, bin.size),
    y=rev(seq(coord.center, coord.end - bin.size, bin.size)))
  idf$chr <- chromosome
  idf <- idf[,c(3,1,2)]
  # fetch values where available
  idf$value <- apply(idf, 1, function(x){
    l <- hic.pruned[(hic.pruned$y == x[2] & hic.pruned$x == x[3]),]
    if(nrow(l) == 0){
      v <- 0
    }else{
      v <- l$value
    }
    return(v)
  })

  return(list(hicp=hic.pruned, idf=idf))
}



#' getVerticalWindowInteractions
#'
#' Extract desired vertical center pole of interactions over the center of every window of window.size.
#' Older non zig zag encoding.
#'
#' @param hic.obj hic object including bin coordinates and an interaction matrix
#' @param bin.df Binned genome / genomic segment data frame as from getBinnedChrom
#' @param window window size in bp of interest usually 1 Mb + 1x bin size
#' @param bin Bin size in bp
#' defaut = paste0(gsub("\\s+", "_", Sys.time()), "_", runif(n=1))
#' @return Dataframe/tibble of interactions encoded in a vertical pole over the center of the window.
#' @examples
#' \dontrun{
#' getVerticalWindowInteractions(hic, binned.genome, window=1e6+5000, bin = 5000)
#' }
#' @export
getVerticalWindowInteractions <- function(hic.obj, bin.df, window, bin){

  expected.values <- window / (2 * bin)

  cdf <- pbapply(bin.df, 1, function(x){
    t <- getBinSubMatrixfunction(hic.obj, as.character(x[1]), as.numeric(x[2]), as.numeric(x[3]), window, bin)
    t$idf <- t$idf[c(nrow(t$idf):1),]  # invert row wise
    d <- data.frame(
      chr=as.character(x[1]),
      start=as.numeric(x[2]),
      end=as.numeric(x[3])
    )
    d[,c((4):(expected.values+3))] <- t$idf$value[c(1:expected.values)]
    return(d)
  })

  cdf <- do.call("rbind", cdf)
  names(cdf)[c(4:(expected.values+3))] <- c(1:expected.values)

  return(as_tibble(cdf))

}


#' getZigZagWindowInteractionsPerl
#'
#' Extract desired vertical zigzag pole of interactions over the center of every window of window.size.
#' Uses helper perl script for faster processing. First creates a query data frame, then saves
#' the dquery and the hic data frame in temporary files in the working directory. Then uses the helper perl script
#' to query and match and read the output file into R again.
#' REQUIRES: perl installed and available  in your default shell. #'
#' @param hic.obj hic object including bin coordinates and an interaction matrix
#' @param bin.df Binned genome / genomic segment data frame as from getBinnedChrom
#' @param window window size in bp of interest usually 1 Mb + 1x bin size
#' @param bin Bin size in bp
#' @param query.pl path to helper perl script as distributed with the deepC repository "match_query_tabel.pl"
#' @param timestamp automatically produces a timestamp to uniquely identify temporary files
#' defaut = paste0(gsub("\\s+", "_", Sys.time()), "_", runif(n=1))
#' @return Dataframe/tibble of interactions encoded in a vertical pole over the center of the window.
#' @examples
#' \dontrun{
#' getZigZagWindowInteractionsPerl(hic, binned.genome, window=1e6+5000, bin = 5000, query.pl = '~/myscripts/deepC/match_query_table.pl')
#' }
#' @export
getZigZagWindowInteractionsPerl <- function(hic.obj,
                                            bin.df,
                                            window,
                                            bin,
                                            query.pl = "~/fusessh/scripts/machine_learning/epigenome_nets/deepC/match_query_table.pl",
                                            timestamp = paste0(gsub("\\s+", "_", Sys.time()), "_", runif(n=1))){

  temp.hic.file <- paste0("temp_hic_file_", timestamp, ".txt")
  temp.query.table <- paste0("temp_query_table_", timestamp, ".txt")
  temp.query.outfile <- paste0("temp_query_outfile_", timestamp, ".txt")

  write.table(hic.obj$matrix.df, file = temp.hic.file, row.names = F, col.names = F, sep="\t", quote=F)
  print("saved temp_hic_file")

  # Make Query Data frame
  expected.values <- window/bin
  q <- apply(bin.df, 1, function(x){
    # get coords
    chromosome <- as.character(x[1])
    coord.start <- as.numeric(x[2])
    coord.end <- as.numeric(x[3])
    coord.center <- (coord.end - coord.start)/2 + coord.start

    # lay out queries
    x <- rep(seq(coord.start, coord.center - bin/2, bin), each = 2)
    x <- x[-length(x)]  # trim last pos

    y <- rep(seq(coord.center - bin/2, coord.end - bin, bin), each = 2)
    y <- y[-length(y)]  # trim last pos
    y <- rev(y)

    idf <- data.frame(x=x, y=y)

    rm(x, y)

    idf$tag <- paste0(idf$x, ':', idf$y)
    idf <- idf[c(nrow(idf):1),]  # invert row wise
    # make new dataframe
    d <- data.frame(
      chr = chromosome,
      start = coord.start,
      end = coord.end)
    # add tags
    d[,c(4:(expected.values+3))] <- idf$tag

    return(d)
  })

  qdf <- do.call("rbind", q)
  # return(qdf)

  # save query table as temp file
  write.table(qdf, file = temp.query.table, col.names = F, row.names = F, quote = F, sep="\t")
  print("saved temp_query_file")

  # match query using perl script
  command <- paste0("perl ", query.pl, " ", temp.hic.file, " ", temp.query.table, " ", temp.query.outfile, " ", bin)
  print("matching ...")
  system(command)
  # read back in match queried data tale as data frame
  mdf <- as.data.frame(read.table(temp.query.outfile, header=F, colClasses = c("character", rep("numeric",2+expected.values))))
  # add names
  names(mdf) <- c("chr", "start", "end", c(1:expected.values))

  # clean temp files
  command <- paste0("rm -f ", temp.hic.file, " ", temp.query.table, " ", temp.query.outfile)
  system(command)

  return(mdf)
}

#' getZigZagWindowInteractionsPerlMemoryFriendly
#'
#' Extract desired vertical zigzag pole of interactions over the center of every window of window.size.
#' Uses helper perl script for faster processing. First creates a query data frame, then saves
#' the dquery and the hic data frame in temporary files in the working directory. Then uses the helper perl script
#' to query and match and read the output file into R again.
#' REQUIRES: perl installed and available  in your default shell.
#' Memory friendly version that uses a second perl helper script for query preparation. Slower then default version.
#' @param hic.obj hic object including bin coordinates and an interaction matrix
#' @param bin.df Binned genome / genomic segment data frame as from getBinnedChrom
#' @param window window size in bp of interest usually 1 Mb + 1x bin size
#' @param bin Bin size in bp
#' @param query.pl path to helper perl script as distributed with the deepC repository "match_query_tabel.pl"
#' @param prepare.pl path to helper perl script as distributed with the deepC repository "prepare_query_tabel.pl"
#' @param timestamp automatically produces a timestamp to uniquely identify temporary files
#' defaut = paste0(gsub("\\s+", "_", Sys.time()), "_", runif(n=1))
#' @return Dataframe/tibble of interactions encoded in a vertical pole over the center of the window.
#' @examples
#' \dontrun{
#' getZigZagWindowInteractionsPerlMemoryFriendly(hic,
#' binned.genome,
#' window=1e6+5000,
#' bin = 5000,
#' query.pl = '~/myscripts/deepC/match_query_table.pl',
#' prepare.pl = '~/myscripts/deepC/prepare_query_tabel.pl')
#' }
#' @export
getZigZagWindowInteractionsPerlMemoryFriendly <- function(hic.obj,
                                                          bin.df,
                                                          window,
                                                          bin,
                                                          query.pl = "~/fusessh/scripts/machine_learning/epigenome_nets/deepC/match_query_table.pl",
                                                          prepare.pl = "~/fusessh/scripts/machine_learning/epigenome_nets/deepC/prepare_query_table.pl",
                                                          timestamp = paste0(gsub("\\s+", "_", Sys.time()), "_", runif(n=1))){

  #with timestamp and random for tempfiles

  temp.hic.file <- paste0("temp_hic_file_", timestamp, ".txt")
  temp.binned.genome.file <- paste0("temp_binned_genome_file_", timestamp, ".txt")
  temp.query.table <- paste0("temp_query_table_", timestamp, ".txt")
  temp.query.outfile <- paste0("temp_query_outfile_", timestamp, ".txt")

  write.table(hic.obj$matrix.df, file = temp.hic.file, row.names = F, col.names = F, sep="\t", quote=F)
  print("saved temp_hic_file")

  write.table(bin.df, file = temp.binned.genome.file, row.names = F, col.names = F, sep="\t", quote=F)
  print("saved temp_binned_genome_file")

  # # Make Query Data frame
  # system("rm -f temp_query_table.txt", ignore.stdout = TRUE) # remove existing temp query table

  expected.values <- window/ bin

  # make query table in perl script
  command <- paste0("perl ", prepare.pl, " ", temp.binned.genome.file, " ", bin, " ", temp.query.table)
  print(command)
  system(command, ignore.stdout = TRUE)
  print("prepared query table")

  # match query using perl script
  command <- paste0("perl ", query.pl, " ", temp.hic.file, " ", temp.query.table, " ", temp.query.outfile, " ", bin)
  print("matching ...")
  system(command, ignore.stdout = TRUE)
  # read back in match queried data tale as data frame
  mdf <- as_tibble(read.table(temp.query.outfile, header=F, colClasses = c("character", rep("numeric",2+expected.values))))
  # add names
  names(mdf) <- c("chr", "start", "end", c(1:expected.values))

  # clean temp files
  command <- paste0("rm -f ", temp.hic.file, " ", temp.query.table, " ", temp.query.outfile, " ", temp.binned.genome.file)
  print(command)
  system(command)

  return(mdf)
}


#' medianImputeZerosDataFrame
#'
#' Sets zero values in an interaction matrix of a hic object in the 4th column
#' onwards to the median of a 2xk x 2xk window
#' padded with median of total matrix
#' Smooth replace according to neighborhood (median of neighborhood)
#' With k being the neighborhood size +- k around central position.
#' @param d hic object including bin coordinates and an interaction matrix
#' @param k Neighborhood size (kxk)
#' @return Hic object with the interaciton matrix median imputed
#' @examples
#' \dontrun{
#' medianImputeZerosDataFrame(hic, k=5)
#' }
#' @export
medianImputeZerosDataFrame <- function(d, k){

  id <- as.matrix(d[,c(4:ncol(d))])

  # padd with median value
  padded.id <- matrix(data=median(id), nrow=2*k+nrow(id), ncol=2*k+ncol(id))
  padded.id[c((k+1):(nrow(padded.id)-k)), c((k+1):(ncol(padded.id)-k))] <- id

  # impure zero values with median of kxk window
  for(r in c((1+k):(nrow(padded.id)-k))){
    for(c in c((1+k):(ncol(padded.id)-k))){
      if(padded.id[r,c] == 0){
        padded.id[r,c] <- median(padded.id[c((r-k):(r+k)), c((c-k):(c+k))])
      }
    }
  }

  id <- cbind(d[,c(1:3)], padded.id[c((k+1):(nrow(padded.id)-k)), c((k+1):(ncol(padded.id)-k))])

  return(id)

}


#' pyramidBin
#'
#' Percentile normalize a deepC matrix into a pyramid scheme of unequal percentiles.
#' First quantile normalize per column with 5 % percentiles and group after
#' into 2x20% | 4x10% | 4x5% percentiles.
#' Assumes data frame with with 4th + columns encoding the interactions.
#' Apply quantile normailsation per column of that matrix.
#' @param df deepC interaction data frame in pos, bin encoding
#' @param return.means TRUE/FALSE if to return a data frame with the mean interaction values per matrix column. default = FALSE
#' @return data frame with interaction in pyramid style percentile bins
#' @examples
#' \dontrun{
#' pyramidBin(hic, return.means=FALSE)
#' }
#' @export
pyramidBin <- function(df, return.means = FALSE){


  # bin into X classes stratified over distance bins
  mat <- as.matrix(df[,c(4:dim(df)[2])])
  # quantile normalize per column 5 % percentiles and group after
  matq <- apply(mat, 2, function(x){
    a <- ntile(as.vector(x), 20)
    return(a)
  })
  matq[matq %in% c(1:4)] <- 1
  matq[matq %in% c(5:8)] <- 2
  matq[matq %in% c(9:10)] <- 3
  matq[matq %in% c(11:12)] <- 4
  matq[matq %in% c(13:14)] <- 5
  matq[matq %in% c(15:16)] <- 6
  matq[matq == 17] <- 7
  matq[matq == 18] <- 8
  matq[matq == 19] <- 9
  matq[matq == 20] <- 10

  # for every distance bin get mean values per quantile
  df.mat <- as_tibble(mat)
  df.mat <- df.mat %>% gather(bin, value)
  df.matq <- as_tibble(matq)
  df.matq <- df.matq %>% gather(bin, quant)
  df.mat$quant <- df.matq$quant  # bind

  # calc means
  df.mean <- df.mat %>%
    group_by(bin, quant) %>%
    summarize(mean = mean(value))
  df.mean <- df.mean %>%
    ungroup() %>%
    mutate(bin = as.integer(bin))

  df[,c(4:dim(df)[2])] <- matq

  if(return.means == TRUE){

    newlist <- list(df = df, df.mean=df.mean)
    return(newlist)

  }else{
    return(df)
  }

}


# Helpers for Deletion Processing and Scoring ================================
# Many functions are experimental and have not been used in the manuscript

#' makeGaussianMask
#'
#' Make Gaussian Mask that weights interactions or interaction differences in the
#' center and the diagonal of a given window.
#' Use plotGaussianMask to visualise it.
#' @param window.size Window size in bp normally 1 MB + 1x bin.size
#' @param bin.size Bin size in bp default = 10000
#' @param diag.sigma.pos Sigma parameter for gaussian for the diagonals along postion dimension (e.g. 75000 for a 1010000 window size and diagonal from the center)
#' @param diag.sigma.bin Sigma parameter for gaussian for the diagonals along bin/genomic distance dimension (e.g. 25 for a 1010000 window size and diagonal from the center)
#' @param cen.sigma.pos Sigma parameter for gaussian for the center along postion dimension (e.g. 200000 for a 1010000 window size and diagonal from the center)
#' @param cen.sigma.bin igma parameter for gaussian for the diagonals along bin/genomic distance dimension (e.g. 40 for a 1010000 window size and diagonal from the center)
#' @return Mask for weighting.
#' @examples
#' \dontrun{
#' makeGaussianMask(window.size = 1010000,
#' bin.size = 10000,
#' diag.sigma.pos = 75000,
#' diag.sigma.bin = 25,
#' cen.sigma.pos = 200000,
#' cen.sigma.bin = 40)
#' }
#' @export
makeGaussianMask <- function(window.size = 1010000,
                             bin.size = 10000,
                             diag.sigma.pos = 75000,
                             diag.sigma.bin = 25,
                             cen.sigma.pos = 200000,
                             cen.sigma.bin = 40){

  halfwindow <- window.size/2
  halfbin <- bin.size/2

  # make grid
  x <- c(
    rep(c(
      seq((-1 * halfwindow + halfbin), halfwindow, by = bin.size),
      seq((-1 * halfwindow), halfwindow - halfbin, by = bin.size)
    ), (prediction.bins - 1 )/2),
    seq((-1 * halfwindow + halfbin), halfwindow, by = bin.size)
  )
  y <- rep(c(1:prediction.bins), each = window.size/bin.size)

  grid <- tibble(
    pos = x,
    bin = y,
    chr = "chrMask"
  )

  mask.df <- grid %>%
    mutate(weight.r = dnorm((bin * bin.size) - pos, mean = pos, sd = diag.sigma.pos) * dnorm(bin, mean = bin, sd = diag.sigma.bin)) %>%
    mutate(weight.l = dnorm((-1 * bin * bin.size) - pos, mean = pos, sd = diag.sigma.pos) * dnorm(bin, mean = bin, sd = diag.sigma.bin)) %>%
    mutate(weight.c = dnorm(pos, mean = 0, sd = cen.sigma.pos) * dnorm(bin, mean = 0, sd = cen.sigma.bin)) %>%
    mutate(weight.s = if_else(weight.l > weight.r, weight.l, weight.r))

  return(mask.df)
}

#' plotGaussianMask
#'
#' Plot a Gaussian Mask that weights interactions.
#' @param mask Gaussian amsk as produced by makeGaussianMask
#' @param bin.size bin size in bp
#' @return ggplot2 plot of Gaussian mask
#' @examples
#' \dontrun{
#' plotGaussianMask(mask, 5000)
#' }
#' @export
plotGaussianMask <- function(mask, bin.size){

  slope <- 1/(bin.size/2)

  m1 <- ggplot(triangularize(mask, bin.size, extra="weight.c"), aes(x = pos, y = bin, group = polyid, fill = weight.c)) +
    geom_polygon() +
    scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = slope, intercept = 0, linetype = "dashed") +
    geom_abline(slope = -slope, intercept = 0, linetype = "dashed")

  m2 <- ggplot(triangularize(mask, bin.size, extra="weight.s"), aes(x = pos, y = bin, group = polyid, fill = weight.s)) +
    geom_polygon() +
    scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = slope, intercept = 0, linetype = "dashed") +
    geom_abline(slope = -slope, intercept = 0, linetype = "dashed")

  newlist <- list(m1 = m1, m2 = m2)
  return(newlist)

}

# Calling TAD Boundaries -----------------------------

#' callInsulationBoundaries
#'
#' Call insulation score based TAD boundaries on a deepC data frame.
#' Calculate the insulation score and approximations of the first and second order derivative.
#' Select minima in the insulation score by looking for zero crossings in teh first order derivative "delta".
#' Choose minima based on neighborhood and filter them on based on their sharpness "delta2".
#' Requires a deepC like data frame (HiC, Skeleton, Prediction)
#' ins.distance - distance from sliding center position in which to calculate the insulation score.
#' For predictions 25 kb works great for HiC need larger distance.
#' Delta2 threshold for local minima to be classed as boundaries.
#' @param d deepc interaction data frame
#' @param ins.dist Genomic distance window based on which which to calculate the insulation. Default 25000
#' Default in Crane et al. is 250 k.
#' @param delta2.thresh Threshold of the delta2 score, the sharpness, (second order derivative approximation)
#' based on which to filter boundaries. Keep above threshold default = 5
#' @return Returns a dataframe based on the HiC/deepC position (bin == 1) with new columns for
#' insulationscore, delta and delta2 (frst and second derivative approcimations and a 0/1 flag if called as minimum
#' also a vector of local minima positions. List("df.ins" = df.ins, "boundaries")
#' @examples
#' \dontrun{
#' callInsulationBoundaries(df, ins.dist = 25000, delta2.thresh = 5)
#' }
#' @export
callInsulationBoundaries <- function(d, ins.dist = 25000, delta2.thresh = 5){

  # calculate running insulation score
  df.ins <- tibble(
    pos = pull(d %>% filter(bin == 1) %>% select(pos)),
    ins = 0
  )

  df.ins$ins <- pbsapply(df.ins$pos, function(x){
    s <- mean(getDistanceWindow(d, x, ins.dist)$value)
    # s <- sum(getDistanceWindow(d, x, ins.dist)$value)
    return(s)
  })

  # calc change in insulation score 1D Sobel
  df.ins$delta <- Sobeln(df.ins$ins)
  # change in change (second derivative of ins score)
  df.ins$delta2 <- Sobeln(df.ins$delta)
  # get zero crossings > minima
  df.ins$minimum <- 0
  df.ins[which(diff(sign(df.ins$delta)) == 2), "minimum"] <- 1 # one before zero crossing
  df.ins[which(df.ins$delta == 0), "minimum"] <- 1 # when zero corssing is exact

  # filter minima with second derivative < 5 (threshold)
  df.ins <- df.ins %>%
    mutate(minimum = if_else(minimum == 1 & delta2 < delta2.thresh, 0, minimum))

  local.minima <- pull(df.ins %>% filter(minimum == 1) %>% select(pos))

  newlist <- list("df.ins" = df.ins, "boundaries" = local.minima)
  return(newlist)

}


#' Sobeln
#'
#' Calculate 1st derivative approximation of profile by applying 1D sobel filter
#' @param profile vector of values / signal / profile
#' @return 1st derivative approximation of profile (numeric vector)
#' @examples
#' \dontrun{
#' Sobeln(in)
#' }
#' @export
Sobeln <- function(profile){


  b=rep(0,length(profile))
  for(i in c(2:(length(profile)-1))){
    b[i] = ( profile[i-1] * -1 ) + ( profile[i] * 0 ) + ( profile[i+1] * 1 )
  }
  b[1] = b[2]
  b[length(profile)] <- b[length(profile)-1]
  return(b)
}


# TEMPORARY DISABLED should be obsolete now -----------------
# DO NOT RUN ? USE!
# # temporary for bug fixing
# pyramidBinBug <- function(df, return.means = FALSE){
#
#   # bin into X classes stratified over distance bins
#   mat <- as.matrix(df[,c(4:dim(df)[2])])
#   # quantile normalize per column 5 % percentiles and group after
#   matq <- apply(mat, 2, function(x){
#     a <- ntile(as.vector(x), 20)
#     return(a)
#   })
#   matq[matq %in% c(1:4)] <- 1
#   matq[matq %in% c(5:8)] <- 2
#   matq[matq %in% c(1:10)] <- 3
#   matq[matq %in% c(11:12)] <- 4
#   matq[matq %in% c(13:14)] <- 5
#   matq[matq %in% c(15:16)] <- 6
#   matq[matq == 17] <- 7
#   matq[matq == 18] <- 8
#   matq[matq == 19] <- 9
#   matq[matq == 20] <- 10
#
#   # for every distance bin get mean values per quantile
#   df.mat <- as_tibble(mat)
#   df.mat <- df.mat %>% gather(bin, value)
#   df.matq <- as_tibble(matq)
#   df.matq <- df.matq %>% gather(bin, quant)
#   df.mat$quant <- df.matq$quant  # bind
#
#   # calc means
#   df.mean <- df.mat %>%
#     group_by(bin, quant) %>%
#     summarize(mean = mean(value))
#   df.mean <- df.mean %>%
#     ungroup() %>%
#     mutate(bin = as.integer(bin))
#
#   df[,c(4:dim(df)[2])] <- matq
#
#   if(return.means == TRUE){
#
#     newlist <- list(df = df, df.mean=df.mean)
#     return(newlist)
#
#   }else{
#     return(df)
#   }
#
# }

# # helper fucntion to read a bedtools deployment prediction file
# readDeepcBedtoolFile <- function(file, chrom="chr16", prediction.bins = 101, gather = TRUE){
#
#   df <- as_tibble(read.table(file, header = F))
#   names(df) <- c("chr", "start", "end", "tag", c(1:prediction.bins))
#
#   if(gather == TRUE){
#     df <- df %>%
#       mutate(pos = start + (end - start)/2) %>%
#       select(chr, pos, c(5:(prediction.bins+5))) %>%
#       gather(bin, value, -chr, -pos) %>%
#       mutate(bin = as.numeric(bin)) %>%
#       filter(chr == chrom)
#   }
#   return(df)
# }
#
# # helper to read skeleton coordinate to bin files as data frame (tibble)
# readSkeletonRegrBins <- function(file, prediction.bins = 50, gather=TRUE){
#
#   sdf <- as_tibble(read.table(file))
#   sdf$V4 <- as.character(sdf$V4)
#   names(sdf) <- c("chr", "start", "end", "qbins")
#
#   df <- sdf[,c(1:3)]
#   names(df) <- c("chr", "start", "end")
#   df <- df %>%
#     mutate(pos = (end - start)/2 + start) %>%
#     select(c(chr, pos))
#
#   # split and fill classes into numeric matrix
#   mat <- matrix(data=0, nrow = dim(sdf)[1], ncol = prediction.bins)
#   b <- as.character(sdf$qbins)
#   for(i in c(1:length(b))){
#     s <- as.numeric(strsplit(b[i], ",")[[1]])
#     mat[i,] <- s
#   }
#
#   # combine and gather for plotting
#   sdf <- as_tibble(cbind(df, mat))
#
#   #gather if specified
#   if(gather == TRUE){
#     sdf <- sdf %>%
#       gather(bin, value, -chr, -pos) %>%
#       mutate(bin = as.numeric(bin))
#   }
#
#   return(sdf)
# }

# NON zig zag encoding version for virtual 4C legasy ================
# # save the non zig zag encoding version for now
# virtual4Cold <- function(idf, chr, pos, bin.size, window.size){
#   downstream.pred.pos <- custom_ceil(pos, bin.size)
#   upstream.pred.pos <- custom_floor(pos, bin.size)
#
#   down.df <- tibble(chr = chr, pos = seq(downstream.pred.pos, downstream.pred.pos + (window.size/4) - bin.size, bin.size), value = 0)
#   up.df <- tibble(chr = chr, pos = rev(seq(upstream.pred.pos - (window.size/4) - bin.size, upstream.pred.pos, bin.size)), value = 0)
#   # fill data frames
#   for(i in c(1:nrow(down.df))){
#     to.put <- idf %>% filter(pos == down.df$pos[i] & bin == i)
#     if(nrow(to.put) > 0){
#       down.df[i,"value"] <- pull(to.put %>% select(value))
#     }else{
#       down.df[i,"value"] <- 0
#     }
#   }
#   for(i in c(1:nrow(up.df))){
#     if(up.df$pos[i] >= 0){
#       to.put <- idf %>% filter(pos == up.df$pos[i] & bin == i)
#       if(nrow(to.put) > 0){
#         up.df[i,"value"] <- pull(to.put %>% select(value))
#       }else{
#         up.df[i,"value"] <- 0
#       }
#     }
#   }
#
#   # adjust pos to represent floor of the interacting bin
#   for(i in c(1:nrow(down.df))){
#     down.df[i,"pos"] <- down.df[i,"pos"] + i * bin.size
#   }
#   for(i in c(1:nrow(up.df))){
#     up.df[i,"pos"] <- up.df[i,"pos"] - i * bin.size
#   }
#
#   # combine and remove negative positions (or above chrom_sizes)
#   v4c.df <- rbind(down.df %>% mutate(direction = "ds"), up.df %>% mutate(direction = "us"))
#   v4c.df <- v4c.df %>%
#     filter(pos >= 0) %>%
#     arrange(pos)
#
#   return(v4c.df)
# }

# old bedtools reliant getBinnedChrom ===========
#' getBinnedChrom
#' 
#' helper function to get a binned genome using bedtools (requires bedtools
#' available) Will Run bedtools on your default shell saving temporary files in
#' your working directory. Create a binned genome data frame on the specified
#' chromosome from start to end positions in window.sized bins and step sized
#' steps. See "bedtools makewindows" which is the fucntion invoked.
# @param chr chromosome
# @param start start position
# @param end end postion
# @param window size of the window / bins
# @param step The step size between windows/bins. The increment.
# @param timestamp automatically produces a timestamp to uniquely identify
#   temporary files defaut = paste0(gsub("\\s+", "_", Sys.time()), "_",
#   runif(n=1))
# @return data frame with genomic bins
# @examples
# \dontrun{
# getBinnedChrom(chr='chr1', start=0, end=10000, window=1000, step=500)
# }
# @export
# getBinnedChrom <- function(chr='chr1', start=0, end=10000, window=1000, step=500, timestamp = paste0(gsub("\\s+", "_", Sys.time()), "_", runif(n=1))){
#   # HELPER FUCNTION: get binned genome using bedtools (requires module load bedtools)
#   
#   tempfile.for.bedtools.call <- paste0("tempfile_for_bedtools_call_", timestamp, ".bed")
#   tempfile.from.bedtools.call <- paste0("tempfile_from_bedtools_call_", timestamp, ".bed")
#   
#   # calc effective end point
#   effective.end <- start + window * (floor((end - start)/window) + 1)
#   # write tempfile for bedtools
#   tdf <- data.frame(chr=chr, start=start, end=effective.end)
#   write.table(tdf, file=tempfile.for.bedtools.call, col.names=F, row.names = F, quote=F, sep="\t")
#   # use bedtools makewindows to make the bins
#   command.string <- paste0("bedtools makewindows -b ", tempfile.for.bedtools.call, "  -w ", format(window, scientific = F), "  -s ", format(step, scientific = F), " >", tempfile.from.bedtools.call)
#   print(command.string)
#   system(command.string)
#   # read back in
#   bg <- as.data.frame(read.table(tempfile.from.bedtools.call, header=F, colClasses=c("character", rep("numeric",2))))
#   # remove temp files
#   system(paste0("rm -f ", tempfile.for.bedtools.call, " ", tempfile.from.bedtools.call))
#   
#   return(bg)
# }

