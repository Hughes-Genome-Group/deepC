#!/usr/bin/env Rscript
# Wrapper script for plotting deepC predictions.
# Required deepC prediction output(s) and optional 1D signals as bigwig tracks.
# Creates plots over the specified genomic range of all data inputs specified.
# Can plot reference and varian prediction
# Can overlap with HiC-skeleton or raw Hi-C data
# Can overlap with to three bigwig tracks

# Script parameters and help messages --------------------------
library(optparse)
option_list = list(
  
  make_option("--sample", type="character", default="my_sample",
              help="Sample name tag to use for naming output files", metavar="character"),
  
  make_option("--out.dir", type="character", default=".",
              help="Path to output directory. Default is current working directory.", metavar="character"),
  
  make_option("--bin.size", type="integer", default="5000",
              help="Genomic bin size of the HiC data.", metavar="integer"),
  
  make_option("--window.size", type="integer", default="1005000",
              help="Window size of the DNA sequence window to feed to deepC.", metavar="integer"),
  
  make_option("--chrom", type="character", default="chrG",
              help="Select a chromosome of plot reigon of interest. Must match 
              with any single chromosome inputs (predictions, hic, skeleton).", metavar="character"),
  
  make_option("--plot.start", type="integer", default=0,
              help="Start position of region to plot.", metavar="integer"),
  
  make_option("--plot.end", type="integer", default=3e+06,
              help="End position of region to plot.", metavar="integer"),
  
  make_option("--plot.width", type="integer", default=10,
              help="Relative width of plot (suggested 10).", metavar="integer"),
  
  make_option("--rel.heights", type="character", default='1',
              help="Relative heigths of sub plots to combine. 
              Supply as comma separated string e.g. \'1,1,0.5,0.5\'.
              Will set to 1 for all if not supplied or if supplied numbers 
              don't match the number of subplots specified.", metavar="character"),
  
  make_option("--plot.height", type="integer", default=15,
              help="Relative height of plot (suggested 15 when plotting all possible data).", metavar="integer"),
  
  make_option("--plot.hic", action="store_true", default=FALSE,
              help="Set TRUE to plot Hi-C data in region.", metavar="boolean"),
  
  make_option("--plot.skeleton", action="store_true", default=FALSE,
              help="Set TRUE to plot Hi-C skeleton in region.", metavar="boolean"),
  
  make_option("--plot.deepc.ref", action="store_true", default=FALSE,
              help="Set TRUE to plot deepC prediction of reference region 
              from provided prediction results. (Can be a different variant to compare.)", metavar="boolean"),
  
  make_option("--plot.deepc.var", action="store_true", default=FALSE,
              help="Set TRUE to plot deepC prediction of variant region 
              from provided prediction results.", metavar="boolean"),
  
  make_option("--fill.deepc.var", action="store_true", default=FALSE,
              help="Set TRUE to fill positions missing in the variant prediction
              with data from the reference prediction. Useful to create large 
              variant plots where only change in a subregion are expected and predicted.", metavar="boolean"),
  
  make_option("--plot.deepc.diff", action="store_true", default=FALSE,
              help="Set TRUE to plot difference ref - var of deepC predictions 
              from provided prediction results.", metavar="boolean"),
  
  make_option("--calc.deepc.diff", action="store_true", default=FALSE,
              help="Set TRUE to compute and report the mean absolute 
              interaction difference ref - var.", metavar="boolean"),
  
  make_option("--plot.tracks", action="store_true", default=FALSE,
              help="Set TRUE to plot supplied bigwig tracks in the region of interest. 
              Up to 3 bigwig files can be supplied.
              Supply in order via (--track.input.1, --track.input.2, --track.input.3);
              Set labels via (--track.name.1, --track.name.2, --track.name.3);
              Set colours via (--track.colour.1, --track.colour.2, --track.colour.3)", metavar="boolean"),
  
  make_option("--track.input.1", type="character", default="None",
              help="Path to bigwig file for track 1.", metavar="character"),
  make_option("--track.input.2", type="character", default="None",
              help="Path to bigwig file for track 2.", metavar="character"),
  make_option("--track.input.3", type="character", default="None",
              help="Path to bigwig file for track 3.", metavar="character"),
  
  make_option("--track.name.1", type="character", default="track 1",
              help="Path to bigwig file for track 1.", metavar="character"),
  make_option("--track.name.2", type="character", default="track 2",
              help="Path to bigwig file for track 2.", metavar="character"),
  make_option("--track.name.3", type="character", default="track 3",
              help="Path to bigwig file for track 3.", metavar="character"),
  
  make_option("--track.colour.1", type="character", default="#E41A1C",
              help="Colour for track 1 in HEX string e.g. #E41A1C", metavar="character"),
  make_option("--track.colour.2", type="character", default="#377EB8",
              help="Colour for track 2 in HEX string e.g. #E377EB8", metavar="character"),
  make_option("--track.colour.3", type="character", default= "#4DAF4A",
              help="Colour for track 3 in HEX string e.g. #4DAF4A", metavar="character"),
  
  make_option("--plot.title.hic", type="character", default='None',
              help="Optional title for hic plot.", metavar="character"),
  make_option("--plot.title.skel", type="character", default='None',
              help="Optional title for hic skeleton plot.", metavar="character"),
  make_option("--plot.title.ref", type="character", default='None',
              help="Optional title for deepC reference plot.", metavar="character"),
  make_option("--plot.title.var", type="character", default='None',
              help="Optional title for deepC variant plot.", metavar="character"),
  make_option("--plot.title.diff", type="character", default='None',
              help="Optional title for deepC diff plot.", metavar="character"),
  make_option("--plot.title.tracks", type="character", default='None',
              help="Optional title for tracks plot.", metavar="character"),

  make_option("--deepc.ref.input", type="character", default='None',
              help="Input for deepC reference prediction, as procuded 
              by run_deploy_net.py etc. tab separate text file with columns: 
              chr start_seq_window end_seq_window predicted_interactions_bin1 predicted_interactions_bin2
              See example in \"test_predict_out\" in tutorials for an example.", metavar="character"),
  
  make_option("--deepc.var.input", type="character", default='None',
              help="Input for deepC variant prediction, as procuded 
              by run_deploy_net.py etc. tab separate text file with columns: 
              chr start_seq_window end_seq_window predicted_interactions_bin1 predicted_interactions_bin2
              See example in \"test_variant_out\" in tutorials for an example.", metavar="character"),

  make_option("--hic.matrix", type="character", default='None',
              help='Single chromosome, sparse interaction matrix.
                Supports two input versions:
                A) A three column matrix with the positions specified
                in column one and two reflecting the leftmost position
                of the genomic bin window and the third column indicating
                the interaction frequency (raw or noarmalised).
                B) A matrix file with ids for the first two columns and interaction
                frequencies in the thirds column AND a coordinate file supplied to
                --hic.coords that links chr and genomic positions to genomic
                window ids. For example the output of HiC-Pro.', metavar="character"),

  make_option("--hic.coords", type="character", default='None',
              help="Only required for two input file input variant. \"None\" otherwise.
                Requires a 4 column format: chr start end interaction_frequency.
                Assumes 0-based coordinates.", metavar="character"),

  make_option("--chromosome.sizes", type="character", default='None',
              help="Chromosome sizes file (UCSC style) two column format 
              (chr chromosome.size). Only required if Hi-C data provided in matrix format.", metavar="character"),

  make_option("--hic.preprocessed ", type="character", default='None',
              help="Preprocessed Hi-C data in deepC format. Quicker to load 
              Hi-C preprocessed especially if plotting multiple times. 
              Otherwise, format raw Hi-C data from sparse contact matrix (slower).
              To get preprocessed Hi-C data use preprocess wrapper with --no.transform flag.", metavar="character"),

  make_option("--skeleton.input", type="character", default='None',
              help="Preprocessed Hi-C skeleton per chromosome. Produced by preprocessing wrapper. 
              If not supplied will attempt to calculate from Hi-C input.", metavar="character"),
  
  make_option("--impute.zeros", action="store_true", default=FALSE,
              help="Set TRUE to attempt to impute interaction frequency 
              zeros with the median of a 5x5 neighbourhood. Only checked 
              if skeleton is calculated here and not provided from preprocessed input.", metavar="boolean"),

  make_option("--keep.median.zero", action="store_true", default=FALSE,
              help="Sets TRUE to keep genomic windows over which the ZigZag 
              interaction pole has median 0 interaction frequency
              (default remove those). Only checked if skeleton is 
              calculated here and not provided from preprocessed input.", metavar="boolean"),

  make_option("--helper", type="character", default="./deepC/helper_for_preprocessing_and_analysis/",
              help="Path to helper directory of deepC repository. 
              Will attempt to read helper scripts and source functions from within here.", metavar="character")
  
  );

# parse options 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Derived parameters
path.to.helper.functions <- opt$helper
bin.size <- opt$bin.size
window.size <- opt$window.size
# number of predictions bins (also called number of classes later).
prediction.bins <- window.size/bin.size

# calculate number of subplots to produce
num.sub.plots <- sum(c(opt$plot.hic, opt$plot.skeleton, opt$plot.deepc.ref, opt$plot.deepc.var, opt$plot.deepc.diff, opt$plot.tracks))
# parse relative sub plot heights
rel.heights <- as.numeric(unlist(strsplit(opt$rel.heights, ',')))
# split and check that relative plot heights work out
if(length(rel.heights) == 1){
  # expand relative heights to number of sub plots
  rel.heights <- rep(rel.heights, num.sub.plots)
}
# check if relative heigths matches 
if(length(rel.heights) != num.sub.plots){
  warning('Number of relative plot heights does not match up with the specified number of sub plots to produce. Seeting all realtive to 1')
  rel.heights <- rep(rel.heights[1], num.sub.plots)
}

# Load required libraries -----------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
theme_set(theme_cowplot())
# load rtracklayer if bigwig plots desired
if(opt$plot.tracks){
  suppressPackageStartupMessages(library(rtracklayer))
}
# source helper functions -----------------------------------
source(paste0(path.to.helper.functions, "/functions_for_HiC.R"))
source(paste0(path.to.helper.functions, "/functions_for_deepC.R"))

# Basic sanity checks --------------------------------------
# check that window size / bin.size is an odd number for
if(prediction.bins %% 2 != 1){
  stop(paste('Genomic window to process divided by the bin.size',
             'should be an odd number to ensure symmetrie over the',
             'central genomic bin. Supplied:', window.size, 'and', bin.size))
}

# 1.1) IF plotting raw Hi-C data is required gather HiC data ------------------
if(opt$plot.hic | (opt$skeleton.input == 'None' & opt$plot.skeleton)){
  
  # If already have a deepC formated 'raw' HiC input file load this here (much quicker)
  if(opt$hic.preprocessed != 'None'){
    print('Loading deepC formated Hi-C data.')
    tryCatch(hic.df <- readDeepcInputHicBed(opt$hic.preprocessed, prediction.bins = prediction.bins, bin.size = bin.size, gather=FALSE),
             error = function(cond){message('Reading the preprocessed Hi-C data failed. Please check the input!')},
             warnings = function(cond){message('Reading the preprocessed HiC data produced warnings. Please check the input!')})
  }else{
    # Otherwise process Hi-C data from   
    print('Procseeing Hi-C data from sparse matrix. If plotting multple times consider pre-processing your Hi-C data!')

    if(opt$chromosome.sizes == 'None'){
      stop('Chromosome sizes file is required when trying to plot hic data from sparse matrix files.')
    }
    # read in chromosome sizes file
    genome.sizes <- tryCatch({read_tsv(opt$chromosome.sizes, col_names = F, progress = F, show_col_types = F)}, 
                             error = function(cond){message('Genome sizes file could not be read!')},
                             warnings = function(cond){message('Reading the genome sizes file procuded warnings!')})
    names(genome.sizes) <- c('chr', 'size')
    
    # get the respective chromosome size 
    if(genome.sizes %>% filter(chr == opt$chrom) %>% nrow() != 1){
      stop(paste('Multiple or no entries for chromosome size found in the supllied sizes file for', opt$chrom))
    }
    chrom.size <- genome.sizes %>% filter(chr == opt$chrom) %>% pull(size)
    # Read in Hi-C data
    # Check with reading mode to use
    # check hic.coord.input if that is the default 'None' string assume matrix only input
    # other wise assume amtrix and coordinate input (e.g. HiC-Pro output)
    if(opt$hic.coords == 'None'){
      # Matrix only input 
      tryCatch(hic <- ImportHicproMatrix(opt$hic.matrix, chr=opt$chrom, bin.size=bin.size),
               error = function(cond){message('Could not read HiC matrix!')},
               warnings = function(cond){message('Reading the hic matrix prodcued errors!')})
    }else{
      # Matrix and genomic window coordinates file
      tryCatch(hic <- ImportHicproMatrix(opt$hic.matrix, coords=opt$hic.coords, bin.size=bin.size),
               error = function(cond){message('Could not read HiC matrix or coordinate file!')},
               warnings = function(cond){message('Reading the hic matrix and.or coordinate file prodcued errors!')})
    }
    
    # remove all interactions that are more then bp_context apart
    hic <- trimHicRange(hic, range = window.size + bin.size)
    
    # get the the first covered position in the HiC data (corrected for centrally encoded bins)
    start.pos <- hic$start.pos - bin.size/2
    
    # to map the hic data into consistent bins for deepC a binned genome template
    binned.genome <- getBinnedChrom(chr=opt$chrom, start=start.pos, end=chrom.size, window=window.size, step=bin.size)
    # filter out ibcombplete megabase bins
    binned.genome <- binned.genome %>%
      mutate(width = end - start) %>%
      filter(width == window.size) %>%
      select(-width)
    
    print('Converting interactions to ZigZag pole coordinates ...')
    hic.df <- getZigZagWindowInteractionsPerlMemoryFriendly(hic,
                                                         binned.genome,
                                                         window.size,
                                                         bin.size,
                                                         prepare.pl = paste0(path.to.helper.functions, "/prepare_query_table.pl"),
                                                         query.pl= paste0(path.to.helper.functions, "/match_query_table.pl"))
    
  }
}

# 2) IF plotting skeleton Hi-C data is required either load skeleton or process from Hi-C input ------------------
if(opt$plot.skeleton){
  
  # If already skeleton in deepC format this here (much quicker)
  if(opt$skeleton.input != 'None'){
    print('Loading preprocessed skeleton Hi-C data.')
    tryCatch(skel.df <-  skel.df <- readDeepcInputHicBed(opt$skeleton.input, prediction.bins = prediction.bins, bin.size = bin.size, gather=FALSE),
             error = function(cond){message('Reading the preprocessed skeleton data failed. Please check the input!')},
             warnings = function(cond){message('Reading the preprocessed skeleton data produced warnings. Please check the input!')})
  }else{
    print('No preprocessed skeleton provided, calculating it from Hi-C input. To speed up consider providing the skeleton data preprocessed.')
    # check if hic.df exists
    if(!exists("hic.df")){
      stop('Provide Hi-C input to compute the skeleton if no preprocessed skeleton is provided! Stopping.')
    }
    # Impute zero values in interaction matrix with Median of a 5x5 neighborhood (optional)
    # padd with median of respective interaction distance for calculating 
    # the median on edge cases and filter median 0 genomic windows from data set
    # first add median column
    skel.df <- hic.df
    skel.df$median.value <- apply(hic.df, 1, function(x){
      m <- median(as.numeric(x[c(4:length(x))]))
      return(m)
    })
    # filter out genomic windows with median 0 interaction values
    if(!opt$keep.median.zero){
      skel.df <- skel.df[skel.df$median.value > 0,]  
    }
    # impute zeros with 5x5 median value
    if(opt$impute.zeros == T){
      print('Trying to impute zeros with 5x5 neighbourhood median values.')  
      tryCatch(skel.df <- medianImputeZerosDataFrame(skel.df, k=5),
               error = function(cond){message('Imputation failed! If matrices are to sparse try without median imputation!')},
               warnings = function(cond){message('Imputation prdocued warnings! If matrices are to sparse try without median imputation!')})
    }
    # remove padded column
    skel.df <- skel.df[,-ncol(skel.df)]
    # convert to tibble
    skel.df <- as_tibble(skel.df)
    # Transform to Hi-C skeleton (pyramid percentile binning)
    print('Transforming interactions to skeleton ...')
    skel.df <- pyramidBin(skel.df)
  }
  # gather and filter for plot.region of interest
  skel.df <- skel.df %>%
    mutate(pos = start+((end-start)/2)) %>%
    filter(pos >= opt$plot.start, pos <= opt$plot.end) %>%
    gather(bin, value, -chr, -start, -end, -pos, factor_key = T) %>%  # gather
    mutate(bin = as.numeric(bin)) %>%
    mutate(pos = if_else(bin %% 2 == 0, pos - bin.size/2, pos)) # adjust positions for zig-zag layout

  # make the sub plot
  # convert to diamond ploygon map for plotting
  tpdf.skel <- triangularize(skel.df, bin = bin.size)
  # make the plot
  p.skel <- ggplot(tpdf.skel, aes(x=pos, y=(bin*bin.size)/1000, fill=value, group = polyid)) + 
    geom_polygon() +
    labs(y = "Genomic Distance [kb]") +
    coord_cartesian(xlim=c(opt$plot.start, opt$plot.end)) +
    scale_fill_gradientn(colours = brewer.pal(9, 'YlOrRd'))
  if(opt$plot.title.skel != 'None'){
    p.skel <- p.skel + ggtitle(paste0(opt$plot.title.skel)) + theme(plot.title = element_text(hjust = 0.5))
  }
}

# 1.2) After skeleton is ready, hic data can be subset and gathered as well ------
if(opt$plot.hic){
  
  hic.df <- hic.df %>%
    mutate(pos = start+((end-start)/2)) %>%
    filter(pos >= opt$plot.start, pos <= opt$plot.end) %>%
    gather(bin, value, -chr, -start, -end, -pos, factor_key = T) %>%  # gather
    mutate(bin = as.numeric(bin)) %>%
    mutate(pos = if_else(bin %% 2 == 0, pos - bin.size/2, pos))  # adjust positions for 
  # zig-zag layout
  
  # To visualize HiC data we uantile squeeze the very low and very high interaction values to increase contrast in the middle value range
  hic.df$value <- SetValueRange(hic.df$value, min=as.numeric(quantile(hic.df$value, .05)), max=as.numeric(quantile(hic.df$value, .95)))
  # ^^^ this might be better to run before subsetting to ROI 
  
  # make the sub plot
  # convert to diamond ploygon map for plotting
  tpdf.hic <- triangularize(hic.df, bin = bin.size)
  # make the plot
  p.hic <- ggplot(tpdf.hic, aes(x=pos, y=(bin*bin.size)/1000, fill=value, group = polyid)) + 
    geom_polygon() +
    labs(y = "Genomic Distance [kb]") +
    coord_cartesian(xlim=c(opt$plot.start, opt$plot.end)) +
    scale_fill_gradientn(colours = brewer.pal(9, 'YlOrRd'))
  if(opt$plot.title.hic != 'None'){
    p.hic <- p.hic + ggtitle(paste0(opt$plot.title.hic)) + theme(plot.title = element_text(hjust = 0.5))
  }
}

# 3.1) Read in deepC (reference) prediction if ref and/or differential plot desired ----------------------
if(opt$plot.deepc.ref | opt$plot.deepc.diff){
  if(opt$deepc.ref.input == 'None'){
    stop('deepc.ref.input required to plot reference prediciton and.or differnetial plot. Provide output of deepC prediction.')
  }
  print('Creating deepC reference predictions.')
  tryCatch(ref.pred <- readDeepcVariantFile(opt$deepc.ref.input, prediction.bins = prediction.bins, bin.size = bin.size, gather = T),
           error = function(cond){message('Reading the reference deepC prediction failed, please check input!')},
           warnings = function(cond){message('Reading the reference deepC prediction produced errrors, please check input!')})
  # subset reference prediction to desired plot region of interest
  ref.pred$df <- ref.pred$df %>% filter(pos >= opt$plot.start, pos <= opt$plot.end)
  # check data frame not empty (close to empty in region of interest)
  if(nrow(ref.pred$df) <= 3){
    stop('Practically no data in prediction (ref) input in the plot region. Check prediction and consider adjusting plot window.')  
  }
  # Makre plot if selcted
  if(opt$plot.deepc.ref){
    trdf <- triangularize(ref.pred$df, bin = bin.size)
    p.ref <- trdf %>%
      ggplot(aes(x = pos, y = (bin*bin.size)/1000, fill = value, group = polyid)) + 
      geom_polygon() +
      labs(y = "Genomic Distance [kb]") +
      coord_cartesian(xlim=c(opt$plot.start, opt$plot.end)) +
      scale_fill_gradientn(colours = brewer.pal(9, 'YlOrRd'))
    
    if(opt$plot.title.ref != 'None'){
      p.ref <- p.ref + ggtitle(paste0(opt$plot.title.ref)) + theme(plot.title = element_text(hjust = 0.5))
    }
  }
}

# 3.2) Read in deepC (variant) prediction if variant and/or differential plot desired ----------------------
if(opt$plot.deepc.var | opt$plot.deepc.diff){
  if(opt$deepc.var.input == 'None'){
    stop('deepc.var.input required to plot variant prediciton and.or differnetial plot. Provide output of deepC prediction.')
  }
  print('Creating deepC variant predictions.')
  tryCatch(var.pred <- readDeepcVariantFile(opt$deepc.var.input, prediction.bins = prediction.bins, bin.size = bin.size, gather = T),
           error = function(cond){message('Reading the variant deepC prediction failed, please check input!')},
           warnings = function(cond){message('Reading the variant deepC prediction produced errrors, please check input!')})
  # subset reference prediction to desired plot region of interest
  var.pred$df <- var.pred$df %>% filter(pos >= opt$plot.start, pos <= opt$plot.end)
  # check data frame not empty (close to empty in region of interest)
  if(nrow(var.pred$df) <= 3){
    stop('Practically no data in prediction (var) input in the plot region. Check prediction and consider adjusting plot window.')  
  }
}

# 4) If selected calculate the difference between the supplied reference and variant -------
if(opt$calc.deepc.diff){
  # calculate the total number of HiC interaction tiles (bin positions) affected by the variant
  # used for normalizing the damage sscore
  covered.tiles <- dim(var.pred$df)[1]/prediction.bins * (window.size/bin.size) 
  # calculate the differntial data frame 
  diff.df <- getDiffDeepC(ref.pred$df, var.pred$df)
  # calculate the sum of the absolut difference and normalize for the number of covered tiles
  mean.abs.diff <- sum(abs(diff.df$diff)) / covered.tiles
  print('Calculating mean of the absolute predicted interaction difference per covered interaction bin (rounded to 5 decimal digits)')
  print(paste0('Mean abs diff (ref - var): ', round(mean.abs.diff, digits = 5)))
}

# 3.3 after optional difference calc, make variat plot and fill variant data frame if specified -------
if(opt$plot.deepc.var){
  
  # fill missing positions from ref if selected
  if(opt$fill.deepc.var){
    if(!exists("ref.pred")){
      stop('Needs reference prediction data loaded to fill variant prediction.')
    }  
    print('Filling missing positions in variant prediction with reference data.')
    to.add <- ref.pred$df %>% filter(!pos %in% var.pred$df$pos)
    var.pred$df <- rbind(to.add, var.pred$df) %>% arrange(pos)
  }

  tvdf <- triangularize(var.pred$df, bin = bin.size)
  p.var <- tvdf %>%
    ggplot(aes(x = pos, y = (bin*bin.size)/1000, fill = value, group = polyid)) + 
    geom_polygon() +
    labs(y = "Genomic Distance [kb]") +
    coord_cartesian(xlim=c(opt$plot.start, opt$plot.end)) +
    scale_fill_gradientn(colours = brewer.pal(9, 'YlOrRd'))
  
  if(opt$plot.title.var != 'None'){
    p.var <- p.var + ggtitle(paste0(opt$plot.title.var)) + theme(plot.title = element_text(hjust = 0.5))
  }
}

# 5) If selected make a differential plot -------------
if(opt$plot.deepc.diff){
  p.diff <- plotDiffDeepC(ref.pred$df, var.pred$df, bin.size = 5000) +
    geom_vline(xintercept = var.pred$variant.start, linetype = "dotted") +
    geom_abline(slope = 1/(bin.size/(2*bin.size/1000)), intercept = -var.pred$variant.start/(bin.size/(2*bin.size/1000)), linetype = "dashed") +
    geom_abline(slope = -1/(bin.size/(2*bin.size/1000)), intercept = var.pred$variant.start/(bin.size/(2*bin.size/1000)), linetype = "dashed") +
    coord_cartesian(xlim=c(opt$plot.start, opt$plot.end))
  
  if(opt$plot.title.diff != 'None'){
    p.diff <- p.diff + ggtitle(paste0(opt$plot.title.diff)) + theme(plot.title = element_text(hjust = 0.5))
  }
  
}

# 6) If selected to plot load 1D signal tracks -----------------------------
if(opt$plot.tracks){
  
  print('Loading bigwig tracks for track plots.')
  suppressMessages(library(rtracklayer))

  # check if tracks supplied
  if(opt$track.input.1 == 'None'){
    stop('Need to supply at least one 1D bigwig track if selected --plot.tracks. Supply input to opt$track.input.1 at least.')
  }
  
  # initialize GRanges object within the plot limits 
  groi <- makeGRangesFromDataFrame(tibble(chrom = "chr17", start = opt$plot.start, end = opt$plot.end))
  
  # For each supplied track, import data over grange of roi and convert to data frame for ggplot2
  selected.track.colours <- c()
  if(opt$track.input.1 != 'None'){
    gr.track.1 <- import(opt$track.input.1, which = groi)  # Import BigWig Signal
    df.track.1 <- as_tibble(gr.track.1)
    df.track.1 <- df.track.1 %>%
      mutate(pos = start + (end - start)/2) %>%
      select(c(pos, width, score))
    df.tracks <- df.track.1 %>% mutate(source = opt$track.name.1)
    selected.track.colours <- c(selected.track.colours, opt$track.colour.1)
  }
  if(opt$track.input.2 != 'None'){
    gr.track.2 <- import(opt$track.input.2, which = groi)  # Import BigWig Signal
    df.track.2 <- as_tibble(gr.track.2)
    df.track.2 <- df.track.2 %>%
      mutate(pos = start + (end - start)/2) %>%
      select(c(pos, width, score))
    # check if df.tracks exists yet, if not no track 1 specified!
    if(!exists('df.tracks')){
      stop('Please supply a bigwig track 1 if attempting to supply track 2 or switch.')
    }
    df.tracks <- rbind(df.tracks, df.track.2 %>% mutate(source = opt$track.name.2))
    selected.track.colours <- c(selected.track.colours, opt$track.colour.2)
  }
  if(opt$track.input.3 != 'None'){
    gr.track.3 <- import(opt$track.input.3, which = groi)  # Import BigWig Signal
    df.track.3 <- as_tibble(gr.track.3)
    df.track.3 <- df.track.3 %>%
      mutate(pos = start + (end - start)/2) %>%
      select(c(pos, width, score))
    if(!exists('df.tracks')){
      stop('Please supply a bigwig track 1 if attempting to supply track 3 or switch.')
    }
    df.tracks <- rbind(df.tracks, df.track.3 %>% mutate(source = opt$track.name.3))
    selected.track.colours <- c(selected.track.colours, opt$track.colour.3)
  }
  
  # Make plot of 1D signal
  p.tracks <- df.tracks %>%
    ggplot(aes(x = pos, y = score, col = source)) + 
    geom_area() +
    scale_color_manual(values = selected.track.colours) +
    xlim(opt$plot.start, opt$plot.end) +
    labs(y='Coverage') +
    facet_wrap(~source, nrow = length(selected.track.colours), strip.position = "right", scales = 'free_y') + 
    theme(strip.background = element_blank(), legend.position = "None")
 
  if(opt$plot.title.tracks != 'None'){
    p.tracks <- p.tracks + ggtitle(paste0(opt$plot.title.tracks)) + theme(plot.title = element_text(hjust = 0.5))
  }
   
}

# 7) Assemble combined plot from everything specified ---------------------
print(paste('Combining', num.sub.plots, 'sub plots ...'))

prio.list <- list()
if(opt$plot.hic){ prio.list <- c(prio.list, list(p.hic + lower_overlay_theme)) }
if(opt$plot.skeleton){ prio.list <- c(prio.list, list(p.skel + lower_overlay_theme)) }
if(opt$plot.deepc.ref){ prio.list <- c(prio.list, list(p.ref + lower_overlay_theme)) }
if(opt$plot.deepc.var){ prio.list <- c(prio.list, list(p.var + lower_overlay_theme)) }
if(opt$plot.deepc.diff){ prio.list <- c(prio.list, list(p.diff + lower_overlay_theme)) }
if(opt$plot.tracks){ prio.list <- c(prio.list, list(p.tracks + lower_overlay_theme)) }
## set lower theme for bottom plot
##prio.list[[length(prio.list)]] <- prio.list[[length(prio.list)]] + lower_overlay_theme

p.combined <- plot_grid(plotlist=prio.list,
                        nrow = num.sub.plots, 
                        align = 'v', axis = 'tblr',
                        rel_heights = rel.heights)

ggsave(p.combined, filename = paste0(opt$out.dir, '/plot_', as.character(bin.size/1000),"kb_", opt$sample, '_', 
                                     opt$chrom, '_', opt$plot.start, '_', opt$plot.end, '.png'), 
                                     width = opt$plot.width, height = opt$plot.height, bg = 'white')


