#!/usr/bin/env Rscript
# Wrapper script for formatting Hi-C data to skeleton transformed input for deepC model training

# Script parameters and help messages --------------------------
library(optparse)
option_list = list(
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

  make_option("--chromosome.sizes", type="character", default='',
              help="Chromosome sizes file (UCSC style) two column format (chr chromosome.size", metavar="character"),

  make_option("--sample", type="character", default="my_sample",
              help="Sample name tag to use for naming output files", metavar="character"),

  make_option("--out.dir", type="character", default=".",
              help="Path to output directory. Default is current working directory.", metavar="character"),

  make_option("--bin.size", type="integer", default="5000",
              help="Genomic bin size of the HiC data.", metavar="integer"),

  make_option("--window.size", type="integer", default="1005000",
              help="Window size of the DNA sequence window to feed to deepC.", metavar="integer"),

  make_option("--chrom", type="character", default="chrG",
              help="Select a single chromosome matching the chromosome sizes names and that is present in your hic matrix (and coord file).", metavar="character"),

  make_option("--helper", type="character", default="./deepC/helper_for_preprocessing_and_analysis/",
              help="Path to helper directory of deepC repository. Will attempt to read helper scripts and source functions from within here.", metavar="character"),

  make_option("--impute.zeros", action="store_true", default=FALSE,
              help="Set TRUE to attempt to impute interaction frequency zeros with the median of a 5x5 neighbourhood.", metavar="boolean"),
  
  make_option("--keep.median.zero", action="store_true", default=FALSE,
              help="Sets TRUE to keep genomic windows over which the ZigZag interaction pole has median 0 interaction frequency (default remove those).", metavar="boolean"),
  
  make_option("--no.transform", action="store_true", default=FALSE,
              help="Sets TRUE to NOT apply the skeleton percentile binning transformation and store Hi-C frequencies as input from the sparse matrix.", metavar="boolean"),

  make_option("--plot.hic", action="store_true", default=FALSE,
              help="Set TRUE to plot an example region plot using the Hi-C data and the region specifcied by --plot.start and --plot.end.", metavar="boolean"),

  make_option("--plot.skeleton", action="store_true", default=FALSE,
              help="Set TRUE to plot an example region plot using the Hi-C skeleton and the region specifcied by --plot.start and --plot.end.", metavar="boolean"),

  make_option("--plot.start", type="integer", default=0,
              help="Start position for the example plots.", metavar="integer"),

  make_option("--plot.end", type="integer", default=3e+06,
              help="End position for the example plots.", metavar="integer"),
  
  make_option("--plot.width", type="integer", default=10,
              help="Relative width of plot (suggested 10).", metavar="integer"),
  
  make_option("--plot.height", type="integer", default=6,
              help="Relative height of plot (suggested 6 for only hic or skeleton and 8 for double plots).", metavar="integer")

  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Derived parameters 
path.to.helper.functions <- opt$helper
bin.size <- opt$bin.size 
window.size <- opt$window.size
# number of predictions bins (also called number of classes later).
prediction.bins <- window.size/bin.size 

# Load required libraries -----------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
theme_set(theme_cowplot())

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

# read in chromosome sizes file ----------------------
genome.sizes <- tryCatch({read_tsv(opt$chromosome.sizes, col_names = F, progress = F, show_col_types = F)}, 
                         error = function(cond){message('Genome sizes file could not be read!')},
                         warnings = function(cond){message('Reading the genome sizes file procuded warnings!')})
names(genome.sizes) <- c('chr', 'size')

# Read in Hi-C data ----------------------------------
# get the respective chromosome size 
if(genome.sizes %>% filter(chr == opt$chrom) %>% nrow() != 1){
  stop(paste('Multiple or no entries for chromosome size found in the supllied sizes file for', opt$chrom))
}

chrom.size <- genome.sizes %>% filter(chr == opt$chrom) %>% pull(size)

# Check with reading mode to use --------------------
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

# remove all interactions that are more then bp_context apart ------
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

# Convert to deepC format ------------------------------------
# Convert HiC interaction matrix into genomic windows of bp_context size with the vertical 
# "zig-zag" pole of associated chromatin interactions
# Get ZigZiag Pole ZickZack Pole of interactions per window
#    /
#    \ 
#    /
#    \ 
#    |
# """""""
#     | 4 |     | 4 |
#       | 2 || 2 |
# | 5 | 3 | 1 | 3 | 5 | 
# uses helper perl scripts from the repository
print('Converting interactions to ZigZag pole coordinates ...')
tdf <- getZigZagWindowInteractionsPerlMemoryFriendly(hic,
                                                    binned.genome,
                                                    window.size,
                                                    bin.size,
                                                    prepare.pl = paste0(path.to.helper.functions, "/prepare_query_table.pl"),
                                                    query.pl= paste0(path.to.helper.functions, "/match_query_table.pl"))

# make a copy of Hi-C data if choose to save Hi-C data plot 
if(opt$plot.hic){
  hdf <- tdf
}

# Impute zero values in interaction matrix with Median of a 5x5 neighborhood (optional)
# padd with median of respective interaction distance for calculating 
# the median on edge cases and filter median 0 genomic windows from data set
# first add median column
tdf$median.value <- apply(tdf, 1, function(x){
  m <- median(as.numeric(x[c(4:length(x))]))
  return(m)
})

# filter out genomic windows with median 0 interaction values
if(!opt$keep.median.zero){
  tdf <- tdf[tdf$median.value > 0,]  
}


# impute zeros with 5x5 median value
if(opt$impute.zeros == T){
  print('Trying to impute zeros with 5x5 neighbourhood median values.')  
  tryCatch(tdf <- medianImputeZerosDataFrame(tdf, k=5),
           error = function(cond){message('Imputation failed! If matrices are to sparse try without median imputation!')},
           warnings = function(cond){message('Imputation prdocued warnings! If matrices are to sparse try without median imputation!')})
}

# remove padded column
tdf <- tdf[,-ncol(tdf)]

# convert to tibble
tdf <- as_tibble(tdf)

# Transform to Hi-C skeleton (pyramid percentile binning)
if(!opt$no.transform){
  print('Transforming interactions to skeleton ...')
  tdf <- pyramidBin(tdf)
}

# Prepare and save output file -------------------------
# collate to single comma separated column
cdf <- tdf
cdf$class <- cdf %>%
  select(-chr, -start, -end) %>%
  unite(col = class, sep = ",") %>%
  pull()

# and trim output
cdf <- cdf %>% select(chr, start, end, class)

# and save in plain txt
print('Saving deepC file ...')
write.table(cdf, file = paste0(opt$out.dir, "/coords_and_hic_skeleton_", as.character(bin.size/1000),"kb_", opt$chrom,"_", opt$sample,".bed"),
            col.names = F, row.names = F, quote = F, sep = "\t")

# Optional: Make example plot(s) ------------------------
if(opt$plot.hic | opt$plot.skeleton){
  print('creating example plot ...')
  # if selected to plot hic data example
  if(opt$plot.hic){
    print('Creating Hi-C plot ...')
    # Make a position column, gather and tranform to numeirc interactions bins
    pdf.hic <- hdf %>%
      mutate(pos = start+((end-start)/2)) %>%
      filter(pos >= opt$plot.start, pos <= opt$plot.end) %>%
      gather(bin, value, -chr, -start, -end, -pos, factor_key = T) %>%  # gather
      mutate(bin = as.numeric(bin)) %>%
      mutate(pos = if_else(bin %% 2 == 0, pos - bin.size/2, pos))  # adjust positions for 
    # zig-zag layout
    
    # To visualize HiC data we uantile squeeze the very low and very high interaction values to increase contrast in the middle value range
    pdf.hic$value <- SetValueRange(pdf.hic$value, min=as.numeric(quantile(pdf.hic$value, .05)), max=as.numeric(quantile(pdf.hic$value, .95)))
    
    # convert to diamond ploygon map for plotting
    tpdf.hic <- triangularize(pdf.hic, bin = bin.size)
    # make the plot
    p.hic <- ggplot(tpdf.hic, aes(x=pos, y=(bin*bin.size)/1000, fill=value, group = polyid)) + 
      geom_polygon() +
      labs(y = "Genomic Distance [kb]") +
      scale_fill_gradientn(colours = brewer.pal(9, 'YlOrRd'))
  }
  # if selected to plot skeleton data
  if(opt$plot.skeleton){
    print('Creating skeleton plot ...')
    
    # Make a position column, gather and tranform to numeirc interactions bins
    # take first 1000 entries
    pdf.skel <- tdf %>%
      mutate(pos = start+((end-start)/2)) %>%
      filter(pos >= opt$plot.start, pos <= opt$plot.end) %>%
      gather(bin, value, -chr, -start, -end, -pos, factor_key = T) %>%  # gather
      mutate(bin = as.numeric(bin)) %>%
      mutate(pos = if_else(bin %% 2 == 0, pos - bin.size/2, pos))  # adjust positions for 
    # zig-zag layout
    if(opt$no.transform){
      # To visualize untransformed data quantile squeeze the very low and very high interaction values
      pdf.skel$value <- SetValueRange(pdf.skel$value, min=as.numeric(quantile(pdf.skel$value, .05)), max=as.numeric(quantile(pdf.skel$value, .95)))
    }
    
    # Plot Coverage of interaction windows
    # convert to diamond ploygon map for plotting
    tpdf.skel <- triangularize(pdf.skel, bin = bin.size)
    # make the plot
    p.skel <- ggplot(tpdf.skel, aes(x=pos, y=(bin*bin.size)/1000, fill=value, group = polyid)) + 
      geom_polygon() +
      labs(y = "Genomic Distance [kb]") +
      scale_fill_gradientn(colours = brewer.pal(9, 'YlOrRd'))
  }
  
  # combine and save plots 
  if(opt$plot.hic & opt$plot.skeleton){
    print('Saving combined plot ...')
    p.combined <- plot_grid(p.hic + upper_overlay_theme,
                            p.skel + lower_overlay_theme,
                            nrow = 2, align = 'v')
    ggsave(p.combined, filename = paste0(opt$out.dir, '/example_plot_', opt$chrom, '_', as.character(bin.size/1000),"kb_", opt$sample, '.png'), width = opt$plot.width, height = opt$plot.height, bg = 'white')
  }else if(opt$plot.hic){
    print('Saving Hi-C plot ...')
    ggsave(p.hic + lower_overlay_theme, filename = paste0(opt$out.dir, '/example_plot_', opt$chrom, '_', as.character(bin.size/1000),"kb_", opt$sample, '.png'), width = opt$plot.width, height = opt$plot.height, bg = 'white')
  }else if(opt$plot.skeleton){
    print('Saving skeleton plot ...')
    ggsave(p.skel + lower_overlay_theme, filename = paste0(opt$out.dir, '/example_plot_', opt$chrom, '_', as.character(bin.size/1000),"kb_", opt$sample, '.png'), width = opt$plot.width, height = opt$plot.height, bg = 'white')
  }
}

