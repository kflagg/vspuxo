# North-South Transect sampling from the easy site,
# with different window sizes.
require(tcltk)
require(spatstat)
source('rfns/data_functions.r')
source('rfns/sampling_functions.r')
source('rfns/spatial_functions.r')

# Paths to command-line versions of GAMV and KT3D
gamv_exe <- 'gamv'
kt3d_exe <- 'kt3d'

# Paths to input/output
fulldir <- 'datasets/easy/full'
sampdir <- 'datasets/easy/sample'
outdir <- 'datasets/easy/exp1'
# Note: GAM/GAMV only support 40 character file paths.


## EXPERIMENT PARAMETERS

nreps <- 3000
nsamp <- 200 # 200 replicates take about 15 hours on my machine

# Sampling plan parameters, 0.99 prob of detecting the smaller TA
spacing <- 225
width <- 6
bg.dens <- 100
fg.dens <- 200

# Spacing is 225, smaller TA has minor axis diameter 800
window.sizes <- c(150, 228, 516, 798, 1500)
ncells <- length(window.sizes)
nobs <- nsamp * ncells

# True TAs and site
TA1 <- ellipse(800/2, 1200/2, c(1558400, 540000), pi/6)
TA2 <- ellipse(2000/2, 900/2, c(1562000, 537000), 0)
TAs <- union.owin(TA1, TA2)
sitewindow <- owin(poly = cbind(x = c(1564294, 1564495, 1556870, 1557126),
                                y = c(535421, 541130, 541085, 535576)))

# Corners of site and number of discretized rows and columns
corners <- vertices(Frame(sitewindow))
xmin <- min(corners$x) + window.sizes / 12
ymin <- min(corners$y) + window.sizes / 12
nx <- ceiling(6 * (max(corners$x) - min(corners$x)) / window.sizes)
ny <- ceiling(6 * (max(corners$y) - min(corners$y)) / window.sizes)


# Starting values for numerically estimating semivariogram parameters:

# There should not be a nugget because the simulation has no measurement
# error or microscale variation.
nug.start <- 0

# The number of anomlies in a window follows a Poisson distribution, and
# sites have little area occupied by TAs, so the expected number of
# background anomalies over the area squared (=density over area)
# is a natural starting point for the sill of the local density estimate.
sill.start <- bg.dens / (width*window.sizes/43560)

# The only locations that should be correlated are locations in the same TA,
# so set an initial range on the same order of magnitude as the TA sizes.
range.start <- 1000

# Basic starting point for power model:
slope.start <- 1
power.start <- 0.5


# SELECT THE SAMPLE
set.seed(37478)
seeds <- sample(2^32-1, nreps)-2^31
samp <- sample(nreps, nsamp)


## RESULT STORAGE

# Matrix to store all responses that are not vectors
results <- data.frame(expand.grid('win' = window.sizes,
                                  'Realization' = samp),
                      'length' = numeric(ncells),
                      'detect' = numeric(ncells),
                      'detect1' = numeric(ncells),
                      'detect2' = numeric(ncells),
                      'dens' = numeric(ncells),
                      'detectarea' = numeric(ncells),
                      'detectarea1' = numeric(ncells),
                      'detectarea2' = numeric(ncells),
                      'identarea' = numeric(ncells),
                      'identcount' = numeric(ncells))

# List of lists to store vectors of distances of false negatives to nearest
# delineated regions
ndist <- array(list(), dim = c(nsamp, length(window.sizes)),
               dimnames = list('Realization' = paste0('r', samp),
                               'win' = paste0('w', window.sizes)))

# List of lists to store vectors of areas of disjoint regions
areas <- array(list(), dim = c(nsamp, length(window.sizes)),
               dimnames = list('Realization' = paste0('r', samp),
                               'win' = paste0('w', window.sizes)))


# Loop for each replicate
pb <- tkProgressBar(max = nsamp+1, min = 1, initial = 1,
                     title = 'Sampling and Kriging',
                     label = 'Sampling and Kriging')
timing <- system.time(for(repl in samp){
  itr <- which(samp==repl)
  r <- (itr - 1) * ncells
  set.seed(seeds[repl])
  setTkProgressBar(pb, itr, label = paste0('Iteration ', which(samp==repl),
                                            ': Sampling rep ', repl))

  ## SAMPLING

  # Read the ground truth file
  load(file = sprintf('%s/easy_full_bg%03d_fg%03d_rep%04d.RData',
                      fulldir, bg.dens, fg.dens, repl))

  # Sample along the transects, starting at a random horizontal coordinate
  sample <- sample.transects.NS(site, width, spacing,
                                offset = runif(1, 0, spacing + width/2))

  # Save the sample
  filepath <- sprintf('%s/easy_sample_sp%04d_bg%03d_fg%03d_rep%04d',
                      sampdir, spacing, bg.dens, fg.dens, repl)
  write.anomaly(sample$anomaly, paste0(filepath, '.anomaly'))
  write.cog(sample$cog, paste0(filepath, '.cog'))

  # Loop for each window size
  for(w in 1:length(window.sizes)){
    setTkProgressBar(pb, itr+w/length(window.sizes),
                      label = paste0('Iteration ', itr, ': Analyzing rep ',
                                  repl, ' with window size ', window.sizes[w]))
    results$length[r+w] <- sample$length

    ## KRIGING

    # Evaluate local density in each window
    datfile <- sprintf('%s/rep%04d_w%04d.dat', outdir, repl, window.sizes[w])
    ldens <- windowed.density.NS(sample, window.sizes[w])
    write.geoeas(ldens, datfile, title = 'Data exported from R')

    # Create GAMV parameter file
    gpar <- sprintf('%s/rep%04d_w%04d_g.par', outdir, repl, window.sizes[w])
    gout <- sprintf('%s/rep%04d_w%04d_g.out', outdir, repl, window.sizes[w])
    cat(gamv_par(datfile, gout, window.sizes[w]), file = gpar)

    # Run GAMV to compute empirical semivariogram
    system2(gamv_exe, input = gpar, wait = TRUE)

    # Read semivariogram and discard lags that were not estimated
    svario <- read.table(gout, row.names = 1,
                         col.names = c('l', 'lag.dist', 'semivariogram',
                                       'n', 'tail.mean', 'head.mean'),
                         header = FALSE, skip = 3)
    svario <- svario[svario$n>0,]

    # Fit parametric semivariograms
    params <- list('sphere' = optim(c(nug.start, sill.start[w], range.start),
                                    sv.wss, lags = svario$lag.dist,
                                    n = svario$n, ghat = svario$semivariogram,
                                    model = sv.sphere, lower = c(0, 0.0001, 0),
                                    upper = c(Inf, Inf, Inf),
                                    method = 'L-BFGS-B'),
                   'expon' = optim(c(nug.start, sill.start[w], range.start),
                                   sv.wss, lags = svario$lag.dist,
                                   n = svario$n, ghat = svario$semivariogram,
                                   model = sv.expon, lower = c(0, 0.0001, 0),
                                   upper = c(Inf, Inf, Inf),
                                   method = 'L-BFGS-B'),
                   'gauss' = optim(c(nug.start, sill.start[w], range.start),
                                   sv.wss, lags = svario$lag.dist,
                                   n = svario$n, ghat = svario$semivariogram,
                                   model = sv.gauss, lower = c(0, 0.0001, 0),
                                   upper = c(Inf, Inf, Inf),
                                   method = 'L-BFGS-B'),
                   'power' = optim(c(nug.start, slope.start, power.start),
                                   sv.wss, lags = svario$lag.dist,
                                   n = svario$n, ghat = svario$semivariogram,
                                   model = sv.power, lower = c(0, 0.0001, 0),
                                   upper = c(Inf, Inf, 2),
                                   method = 'L-BFGS-B'))

    # Find the model with the smallest sum of squares
    # The indices match GSLIB's model type numbers
    type <- order(sapply(params, function(x){return(x$value)}))[1]

    # Create KT3D parameter file
    kpar <- sprintf('%s/rep%04d_w%04d_k.par', outdir, repl, window.sizes[w])
    kout <- sprintf('%s/rep%04d_w%04d_k.out', outdir, repl, window.sizes[w])
    kdbg <- sprintf('%s/rep%04d_w%04d.dbg', outdir, repl, window.sizes[w])
    cat(kt3d_par(datfile, kdbg, kout, nx[w], ny[w], xmin[w], ymin[w],
                 window.sizes[w], params[[type]]$par[1],
                 params[[type]]$par[2], params[[type]]$par[3], type),
        file = kpar)

    # Run KT3D to do the kriging
    # Note: Value of -999 indicates that the value that was not computed
    system2(kt3d_exe, input = kpar, wait = TRUE)


    ## ANALYSIS

    ## Read and clean KT3D output
    krige.out <- read.geoeas(kout)
    krige.out[krige.out$Estimate == -999,] <- rep(NA, 2)
    krige.out$EstimationVariance[krige.out$EstimationVariance < 0] <- 0
    kest <- im(matrix(krige.out$Estimate, nrow = ny[w], byrow = TRUE),
               seq(xmin[w], length.out = nx[w], by = window.sizes[w]/6),
               seq(ymin[w], length.out = ny[w], by = window.sizes[w]/6),
               unitname = c('foot', 'feet'))
    kvar <- im(matrix(krige.out$EstimationVariance, nrow = ny[w], byrow = TRUE),
               seq(xmin[w], length.out = nx[w], by = window.sizes[w]/6),
               seq(ymin[w], length.out = ny[w], by = window.sizes[w]/6),
               unitname = c('foot', 'feet'))

    # Get the delineated regions, if there are any
    highdens <- kest > bg.dens + qnorm(0.95) * sqrt(kvar)
    if(sum(highdens) > 0){
      identified <- connected(highdens, background = FALSE)
      areas[[itr, w]] <- sapply(levels(identified$v), function(x){
          return(area(Window(connected(identified==x,
                                       background = FALSE))[sitewindow]))
      })

      # Ignore regions with non-positive area
      ignore <- which(areas[[itr, w]] <= 0)
      areas[[itr, w]][ignore] <- NA
      results$identcount[r+w] <- sum(!is.na(areas[[itr, w]]))
      for(i in ignore){
        identified$v[identified$v==i] <- NA
      }
    }
    if(results$identcount[r+w] > 0){
      idboundary <- as.polygonal(Window(identified))[sitewindow]
      idpoints <- site
      Window(idpoints) <- idboundary
      missedpoints <- site
      Window(missedpoints) <- complement.owin(idboundary,
        frame = dilation(Frame(sitewindow), window.sizes[w]))

      results$detect[r+w] <- sum(marks(idpoints)>0) / sum(marks(site)>0)
      results$detect1[r+w] <- sum(marks(idpoints)==1) / sum(marks(site)==1)
      results$detect2[r+w] <- sum(marks(idpoints)==2) / sum(marks(site)==2)
      results$dens[r+w] <- sum(marks(idpoints)>0) / area(idpoints) * 43560
      a <- intersect.owin(Window(idpoints), TAs, fatal = FALSE)
      results$detectarea[r+w] <- ifelse(is.null(a), 0, area(a))
      a <- intersect.owin(Window(idpoints), TA1, fatal = FALSE)
      results$detectarea1[r+w] <- ifelse(is.null(a), 0, area(a))
      a <- intersect.owin(Window(idpoints), TA2, fatal = FALSE)
      results$detectarea2[r+w] <- ifelse(is.null(a), 0, area(a))
      results$identarea[r+w] <- area(idpoints)
      ndist[[itr, w]] <- nncross(missedpoints[marks(missedpoints)>0],
                                 edges(idboundary), what = 'dist')
    }
  }
})
setTkProgressBar(pb, nsamp+1, label = 'Done')
invisible(close(pb))
print(timing)

save(results, file = paste0('datasets/easy/results/easy_winresults_sp', spacing,
                            '_fg', fg.dens, '.RData'))
save(ndist, file = paste0('datasets/easy/results/easy_winndist_sp', spacing,
                          '_fg', fg.dens, '.RData'))
save(areas, file = paste0('datasets/easy/results/easy_winareas_sp', spacing,
                          '_fg', fg.dens, '.RData'))
