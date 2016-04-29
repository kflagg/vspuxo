# North-South Transect sampling from the medium site,
# with twelve different sampling plans.
require(tcltk)
require(spatstat)
source('rfns/data_functions.r')
source('rfns/sampling_functions.r')
source('rfns/spatial_functions.r')

# Paths to command-line versions of GAMV and KT3D
gamv_exe <- 'gamv'
kt3d_exe <- 'kt3d'

# Paths to input/output
fulldir <- 'datasets/medium/full'
sampdir <- 'datasets/medium/sample'
outdir <- 'datasets/medium/exp2'


## EXPERIMENT PARAMETERS

nreps <- 3000
nsamp <- 100 # Takes about 60 hours

# Site parameters
width <- 6
bg.dens <- 100
r.dens <- 100
t.dens <- 200
a.dens <- 200

# Sampling plan parameters (treatments)
ta.prior <- c('Small', 'T1', 'A', 'Large')
fg.prior <- c(100, 200, 400)
spacings <- matrix(c(30,  70, 175,  320,
                     65, 135, 400,  780,
                    130, 315, 785, 1145), ncol = 3)
window.sizes <- 0.9 * c(424, 600, 1500, 2121)
ncells <- length(ta.prior) * length(fg.prior)
nobs <- nsamp * ncells

# True TAs and site
T1 <- ellipse(1000/2, 600/2, c(1558000, 540000), -pi/9)
T2 <- disc(800/2, c(1558300, 537500))
A  <- disc(1500/2, c(1561200, 539200))
TAs <- union.owin(T1, T2, A)
sitewindow <- owin(poly = cbind(x = c(1564294, 1564495, 1556870, 1557126),
                                y = c(535421, 541130, 541085, 535576)))
roadwindow <- intersect.owin(
  dilation(psp(x0 = c(1559750, 1560000, 1560250, 1560700,
                      1560700, 1560000, 1558050),
               y0 = c(535421, 536400, 536750, 537000,
                      537000, 537850, 538500),
               x1 = c(1560000, 1560250, 1560700, 1564495,
                      1560000, 1558050, 1557550),
               y1 = c(536400, 536750, 537000, 538000,
                      537850, 538500, 538900),
               window = boundingbox(sitewindow)),
           25), sitewindow)

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
set.seed(73578)
seeds <- sample(2^32-1, nreps)-2^31
samp <- sample(nreps, nsamp)


## RESULT STORAGE

# Matrix to store all responses that are not vectors
results2m <- data.frame(expand.grid('Target' = ta.prior,
                                    'fg' = fg.prior,
                                    'Realization' = samp),
                        'length' = numeric(ncells),
                        'detect' = numeric(ncells),
                        'detect1' = numeric(ncells),
                        'detect3' = numeric(ncells),
                        'detect4' = numeric(ncells),
                        'detect5' = numeric(ncells),
                        'dens' = numeric(ncells),
                        'detectarea' = numeric(ncells),
                        'detectarea1' = numeric(ncells),
                        'detectarea3' = numeric(ncells),
                        'detectarea4' = numeric(ncells),
                        'detectarea5' = numeric(ncells),
                        'identarea' = numeric(ncells),
                        'identcount' = numeric(ncells))

# List of lists to store vectors of distances of false negatives to nearest
# delineated regions
ndist2m <- array(list(), dim = c(nsamp, length(ta.prior), length(fg.prior)),
                 dimnames = list('Realization' = paste0('r', samp),
                                 'Target' = paste0('Target', ta.prior),
                                 'fg' = paste0('fg', fg.prior)))

# List of lists to store vectors of areas of disjoint regions
areas2m <- array(list(), dim = c(nsamp, length(ta.prior), length(fg.prior)),
                 dimnames = list('Realization' = paste0('r', samp),
                                 'Target' = paste0('Target', ta.prior),
                                 'fg' = paste0('fg', fg.prior)))


# Loop for each realization
pb <- TkProgressBar(max = nsamp+1, min = 1, initial = 1,
                     title = 'Sampling and Kriging',
                     label = 'Sampling and Kriging')
timing <- system.time(for(itr in seq_along(samp)){
  r <- (itr-1) * ncells
  repl <- results2m$Realization[r+1]
  set.seed(seeds[repl])
  setTkProgressBar(pb, r, label = paste0('Iteration ', itr,
                                          ': Loading rep ', repl))

  # Read the ground truth file
  load(file = sprintf('%s/medium_full_bg%03d_ro%03d_t%03d_a%03d_rep%04d.RData',
                      fulldir, bg.dens, r.dens, t.dens, a.dens, repl))

  # Loop for each treatment combination
  for(trt in seq_len(ncells)){
    ta <- which(ta.prior==results2m$Target[r+trt])
    fg <- which(fg.prior==results2m$fg[r+trt])
    setTkProgressBar(pb, itr+trt/ncells,
                      label = paste0('Iteration ', itr,
                                     ': Analyzing rep ', repl,
                                     ' with spacing ', spacings[ta, fg]))


    ## SAMPLING

    # Sample along the transects, starting at a random horzontal coordinate
    sample <- sample.transects.NS(site, width, spacings[ta, fg],
                                  offset = runif(1, 0, spacings[ta, fg] + width/2))

    # Save the sample
    filepath <- sprintf('%s/medium_sample_t%s_p%03d_bg%03d_rep%04d',
                        sampdir, ta.prior[ta], fg.prior[fg], bg.dens, repl)
    write.anomaly(sample$anomaly, paste0(filepath, '.anomaly'))
    write.cog(sample$cog, paste0(filepath, '.cog'))
    results2m$length[r+trt] <- sample$length


    ## KRIGING

    # Evaluate local density in each window
    datfile <- sprintf('%s/rep%04d_s%04d.dat', outdir, repl, spacings[ta, fg])
    ldens <- windowed.density.NS(sample, window.sizes[ta])
    write.geoeas(ldens, datfile, title = 'Data exported from R')

    # Create GAMV parameter file
    gpar <- sprintf('%s/rep%04d_s%04d_g.par', outdir, repl, spacings[ta, fg])
    gout <- sprintf('%s/rep%04d_s%04d_g.out', outdir, repl, spacings[ta, fg])
    cat(gamv_par(datfile, gout, window.sizes[ta]), file = gpar)

    # Run GAMV to compute empirical semivariogram
    system2(gamv_exe, input = gpar, wait = TRUE)

    # Read semivariogram and discard lags that were not estimated
    svario <- read.table(gout, row.names = 1,
                         col.names = c('l', 'lag.dist', 'semivariogram',
                                       'n', 'tail.mean', 'head.mean'),
                         header = FALSE, skip = 3)
    svario <- svario[svario$n>0,]

    # Fit parametric semivariograms
    params <- list('sphere' = optim(c(nug.start, sill.start[ta], range.start),
                                    sv.wss, lags = svario$lag.dist,
                                    n = svario$n, ghat = svario$semivariogram,
                                    model = sv.sphere, lower = c(0, 0.0001, 0),
                                    upper = c(Inf, Inf, Inf),
                                    method = 'L-BFGS-B'),
                    'expon' = optim(c(nug.start, sill.start[ta], range.start),
                                    sv.wss, lags = svario$lag.dist,
                                    n = svario$n, ghat = svario$semivariogram,
                                    model = sv.expon, lower = c(0, 0.0001, 0),
                                    upper = c(Inf, Inf, Inf),
                                    method = 'L-BFGS-B'),
                    'gauss' = optim(c(nug.start, sill.start[ta], range.start),
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
    kpar <- sprintf('%s/rep%04d_s%04d_k.par', outdir, repl, spacings[ta, fg])
    kout <- sprintf('%s/rep%04d_s%04d_k.out', outdir, repl, spacings[ta, fg])
    kdbg <- sprintf('%s/rep%04d_s%04d.dbg', outdir, repl, spacings[ta, fg])
    cat(kt3d_par(datfile, kdbg, kout, nx[ta], ny[ta], xmin[ta], ymin[ta],
                 window.sizes[ta], params[[type]]$par[1],
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
    kest <- im(matrix(krige.out$Estimate, nrow = ny[ta], byrow = TRUE),
               seq(xmin[ta], length.out = nx[ta], by = window.sizes[ta]/6),
               seq(ymin[ta], length.out = ny[ta], by = window.sizes[ta]/6),
               unitname = c('foot', 'feet'))
    kvar <- im(matrix(krige.out$EstimationVariance, nrow = ny[ta], byrow = TRUE),
               seq(xmin[ta], length.out = nx[ta], by = window.sizes[ta]/6),
               seq(ymin[ta], length.out = ny[ta], by = window.sizes[ta]/6),
               unitname = c('foot', 'feet'))

    # Get the delineated regions, if there are any
    highdens <- kest > bg.dens + qnorm(0.95) * sqrt(kvar)
    if(sum(highdens) > 0){
      identified <- connected(highdens, background = FALSE)
      areas2m[[itr, ta, fg]] <- sapply(levels(identified$v), function(x){
          return(area(Window(connected(identified==x,
                                       background = FALSE))[sitewindow]))
      })

      # Ignore regions less than 3 acres
      ignore <- which(areas2m[[itr, ta, fg]] < 3*43560)
      areas2m[[itr, ta, fg]][ignore] <- NA
      results2m$identcount[r+trt] <- sum(!is.na(areas2m[[itr, ta, fg]]))
      for(i in ignore){
        identified$v[identified$v==i] <- NA
      }
    }
    if(results2m$identcount[r+trt] > 0){
      idboundary <- as.polygonal(Window(identified))[sitewindow]
      idpoints <- site
      Window(idpoints) <- idboundary
      missedpoints <- site
      Window(missedpoints) <- complement.owin(idboundary,
        frame = dilation(Frame(sitewindow), window.sizes[ta]))

      results2m$detect[r+trt] <- sum(marks(idpoints)>2) / sum(marks(site)>2)
      results2m$detect1[r+trt] <- sum(marks(idpoints)==1) / sum(marks(site)==1)
      results2m$detect3[r+trt] <- sum(marks(idpoints)==3) / sum(marks(site)==3)
      results2m$detect4[r+trt] <- sum(marks(idpoints)==4) / sum(marks(site)==4)
      results2m$detect5[r+trt] <- sum(marks(idpoints)==5) / sum(marks(site)==5)
      results2m$dens[r+trt] <- sum(marks(idpoints)>0) / area(idpoints) * 43560
      a <- intersect.owin(Window(idpoints), TAs, fatal = FALSE)
      results2m$detectarea[r+trt] <- ifelse(is.null(a), 0, area(a))
      a <- intersect.owin(Window(idpoints), roadwindow, fatal = FALSE)
      results2m$detectarea1[r+trt] <- ifelse(is.null(a), 0, area(a))
      a <- intersect.owin(Window(idpoints), T1, fatal = FALSE)
      results2m$detectarea3[r+trt] <- ifelse(is.null(a), 0, area(a))
      a <- intersect.owin(Window(idpoints), T2, fatal = FALSE)
      results2m$detectarea4[r+trt] <- ifelse(is.null(a), 0, area(a))
      a <- intersect.owin(Window(idpoints), A, fatal = FALSE)
      results2m$detectarea5[r+trt] <- ifelse(is.null(a), 0, area(a))
      results2m$identarea[r+trt] <- area(idpoints)
      ndist2m[[itr, ta, fg]] <- nncross(missedpoints[marks(missedpoints)>0],
                                        edges(idboundary), what = 'dist')
    }
  }
})
setTkProgressBar(pb, nsamp+1, label = 'Done')
invisible(close(pb))
print(timing)

save(results2m, file = 'datasets/medium/results/medium_exp2results.RData')
save(ndist2m, file = 'datasets/medium/results/medium_exp2ndist.RData')
save(areas2m, file = 'datasets/medium/results/medium_exp2areas.RData')
