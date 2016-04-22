# Create several detection probability curves for the small TA at the easy site
require(RDCOMClient)

nreps <- 10
window.sizes <- c(80, 160, 400, 640, 720, 1000, 1500)
min.spacing <- 1
max.spacing <- 600
probs <- data.frame('spacing' = min.spacing:max.spacing,
                    matrix(numeric(),
                           nrow = length(min.spacing:max.spacing),
                           ncol = length(window.sizes) * nreps))
colnames(probs)[-1] <- paste0('w', rep(window.sizes,
                                       rep(nreps, length(window.sizes))),
                                       '.', rep(1:nreps, length(window.sizes)))
major.diam <- 1200
minor.diam <- 800
fg.dens <- 200
bg.dens <- 100
width <- 6

# Start VSP
vsp <- COMCreate('VSample.Document')

pb <- winProgressBar(title = 'Computing Detection Probabilities',
                     label = 'Computing Detection Probabilities',
                     min = 1, max = length(window.sizes)+1, initial = 1)
for(w in window.sizes){
  for(i in 1:nreps){
    setWinProgressBar(pb, which(window.sizes==w)+(i-1)/nreps,
                      label = paste0(w, ' foot window, rep ', i))

    # If this function returns a positive integer, that number is the index of
    # the argument that it didn't like
    vsp$Plan()$UXOPowerCurveS(width, # xsect width
                              0, # map units (feet = 0)
                              0, # xsect pattern (parallel = 0)
                              fg.dens, # fg density
                              3, # density units (per acre = 3?)
                              2, # density at (center = 2?)
                              3, # swath ratio
                              major.diam/2, # target radius
                              minor.diam/major.diam, # target shape
                              TRUE, # random angle
                              0, # angle
                              bg.dens, # bg density
                              1, # 1-false negative rate
                              0.05, # alpha
                              w, # window length
                              0.03, # min precision
                              0.01, # max error
                              min.spacing, # min spacing
                              max.spacing, # max spacing
                              TRUE) # bivariate normal?
    probs[,paste0('w', w, '.', i)] <- sapply(probs$spacing, function(x){
        vsp$Plan()$UXOPowerCurveY(x)
      })
  }
}
setWinProgressBar(pb, length(window.sizes)+1, label = 'Done')
close(pb)

# Close VSP
vsp$Window()$Close('VSampl1')

save(probs, file = 'datasets/plan/detectionprobs.RData')
