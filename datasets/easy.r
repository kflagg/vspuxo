# Generate an easy site with two elliptical TAs
require(spatstat)
source('rfns/data_functions.r')
source('rfns/spatial_functions.r')

nreps <- 3000

bg.dens <- 100 / 43560 # 100 per acre, converted to square feet
fg.dens <- 200 / 43560 # 200 per acre above bg

# 952.375 acres. At 100 bg per acre, we expect 95237.5 bg anomalies.
sitewindow <- owin(poly = cbind(x = c(1564294, 1564495, 1556870, 1557126),
                                y = c(535421, 541130, 541085, 535576)))

# Generate a vector of random seeds and reseed each iteration
# so the whole set of reps doesn't need to be generated at once.
# Valid seeds are 32 bit signed integers.
set.seed(783614)
seeds <- sample(2^32-1, nreps)-2^31

# Loop to generate many
cat(sprintf('Simulating %d Easy Sites\n', nreps))
pb <- txtProgressBar(max = nreps, style = 3)
timing <- system.time(for(repl in 1:nreps){
  set.seed(seeds[repl])

  # Uniform background
  bg.anomalies <- rpoispp(lambda = bg.dens, win = sitewindow)
  marks(bg.anomalies) <- 0

  # Target Area 1
  fg1 <- rpoispp(lambda = gauss.elliptic, win = sitewindow, mu.x = 1558400, mu.y = 540000,
                 s.a = 800/(2*qnorm(0.995)), s.b = 1200/(2*qnorm(0.995)), r = pi/6, maxrate = fg.dens)
  marks(fg1) <- 1

  # Target Area 2
  fg2 <- rpoispp(lambda = gauss.elliptic, win = sitewindow, mu.x = 1562000, mu.y = 537000,
                 s.a = 2000/(2*qnorm(0.995)), s.b = 900/(2*qnorm(0.995)), r = 0, maxrate = fg.dens)
  marks(fg2) <- 2

  site <- superimpose(bg.anomalies, fg1, fg2)
  save(site, file = sprintf('datasets/easy/full/easy_full_bg%03d_fg%03d_rep%04d.RData',
                            bg.dens*43560, fg.dens*43560, repl))

  setTxtProgressBar(pb, repl)
})
close(pb)
print(timing)
