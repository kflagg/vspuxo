# Generate a realistic site with three TAs and some roads
require(spatstat)
source('rfns/data_functions.r')
source('rfns/spatial_functions.r')

nreps <- 3000


bg.dens <- 100 / 43560 # 100 per acre, converted to square feet
road.dens <- 100 / 43560
tank.dens <- 200 / 43560
art.dens <- 200 / 43560

# 952.375 acres. At 100 bg per acre, we expect 95237.5 bg anomalies.
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

# Generate a vector of random seeds and reseed each iteration
# so the whole set of reps doesn't need to be generated at once.
# Valid seeds are 32 bit signed integers.
set.seed(46347)
seeds <- sample(2^32-1, nreps)-2^31

# Loop to generate many
cat(sprintf('Simulating %d Medium Sites\n', nreps))
pb <- txtProgressBar(max = nreps, style = 3)
timing <- system.time(for(repl in 1:nreps){
  set.seed(seeds[repl])

  # Uniform background
  bg.homog <- rpoispp(lambda = bg.dens, win = sitewindow)
  marks(bg.homog) <- 0

  # Uniform background
  bg.road <- rpoispp(lambda = road.dens, win = roadwindow)
  marks(bg.road) <- 1

  # Tank Area 1
  tank1 <- rpoispp(lambda = gauss.elliptic, win = sitewindow, mu.x = 1558000, mu.y = 540000,
                   s.a = 1000/(2*qnorm(0.995)), s.b = 600/(2*qnorm(0.995)), r = -pi/9, maxrate = tank.dens)
  marks(tank1) <- 3

  # Tank Area 2
  tank2 <- rpoispp(lambda = gauss.elliptic, win = sitewindow, mu.x = 1558300, mu.y = 537500,
                   s.a = 800/(2*qnorm(0.995)), s.b = 800/(2*qnorm(0.995)), r = 0, maxrate = tank.dens)
  marks(tank2) <- 4

  # Artillery Area
  art <- rpoispp(lambda = gauss.elliptic, win = sitewindow, mu.x = 1561200, mu.y = 539200,
                 s.a = 1500/(2*qnorm(0.995)), s.b = 1500/(2*qnorm(0.995)), r = 0, maxrate = art.dens)
  marks(art) <- 5

  site <- superimpose(bg.homog, bg.road, tank1, tank2, art)
  save(site, file = sprintf('datasets/medium/full/medium_full_bg%03d_ro%03d_t%03d_a%03d_rep%04d.RData',
                            bg.dens*43560, road.dens*43560, tank.dens*43560, art.dens*43560, repl))

  setTxtProgressBar(pb, repl)
})
close(pb)
print(timing)
