# Take a sample along evenly-spaced transects going north-south
sample.transects.NS <- function(data, width, spacing, offset = 0){
  first <- min(vertices(Frame(data))$x) + width/2 + offset
  last <- max(vertices(Frame(data))$x) - width/2
  bottom <- min(vertices(Frame(data))$y)
  top <- max(vertices(Frame(data))$y)
  x <- seq(first, last, spacing+width)
  xsect <- do.call(union.owin, lapply(x, function(i){
      owin(c(i-width/2, i+width/2), c(bottom, top))
    }))
  sample.cog <- data.frame('x' = sort(rep(x, 2)),
                           'y' = rep(c(bottom-width, top+width,
                                       top+width, bottom-width),
                                     ceiling(length(x)/2))[1:(2*length(x))])
  sample.data <- data
  Window(sample.data) <- intersect.owin(xsect, Window(data))
  ybd <- sapply(x, function(i){
      return(range(crossing.psp(
        edges(sample.data),
        psp(x0 = i, x1 = i, y0 = min(sample.cog$y), y1 = max(sample.cog$y),
            window = owin(c(i-1, i+1), range(sample.cog$y))))$y))
    })
  sample.length <- sum(ybd[2,]-ybd[1,])
  return(list('anomaly' = sample.data,
              'cog' = sample.cog,
              'length' = sample.length))
}

# Compute the density in one window
local.density <- function(center, data, winsize){
  new.win <- intersect.owin(Window(data),
                            disc(radius = winsize/2, centre = center))
  Window(data) <- new.win
  return(npoints(data)/area(new.win))
}

# Compute the density in evenly-spaced windows centered along transects
windowed.density.NS <- function(data, winsize){
  x <- unique(data$cog$x)
  ybd <- sapply(x, function(i){
      return(range(crossing.psp(
        edges(data$anomaly),
        psp(x0 = i, x1 = i, y0 = min(data$cog$y), y1 = max(data$cog$y),
            window = owin(c(i-1, i+1), range(data$cog$y))))$y))
    })
  y <- lapply(seq_len(ncol(ybd)), function(i){
      return(seq(min(ybd[,i]), max(ybd[,i]), winsize/6))
  })
  centers <- cbind('X' = rep(x, sapply(y, length)), 'Y' = unlist(y))
  return(data.frame(centers, 'Z' = 0,
    'density' = 43560 * apply(centers, 1, local.density,
                              data = data$anomaly, winsize = winsize)))
}

# Par file templates
gamv_par <- function(dat, out, win){
  return(paste0('                  Parameters for GAMV
                  *******************

START OF PARAMETERS:
', dat, '    \\file with data
1  2  0    \\  columns for X, Y, Z coordinates
1  4    \\  number of variables,col numbers
-1.0e+21	1.0e+21    \\  trimming limits
', out, '    \\file for variogram output
36    \\number of lags
', win/6, '    \\lag separation distance
', win/12, '    \\lag tolerance
1    \\number of directions
0.00  90.00 7625.00 0.00  90.00  50.00  \\azm,atol,bandh,dip,dtol,bandv
0    \\standardize sills? (0=no, 1=yes)
1    \\number of variograms
1  1  1  \\tail var., head var., variogram type
')
  )
}
kt3d_par <- function(dat, dbg, out, nx, ny, xmin, ymin, win,
                     nug, sill, range, type){
  return(paste0('                 Parameters for KT3D
                 *******************

START OF PARAMETERS:
', dat, '                \\file with data
0  1  2  0  4  0                 \\   columns for DH,X,Y,Z,var,sec var
-1.0e21   1.0e21                 \\   trimming limits
0                                \\option: 0=grid, 1=cross, 2=jackknife
nodata                           \\file with jackknife data
1   2   0    4    0              \\   columns for X,Y,Z,vr and sec var
1                                \\debugging level: 0,1,2,3
', dbg, '                \\file for debugging output
', out, '              \\file for kriged output
', nx, '  ', xmin, '  ', win/6, '             \\nx,xmn,xsiz
', ny, '  ', ymin, '  ', win/6, '             \\ny,ymn,ysiz
1    0.5    1.0                  \\nz,zmn,zsiz
1    1      1                    \\x,y and z block discretization
2    50                          \\min, max data for kriging
8                                \\max per octant (0-> not used)
', win*2, '  ', win*2, '  0.0                \\maximum search radii
 0.0   0.0   0.0                 \\angles for search ellipsoid
1     0                          \\0=SK,1=OK,2=non-st SK,3=exdrift
0 0 0 0 0 0 0 0 0                \\drift: x,y,z,xx,yy,zz,xy,xz,zy
0                                \\0, variable; 1, estimate trend
nodata                           \\gridded file with drift/mean
4                                \\  column number in gridded file
1  ', nug,'                         \\nst, nugget effect
', type, '  ', sill, '  0.0   0.0   0.0     \\it,cc,ang1,ang2,ang3
         ', range, '  ', range, '  0.0     \\a_hmax, a_hmin, a_vert
')
  )
}
