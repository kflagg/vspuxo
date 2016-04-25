# Redefine the histogram plotting method to include a vertical offset
plot.histogram <- function(x, freq = equidist, density = NULL, angle = 45,
                           col = NULL, border = par("fg"), lty = NULL,
                           main = paste("Histogram of", paste(x$xname, collapse = "\n")),
                           sub = NULL, xlab = x$xname, ylab, xlim = range(x$breaks), ylim = NULL,
                           axes = TRUE, labels = FALSE, add = FALSE, ann = TRUE, voffset, ...){
  if(missing(voffset)) voffset <- 0
  equidist <- if (is.logical(x$equidist))
    x$equidist
  else{
    h <- diff(x$breaks)
    diff(range(h)) < 1e-07 * mean(h)
  }
  if(freq && !equidist)
    warning("the AREAS in the plot are wrong -- rather use 'freq = FALSE'")
  y <- if(freq)
    x$counts
  else x$density
  nB <- length(x$breaks)
  if(is.null(y) || 0L == nB)
    stop("'x' is wrongly structured")
  dev.hold()
  on.exit(dev.flush())
  if(!add){
    if (is.null(ylim))
      ylim <- range(y, 0)
    if(missing(ylab))
      ylab <- if (!freq)
        "Density"
    else "Frequency"
    plot.new()
    plot.window(xlim, ylim, "", ...)
    if(ann)
      title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
    if(axes){
      axis(1, ...)
      axis(2, ...)
    }
  }
  rect(x$breaks[-nB], voffset, x$breaks[-1L], y + voffset, col = col,
       border = border, angle = angle, density = density, lty = lty)
  if((logl <- is.logical(labels) && labels) || is.character(labels))
    text(x$mids, y, labels = if (logl){
      if (freq)
        x$counts
      else round(x$density, 3)
    }
    else labels, adj = c(0.5, -0.5))
  invisible()
}

# Function to plot overlaid histograms of a quantitative response by
# one or two categorical explanatory variables.
histby <- function(y, ...)
  UseMethod('histby')

histby.default <- function(y, x1, x2 = NULL, breaks = 'Sturges',
    x1levels = NULL, x2levels = NULL, reversex1 = FALSE, reversex2 = FALSE,
    space1 = 0.05, space2 = 0.05, themax, caterpillar = TRUE,
    quantiles = if(caterpillar == TRUE) c(0, 0.5, 0.95) else c(0.025, 0.25, 0.5, 0.75, 0.975),
    qcol = NULL, lwd = 2, col = NULL, border = par('fg'),
    freq = FALSE, at = NULL, xlim = range(y, na.rm = TRUE), ylim = NULL,
    labels = NULL, main = paste("Histogram of", yname), xlab = yname, ylab = NA,
    yaxt = par('yaxt'), xaxt = par('xaxt'), log = '', prettyx = FALSE, ...){
  yname <- paste(deparse(substitute(y), 500), collapse = '\n')
  if(missing(x1))
    x1 <- rep(TRUE, length(y))
  if(is.null(x2))
    x2 <- rep(TRUE, length(y))
  if(length(y) != length(x1) | length(y) != length(x2))
    stop('Variable lengths differ.')
  if(is.null(x1levels))
    x1levels <- unique(x1)
  if(is.null(x2levels))
    x2levels <- unique(x2)
  if(!(is.null(col) | length(col)==1))
    col <- rep(col, length(x1levels))
  if(length(border)==1)
    border <- rep(border, length(x1levels))
  if(length(lwd)==1)
    lwd <- lwd + seq_along(quantiles)-1
  if(caterpillar == TRUE)
    quantiles <- sort(quantiles, decreasing = TRUE)
  if(is.null(qcol))
    if(caterpillar == TRUE)
      qcol <- replicate(length(quantiles), rep(border, length(x1levels)))
  else
    qcol <- rep(border, length(x1levels))

  hists <- array(list(), dim = c(length(x1levels), length(x2levels)),
                 dimnames = list(NULL, NULL))
  for(l2 in seq_along(x2levels)){
    for(l1 in seq_along(x1levels)){
      if(log == 'x')
      {
        hists[[l1, l2]] <- hist(log(y[x1==x1levels[l1] & x2==x2levels[l2]]),
          plot = FALSE, breaks = breaks, warn.unused = FALSE, ...)
        hists[[l1, l2]]$breaks <- exp(hists[[l1, l2]]$breaks)
        hists[[l1, l2]]$mids <- exp(hists[[l1, l2]]$mids)
      }
      else
        hists[[l1, l2]] <- hist(y[x1==x1levels[l1] & x2==x2levels[l2]],
          plot = FALSE, breaks = breaks, warn.unused = FALSE, ...)
    }
  }

  if(missing(themax))
    themax <- max(sapply(hists, function(h){
      return(ifelse(freq, max(h$counts), max(h$density)))
    }))
  offset1 <- themax*(seq_along(x1levels)-1)*space1
  if(reversex1)
    offset1 <- offset1[length(offset1):1]
  offset2 <- (themax+max(offset1))*(seq_along(x2levels)-1)*(1+space2)
  if(reversex2)
    offset2 <- offset2[length(offset2):1]

  if(is.null(ylim))
    ylim <- c(0, themax+max(offset1)+ifelse(length(offset2>0), max(offset2), 0))

  plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE, log = log,
       main = main, ylab = ylab, xlab = xlab, ...)
  for(l2 in seq_along(x2levels)){
    for(l1 in seq_along(x1levels)){
      plot(hists[[l1, l2]], col = col[l1], border = border[l1], freq = freq,
           voffset = offset1[l1]+offset2[l2], add = TRUE, ...)
      if(quantiles[1] != FALSE)
        if(caterpillar == TRUE)
          segments(x0 = quantile(y[x1==x1levels[l1] & x2==x2levels[l2]], (1-quantiles)/2),
                   x1 = quantile(y[x1==x1levels[l1] & x2==x2levels[l2]], 1-(1-quantiles)/2),
                   y0 = offset1[l1]+offset2[l2], col = qcol[l1,], lty = 1, lwd = lwd, ...)
        else
          segments(x0 = quantile(y[x1==x1levels[l1] & x2==x2levels[l2]], quantiles),
                   y0 = rep(offset1[l1], length(quantiles)),
                   y1 = rep(offset1[l1] + offset2[l2], length(quantiles)) + themax*space1,
                   col = qcol[l1], ...)
    }
  }
  if(is.null(labels)){
    if(length(x2levels) > 1)
      labels <- x2levels
    else
      labels = FALSE
  }
  if(xaxt != 'n'){
    if(prettyx)
      axis(1, at = axTicks(1), labels = prettyNum(axTicks(1), big.mark = ','))
    else
      axis(1)
  }
  if(yaxt != 'n')
    axis(2, at = offset2+(themax+max(offset1))/2,
         tick = FALSE, labels = labels)

  invisible(hists)
}

histby.formula <- function(formula, data = NULL, ...){
  vars <- eval(attr(terms(formula), 'variables'), data, enclos = parent.frame())
  resp <- eval(attr(terms(formula), 'response'), data, enclos = parent.frame())
  if(resp==0)
    stop('No response variable specified.')
  y <- vars[[resp]]
  x <- vars[-resp]
  x1 <- x[[1]]
  if(length(x)>2)
    warning('Only the first variables are used.')
  if(length(x)<2){
    x2 <- NULL
  }
  else
    x2 <- x[[2]]
  invisible(histby(y, x1, x2, ...))
}

histby.data.frame <- function(data, y, x1, x2, ...){
  if(missing(x2))
    invisible(histby(data[,y], data[,x1], ...))
  else
    invisible(histby(data[,y], data[,x1], data[,x2], ...))
}

histby.array <- function(data, breaks = 'Sturges',
    x1levels = NULL, x2levels = NULL, reversex1 = FALSE, reversex2 = FALSE,
    space1 = 0.05, space2 = 0.05, themax,
    quantiles = c(0.25, 0.5, 0.75), qcol = border,
    col = NULL, border = par('fg'), freq = FALSE, at = NULL,
    xlim = range(data, na.rm = TRUE), ylim = NULL, labels = NULL,
    yaxt = par('yaxt'), xaxt = par('xaxt'), prettyx = FALSE, ...){
  dims <- length(dim(data)) - 1
  if(dims < 1 | dims > 2)
    stop('Must have one or two predictors.')
  if(is.null(x1levels))
    x1levels <- dimnames(data)[[2]]
  if(is.null(x2levels) & dims == 2)
    x2levels <- dimnames(data)[[3]]
  if(!(is.null(col) | length(col) == 1))
    col <- rep(col, length(x1levels))
  if(length(border) == 1)
    border <- rep(border, length(x1levels))

  if(dims == 1){
    hists <- array(list(), dim = length(x1levels), dimnames = list(NULL))
    for(l1 in seq_along(x1levels)){
      hists[[l1]] <- hist(unlist(data[, x1levels[l1]]),
                          plot = FALSE, breaks = breaks, warn.unused = FALSE, ...)
    }
  }else{
    hists <- array(list(), dim = c(length(x1levels), length(x2levels)),
                   dimnames = list(NULL, NULL))
    for(l2 in seq_along(x2levels)){
      for(l1 in seq_along(x1levels)){
        hists[[l1, l2]] <- hist(unlist(data[, x1levels[l1], x2levels[l2]]),
                                plot = FALSE, breaks = breaks, warn.unused = FALSE, ...)
      }
    }
  }

  if(missing(themax))
    themax <- max(unlist(sapply(hists, function(h){
      return(ifelse(freq, max(h$counts), max(h$density)))
    })))
  offset1 <- themax*(seq_along(x1levels)-1)*space1
  if(reversex1)
    offset1 <- offset1[length(offset1):1]
  offset2 <- (themax+max(offset1))*(seq_along(x2levels)-1)*(1+space2)
  if(reversex2)
    offset2 <- offset2[length(offset2):1]

  if(is.null(ylim))
    ylim <- c(0, themax+max(offset1)+ifelse(length(offset2>0), max(offset2), 0))

  plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE, ...)
  if(dims == 1){
    for(l1 in seq_along(x1levels)){
      plot(hists[[l1]], col = col[l1], border = border[l1], freq = freq,
           voffset = offset1[l1], add = TRUE, ...)
      if(quantiles[1] != FALSE)
        segments(x0 = quantile(unlist(data[, x1levels[l1]]), quantiles, na.rm = TRUE),
                 y0 = rep(offset1[l1], length(quantiles)),
                 y1 = rep(offset1[l1], length(quantiles)) + themax*space1,
                 col = qcol[l1], ...)
    }
  }else{
    for(l2 in seq_along(x2levels)){
      for(l1 in seq_along(x1levels)){
        plot(hists[[l1, l2]], col = col[l1], border = border[l1], freq = freq,
             voffset = offset1[l1] + offset2[l2], add = TRUE, ...)
        if(quantiles[1] != FALSE)
          segments(x0 = quantile(unlist(data[, x1levels[l1]]), quantiles, na.rm = TRUE),
                   y0 = rep(offset1[l1] + offset2[l2], length(quantiles)),
                   y1 = rep(offset1[l1] + offset2[l2], length(quantiles)) + themax*space1,
                   col = qcol[l1], ...)
      }
    }
  }
  if(is.null(labels)){
    if(length(x2levels) > 1)
      labels <- x2levels
    else
      labels = FALSE
  }
  if(xaxt != 'n'){
    if(prettyx)
      axis(1, at = axTicks(1), labels = prettyNum(axTicks(1), big.mark = ','))
    else
      axis(1)
  }
  if(yaxt != 'n')
    axis(2, at = offset2+(themax+max(offset1))/2,
         tick = FALSE, labels = labels)

  invisible(hists)
}

histby.list <- function(data, breaks = 'Sturges',
    x1levels = NULL, x2levels = NULL, reversex1 = FALSE, reversex2 = FALSE,
    space1 = 0.05, space2 = 0.05, themax,
    quantiles = c(0.25, 0.5, 0.75), qcol = border, qpch = '|',
    col = NULL, border = par('fg'), freq = FALSE, at = NULL,
    xlim = range(data, na.rm = TRUE), ylim = NULL, labels = NULL,
    yaxt = par('yaxt'), xaxt = par('xaxt'), prettyx = FALSE, ...){
  if(is.null(dim(data))){
    dims <- 1
    if(is.null(x1levels))
      x1levels <- names(data)
  }else{
    dims <- length(dim(data))
    if(is.null(x1levels))
      x1levels <- dimnames(data)[[1]]
  }
  if(dims < 1 | dims > 2)
    stop('Must have one or two predictors.')
  if(is.null(x2levels) & dims == 2)
    x2levels <- dimnames(data)[[1]]
  if(!(is.null(col) | length(col) == 1))
    col <- rep(col, length(x1levels))
  if(length(border) == 1)
    border <- rep(border, length(x1levels))

  if(dims == 1){
    hists <- array(list(), dim = length(x1levels), dimnames = list(NULL))
    for(l1 in seq_along(x1levels)){
      hists[[l1]] <- hist(unlist(data[x1levels[l1]]),
        plot = FALSE, breaks = breaks, warn.unused = FALSE, ...)
    }
  }else{
    hists <- array(list(), dim = c(length(x1levels), length(x2levels)),
                   dimnames = list(NULL, NULL))
    for(l2 in seq_along(x2levels)){
      for(l1 in seq_along(x1levels)){
        hists[[l1, l2]] <- hist(unlist(data[x1levels[l1], x2levels[l2]]),
          plot = FALSE, breaks = breaks, warn.unused = FALSE, ...)
      }
    }
  }

  if(missing(themax))
    themax <- max(unlist(sapply(hists, function(h){
      return(ifelse(freq, max(h$counts), max(h$density)))
    })))
  offset1 <- themax*(seq_along(x1levels)-1)*space1
  if(reversex1)
    offset1 <- offset1[length(offset1):1]
  offset2 <- (themax+max(offset1))*(seq_along(x2levels)-1)*(1+space2)
  if(reversex2)
    offset2 <- offset2[length(offset2):1]

  if(is.null(ylim))
    ylim <- c(0, themax+max(offset1)+ifelse(length(offset2>0), max(offset2), 0))

  plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE, ...)
  if(dims == 1){
    for(l1 in seq_along(x1levels)){
      plot(hists[[l1]], col = col[l1], border = border[l1], freq = freq,
        voffset = offset1[l1], add = TRUE, ...)
      if(quantiles[1] != FALSE)
        segments(x0 = quantile(unlist(data[x1levels[l1]]), quantiles, na.rm = TRUE),
                 y0 = rep(offset1[l1], length(quantiles)),
                 y1 = rep(offset1[l1], length(quantiles)) + themax*space1,
                 col = qcol[l1], ...)
    }
  }else{
    for(l2 in seq_along(x2levels)){
      for(l1 in seq_along(x1levels)){
        plot(hists[[l1, l2]], col = col[l1], border = border[l1], freq = freq,
          voffset = offset1[l1] + offset2[l2], add = TRUE, ...)
        if(quantiles[1] != FALSE)
          segments(x0 = quantile(unlist(data[, x1levels[l1]]), quantiles, na.rm = TRUE),
                   y0 = rep(offset1[l1] + offset2[l2], length(quantiles)),
                   y1 = rep(offset1[l1] + offset2[l2], length(quantiles)) + themax*space1,
                   col = qcol[l1], ...)
      }
    }
  }
  if(is.null(labels)){
    if(length(x2levels) > 1)
      labels <- x2levels
    else
      labels = FALSE
  }
  if(xaxt != 'n'){
    if(prettyx)
      axis(1, at = axTicks(1), labels = prettyNum(axTicks(1), big.mark = ','))
    else
      axis(1)
  }
  if(yaxt != 'n')
    axis(2, at = offset2+(themax+max(offset1))/2,
         tick = FALSE, labels = labels)

  invisible(hists)
}

# Function to plot caterpillar plots of a quantitative response by
# one or two categorical explanatory variables.
caterplot <- function(y, ...)
  UseMethod('caterplot')

caterplot.default <- function(y, x1, x2 = NULL,
    x1levels = NULL, x2levels = NULL, reversex1 = FALSE, reversex2 = FALSE,
    log = FALSE, vertical = FALSE,
    proportions = c(0.5, 0.95), type = 7, median = TRUE, mean = TRUE,
    out = TRUE, pch = 21, out.pch = pch, median.pch = 19, mean.pch = NULL,
    col = NULL, bg = par('bg'), out.col = col, out.bg = bg,
    cex = par('cex'), out.cex = 0.83*cex, lwd = par('lwd'), at = NULL,
    xlim = range(y, na.rm = TRUE), ylim = NULL, labels = NULL,
    main = paste("Caterpillar Plot of", yname), xlab = yname, ylab = NA,
    yaxt = par('yaxt'), xaxt = par('xaxt'), prettyx = FALSE, ...){
  yname <- paste(deparse(substitute(y), 500), collapse = '\n')
  if(missing(x1))
    x1 <- rep(TRUE, length(y))
  if(is.null(x2))
    x2 <- rep(TRUE, length(y))
  if(length(y) != length(x1) | length(y) != length(x2))
    stop('Variable lengths differ.')
  if(is.null(x1levels))
    x1levels <- unique(x1)
  if(is.null(x2levels))
    x2levels <- unique(x2)
  if(is.null(at))
    at <- seq_len(length(x1levels)*length(x2levels))
  if(is.null(col))
    col <- replicate(length(quantiles), rep(par('col'), length(x1levels)))
  if(length(bg)==1)
    bg <- rep(bg, length(x1levels))
  if(length(out.col)==1)
    out.col <- rep(out.col, length(x1levels))
  if(length(out.bg)==1)
    out.bg <- rep(out.bg, length(x1levels))
  if(length(lwd)==1)
    lwd <- lwd + 2*(seq_along(proportions)-1)
  if(vertical == TRUE){
    if(log == TRUE)
      log <- 'y'
    else
      log <- ''
    if(is.null(mean.pch))
      mean.pch <- '-'
  }else{
    if(log == TRUE)
      log <- 'x'
    else
      log <- ''
    if(is.null(mean.pch))
        mean.pch <- '|'
  }
  if(mean.pch == '-')
    mean.cex <- 2*cex
  else
    mean.cex <- cex
  proportions <- sort(proportions, decreasing = TRUE)

  cats <- array(numeric(), dim = c(length(x1levels), length(x2levels),
                                   length(proportions), 2),
                dimnames = list(NULL, NULL, NULL, c('Lower', 'Upper')))
  for(l2 in seq_along(x2levels)){
    for(l1 in seq_along(x1levels)){
      cats[l1, l2, , 'Lower'] <- quantile(y[x1==x1levels[l1] & x2==x2levels[l2]],
        (1-proportions)/2, na.rm = TRUE, type = type)
      cats[l1, l2, , 'Upper'] <- quantile(y[x1==x1levels[l1] & x2==x2levels[l2]],
        1-(1-proportions)/2, na.rm = TRUE, type = type)
    }
  }

  if(is.null(ylim))
    ylim <- c(0.5, length(x1levels) * length(x2levels) + 0.5)
  if(is.null(labels)){
    if(length(x2levels) > 1)
      labels <- x2levels
    else
      labels = FALSE
  }

  if(vertical == TRUE){
    plot(x = NA, y = NA, xlim = ylim, ylim = xlim, axes = FALSE,
         main = main, ylab = xlab, xlab = ylab, log = log, ...)
    for(l2 in seq_along(x2levels)){
      for(l1 in seq_along(x1levels)){
        segments(y0 = cats[l1, l2, , 1], y1 = cats[l1, l2, , 2],
                 x0 = at[l1 + length(x1levels) * (l2 - 1)],
                 col = col[l1,], bg = bg[l1], lty = 1, lwd = lwd, ...)
        if(median != FALSE)
          points(y = quantile(y[x1==x1levels[l1] & x2==x2levels[l2]], 0.5, na.rm = TRUE, type = type),
                 x = at[l1 + length(x1levels) * (l2 - 1)], cex = cex,
                 pch = median.pch, col = col[l1, 1], bg = bg[l1], lty = 1, ...)
        if(mean != FALSE)
          points(y = mean(y[x1==x1levels[l1] & x2==x2levels[l2]], na.rm = TRUE),
                 x = at[l1 + length(x1levels) * (l2 - 1)], cex = mean.cex,
                 pch = mean.pch, col = col[l1, 1], bg = bg[l1], ...)
        if(out != FALSE){
          outliers <- y[x1==x1levels[l1] & x2==x2levels[l2]]
          outliers <- outliers[outliers < min(cats[l1, l2,,]) | outliers > max(cats[l1, l2,,])]
          points(y = outliers, x = rep(at[l1 + length(x1levels) * (l2 - 1)], length(outliers)),
                 cex = out.cex, pch = out.pch, col = out.col[l1], bg = out.bg[l1], ...)
        }
      }
    }
    if(xaxt != 'n'){
      if(prettyx)
        axis(2, at = axTicks(2), labels = prettyNum(axTicks(2), big.mark = ','))
      else
        axis(2)
    }
    if(yaxt != 'n')
      axis(1, at = mean(seq_along(x1levels))+length(x1levels)*(seq_along(x2levels)-1),
           tick = FALSE, labels = labels)
  }else{
    plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE,
         main = main, ylab = ylab, xlab = xlab, log = log, ...)
    for(l2 in seq_along(x2levels)){
      for(l1 in seq_along(x1levels)){
        segments(x0 = cats[l1, l2, , 1], x1 = cats[l1, l2, , 2],
                 y0 = at[l1 + length(x1levels) * (l2 - 1)],
                 col = col[l1], bg = bg[l1], lty = 1, lwd = lwd, ...)
        if(median != FALSE)
          points(x = quantile(y[x1==x1levels[l1] & x2==x2levels[l2]], 0.5, na.rm = TRUE, type = type),
                 y = at[l1 + length(x1levels) * (l2 - 1)], cex = cex,
                 pch = median.pch, col = col[l1], bg = bg[l1], lty = 1, ...)
        if(mean != FALSE)
          points(x = mean(y[x1==x1levels[l1] & x2==x2levels[l2]], na.rm = TRUE),
                 y = at[l1 + length(x1levels) * (l2 - 1)], cex = cex,
                 pch = mean.pch, col = col[l1], bg = bg[l1], ...)
        if(out != FALSE){
          outliers <- y[x1==x1levels[l1] & x2==x2levels[l2]]
          outliers <- outliers[outliers < min(cats[l1, l2,,]) | outliers > max(cats[l1, l2,,])]
          points(x = outliers, y = rep(at[l1 + length(x1levels) * (l2 - 1)], length(outliers)),
                 cex = out.cex, pch = out.pch, col = out.col[l1], bg = out.bg[l1], ...)
        }
      }
    }
    if(xaxt != 'n'){
      if(prettyx)
        axis(1, at = axTicks(1), labels = prettyNum(axTicks(1), big.mark = ','))
      else
        axis(1)
    }
    if(yaxt != 'n')
      axis(2, at = mean(seq_along(x1levels))+length(x1levels)*(seq_along(x2levels)-1),
           tick = FALSE, labels = labels)
  }

  invisible(cats)
}

caterplot.formula <- function(formula, data = NULL, ...){
  vars <- eval(attr(terms(formula), 'variables'), data, enclos = parent.frame())
  resp <- eval(attr(terms(formula), 'response'), data, enclos = parent.frame())
  if(resp==0)
    stop('No response variable specified.')
  y <- vars[[resp]]
  x <- vars[-resp]
  x1 <- x[[1]]
  if(length(x)>2)
    warning('Only the first variables are used.')
  if(length(x)<2){
    x2 <- NULL
  }
  else
    x2 <- x[[2]]
  invisible(caterplot(y, x1, x2, ...))
}
