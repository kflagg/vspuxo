# Write VSP anomaly files
write.anomaly <- function(data, file, col.names = TRUE, row.names = FALSE){
  return(write.table(data, file, row.names = row.names, col.names = col.names,
                     quote = FALSE, sep = ','))
}

# Read VSP anomaly files
read.anomaly <- function(file, header = TRUE, sep = ',', ...){
  return(read.table(file = file, header = header, sep = sep, ...))
}

# Write VSP cog files
write.cog <- function(data, file, col.names = TRUE, row.names = FALSE){
  return(write.table(data, file, row.names = row.names, col.names = col.names,
                     quote = FALSE, sep = ','))
}

# Read and write GeoEAS files such as used by GAM/GAMV and KT3D
read.geoeas <- function(file){
  header <- readLines(file, n = 2)
  headstr <- strsplit(header[2], ' ')[[1]]
  headstr <- headstr[headstr!='']
  ncol <- as.numeric(headstr[1])
  header2 <- readLines(file, n = 2 + ncol)
  cnames <- make.names(header2[2+(1:ncol)])
  return(read.table(file = file, col.names = cnames, skip = 2 + ncol))
}
write.geoeas <- function(data, file, title = file){
  cat(c(title, length(colnames(data)), colnames(data)), sep = '\n', file = file)
  return(write.table(data, file = file,
                     row.names = FALSE, col.names = FALSE, append = TRUE))
}
