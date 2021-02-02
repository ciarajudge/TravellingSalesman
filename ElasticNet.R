set.seed(123)
library(dplyr)
library(stringr)
library(argparser)
library(progress)
library(RColorBrewer)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]


###########################Examine files###########################
files <- list.files("./TSP_datasets")
tsp_files <- c()
for (i in 1:length(files)) {
  files[i] <- unlist(str_split(files[i], "\\."))[1]
  tsp_files[i] <- paste0(c(files[i], ".tsp"), collapse = "") 
}
files <- unique(files)
tsp_files <- unique(tsp_files)

filelengths <- c()
finalfiles <- c()
for (file in tsp_files) {
  lines <- readLines(paste0(c("TSP_datasets/", file), collapse = ""))
  test <- T
  for (line in lines){
    if("EDGE_WEIGHT_SECTION" %in% line)
      test <- F
  }
  if (length(lines) > 70 ) {
    test <- F
  }
  if (test) {
    filelengths <- append(filelengths, length(lines))
    finalfiles <- append(finalfiles, unlist(str_split(file, "\\."))[1])
  }
}

hist(filelengths, breaks = (30))


###########################Create Maps for General Tuning####################
num_cities <- 20
map <- data.frame(ID = 1:num_cities, V2 = runif(num_cities), V3 = runif(num_cities))
map2 <- data.frame(ID = 1:num_cities, V2 = runif(num_cities), V3 = runif(num_cities))
map3 <- data.frame(ID = 1:num_cities, V2 = runif(num_cities), V3 = runif(num_cities))


###########################Import Map and Solution###########################
problemparser <- function(filename) {
  rev1 <- data.frame(read.table(paste0(c("TSP_datasets/", filename, ".tsp"), collapse = ""), header = F, sep = "", fill = TRUE, stringsAsFactors = F))
  numskip <- 0
  for (i in 1:nrow(rev1)) {
    if ("EDGE_WEIGHT_SECTION" %in% rev1[i,]) {
      stop("This problem is not in a suitable format")
    }
  }
  for (i in 1:nrow(rev1)) {
    numskip <- numskip + 1
    if ("NODE_COORD_SECTION" %in% rev1[i,]) {
      break
    }
  }
  problem <- data.frame(read.table(paste0(c("TSP_datasets/", filename, ".tsp"), collapse = ""), skip = numskip, header = F, sep = "", fill = TRUE, stringsAsFactors = F))
  problem <- problem[1:(nrow(problem)-1),]
  return(problem)
}

solutionparser <- function(filename, problem) {
  rev1 <- data.frame(read.table(paste0(c("TSP_datasets/", filename, ".opt.tour"), collapse = ""), header = F, sep = "", fill = TRUE, stringsAsFactors = F))
  numskip <- 0
  for (i in 1:nrow(rev1)) {
    numskip <- numskip + 1
    if ("TOUR_SECTION" %in% rev1[i,]) {
      break
    }
  }
  solution <- data.frame(read.table(paste0(c("TSP_datasets/", filename, ".opt.tour"), collapse = ""), skip = numskip, header = F, sep = "", fill = TRUE, stringsAsFactors = F))
  solution <- data.frame(V1 = solution[1:(nrow(problem)),])
  return(solution)
}






#Distance Function: Takes a vector with the path of the cities and a table 
#with the city coordinates, calculates distance of the path
distancefunction <- function(vector, table) {
  distances <- c()
  startx <- table[vector[1], 2]
  endx <- table[vector[length(vector)], 2]
  starty <- table[vector[1], 3]
  endy <- table[vector[length(vector)], 3]
  returntostart <- (((endx-startx)^2) + ((endy-starty)^2))^(1/2)
  distances <- append(distances, returntostart)
  for (i in 1:(length(vector)-1)){
    x1 <- table[vector[i], 2]
    x2 <- table[vector[i+1], 2]
    y1 <- table[vector[i], 3]
    y2 <- table[vector[i+1], 3]
    abdistance <- (((x2-x1)^2) + ((y2-y1)^2))^(1/2)
    distances <- append(distances, abdistance)
  }
  distance <- sum(unlist(distances))
  return(distance)
}

#Adjustpath takes an input vector with a path and randomly adjusts it. Used in SA.
adjustpath <- function(vector) {
  points <- sample(1:length(vector), 2, replace = FALSE)
  vector[points[1]:points[2]] <- vector[points[2]:points[1]]
  return(vector)
}

#swappath takes an input vector with a path and randomly swaps two cities. Used in SA.
swappath <- function(vector, length) {
  points <- sample(1:length(vector), 1, replace = FALSE)
  vector <- paste0(c(vector[points[1]:length(vector)], vector[1:(points[1]-1)]))
  vector[1:length] <- vector[length:1]
  return(vector)
}

#mutate takes a vector with a path and randomly swaps two cities n times where n = mutation rate. Used in GA.
mutate <- function(vector, rate) {
  for (i in 1:rate) {
    points <- sample(1:length(vector), 2, replace = FALSE)
    a <- vector[points[1]]
    b <- vector[points[2]]
    vector[points[1]] <- b
    vector[points[2]] <- a
    return(vector)
  }
}

#crossover takes two vectors and produces two new vectors descended 
#from a crossover event of the inputs. Used in GA.
crossover <- function(vector1, vector2) {
  point <- sample(1:length(vector1), 1)
  newvec1 <- c(vector1[0:point])
  newvec2 <- c()
  for (x in vector2){
    if (!(x %in% newvec1)) {
      newvec1 <- append(newvec1, x)
    }
    else {
      newvec2 <- append(newvec2, x)
    }
  }
  for (x in vector1) {
    if (!(x %in% newvec2)) {
      newvec2 <- append(newvec2, x)
    }
  }
  return(list(newvec1, newvec2))
}

#newpopulation takes a population of paths, their fitness and a mutation rate, and uses
#the functions crossover and mutate to create new, theoretically fitter, population. Used in GA.
newpopulation <- function(population, fitness, mutationrate) {
  newpop <- list()
  sumf <- sum(fitness)
  for (i in 1:length(population)) {
    fit <- fitness[i] / sumf
    fitness[i] <- fit 
  }
  for (i in seq(1,length(fitness), 2)) {
    rand1 <- runif(1)
    rand2 <- runif(1)
    index1 <- 0
    index2 <- 0
    while (isTRUE(rand1 > 0)) {
      index1 <- index1 + 1
      rand1 <- rand1 - fitness[index1]
    }
    while (isTRUE(rand2 > 0)) {
      index2 <- index2 + 1
      rand2 <- rand2 - fitness[index2]
    }
    newmembers <- crossover(population[[index1]], population[[index2]])
    newpop[[i]] <- mutate(newmembers[[1]], mutationrate)
    newpop[[i+1]] <- mutate(newmembers[[2]], mutationrate)
  }
  return(newpop)
}

#plottable takes an order of cities and a table with their coordinates and maps
#the path through all the cities. Used in SA, GA and EN.
plottable <- function(vector, table) {
  startx <- table[vector[1], 2]
  endx <- table[vector[length(vector)], 2]
  starty <- table[vector[1], 3]
  endy <- table[vector[length(vector)], 3]
  plot(c(startx, endx), c(starty, endy), type = "l", ylim = c(min(table$V3), max(table$V3)), 
       xlim = c(min(table$V2), max(table$V2)), xlab = "X", ylab = "Y")
  points(table$V2, table$V3, pch = 19)
  for (i in 1:(length(vector)-1)){
    x1 <- table[vector[i], 2]
    x2 <- table[vector[i+1], 2]
    y1 <- table[vector[i], 3]
    y2 <- table[vector[i+1], 3]
    lines(c(x1,x2), c(y1,y2))
  }
}

#Phi is used to create the weights table for Elastic Net using d and K.
phi <- function(d, K) {
  exp(-(d^2)/(2*(K^2)))
}

#euclidiandist takes in two vectors of xy coordinates and returns the 
#euclidian distance between them. Used in EN.
euclidiandist <- function(vec1, vec2) {
  dx <- (vec1[1] - vec2[1])^2
  dy <- (vec1[2] - vec2[2])^2
  dist <- (dx+dy)^0.5
  return(dist)
}

#Startingcircle takes the number of points (M), the angle between them, the radius of
#the starting circle and it's centerpoint and returns a dataframe with the initial coordinates
#of the net. Used in Elastic Net
startingcircle <- function(numpoints, angle, r, center_coord) {
  xj <- c()
  yj <- c()
  for (i in 1:numpoints) {
    a <- (angle*(i-1))
    xj[i] <- center_coord[1] + ((cos(a))*r)
    xj[(2*numpoints) - (i-1)] <- center_coord[1] - ((cos(a))*r)
    xj[(2*numpoints) + i] <- center_coord[1] - ((cos(a))*r)
    xj[(4*numpoints) - (i-1)] <- center_coord[1] + ((cos(a))*r)
    yj[i] <- center_coord[2] + ((sin(a))*r)
    yj[(2*numpoints) - (i-1)] <- center_coord[2] + ((sin(a))*r)
    yj[(2*numpoints) + i] <- center_coord[2] - ((sin(a))*r)
    yj[(4*numpoints) - (i-1)] <- center_coord[2] - ((sin(a))*r)
  }
  Y <- data.frame(cbind(xj, yj))
  return(Y)
}

#datasetnormaliser takes a two column table with x and y coordinates for a
#list of points and normalises the coordinates. Used in GA, SA and EN.
datasetnormaliser <- function(table) {
  for (i in 1:ncol(table)) {
    max <- max(table[,i])
    min <- min(table[,i])
    for (j in 1:nrow(table)) {
      table[j,i] <- (table[j,i] - min) / (max - min)
    }
  }
  return(table)
}

#enetdistance takes in a table with the coordinates of the net and computes
#the net length. Used in Elastic Net.
enetdistance <- function(table) {
  distances <- c()
  startx <- table[1, 1]
  endx <- table[nrow(table), 1]
  starty <- table[1, 2]
  endy <- table[nrow(table), 2]
  returntostart <- (((endx-startx)^2) + ((endy-starty)^2))^(1/2)
  distances <- append(distances, returntostart)
  for (i in 1:(nrow(table)-1)){
    x1 <- table[i,1]
    x2 <- table[i+1, 1]
    y1 <- table[i,2]
    y2 <- table[i+1,2]
    abdistance <- (((x2-x1)^2) + ((y2-y1)^2))^(1/2)
    distances <- append(distances, abdistance)
  }
  distance <- sum(unlist(distances))
  return(distance)
}

#pathfromnetgenerator takes a table with the net coordinates and a table with
#the city coordinates, and can be used to find the path order upon completion
#of net expansion. Used in EN.
pathfromnetgenerator <- function(table1, table2) {
  path <- c()
  for (p in 1:nrow(table1)) {
    eudist <- Inf
    city <- 0
    point <- table1[p,]
    for (c in 1:nrow(table2)) {
      if (euclidiandist(point, table2[c,]) < eudist) {
        eudist <- euclidiandist(point, table2[c,])
        city <- c
      }
    }
    path[p] <- city
  }
  path <- unique(path)
  return(path)
}

#The remaining functions are all used in Elastic Net with sapply to optimise the process
dxyfinder <- function(point, num_points, num_cities, P, Q, Wpq, alpha, beta, K) {
  yj <- as.numeric(as.vector(Q[point,]))
  if (point == 1){
    ybefore <- as.numeric(as.vector(Q[(num_points),]))
    yafter <- as.numeric(as.vector(Q[(point+1),]))
  }
  else if (point == num_points){
    ybefore <- as.numeric(as.vector(Q[(point-1),]))
    yafter <- as.numeric(as.vector(Q[(1),]))
  }
  else {
    ybefore <- as.numeric(as.vector(Q[(point-1),]))
    yafter <- as.numeric(as.vector(Q[(point+1),]))
  }
  sumd <- c(0,0)
  for (city in 1:num_cities) {
    xi <- as.numeric(as.vector(P[city,1:2]))
    wij <- Wpq[point, city]
    sumd <- sumd + (wij*(xi - yj))
  }
  delta <- (alpha*sumd)+(beta*K*(yafter-(2*yj)+ybefore))
}

deltafinder <- function(city, point, yj) {
  xi <- as.numeric(as.vector(P[city,1:2]))
  wij <- Wpq[point, city]
  sumd <- (wij*(xi - yj))
}

Dpqfinder <- function(city, num_points, P, Q) {
  dist <- c()
  for (point in 1:num_points) {
    dist[point] <-  euclidiandist((as.numeric(unlist(P[city,]))), (as.numeric(unlist(Q[point,]))))
  }
  return(dist)
}

Wpqfinder <- function(city, num_points, Dpq, K) {
  W <- c()
  for (point in 1:num_points){
    W[point] <- phi(Dpq[point, city], K) / sum(phi(Dpq[1:num_points,city], K))
  }
  return(W)
}


#########################Elastic Net Function###################
#The Elastic Net Function takes in a table with city coordinates, and returns a table of xy coordinates
#for the net every 500 iterations. This allows the net movement to be examined, and the final path to
#be calculated using the last two columns of the output from the elastic net function and the 
#"pathfromnetgenerator"
elasticnet <- function(table) {
  table <- cbind(table$V2, table$V3)
  table <- datasetnormaliser(table)
  num_cities <- nrow(table)
  num_points <- 2*num_cities
  alpha <- 0.3
  beta <- 2.5
  K <- 0.2
  
  quarter <- round(num_points/4)
  num_points <- 4*quarter
  center_coord <- c(mean(table[,1]), mean(table[,2]))
  angle <- 1.5708/quarter
  r <- min(c((max(table[,1])- mean(table[,1]))/3, (max(table[,2])- mean(table[,2]))/3))
  r <- 0.1
  
  P <- table
  Q <- startingcircle(quarter, angle, r, center_coord)
  Qtotal <- Q
  
  i <- 0
  pathlength <- length(pathfromnetgenerator(Q,P))
  while (pathlength != num_cities) {
    i <- i+1
    Dpq <- sapply(1:num_cities, Dpqfinder, num_points, P, Q)
    Wpq <- sapply(1:num_cities, Wpqfinder, num_points, Dpq, K)
    Qdxy <- (sapply(1:num_points, dxyfinder, num_points, num_cities, P, Q, Wpq, alpha, beta, K))
    Qdxy <- t(Qdxy)
    Q <- Q+Qdxy
    if (i %% 500 == 0 ) {
      Qtotal <- cbind(Qtotal, Q)
      print(pathlength)
    }
    if (i %% 50 == 0) {
      K <- 0.99*K
    }
    pathlength <- length(pathfromnetgenerator(Q,P))
  }
  Qtotal <- cbind(Qtotal, Q)
  print(enetdistance(Q))
  return(Qtotal)
}


