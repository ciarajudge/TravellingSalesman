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

#####################Simulated Annealing Function################
#The Simulated Annealing Function takes a table with city coordinates, and values for T and the cooling rate of T.
#It returns a table with the order of cities every niterations/10 iterations. This allows the path to be viewed as it
#developed over iterations.
simulatedannealing <- function(table, tee, cooling) {
  num_cities <- nrow(table)
  niterations <- 10000
  litmus <- c(seq(1, niterations, (niterations/10)), niterations)
  T <- tee
  pb <- progress_bar$new(total = niterations)
  pb$tick(0)
  initialorder <- sample(1:num_cities)
  order <- initialorder
  outputtable <- data.frame(matrix(ncol = 0, nrow = num_cities))
  for (i in 1:niterations) {
    pb$tick()
    T <- T*cooling
    neighbour <- adjustpath(order)
    if (distancefunction(neighbour, table) < distancefunction(order, table)) {
      order <- neighbour
    }
    else {
      if (exp(((distancefunction(neighbour, table)-distancefunction(order, table))/T)) < runif(1)) {
        order <- neighbour
      }
      else {
        order <- order
      }
    }
    if (any(litmus == i)) {
      outputtable <- cbind(outputtable, order)
    }
  }
  colnames(outputtable) <- as.character(litmus)
  return(outputtable)
}
