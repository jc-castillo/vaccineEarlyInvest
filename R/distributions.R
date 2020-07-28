#' Net benefits for a country
#'
#' Expected net benefits for a country that invests in a certain portfolio
#'
#' @param capacities Vector of capacities for each candidate
#' @param dcandidate `data.table` with candidate information
#' @param targetPermutations `data.table` with permutations of targets
#' @param dplatforms `data.table` with platform information
#' @param par `Parameters` object with model parameters
#' @param grid Size of the capacity grid
#' @param price Price per vaccine / year
#' @param lambda Lagrange multiplier. Should be one, unless trying to solve problem with constrained budget.
#'
#' @return Expected benefits for the country
#' @export
countryNetBenefits <- function(capacities, dcandidate, targetPermutations, dplatforms, grid, price, par, lambda=1) {
  netBenefits <- countryExpectedBenefits(capacities, dcandidate, targetPermutations, dplatforms, par, grid=grid) -
    lambda * priceTakerCost(capacities, price)

  return(netBenefits)
}

#' Expected benefits for one country
#'
#' Compute expected benefits from a portfolio for one country
#'
#' @param capacities Vector of capacities for each candidate
#' @param dcandidate `data.table` with candidate information
#' @param targetPermutations `data.table` with permutations of targets
#' @param dplatforms `data.table` with platform information
#' @param par `Parameters` object with model parameters
#' @param grid Size of the capacity grid
#'
#' @return Expected benefits from vaccination
#' @export
countryExpectedBenefits <- function(capacities, dcandidate, targetPermutations, dplatforms, par, grid=1) {
  dcandidate[, capacity := capacities]
  distribution <- overallDistribution(dcandidate, targetPermutations, dplatforms,
                                      poverall=par$poverall, psubcat=par$psubcat, grid=grid)

  distribution[, progBen := benefits(par$totmonthben, list(capacity/1000, par$afterCapacity/1000),
                                     c(par$TT, par$tau), par)]
  distribution[, noProgBen := benefits(par$totmonthben, list(1e-10, par$counterCapacity/1000),
                                       c(par$TT, par$tau), par)]
  distribution[, socialBenefit := progBen - noProgBen]
  distribution[capacity==0, socialBenefit := 0]

  return(sum(distribution[, prob*socialBenefit]))
}

#' Expected benefits for one country
#'
#' Compute expected benefits from a portfolio for one country
#'
#' @param capacities Vector of capacities for each candidate
#' @param dcandidate `data.table` with candidate information
#' @param targetPermutations `data.table` with permutations of targets
#' @param dplatforms `data.table` with platform information
#' @param par `Parameters` object with model parameters
#' @param grid Size of the capacity grid
#'
#' @return Expected benefits from vaccination
#' @export
countryDistribution <- function(capacities, dcandidate, targetPermutations, dplatforms, par, grid=1) {
  dcandidate[, capacity := capacities]
  distribution <- overallDistribution(dcandidate, targetPermutations, dplatforms,
                                      poverall=par$poverall, psubcat=par$psubcat, grid=grid)

  distribution[, progBen := benefits(par$totmonthben, list(capacity/1000, par$afterCapacity/1000),
                                     c(par$TT, par$tau), par)]
  distribution[, noProgBen := benefits(par$totmonthben, list(1e-10, par$counterCapacity/1000),
                                       c(par$TT, par$tau), par)]
  distribution[, socialBenefit := progBen - noProgBen]
  distribution[capacity==0, socialBenefit := 0]

  return(distribution)
}

#' Cost for a price taker
#'
#' Computes the total cost of a portfolio for a player that is a price taker
#'
#' @param capacities Vector of capacities for all candidates
#' @param price Price per vaccine / year
#'
#' @return Total cost of the portfolio
#' @export
priceTakerCost <- function(capacities, price) {
  baseMgCost <- price * 12 / 1000
  cost <- baseMgCost * sum(capacities)
  return(cost)
}

#' Title
#'
#' @param capacities
#' @param distribution
#' @param par
#'
#' @return
#' @export
#'
#' @examples
socialCost <- function(capacities, distribution, par) {

  totcap <- sum(capacities)
  baseMgCost <- par$c * 12 / 1000
  cost <- if_else(totcap <= par$capkink,
                  baseMgCost * totcap,
                  baseMgCost * (totcap^(par$mgcostelast + 1) / par$capkink^(par$mgcostelast) +
                                  par$mgcostelast * par$capkink) / (par$mgcostelast + 1)
  )

  recover <- par$fracScrap * baseMgCost * (totcap - sum(distribution[, capacity * prob]))

  return(cost - recover)
}

#' Candidate success draws
#'
#' Generate Monte Carlo draws for the success and failure of candidates
#'
#' @param dcandidate `data.table` with information about candidates
#' @param par `Parameters` object with model parameters
#' @param seed Randome seed
#'
#' @return `data.table` with a summary of the candidates
#' @export
candidateDraws <- function(dcandidate, par, seed=30) {
  dsubcat <- dcandidate[, c("Platform", "Subcategory")]
  dsubcat <- unique(dsubcat)
  dplat <- dcandidate[, c("Platform", "pplat")]
  dplat <- unique(dplat)
  dtarget <- dcandidate[, c("Target", "ptarget")]
  dtarget <- unique(dtarget)

  # Dataset with all replications. Cartesian product of that dataset with datesets by candidate, subcategory, and platform
  ddraws <- data.table(r=1:par$replications)
  setkey(ddraws, r)
  dplatdraws <- crossJoin(ddraws, dplat)
  setkey(dplatdraws, r, Platform)
  dsubcatdraws <- crossJoin(ddraws, dsubcat)
  setkey(dsubcatdraws, r, Platform, Subcategory)
  dtargetdraws <- crossJoin(ddraws, dtarget)
  setkey(dtargetdraws, r, Target)
  dcandidate[, candInd := .I]
  dcanddraws <- crossJoin(ddraws, dcandidate)
  setkey(dcanddraws, r, Platform, Subcategory)

  # Draw main random variables
  set.seed(seed)
  ddraws[, y := rbernoulli(.N, par$poverall)]
  dplatdraws[, y := rbernoulli(.N, pplat)]
  dsubcatdraws[, y := rbernoulli(.N, par$psubcat)]
  dtargetdraws[, y := rbernoulli(.N, ptarget)]
  dcanddraws[, y := rbernoulli(.N, pcand)]

  # Input random variables into main dataset with all candidates and all draws
  dcanddraws[, ysubcand := dsubcatdraws[.(dcanddraws$r, dcanddraws$Platform, dcanddraws$Subcategory), y]]
  dcanddraws[, yplat := dplatdraws[.(dcanddraws$r, dcanddraws$Platform), y]]
  dcanddraws[, ytarget := dtargetdraws[.(dcanddraws$r, dcanddraws$Target), y]]
  dcanddraws[, yoverall := ddraws[.(dcanddraws$r), y]]
  dcanddraws[, success := yoverall * yplat * ysubcand * ytarget * y]

  return(dcanddraws)
}

#' Maximize objective function staying within a grid
#'
#' Implements an algorithm to find the optimum of the objective function without ever moving outside of a grid with
#' a fixed step size.
#'
#' @param capacities Starting point
#' @param objectiveFun Function to optimize
#' @param step Step size of the grid
#' @param verbose Whether to print extensive output (2), some output (1), or no output (0)
#' @param ... Other parameters to pass on to the objective function
#'
#' @return Optimal capacities
#' @export
optimizeGrid <- function(capacities, objectiveFun, step=100, verbose=1, ...) {
  # Setting up values for main optimization
  max <- objectiveFun(capacities, step, ...)
  end <- F
  ncand <- length(capacities)

  improvement <- rep(1e12, ncand) # Vector that keeps track of last improvement obtained when moving in each dimension.
                                  # Initialize assuming every dimension can potentially lead to a large improvement
  l <- 100000

  if (verbose==1) {
    print("Objective function: ")
  }

  # Main loop that moves capacities towards optimum
  while (!end) {
    bestdir <- -1
    best <- max

    improved <- F
    j <- 1

    # Within each step, loop over all dimensions until one of them leads to an improvement
    while (j <= ncand & !improved) {
      i <- which.max(improvement) # Choose dimension that had the largest previous improvement

      # Try increasing dimension i by one grid step
      captemp <- capacities
      captemp[i] <- captemp[i] + step
      nval <- objectiveFun(captemp, grid=step, ...)

      if (nval > best) { # Mark down if there's an improvement
        improvement[i] <- (nval - best)/step
        best <- nval
        bestdir <- i
        positive <- T
        improved <- T
      }

      if (capacities[i] - step >= 0 & !improved) {
        # If no improvement, try decreasing that dimension by one grid step.
        # Constrains analysis to positive capacities
        captemp <- capacities
        captemp[i] <- captemp[i] - step
        nval <- objectiveFun(captemp, grid=step, ...)

        if (nval > best) { # Mark down if there's an improvement
          improvement[i] <- (nval - best)/step
          best <- nval
          bestdir <- i
          positive <- F
          improved <- T
        }
      }

      if (!improved) {
        # If no improvement, assign a very small value to the improvement vector. The very small value gets smaller
        # in each step so that future steps give priority to dimensions that haven't been considered in a while.
        improvement[i] <- exp(-l / 5000)
        l <- l+1
      }

      j <- j+1
    }

    if (bestdir == -1) {
      # Stop when no dimension leads to an improvement
      end <- T
    } else {
      # Update capacities when there is an improvement
      if (positive == T) {
        capacities[bestdir] <- capacities[bestdir] + step
      } else {
        capacities[bestdir] <- capacities[bestdir] - step
      }
      max <- best
    }

    if (verbose == 2) {
      print(capacities)
      print(best)
      print(improvement, digits=3)
    } else if (verbose==1) {
      cat(best, ", ")
    }
  }

  if (verbose==1) {
    cat("\n")
  }

  return(capacities)
}


#' Target permutations
#'
#' Create data.table with information about all possible permutations of target success/failure
#'
#' @param targets vector with target names
#' @param probs vector with target probabilities
#'
#' @return data.table with all permutations and their probabilities
#' @export
getTargetPermutations <- function(targets, probs) {
  tperm <- permutations(2,length(targets),v=c(0,1),repeats.allowed=TRUE)

  perm <- vector(mode="numeric", length=nrow(tperm))
  temp <- matrix(0,nrow=nrow(tperm),ncol=ncol(tperm))
  for (i in 1:length(probs)){
    temp[,i]<- if_else(tperm[,i]==1,probs[i], 1-probs[i])
  }
  for (i in 1:length(perm)){
    perm[i]<- prod(temp[i,])
  }
  tperm <- data.table(tperm)
  tperm <- tperm[perm!=0,]
  names(tperm) <- targets
  tperm[, probability := perm[perm!=0]]
  tperm[, perm_index := seq.int(nrow(tperm))]
  tperm <- data.table(tperm)

  targetPermutations <- melt(tperm, id.vars=c("probability", "perm_index"), variable.name="Target", value.name="success")
  setkey(targetPermutations, perm_index, Target)

  return(targetPermutations)
}

#' Overall distribution of capacity
#'
#' Computes the distribution of total capacity based on data by candidate and platform
#'
#' @param dcandidate data.table with information by candidate
#' @param targetPermutations data.table with information about target permutations
#' @param dplatforms data.table with information by platform
#' @param grid Resolution of the grid to compute all distributions
#' @param target Whether to allow for target distributions (otherwise all targets succeed with probability 1)
#' @param poverall Overall probability that a vaccine is feasible
#' @param psubcat Probability that no problem shows up at the subcategory level
#'
#' @return data.table with the overall distribution
#' @export
overallDistribution <- function(dcandidate, targetPermutations, dplatforms, poverall, psubcat, grid=1, target=T) {
  # Map capacities to integers
  dcandidate[, capacity := round(capacity / grid)]

  if (!target) {
    # No need to worry about permutations if all targets are successful
    dcandidate[, pcandperm := pcand]

    overallDist <- permutationDistribution(dcandidate, dplatforms, poverall, psubcat)
  } else {

    perm_indices <- unique(targetPermutations$perm_index)
    dist <- c(0)

    # Loop over all permutations to add the distribution within each permutation
    for (i in perm_indices) {

      dcandidate[, pcandperm := pcand * targetPermutations[.(i, dcandidate$Target), success]]

      ndist <- targetPermutations[.(i), probability][1] *
        permutationDistribution(dcandidate, dplatforms, poverall, psubcat)[, prob]
      olddist <- dist

      if (length(ndist) > length(olddist)) {
        dist <- ndist
        dist[1:length(olddist)] <- dist[1:length(olddist)] + olddist
      } else {
        dist <- olddist
        dist[1:length(ndist)] <- dist[1:length(ndist)] + ndist
      }
    }

    maxcap <- max(c(which(dist>1e-14), 1)) - 1
    overallDist <- data.table(capacity=0:maxcap, prob=dist[1:(maxcap + 1)])
  }

  # Discount by overall probability
  overallDist[, prob := poverall * prob]
  overallDist[1, prob := prob + (1-poverall)]

  # Map distributions back to original space
  dcandidate[, capacity := grid * capacity]
  overallDist[, capacity := grid * capacity]

  return(overallDist)
}

#' Total capacity distribution for one permutation of platforms
#'
#' Compute the distribution of total capacity across candidates assuming one particular outcome for target successes
#'
#' @param dcandidate data.table with information by candidate
#' @param dplatforms data.table with information by platform
#' @param poverall Overall probability that a vaccine is feasible
#' @param psubcat Probability that no problem shows up at the subcategory
#'
#' @return data.table with the distribution of total capacity
#' @import matrixStats
#' @importFrom stats fft mvfft
permutationDistribution <- function(dcandidate, dplatforms, poverall, psubcat) {

  # t0 <- proc.time()
  # Compute capacity by subcategory
  subcatDists <- rbindlist(lapply(unique(dcandidate$Subcategory), subcatDistribution, dcandidate))
  setkey(subcatDists, Platform, Subcategory, capacity)

  # t1 <- proc.time()
  # Compute capacity by platform
  platDists <- rbindlist(lapply(unique(dcandidate$Platform), platformDistribution, subcatDists, psubcat))
  setkey(platDists, Platform, capacity)

  # t2 <- proc.time()
  # Create matrix where each column represents the distribution of capacity in one platform
  probTable <- dcast(platDists, capacity ~ Platform, value.var="prob")
  probTable[, capacity := NULL]
  probTableMat <- as.matrix(probTable)
  probTableMat[is.na(probTableMat)] <- 0

  for (i in 1:ncol(probTableMat)) {
    probTableMat[, i] <- dplatforms$pplat[i] * probTableMat[, i]
  }
  probTableMat[1, ] <- probTableMat[1, ] + (1-dplatforms$pplat)

  veclen <- nrow(platDists)
  nplats <- length(unique(platDists$Platform))
  distributionMat <- matrix(0, veclen, nplats)
  distributionMat[1:nrow(probTableMat), ] <- probTableMat

  # Find the distribution of the sum through the product of the convolution of the fourier transforms
  fourier <- mvfft(distributionMat)
  conv <- rowProds(fourier)
  dist <- Re(fft(conv, inverse=T))/veclen

  # Trim vector to only include nonzero values
  maxcap <- max(c(which(dist>1e-14), 1)) - 1

  # Create data.table summarizing information about distribution
  overallDist <- data.table(capacity=0:maxcap, prob=dist[1:(maxcap + 1)])
  overallDist$prob <- poverall * overallDist$prob
  overallDist[1, prob := prob + (1-poverall)]

  # t3 <- proc.time()

  # print((t1-t0)["elapsed"])
  # print((t2-t1)["elapsed"])
  # print((t3-t2)["elapsed"])

  return(overallDist)
}


#' Distribution of capacity in a platform
#'
#' Finds the distribution of the total capacity in a platform depending on distributions for subcategories in
#' the platform
#'
#' @param plat Character with the subcategory name
#' @param subcatDists Data.table with the distributions for every subcategory in the platform
#' @param psubcat Probability that no problem shows up at the subcategory
#'
#' @return A data.table summarizing the distribution of total capacity in the platfom
#' @importFrom stats fft mvfft
platformDistribution <- function(plat, subcatDists, psubcat) {
  dplat <- subcatDists[.(plat), on=.(Platform)]

  # Creating matrix where each column represents the distribution for one subcategory
  probTable <- dcast(dplat, capacity ~ Subcategory, value.var="prob")
  probTable[, capacity := NULL]
  probTableMat <- as.matrix(probTable)
  probTableMat[is.na(probTableMat)] <- 0

  for (i in 1:ncol(probTableMat)) {
    probTableMat[, i] <- psubcat * probTableMat[, i]
  }
  probTableMat[1, ] <- probTableMat[1, ] + (1-psubcat)

  veclen <- nrow(dplat)
  nsubcats <- length(unique(dplat$Subcategory))
  distributionMat <- matrix(0, veclen, nsubcats)
  distributionMat[1:nrow(probTableMat), ] <- probTableMat

  # Find the distribution of the sum through the product of the convolution of the fourier transforms
  fourier <- mvfft(distributionMat)
  conv <- rowProds(fourier)
  dist <- Re(fft(conv, inverse=T))/veclen

  # Trim vector to only include nonzero values
  maxcap <- max(c(which(dist>1e-14), 1)) - 1

  # Creating data.table with the data on the distribution
  platDist <- data.table(Platform=plat, capacity=0:maxcap, prob=dist[1:(maxcap + 1)])

  return(platDist)
}

#' Distribution of capacity in a subcategory
#'
#' Finds the distribution of the total capacity in a subcategory depending on individual candidate characteristics
#'
#' @param sub Character with the subcategory name
#' @param dcandidate Data.table with information by candidate
#'
#' @return A data.table summarizing the distribution of total capacity in the subcategory
#' @importFrom stats fft mvfft
subcatDistribution <- function(sub, dcandidate) {
  dsub <- dcandidate[Subcategory == sub]

  # Creating matrix where each column represents the distribution for one candidate
  rsub <- nrow(dsub)
  veclen <- sum(dsub$capacity) + 1

  distributionMat <- matrix(0, veclen, rsub)
  distributionMat[1, ] <- 1-dsub$pcandperm
  is <- seq_len(rsub)
  js <- dsub$capacity + 1
  ks <- veclen * (is-1) + js
  distributionMat[ks] <- distributionMat[ks] + dsub$pcandperm

  # Find the distribution of the sum through the product of the convolution of the fourier transforms
  fourier <- mvfft(distributionMat)
  conv <- rowProds(fourier)
  dist <- Re(fft(conv, inverse=T))/veclen

  # Trim vector to only include nonzero values
  maxcap <- max(c(which(dist>1e-14), 1)) - 1

  # Creating data.table with data on the distribution
  subcatDist <- data.table(Subcategory=sub, Platform=dsub$Platform[1], capacity=0:maxcap, prob=dist[1:(maxcap + 1)])

  return(subcatDist)
}

#' Sum of distributions
#'
#' Distribution of the sum of two indepenent discrete random variables
#'
#' @param dist1 First distribution. Vector where position i corresponds to the probability that the probability is i-1
#' @param dist2 Second distribution. Vector where position i corresponds to the probability that the probability is i-1
#'
#' @return DIstribution of the sum. Vector where position i corresponds to the probability that the probability is i-1
#' @importFrom stats fft mvfft
sumdist <- function(dist1, dist2) {
  l <- length(dist1)+length(dist2)-1

  d1 <- rep(0, l)
  d2 <- rep(0, l)
  d1[1:length(dist1)] <- dist1
  d2[1:length(dist2)] <- dist2

  dist <- Re(fft(fft(d1) * fft(d2), inverse=T))/l

  return(dist)
}

#' Sum of distribution with itself
#'
#' Computes the sum of the distribution of N iid discrete random variables
#'
#' @param dist Input distribution. Vector where position i corresponds to the probability that the probability is i-1
#' @param N Number of times to add
#'
#' @return Distribution of the sum. Vector where position i corresponds to the probability that the probability is i-1
#' @importFrom stats fft mvfft
sumdistSelf <- function(dist, N=2) {
  l <- N*length(dist)-(N-1)

  dd <- rep(0, l)
  dd[1:length(dist)] <- dist

  dist <- Re(fft(fft(dd)^N, inverse=T))/l

  return(dist)
}
