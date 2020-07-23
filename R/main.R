portfolioPriceTaker <- function(population, gdp_pc, frac_high_risk, loss2yr, price, steps=c(10,1,0.1), candidateFile=NULL) {
  par <- Parameters$new(global=F, inputfile="Default", maxcand=30, monthben=500, popshare=population/7800,
                         gdpshare=population*gdp_pc/1e3/87.3, fracHighRisk=frac_high_risk, afterCapacity=population/7800*500,
                         counterCapacity=population/7800*500, econlossratio=loss2yr/0.138)

  d <- loadData(par, candidateFile)

  d$Target <- "Other"
  d$Target[1:5]<-"Spike"
  d$Target[10:15]<-"Recombinant"

  dordered <- candidatesFung(d, par)$dordered

  dordered <- dordered[,1:11]
  dplatforms <- unique(dordered[, .(Platform, pplat)])
  setkey(dplatforms, Platform, pplat)
  dcandidate <- copy(dordered)

  # Getting information about permutations
  targets <- c("Spike","Recombinant","Other")
  probs <- c(as.numeric(par$pspike), as.numeric(par$precombinant), as.numeric(par$potherprotein))
  targetPermutations <- getTargetPermutations(targets, probs)

  objectiveFun <- function(capacities, grid, price, par) {
    countryNetBenefits(capacities, dcandidate, targetPermutations, dplatforms, grid, price, par, lambda=1)
  }

  capacities <- rep(0, nrow(dcandidate))

  for (s in steps) {
    capacities <- optimize(capacities, objectiveFun, step=s, price=price, par=par)
  }

  return(capacities)
}
