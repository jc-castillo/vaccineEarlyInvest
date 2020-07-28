#' Title
#'
#' @param capacities
#' @param par
#' @param benefitsTable
#' @param grid
#' @param dcandidate
#'
#' @return
#' @export
#'
#' @examples
globalNetBenefits <- function(capacities, dcandidate, targetPermutations, dplatforms, benefitsTable, par, grid=1) {
  dcandidate[, capacity := capacities]
  distribution <- overallDistribution(dcandidate, targetPermutations, dplatforms, grid=grid,
                                      poverall=par$poverall, psubcat=par$psubcat)

  of <- expectedBenefitsTable(distribution, benefitsTable, grid=grid) -
    socialCost(capacities, distribution, par)

  return(of)
}

#' Title
#'
#' @param distribution
#' @param benefitsTable
#' @param grid
#'
#' @return
#' @export
#'
#' @examples
expectedBenefitsTable <- function(distribution, benefitsTable, grid=1) {
  distribution[, socialBenefit := benefitsTable[.(grid * round(distribution$capacity/grid)), socialBenefit]]

  return(sum(distribution[, prob*socialBenefit]))
}


#' Title
#'
#' @param countryData
#'
#' @return
#' @export
countryParameters <- function(countryData) {


  countryPar <- function(i) {
    pop <- countryData[i, pop]
    gdp <- countryData[i, gdp]
    frac <- countryData[i, frac]
    cumulative_loss <- countryData[i, cumulative_loss]

    npar <- Parameters$new(population=pop, gdp_pc=gdp, frac_high_risk=frac, loss2yr=cumulative_loss)

    return(npar)
  }

  countryPars <- lapply(1:nrow(countryData), countryPar)

  return(countryPars)
}

#' Title
#'
#' @param countryData
#' @param max
#' @param grid
#'
#' @return
#' @export
getBenefitsTable <- function(countryData, max=1000, grid=1) {
  countryPars <- countryParameters(countryData)

  distribution <- data.table(capacity = seq(0, max, grid))

  distribution[, progBen := 0]
  distribution[, noProgBen := 0]

  GDPshares <- sapply(countryPars, function(x) {x$gdpshare})
  GDPshares <- GDPshares / sum(GDPshares)

  for (i in 1:length(countryPars)) {
    npar <- countryPars[[i]]

    distribution[, progBen := progBen +
                   benefits(npar$totmonthben,
                            list(GDPshares[i] * capacity/1000, npar$afterCapacity/1000),
                            c(npar$TT, npar$tau), npar)]
    distribution[, noProgBen := noProgBen +
                   benefits(npar$totmonthben,
                            list(1e-10, npar$counterCapacity/1000),
                            c(npar$TT, npar$tau), npar)]
  }
  distribution[, socialBenefit := progBen - noProgBen]
  distribution[capacity==0, socialBenefit := 0]

  setkey(distribution, capacity)

  return(distribution)
}


#' Load data from countries
#'
#' Loads a .xlsx file with country data, including demographics, GDP, and economic impact due to Covid-19
#'
#' @param filename File name (with path) of the xlsx file with country data
#' @param Gavi Logical, whether to treat Gavi countries as a blcok
#'
#' @return a `data.table` with country data
#' @export
#'
#' @importFrom readxl read_excel
#' @import data.table
loadCountryData <- function(filename, Gavi=F) {
  rawData <- data.table(read_excel("../../Other/highrisk_clean_bt.xlsx"))

  #drop countries without needed data (DISCUSS HOW TO COMPARE TO GLOBAL SCENARIO)
  data <- rawData[!is.na(populationtotal) & !is.na(gdp) & !is.na(frac_highrisk)]

  if (Gavi) {
    data[`GAVI eligibility`=="GAVI eligible", country := "GAVI"]
  }

  data <- data[, .(pop=sum(populationtotal)/1000000, gdp=weighted.mean(gdp, w=populationtotal)/1000,
                   frac=weighted.mean(frac_highrisk, w=populationtotal),
                   monthly_loss=weighted.mean(monthly_loss, w=gdp*populationtotal),
                   cumulative_loss=weighted.mean(cumulative_loss, w=gdp*populationtotal)), by="country"]

  setorderv(data, "gdp", order=-1)
  setnames(data, "country", "name")

  data[frac==0, frac := 1e-6]

  return(data)
}
