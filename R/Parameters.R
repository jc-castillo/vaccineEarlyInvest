library("R6")

#' Model parameters
#'
#' @description Class to handle parameters of the model for investment in early vaccine capacity
#'
#' @import data.table
#' @export
Parameters <- R6Class("Parameters", list(
  #' @field inputfile World or US input file
  inputfile="Default",
  #' @field poverall Probability that no problem at the overall level prevents vaccine feasibility
  poverall=0.9,
  #' @field pvector Probability that there's no problem at the viral vector platform level
  pvector=0.8,
  #' @field psubunit Probability that there's no problem at the protein subunit platform level
  psubunit=0.8,
  #' @field prna Probability that there's no problem at the RNA platform level
  prna=0.6,
  #' @field pdna Probability that there's no problem at the DNA platform level
  pdna=0.4,
  #' @field pattenuated Probability that there's no problem at the live attenuated platform level
  pattenuated=0.8,
  #' @field pinactivated Probability that there's no problem at the inactivated platform level
  pinactivated=0.8,
  #' @field pvlp Probability that there's no problem at the virus-like particle platform level
  pvlp=0.8,
  #' @field ppreclinical Probability that there's no problem at the candidate level when a vaccine is in preclincal trials
  ppreclinical=0.14,
  #' @field plivebac Probability that there's no problem at the live attenuated bacteria platform level
  plivebac=0,
  #' @field paapc Probability that there's no problem at the artificial antigen presenting cells platform level
  paapc=0,
  #' @field pdendritic Probability that there's no problem at the dendritic cells platform level
  pdendritic=0,
  #' @field psav Probability that there's no problem at the self-assembling vaccine platform level
  psav=0,
  #' @field punknown Probability that there's no problem for unknown platform
  punknown=0,
  #' @field prepurposed Probability that there's no problem with a repurposed candidate
  prepurposed=0,
  #' @field pphase1 Probability that there's no problem at the candidate level when a vaccine is in phase 1 trials
  pphase1=0.23,
  #' @field pphase2 Probability that there's no problem at the candidate level when a vaccine is in phase 2 trials
  pphase2=0.32,
  #' @field pphase3 Probability that there's no problem at the candidate level when a vaccine is in phase 2 trials
  pphase3=0.45,
  #' @field psubcat Probability that there's no problem at subcategory level
  psubcat=0.8,
  #' @field pspike Probability that there's no problem at the candidate level when a vaccine targets spike proten
  pspike=1.0,
  #' @field precombinant Probability that there's no problem at the candidate level when a vaccine targets recombinant proten
  precombinant=1.0,
  #' @field potherprotein Probability that there's no problem at the candidate level when a vaccine targets some other proten
  potherprotein=1.0,

  # Country / coalition specific parameters
  #' @field popshare Share of world population accounted for by country / coalition
  popshare=1,
  #' @field gdpshare Share of world GDP accounted for by country /coalition
  gdpshare=1,
  #' @field econlossratio Ratio of economic losses in the country / coalition as a \% relative to world
  econlossratio=1,
  #' @field spillovers How much do countries in the country /coalition care about spillovers
  spillovers=0,
  #' @field outlifefactor Extent to which countries care about international lives
  outlifefactor=0,
  #' @field outlivesratio Fraction of world that lives outside country /coalition
  outlivesratio=0,
  #' @field fracHighRisk Fraction of population in the country /coalition at high risk
  fracHighRisk=0.12841,

  # Parameters from model that converts production capacity into benefits
  #' @field TT Time that the vaccine is brought forward (months)
  TT=6,
  #' @field tau Period of analysis (months)
  tau=24,
  #' @field damageshare Fraction of damage accounted for by some fraction of the population
  damageshare=0.75,
  #' @field vaccshare Fraction of population that accounts for some fraction of damage
  vaccshare=1.5/7.8,
  #' @field monthben Economic benefit for the whole world per month (in billion $)
  monthben=500,
  #' @field worldmortality World mortality
  worldmortality=200000,
  #' @field mortality Monthly internal mortality
  mortality=NULL,
  #' @field statLife Value of statistical life in million $
  statLife=1.18,
  #' @field lifeExp Life expectancy
  lifeExp=71,
  #' @field yearsLost Years lost per covid death
  yearsLost=10,
  #' @field sharm Expected share of harm avoided due to an alternative drug being developed
  sharm=0.5,
  #' @field maxcand Number of candidates to consider
  maxcand=30,
  #' @field hben Monthly health benefits
  hben=NA,
  #' @field damageint Integral that measures fraction of damages
  damageint=NA,
  #' @field totmonthben Total benefit of bringing forward vaccine for one month
  totmonthben=NA,
  #' @field alpha Parameter for heterogeneity of benefits of vaccination with a smooth functional form.
  #'  a=1 corresponds to perfectly homogeneous benefits, a=0 corresponds to all benefits to the first person.
  alpha=NA,
  #' @field benefitdist String describing the form of the benefit distribution
  benefitdist="piecewiseLinearGen",
  #' @field piecewisepar Parameters of a piecewise linear benefit distribution
  piecewisepar=NA,
  #' @field fracneeded Fraction of population needed to reopen economy
  fracneeded=0.7,
  #' @field effpop Effective world population that needs to be vaccinated to get 100\% of benefits
  effpop=0.7*7.8,

  # Parameters for cost of capacity
  #' @field c Cost of capacity per dose / year
  c=2,
  #' @field capacity Constraint on capacity (doses / month)
  capacity=100,
  #' @field capkink Capacity at which marginal cost has a kink (doses / month)
  capkink=100,
  #' @field capkink_nucleic Capacity at which marginal cost has a kink for RNA and DNA (doses / month)
  capkink_nucleic=NULL,
  #' @field capkink_subunits Capacity at which marginal cost has a kink for protein subunit and viral vector (doses / month)
  capkink_subunits=NULL,
  #' @field capkink_live Capacity at which marginal cost has a kink for live attenuated and deactivated (doses / month)
  capkink_live=NULL,
  #' @field mgcostelast Elasticity of marginal cost after kink
  mgcostelast=3,
  #' @field afterCapacity Worldwide capacity after the end of the program (doses / month)
  afterCapacity=500,
  #' @field counterCapacity Worldwide capacity in the counterfactual (doses / month)
  counterCapacity=500,
  #' @field pop World population population (billions)
  pop=7.8,
  #' @field cprod Production cost per vaccine ($)
  cprod=1,
  #' @field replications Number of replications for Monte Carlo
  replications=5000,
  #' @field fcapacity Fraction of capacity cost that must be borne by the firm
  fcapacity=0.15,
  #' @field fcostred Fraction of cost reduction after the second vaccine within a platform
  fcostred=0,
  #' @field outsidefactor How much more likely are firms to be successful in an outside market
  outsidefactor=2,

  # Parameters related to procurement scheme
  #' @field pearly Price of early batch of vaccines
  pearly=35,
  #' @field plate Price of late batch of vaccines
  plate=5,
  #' @field nearly Number or vaccines in early batch
  nearly=1,
  #' @field nlate Number or vaccines in late batch
  nlate=3.5,
  #' @field riskPremium Premium that must given to candidates to participate despite risk (bn $)
  riskPremium=0.65,
  #' @field overallFung Fraction of capacity cost that is fungible across platforms
  overallFung=0,
  #' @field platFung Fraction of capacity cost that is fungible within platforms
  platFung=0,#0.2,
  #' @field subcatFung Fraction of capacity cost that is fungible within subcategory
  subcatFung=0,#0.6,
  #' @field protVectorFung Fraction of capacity cost that is fungible between protein and viral vectors
  protVectorFung=0,#0.4,
  #' @field fracScrap Fraction of cost that can be saved if a candidate is not successful
  fracScrap=0.3,
  #' @field budget Program budget
  budget=100,

  #' @field capcost Total capacity cost per candidate
  capcost=NA,
  #' @field global Whether the the analysis is for the whole wolrd
  global=F,

  # Parameters related to finding firm's incentives to participate
  #' @field outcap Capacity to build in an outside market
  outcap=0.1, # Billions
  #' @field outprob Probability that copy vaccine will be successful, conditional on own success
  outprob=0.2,
  #' @field outquantity Amount of outside vaccines that will be sold at high prices
  outquantity=0.5, # Billions

  #' @description
  #' Create a new `Parameters` object.
  #' @param input Determines how to initialize object. Can be a string telling which default parameters to use. Can also be
  #' the `input` object (of class `reactivevalues`) with the inputs from a shiny app, in which case all inputs are copied into
  #' fields.
  #' @param population Country population (in millions)
  #' @param gdp_pc Country GDP per capita (in thousand $)
  #' @param frac_high_risk Fraction of population that is high risk
  #' @param loss2yr Cumulative percent of GDP lost because of pandemic over two years
  #' @param ... Set fields at non-default values
  #' @return A new `Parameters` object.
  initialize = function(input=NULL, population=NULL, gdp_pc=NULL, frac_high_risk=NULL, loss2yr=NULL, ...) {
    if (class(input) == "reactivevalues") { # Copy all parameters if the input comes from a shiny app interface
      nms <- names(input)
      for (nm in nms) {
        if (nm %in% names(self)) {
          if (class(self[[nm]]) != "function") {
            self[[nm]] <- input[[nm]]
          }
        }
      }
    }

    self$capkink_nucleic <- self$capkink
    self$capkink_subunits <- self$capkink
    self$capkink_live <- self$capkink

    parlist <- list(...)

    for (i in seq_len(length(parlist))) {
      nm <- names(parlist)[i]
      if (nm %in% names(self)) {
        if (class(self[[nm]]) != "function") {
          self[[nm]] <- parlist[[nm]]
        }
      }
    }

    # Setting parameters related to country
    if (!is.null(population) & !is.null(gdp_pc) & !is.null(frac_high_risk)) {
      self$popshare <- population/7800
      self$gdpshare <- population*gdp_pc/1e3/87.3
      self$fracHighRisk <- frac_high_risk
      self$afterCapacity <- population/7800*500
      self$counterCapacity <- population/7800*500
    } else if (!is.null(population) | !is.null(gdp_pc) | !is.null(frac_high_risk)) {
      stop("The arguments for population, gdp_pc, and frac_high_risk must all be provided or none of them should be provided")
    }

    if (!is.null(loss2yr)) {
      self$econlossratio <- loss2yr/0.138
    }



    for (i in seq_len(length(parlist))) {
      nm <- names(parlist)[i]
      if (nm %in% names(self)) {
        if (class(self[[nm]]) != "function") {
          self[[nm]] <- parlist[[nm]]
        }
      }
    }
    
    self$setDerivedParameters()
  },

  #' @description
  #' Compute additional parameters that are functions of the input parameters
  setDerivedParameters = function() {
    if (is.null(self$mortality)) {
      self$mortality <- self$popshare * self$worldmortality
    }

    if (!self$global) {
      popratio <- self$popshare / (0.33/7.8)
      gdpratio <- self$gdpshare / (20.54/87.3)
      gdppcratio <- gdpratio / popratio
      self$statLife=7 * gdppcratio
    }

    self$capcost <- self$c * self$capacity * 12 / 1000
    self$hben <- self$mortality * self$statLife * self$yearsLost / self$lifeExp / 1000
    self$effpop <- self$fracneeded * self$pop

    if (self$global) {
      self$totmonthben <- (self$monthben + self$hben) * (1 - self$sharm)
    } else {
      self$totmonthben <- (self$econlossratio * self$monthben * self$gdpshare + self$hben) * (1 - self$sharm)
      damageratio <- 0.25 / (1-0.25) * (1-0.75) / 0.75
      self$damageshare <- 1 / (1 + (1 - self$fracHighRisk) / self$fracHighRisk * damageratio)
      self$vaccshare <- self$fracHighRisk

      self$benefitdist <- "piecewiseLinearGen"

      piecewisepar <- list(vaccshares=c(self$fracHighRisk * self$popshare / self$fracneeded, self$popshare),
                           damageshares=c(self$damageshare, 1))
      self$piecewisepar <- piecewisepar
    }

    # Compute damage parameters based on damage distribution
    if (self$benefitdist=="pnorm") {
      as <- seq(0,10,0.01)
      ys <- self$vaccshare^as + (1-self$damageshare)^as
      self$alpha <- spline(ys, as, xout=1)$y
      self$damageint <- 1-gamma(1+1/self$alpha)^2/gamma((2+self$alpha)/self$alpha)
    } else if (self$benefitdist=="piecewiseLinear") {
      slope1 <- self$damageshare / self$vaccshare
      int1 <- 0
      slope2 <- (1-self$damageshare) / (1 - self$vaccshare)
      int2 <- self$damageshare - slope2 * self$vaccshare
      self$piecewisepar <- list(slope1=slope1, int1=int1, slope2=slope2, int2=int2)
    }
  }

)
)
