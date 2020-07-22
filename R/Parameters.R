library("R6")

#' Class to handle model parameters
Parameters <- R6Class("Parameters", list(
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
  #' @field psubunit Probability that there's no problem at the inactivated platform level
  pinactivated=0.8,
  #' @field pvlp Probability that there's no problem at the virus-like particle platform level
  pvlp=0.8,
  #' @field ppreclinical Probability that there's no problem at the candidate level when a vaccine is in preclincal trials
  ppreclinical=0.14,
  #' @field pphase1 Probability that there's no problem at the candidate level when a vaccine is in phase 1 trials
  pphase1=0.23,
  #' @field pphase2 Probability that there's no problem at the candidate level when a vaccine is in phase 2 trials
  pphase2=0.32,
  #' @field psubcat Probability that there's no problem at subcategory level
  psubcat=0.8,
  #' @field TT Time that the vaccine is brought forward (months)
  TT=6,
  #' @field tau Period of analysis (months)
  tau=24,
  #' @field damageshare Fraction of damage accounted for by some fraction of the population
  damageshare=0.75,
  #' @field vaccshare Fraction of population that accounts for some fraction of damage
  vaccshare=1.5/7.8,
  #' @field mmonthben Economic benefit for the whole world per month (in billion $)
  monthben=375,
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
  #' @field c Cost of capacity per dose / year
  c=1.5,#2,
  #' @field capacity Constraint on capacity (doses / month)
  capacity=500,
  #' @field capkink Capacity at which marginal cost has a kink (doses / month)
  capkink=100,
  #' @field capkink Capacity at which marginal cost has a kink for RNA and DNA (doses / month)
  capkink_nucleic=NULL,
  #' @field capkink Capacity at which marginal cost has a kink for protein subunit and viral vector (doses / month)
  capkink_subunits=NULL,
  #' @field capkink Capacity at which marginal cost has a kink for live attenuated and deactivated (doses / month)
  capkink_live=NULL,
  #' @field mgcostelast Elasticity of marginal cost after kink
  mgcostelast=1,
  #' @field afterCapacity Capacity after the end of the program (doses / month)
  afterCapacity=500,
  #' @field afterCapacity Capacity in the counterfactual (doses / month)
  counterCapacity=500,
  #' @field pop Size of total population (billions)
  pop=7.8,
  #pop=3,
  #' @field cprod Production cost per vaccine ($)
  cprod=1,
  #' @field replications Number of replications for Monte Carlo
  replications=5000,
  #' @field maxcand Number of candidates to consider
  maxcand=35,
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
  #' @field fcapacity Fraction of capacity cost that must be borne by the firm
  fcapacity=0.15,
  #' @field fcostred Fraction of cost reduction after the second vaccine within a platform
  fcostred=0,
  #' @field outsidefactor How much more likely are firms to be successful in an outside market
  outsidefactor=2,
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
  #overallFung=0.3,
  overallFung=0,#0.2,
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
  #' @field hben Monthly health benefits
  hben=NA,
  #' @field damageint Integral that measures fraction of damages
  damageint=NA,
  #' @field totmonthben Total benefit of bringing forward vaccine for one month
  totmonthben=NA,
  #' @field alpha Parameter for heterogeneity of benefits of vaccination. a=1 corresponds to perfectly
  #'  homogeneous benefits, a=0 corresponds to all benefits to the first person.
  alpha=NA,
  #' @field pspike Probability that there's no problem at the candidate level when a vaccine targets spike proten
  pspike=1.0,
  #' @field precombinant Probability that there's no problem at the candidate level when a vaccine targets recombinant proten
  precombinant=1.0,
  #' @field potherprotein Probability that there's no problem at the candidate level when a vaccine targets some other proten
  potherprotein=1.0,
  
  #' @field benefitdist String describing the form of the benefit distribution
  benefitdist="pnorm",
  #' @field piecewisepar Parameters of a piecewise linear benefit distribution
  piecewisepar=NA,
  
  #' @field popshare Share of world population accounted for by coalition
  popshare=1,
  #' @field gdpshare Share of world GDP accounted for by coalition
  gdpshare=1,
  #' @field econlossratio Ratio of economic losses as a % relative to world
  econlossratio=1,
  #' @field spillovers How much do countries in the coalition care about spillovers
  spillovers=0,
  #' @field outlifefactor Extent to which countries care about international lives
  outlifefactor=0,
  #' @field outlivesratio Fraction of world that lives outside coalition
  outlivesratio=0,
  #' @field fracHighRisk Fraction of population at high risk
  fracHighRisk=0.12841,
  
  #' @field global Whether the coalition is the global coalition
  global=T,
  
  #' @field global String determining the candidate input file to use
  inputfile="Default",
  
  #' @field global Fraction of population needed to reopen economy
  fracneeded=0.7,
  #' @field global Effective world population that needs to be vaccinated to get 100% of benefits
  effpop=0.7*7.8,
  
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
  #' @return A new `Parameters` object.
  initialize = function(input="International", ...) {
    if (class(input) == "reactivevalues") {
      nms <- names(input)
      for (nm in nms) {
        if (nm %in% names(self)) {
          if (class(self[[nm]]) != "function") {
            self[[nm]] <- input[[nm]]
          }
        }
      }
    } else if (input == "US") { # Defaults for the US
      self$global=F
      self$popshare=0.33/7.8
      self$gdpshare=20.54/87.3
      self$fracHighRisk=0.2452

    } else if (input == "EU") {
      self$global=F
      self$popshare=446/7800
      self$gdpshare=16/87.3
      self$fracHighRisk=0.2677
      
    } else if (input == "EU+") { # EU + UK + Japan + S. Korea + Canada + Switzerland
      self$global=F
      self$popshare=737/7800
      self$gdpshare=27.8/87.3
      self$fracHighRisk=0.268
      # GDP: 18.7 + 5.1 + 1.7 + 1.6 + 0.7
      # Pop: 513 + 127 + 51 + 37 + 9
      
    } else if (input == "Rich") { # For US + EU + UK + Japan + Korea + Canada + Australia + Switzerland + Norway + New Zealand
      self$global=F
      self$popshare=1100/7800
      self$gdpshare=50.36/87.3
      self$fracHighRisk=0.2604
      # GDP: 20.5 + 18.7 + 5.1 + 1.7 + 1.6 + 1.43 + 0.7 + 0.43 + 0.2
      # Pop: 328 + 513 + 127 + 51 + 37 + 25 + 9 + 5.6 + 5
      
    } else if (input == "BRIC") { # Defaults for the US
      self$global=F
      self$popshare=3080/7800
      self$gdpshare=19.6/87.3
      self$fracHighRisk=0.12839
      
    } else if (input == "Rest") { # Defaults for the rest of the world
      self$global=F
      self$popshare=3620/7800
      self$gdpshare=19.4/87.3
      self$fracHighRisk=0.0869
      
    } else if (input != "International") { # The defaults are for the international version. Error message for other inputs.
      stop("Wrong name when initializing parameters")
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
    
    # if (!is.null(benefitdist)) {
    #   self$benefitdist <- benefitdist
    # }
    
    self$setCountryParameters(input)
    self$setDerivedParameters()
    
    for (i in seq_len(length(parlist))) {
      nm <- names(parlist)[i]
      if (nm %in% names(self)) {
        if (class(self[[nm]]) != "function") {
          self[[nm]] <- parlist[[nm]]
        }
      }
    }
  },

  #' @description
  #' Compute additional parameters that are functions of the input parameters
  setDerivedParameters = function() {
    if (is.null(self$mortality)) {
      self$mortality <- self$popshare * self$worldmortality
    }
    
    self$capcost <- self$c * self$capacity * 12 / 1000
    self$hben <- self$mortality * self$statLife * self$yearsLost / self$lifeExp / 1000
    self$effpop <- self$fracneeded * self$pop
    
    if (self$global) {
      self$totmonthben <- (self$monthben + self$hben) * (1 - self$sharm)
    } else {
      self$outlivesratio <- (150000 - self$mortality) / self$mortality
      
      # Commented out: for altruistic programs
      # self$totmonthben <- (self$monthben * self$gdpshare * (1+self$spillovers) + 
      #                        self$hben * (1 + self$outlifefactor * self$outlivesratio)) * (1 - self$sharm)
      # self$damageshare <- (self$monthben * self$gdpshare + self$hben) * (1 - self$sharm) / self$totmonthben
      
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
      
      # browser()
      # x <- seq(0,1,0.01)
      # y <- ifelse(x < self$vaccshare, int1 + slope1 * x,  int2 + slope2 * x)
      # ggplot() + geom_line(aes(x, y))
    }
  },
  
  setCountryParameters = function(input) {
    if (length(input) == 1) {
      if (input %in% c("US", "EU", "EU+", "Rich", "BRIC", "Rest")) {
        popratio <- self$popshare / (0.33/7.8)
        gdpratio <- self$gdpshare / (20.54/87.3)
        gdppcratio <- gdpratio / popratio
        self$statLife=7 * gdppcratio
        self$lifeExp=78.5
        self$yearsLost=10
        self$capacity=300
        
        self$inputfile <- "US"
        self$maxcand <- 19
      }
    }
  }
  
)
)

defaults <- Parameters$new()
defaultsUS <- Parameters$new("US")
