#' Bipartite exponential random graph models with nodal random effects
#'
#' Function to fit ERGMs with nodal random effects to bipartite networks.
#' Estimation is carried out combining the stepping algorithm
#' (see \code{\link[ergm]{control.ergm}}) for the fixed parameters and
#' pseudo-likelihood (see \code{\link[mgcv]{bam}}) for the nodal random effects.
#'
#' @param formula formula; an \code{R} formula object, of the form
#' <network> ~ <model terms> where <network> is a \code{\link[network]{network}}
#' object and <model terms> are \code{\link[ergm]{ergm.terms}}. The terms
#' \code{b1sociality} and \code{b2sociality} are added to the formula argument.
#'
#' @param iter count; number of iterations used for iterative estimation
#' of the model coefficients and the random effects, 10 iterations are usually
#' enough.
#'
#' @details This implemented algorithm can only handle undirected bipartite
#' networks.
#'
#' @author
#' Sevag Kevork
#'
#' @references Sevag Kevork, Göran Kaeurmann (2021).
#' Bipartite exponential random graph models with nodal random effects
#'
#' @import ergm
#' @import network
#' @import mgcv
#'
#' @export
#'
bimergm <- function(formula, iter){

  terms.vec <- attr(terms.formula(formula), "term.labels")

  index <- 1:length(terms.vec)

  formula.string <- paste(attr(terms(formula, keep.order = TRUE),
                               "variables")[[2]], "~",
                          paste(terms.vec[index], collapse = " + "),
                          "+ offset(b1sociality(nodes = 1:nb1)) + offset(b2sociality(nodes = 1:nb2))")

  formula <- as.formula(formula.string)

  formula.string.mple <- paste(attr(terms(formula, keep.order = TRUE),
                                    "variables")[[2]], "~",
                               paste(terms.vec[index], collapse = " + "))

  formula.mple <- as.formula(formula.string.mple)

  net <- ergm::ergm.getnetwork(formula)

  nb1 <- net$gal$bipartite
  nb2 <- net$gal$n - nb1
  n <- net$gal$n

  ndyads <- network::network.dyadcount(net)

  dta.array <- ergm::ergmMPLE(formula.mple,
                              output = "array",
                              control = control.ergm(MPLE.max.dyad.types = ndyads*2))

  ncoef <- length(dta.array$predictor[1,2,])

  dta <- matrix(0, nrow = ndyads, ncol = 3+ncoef)

  idx <- 1
  for(tail in 1:nb1){
    for(head in 1:nb2){
      dta[idx, ] <- c(dta.array$response[tail, head],
                      dta.array$predictor[tail,head, ],
                      tail,
                      head)
      idx <- idx + 1
    }
  }

  dta <- data.frame(dta)
  nm <- c("Y", names(dta.array$predictor[tail, head, ]),
          "Sociality1", "Sociality2")
  names(dta) <- nm

  dta$Sociality1 <- as.factor(dta$Sociality1)
  dta$Sociality2 <- as.factor(dta$Sociality2)

  dtabam <- mgcv::bam(Y ~ s(Sociality1, bs = "re") +
                        s(Sociality2, bs = "re"),
                      data = dta,
                      family = "binomial",
                      discrete = TRUE,
                      nthreads = 5)

  refb1 <- dtabam$coefficients[2:(1+nb1)]
  refb2 <- dtabam$coefficients[(nb1+2):(n+1)]

  results.coef.matrix <- matrix(nrow = iter, ncol = ncoef,
                                dimnames = list(c(), c(names(dta.array$predictor[tail, head, ]))))

  results.refb1 <- matrix(nrow = length(unique(dta$Sociality1)),
                          ncol = iter)
  results.refb2 <- matrix(nrow = length(unique(dta$Sociality2)),
                          ncol = iter)


  for(i in 1:iter){

    cat(paste("Iteration..."));cat("\n");cat(paste(i));cat("\n")

    model.ergm <- ergm(formula,
                       offset.coef = c(refb1, refb2),
                       control = control.ergm(main.method = "Stepping",
                                              Step.maxit = 10,
                                              Step.MCMC.samplesize = 300,
                                              Step.gridsize = 200))

    results.coef.matrix[i,] <- model.ergm$coef[1:ncoef] ################# i=1

    ## take the results and build a formula to fit the gam
    off <- rep("", times=length(ncoef))
    for(z in 1:ncoef){
      off[z] <- model.ergm$coef[z] * dta[1+z]
      off[[z]] <- as.matrix(off[[z]])
      z <- z +1
    }
    off <- do.call(cbind, off)
    off <- apply(off, 1, sum)

    model.bam <- mgcv::bam(Y ~ offset(off) + s(Sociality1, bs = "re") +
                             s(Sociality2, bs = "re"),
                           data = dta,
                           family = "binomial",
                           discrete = TRUE,
                           nthreads = 5)

    refb1 <- model.bam$coefficients[2:(1+nb1)]
    refb2 <- model.bam$coefficients[(nb1+2):(n+1)]

    results.refb1[,i] <- refb1
    results.refb2[,i] <- refb2

    i <- i + 1
  }

  out = list(formula = formula,
             iterations = iter,
             model = model.ergm,
             coefs = results.coef.matrix[iter,],
             coefs.all = results.coef.matrix,
             b1random = results.refb1,
             b2random = results.refb2)

  return(out)
}
