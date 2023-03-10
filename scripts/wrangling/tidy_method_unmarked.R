# Generate unmarked tidy method
library(broom)

# Write required tidy and glance methods
tidy.unmarkedEstimate <- 
  function(x, ...){
    submod <- paste0("[",x@short.name,"]")
    nul <- capture.output(out <- summary(x))
    names(out) <- c("estimate", "std.error", "statistic", "p.value")
    out <- cbind(term = paste(submod, rownames(out)), out)
    rownames(out) <- NULL
    tibble::as_tibble(out)
  }

tidy.unmarkedFit <- 
  function(x, submodel, ...){
    if(missing(submodel)) submodel <- names(x)
    stopifnot(all(submodel %in% names(x)))
    out <- lapply(x@estimates@estimates[submodel], tidy)
    do.call("rbind", out)
  }

glance.unmarkedFit <- 
  function(x, ...){
    tibble::tibble(AIC=x@AIC, nobs=numSites(x@data))
  }


