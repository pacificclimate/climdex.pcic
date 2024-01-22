
#' Min and max indices list
#' This constant contains a list of named lists 
#' each list contains two elements: "stat" and "var." 
#' This structure allows you to access the statistical function and variable 
#' associated with each index using the index name as a key.
#' @export
climdex.min.max.idx.list <- list(
  "tnn" = list(stat = "min", var = "tmin"), "tnx" = list(stat = "max", var = "tmin"),
  "txn" = list(stat = "min", var = "tmax"), "txx" = list(stat = "max", var = "tmax"),
  "rx1day" = list(stat = "max", var = "prec"), "rx5day" = list(stat = "max", var = "prec")
)

#' Mean indices list
#' This constant contains a list of named lists 
#' each list contains two elements: "stat" and "var." 
#' This structure allows you to access the statistical function and variable 
#' associated with each index using the index name as a key.
#' @export
climdex.mean.idx.list <- list("dtr" = list(var = c("tmin", "tmax")))
#
#' Boostrap indices list
#' This constant contains a list of named lists 
#' each list contains two elements: "stat" and "var." 
#' This structure allows you to access the statistical function and variable 
#' associated with each index using the index name as a key.
#' @export
climdex.bootstrap.idx.list <- list("tn10p" = list(var = "tmin"),
                                   "tn90p" = list(var = "tmin"),
                                   "tx10p" = list(var = "tmax"),
                                   "tx90p" = list(var = "tmax"))