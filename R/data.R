#' Bids Received by U.S. Firms
#'
#' A dataset with bids received by U.S. firms.
#'
#' @format A data frame with 126 rows and 13 variables: \describe{
#'   \item{docno}{doc no.} \item{weeks}{weeks} \item{numbids}{count}
#'   \item{takeover}{delta(1 if taken over)} \item{bidprem}{bid Premium}
#'   \item{insthold}{institutional holdings} \item{size}{size measured in
#'   billions} \item{leglrest}{legal restructuring} \item{rearest}{real
#'   restructuring} \item{finrest}{financial restructuring}
#'   \item{regulatn}{regulation} \item{whtknght}{white knight} }
#' @source Sanjiv Jaggia and Satish Thosar (1993) "Multiple Bids as a
#'   Consequence of Target Management Resistance", Review of Quantitative
#'   Finance and Accounting, 3(4), pp. 447-457.
#'
#'   A. Colin Cameron and Per Johansson (1997) "Count Data Regression Models using
#'   Series Expansions: with Applications", Journal of Applied Econometrics, 12,
#'   pp. 203-223.
"Bids"

#' Customer profile for a household supplies company
#'
#' An observation corresponds to the census tracts within 10-mile radius around
#' certain store.
#'
#' @format A data frame with 110 rows and 6 variables: \describe{
#'   \item{ncust}{number of customer of the census tracts who visit the store.}
#'   \item{nhu}{number of housing units in the census tracts} \item{aid}{average
#'   income in dollars} \item{aha}{average housing unit in years}
#'   \item{dnc}{distance to the nearest competitor in miles} \item{ds}{distance
#'   to store in miles} }
#' @source
#' http://www.leg.ufpr.br/lib/exe/fetch.php/publications:papercompanions:ptwdataset4.txt
#'
"CustomerProfile"
