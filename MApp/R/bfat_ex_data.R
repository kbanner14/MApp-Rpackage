#' Body fat data for 251 men between the ages of 22 and 81 yrs.
#'
#' A dataset containing percent body fat measurements along with 13 other body
#' size measurements. In, Model Averaging: A Tutorial, Hoeting et. al. (1999) use these data in an example. We have chosen to keep only the variables used in Hoeting et. al. (1999) for this dataset.
#' @format A data frame with 251 rows and 14 variables: \describe{
#'   \item{bf.brozek}{percent body fat, calculated using Brozek's
#'          equation (457/Density - 414.2)}
#'   \item{age}{age, in years}
#'   \item{weight}{weight, in pounds}
#'   \item{height}{height, in inches}
#'   \item{neck}{neck circumference, in cm}
#'   \item{chest}{chest circumference, in cm}
#'   \item{abdomen}{circumference around abdomen "at the umbilicus
#'          and level with the iliac crest", in cm}
#'   \item{hip}{hip circumference, in cm}
#'   \item{thigh}{thigh circumference, in cm}
#'   \item{knee}{knee circumference, in cm}
#'   \item{ankle}{ankle circumference, in cm}
#'   \item{ex.bicep}{extended bicep, in cm}
#'   \item{forearm}{forearm circumference, cm}
#'   \item{wrist}{wrist circumference "distal to the
#'                 styloid processes", cm}}
#' @source Johnson RW, (1996). Fitting percentage of body fat to simple
#             body measurements. Journal of Statistics Education 4(1)
#'   \url{http://www.amstat.org/publications/jse/v4n1/datasets.johnson.html}
#'
"bfat"
