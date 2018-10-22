#' Computation of Health Indicators for dementia in France
#'
#' This function computes many health indicators without risk factor on transitions for a given year with all french data.
#'
#' @param t year of the projections for health indicators.
#' @param nb_people number of people whose trajectory will be simulated for each generation. Default is \code{100}.
#' @param nb_iter number of iterations for the algorithm. Default is \code{0}.
#' @param gender gender for computation. \code{"W"} for women and \code{"M"} for men. Default is \code{"W"}.
#'
#' @return a list containing the health indicators
#'
#' @export
#'
#' @examples
#' estimDementia(t = 2040,
#' nb_people = 10000,
#' nb_iter = 100,
#' gender = "W")
estimDementia <- function (t,
                            nb_people = 100,
                            nb_iter = 0,
                            gender = "W")

{
  t_FR <- t;
  nb_people_FR <- nb_people;
  nb_iter_FR <- nb_iter;
  gender_FR <- gender;

  theta01_values_1 <- theta01_cas_1_6_values
  theta02_values_1 <- theta02_increase_values
  theta12_values_1 <- theta02_increase_values

  theta01_values_1[,2] <- 1
  theta02_values_1[,2] <- 1
  theta12_values_1[,2] <- 1

  theta01_1 <- theta01_cas_1_6
  theta02_1 <- theta02_increase
  theta12_1 <- theta02_increase

  theta01_1[,2] <- 1
  theta02_1[,2] <- 1
  theta12_1[,2] <- 1

  FR_indicators <- estimHI(t = t_FR,
                           intervention = 0,
                           year_intervention = NULL,
                           nb_people = nb_people_FR,
                           nb_iter = nb_people_FR,
                           data_pop = pop,
                           gender = gender_FR,
                           data_a01_values = a01_constant_values,
                           data_a02_values = a02_constant_values,
                           data_theta01_values = theta01_values_1,
                           data_theta02_values = theta02_values_1,
                           data_theta12_values = theta12_values_1,
                           data_prev_values = prevconso_values,
                           data_incid_values = incidconso_values,
                           data_rr_DvsND_values = rr_DvsND_values,
                           data_a01 = a01_constant,
                           data_a02 = a02_constant,
                           data_theta01 = theta01_1,
                           data_theta02 = theta02_1,
                           data_theta12 = theta02_1,
                           data_prev = prevconso,
                           data_incid = incidconso,
                           data_rr_DvsND = rr_DvsND);

  return(FR_indicators)
}
