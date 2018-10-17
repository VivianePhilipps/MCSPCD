#' Computation of Health Indicators for dementia in France
#'
#' This function computes many health indicators under several interventions of intervention in risk factor distribution for a given year with all french data.
#'
#' @param t year of the projections for health indicators.
#' @param intervention 0 = no change; 1 = reduction by two of risk factor distribution; 2 = risk factor distribution considered as null. Default is \code{0}.
#' @param year_intervention year of the intervention in risk factor distribution takes place. Default is \code{NULL}.
#' @param nb_people number of people whose trajectory will be simulated for each generation. Default is \code{100}.
#' @param nb_iter number of iterations for the algorithm. Default is \code{0}.
#' @param gender gender for computation. \code{"W"} for women and \code{"M"} for men. Default is \code{"W"}.
#' @param a01_cst \code{TRUE} without reduction for a01 during time. \code{TRUE} with a reduction for a01 during time. Default is \code{TRUE}.
#' @param theta02_cst \code{TRUE} without increase for theta02. \code{TRUE} with an increase for theta02. Default is \code{FALSE}.
#'
#' @return a list containing the health indicators
#'
#' @export
#'
#' @examples
#' FrenchDementia(t = 2040,
#' intervention = 1,
#' year_intervention = 2020,
#' nb_people = 10000,
#' nb_iter = 100,
#' gender = "W",
#' a01_cst = TRUE,
#' theta02_cst = FALSE)
FrenchDementia <- function (t,
                            intervention = 0,
                            year_intervention = NULL,
                            nb_people = 100,
                            nb_iter = 0,
                            gender = "W",
                            a01_cst = TRUE,
                            theta02_cst = FALSE)

{
  t_FR <- t;
  intervention_FR <- intervention;
  year_intervention_FR <- year_intervention;
  nb_people_FR <- nb_people;
  nb_iter_FR <- nb_iter;
  gender_FR <- gender;

  if (a01_cst == TRUE) {
    data_for_a01_values <- a01_constant_values
    data_for_a01 <- a01_constant
  } else {
    if (a01_cst == FALSE) {
    data_for_a01_values <- a01_reduction_values
    data_for_a01 <- a01_reduction
    }
  }

  if (theta02_cst == TRUE) {
    data_for_theta02_values <- theta02_1_values
    data_for_theta02 <- theta02_1
  } else {
    if (theta02_cst == FALSE) {
      data_for_theta02_values <- theta02_increase_values
      data_for_theta02 <- theta02_increase
    }
  }

  FR_indicators <- estimHI(t = t_FR,
                           intervention = intervention_FR,
                           year_intervention = year_intervention_FR,
                           nb_people = nb_people_FR,
                           nb_iter = nb_iter_FR,
                           data_pop = pop,
                           gender = gender_FR,
                           data_a01_values = data_for_a01_values,
                           data_a02_values = a02_constant_values,
                           data_theta01_values = theta01_cas_1_6_values,
                           data_theta02_values = data_for_theta02_values,
                           data_theta12_values = data_for_theta02_values,
                           data_prev_values = prevconso_values,
                           data_incid_values = incidconso_values,
                           data_rr_DvsND_values = rr_DvsND_values,
                           data_a01 = data_for_a01,
                           data_a02 = a02_constant,
                           data_theta01 = theta01_cas_1_6,
                           data_theta02 = data_for_theta02,
                           data_theta12 = data_for_theta02,
                           data_prev = prevconso,
                           data_incid = incidconso,
                           data_rr_DvsND = rr_DvsND);

  return(FR_indicators)
}
