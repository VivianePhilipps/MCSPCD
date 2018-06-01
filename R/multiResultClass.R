#' Combination of multiple results
#'
#' This function organizes the different results of the \code{varHI} function.
#'
#' @param LE_overall Default is \code{NULL}.
#' @param LE_overall_exp Default is \code{NULL}.
#' @param LE_overall_nonexp Default is \code{NULL}.
#' @param LE_without_dis Default is \code{NULL}.
#' @param LE_without_dis_exp Default is \code{NULL}.
#' @param LE_without_dis_nonexp Default is \code{NULL}.
#' @param LE_dis Default is \code{NULL}.
#' @param LE_dis_exp Default is \code{NULL}.
#' @param LE_dis_nonexp Default is \code{NULL}.
#' @param LE_non_dis Default is \code{NULL}.
#' @param LE_non_dis_exp Default is \code{NULL}.
#' @param LE_non_dis_nonexp Default is \code{NULL}.
#' @param np_age Default is \code{NULL}.
#' @param tp_dis Default is \code{NULL}.
#' @param nsurvival Default is \code{NULL}.
#' @param rsurvival Default is \code{NULL}.
#' @param nb_dis Default is \code{NULL}.
#' @param nb_dis_exp Default is \code{NULL}.
#' @param nb_dis_nonexp Default is \code{NULL}.
#' @param p_dis Default is \code{NULL}.
#' @param a_dis Default is \code{NULL}.
#' @param m_exp Default is \code{NULL}.
#' @param p_exp Default is \code{NULL}.
#' @param mortality_r Default is \code{NULL}.
#'
#' @return a list of the different results of the \code{varHI} function.
#'
#' @export
#'
#' @examples
#' multiResultClass()
multiResultClass <- function(LE_overall = NULL,
                             LE_overall_exp = NULL,
                             LE_overall_nonexp = NULL,
                             LE_without_dis = NULL,
                             LE_without_dis_exp = NULL,
                             LE_without_dis_nonexp = NULL,
                             LE_dis = NULL,
                             LE_dis_exp = NULL,
                             LE_dis_nonexp = NULL,
                             LE_non_dis = NULL,
                             LE_non_dis_exp = NULL,
                             LE_non_dis_nonexp = NULL,
                             np_age = NULL,
                             tp_dis = NULL,
                             nsurvival = NULL,
                             rsurvival = NULL,
                             nb_dis = NULL,
                             nb_dis_exp = NULL,
                             nb_dis_nonexp = NULL,
                             p_dis = NULL,
                             a_dis = NULL,
                             m_exp = NULL,
                             p_exp = NULL,
                             mortality_r = NULL)
{
  me <- list(LE_overall = LE_overall,
             LE_overall_exp = LE_overall_exp,
             LE_overall_nonexp = LE_overall_nonexp,
             LE_without_dis = LE_without_dis,
             LE_without_dis_exp = LE_without_dis_exp,
             LE_without_dis_nonexp = LE_without_dis_nonexp,
             LE_dis = LE_dis,
             LE_dis_exp = LE_dis_exp,
             LE_dis_nonexp = LE_dis_nonexp,
             LE_non_dis = LE_non_dis,
             LE_non_dis_exp = LE_non_dis_exp,
             LE_non_dis_nonexp = LE_non_dis_nonexp,
             np_age = np_age,
             tp_dis = tp_dis,
             nsurvival = nsurvival,
             rsurvival = rsurvival,
             nb_dis = nb_dis,
             nb_dis_exp = nb_dis_exp,
             nb_dis_nonexp = nb_dis_nonexp,
             p_dis = p_dis,
             a_dis = a_dis,
             m_exp = m_exp,
             p_exp = p_exp,
             mortality_r = mortality_r
  )

  return(me)

}
