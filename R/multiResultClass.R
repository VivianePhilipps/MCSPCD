#' Combination of multiple results
#'
#' This function organizes the different results of the \code{varHI} function.
#'
#' @param ev_gen Default is \code{NULL}.
#' @param ev_gen_conso Default is \code{NULL}.
#' @param ev_gen_nonconso Default is \code{NULL}.
#' @param ev_sans_mal Default is \code{NULL}.
#' @param ev_sans_mal_conso Default is \code{NULL}.
#' @param ev_sans_mal_nonconso Default is \code{NULL}.
#' @param ev_mal Default is \code{NULL}.
#' @param ev_mal_conso Default is \code{NULL}.
#' @param ev_mal_nonconso Default is \code{NULL}.
#' @param ev_non_mal Default is \code{NULL}.
#' @param ev_non_mal_conso Default is \code{NULL}.
#' @param ev_non_mal_nonconso Default is \code{NULL}.
#' @param np_age Default is \code{NULL}.
#' @param tp_dem Default is \code{NULL}.
#' @param nsurvie Default is \code{NULL}.
#' @param tsurvie Default is \code{NULL}.
#' @param nm_dem Default is \code{NULL}.
#' @param nm_dem_conso Default is \code{NULL}.
#' @param nm_dem_nonconso Default is \code{NULL}.
#' @param p_dem Default is \code{NULL}.
#' @param a_dem Default is \code{NULL}.
#' @param m_conso Default is \code{NULL}.
#' @param p_conso Default is \code{NULL}.
#' @param q_mortalite Default is \code{NULL}.
#'
#' @return a list of the different results of the \code{chi_mcmc} function.
#'
#' @export
#'
#' @examples
#' multiResultClass()
multiResultClass <- function(ev_gen = NULL,
                             ev_gen_conso = NULL,
                             ev_gen_nonconso = NULL,
                             ev_sans_mal = NULL,
                             ev_sans_mal_conso = NULL,
                             ev_sans_mal_nonconso = NULL,
                             ev_mal = NULL,
                             ev_mal_conso = NULL,
                             ev_mal_nonconso = NULL,
                             ev_non_mal = NULL,
                             ev_non_mal_conso = NULL,
                             ev_non_mal_nonconso = NULL,
                             np_age = NULL,
                             tp_dem = NULL,
                             nsurvie = NULL,
                             tsurvie = NULL,
                             nm_dem = NULL,
                             nm_dem_conso = NULL,
                             nm_dem_nonconso = NULL,
                             p_dem = NULL,
                             a_dem = NULL,
                             m_conso = NULL,
                             p_conso = NULL,
                             q_mortalite = NULL)
{
  me <- list(ev_gen = ev_gen,
             ev_gen_conso = ev_gen_conso,
             ev_gen_nonconso = ev_gen_nonconso,
             ev_sans_mal = ev_sans_mal,
             ev_sans_mal_conso = ev_sans_mal_conso,
             ev_sans_mal_nonconso = ev_sans_mal_nonconso,
             ev_mal = ev_mal,
             ev_mal_conso = ev_mal_conso,
             ev_mal_nonconso = ev_mal_nonconso,
             ev_non_mal = ev_non_mal,
             ev_non_mal_conso = ev_non_mal_conso,
             ev_non_mal_nonconso = ev_non_mal_nonconso,
             np_age = np_age,
             tp_dem = tp_dem,
             nsurvie = nsurvie,
             tsurvie = tsurvie,
             nm_dem = nm_dem,
             nm_dem_conso = nm_dem_conso,
             nm_dem_nonconso = nm_dem_nonconso,
             p_dem = p_dem,
             a_dem = a_dem,
             m_conso = m_conso,
             p_conso = p_conso,
             q_mortalite = q_mortalite
  )

  return(me)

}
