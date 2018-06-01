#' Computation of the variability of the Health Indicators
#'
#' This function computes many iterations of the \code{estimHI} function for compute the variability of the health indicators.
#'
#' @param t year of the projections for health indicators.
#' @param intervention 0 = no change; 1 = reduction by two of risk factor distribution; 2 = risk factor distribution considered as null. Default is \code{0}.
#' @param year_intervention year of the intervention in risk factor distribution takes place.
#' @param nb_people number of people whose trajectory will be simulated for each generation.
#' @param nb_iter number of iterations for the algorithm.
#' @param data_pop data source for demographics data.
#' @param gender gender for computation. "W" for women and "M" for men.
#' @param data_prev data source for the prevalence of the exposition.
#' @param data_incid data source for the incidence of the exposition.
#' @param a010 incidence of disease on non exposed peoples.
#' @param a011 incidence of disease on exposed peoples.
#' @param a01_global global incidence of disease.
#' @param a020 mortality of healthy subjects on non exposed peoples.
#' @param a021 mortality of healthy subjects on exposed peoples.
#' @param a02_global global mortality of healthy subjects.
#' @param a120 mortality of diseased subjects on non exposed peoples.
#' @param a121 mortality of diseased subjects on exposed peoples.
#' @param a12_global global mortality of diseased subjects.
#' @param data_a01 data source for the incidence of disease.
#' @param data_theta01 data source for the relative risks associated with the exposure for disease.
#' @param data_a02 data source for the mortality of healthy subjects.
#' @param data_theta02 data source for the relative risks associated with the exposure for mortality among healthy subjects.
#' @param data_theta12 data source for the relative risks associated with the exposure for mortality among diseased subjects.
#' @param RR relative risks associated with the disease for mortality
#' @param prb_dem life-long probability of disease
#' @param age_dem average age at disease onset
#'
#' @return a list containing the variability of health indicators
#'
#' @export
#'
#' @examples
#' varHI(t = t,
#' intervention = intervention,
#' year_intervention = year_intervention,
#' nb_people = nb_people,
#' nb_iter = nb_iter,
#' data_pop = data_pop,
#' gender = gender,
#' data_prev = data_prev,
#' data_incid = data_incid,
#' a010 = a010,
#' a011 = a011,
#' a01_global = a01_global,
#' a020 = a020,
#' a021 = a021,
#' a02_global = a02_global,
#' a120 = a120,
#' a121 = a121,
#' a12_global = a12_global,
#' data_a01 = data_a01,
#' data_theta01 = data_theta01,
#' data_a02 = data_a02,
#' data_theta02 = data_theta02,
#' data_theta12 = data_theta12,
#' RR = RR,
#' prb_dem = prb_dem,
#' age_dem = age_dem)
varHI <- function(t, intervention, year_intervention, nb_people, nb_iter, data_pop, gender,
                  data_prev, data_incid,
                  a010, a011, a01_global, a020, a021, a02_global, a120, a121, a12_global,
                  data_a01, data_theta01, data_a02, data_theta02, data_theta12,
                  RR, prb_dem, age_dem)

{

  #setup parallel backend to use many processors
  library(doParallel)
  cores <- detectCores()-1
  cl <- makeCluster(cores) #not to overload your computer
  registerDoParallel(cl)

  indicateurs <- foreach(it=1:nb_iter, .combine='cbind', .verbose=T, .export="multiResultClass") %dopar% { # number of iterations for each generation

    result <- multiResultClass()

    ###############################
    ###          STEP 1         ###
    ###############################

    ### Computation of number of exposed peoples at least one time

    pr_conso_benzo <- matrix(c(0),
                             nrow=105-65+1,
                             ncol=2,
                             byrow=T);
    pr_conso_benzo[,1] <- c(65:105);

    conso_benzo <- data_prev[which(data_prev[,3]%in%(gender)),];

    conso_benzo <- conso_benzo[(1+41*(it-1)):(41*it),];

    incid_benzo <- data_incid[which(data_incid[,3]%in%(gender)),];

    incid_benzo <- incid_benzo[(1+40*(it-1)):(40*it),];

    for (age in 65:105) {

      an_naiss <- t-age;

      an0 <- an_naiss + 65;

      annee <- an0;

      etat <- matrix(c(0),
                     nrow=nb_people,
                     ncol=105-65+1,
                     byrow=T);

      colnames(etat) <- c(65:105);

      for (i in 1:nrow(etat)) {

        alea0 <- runif(1, 0, 1);

        if (intervention == 0) {

          donnees_conso <- conso_benzo

        } else {

          if (intervention == 1) {

            if (annee < year_intervention) {

              donnees_conso <- conso_benzo

            } else {

              donnees_conso <- conso_benzo
              donnees_conso[,2] <- conso_benzo[,2] / 2

            }

          } else {

            if (intervention == 2) {

              if (annee < year_intervention) {

                donnees_conso <- conso_benzo

              } else {

                donnees_conso <- conso_benzo
                donnees_conso[,2] <- 0

              }

            }

          }

        }

        if (alea0 <= donnees_conso[which(donnees_conso[,1]%in%(65) & donnees_conso[,3]%in%(gender)),2]) {

          etat[i,1] <- "01"

        } else {

          etat[i,1] <- "00"

        }

      };

      for (i in 1:nrow(etat)) {

        for (j in 2:ncol(etat)) {

          annee <- an0 + (j-1)

          alea <- runif(1, 0, 1);

          alea0 <- runif(1, 0, 1);

          if (intervention == 0) {

            donnees_incid <- incid_benzo

          } else {

            if (intervention == 1) {

              if (annee < year_intervention) {

                donnees_incid <- incid_benzo

              } else {

                donnees_incid <- incid_benzo
                donnees_incid[,2] <- incid_benzo[,2] / 2

              }

            } else {

              if (intervention == 2) {

                if (annee < year_intervention) {

                  donnees_incid <- incid_benzo

                } else {

                  donnees_incid <- incid_benzo
                  donnees_incid[,2] <- 0

                }

              }

            }

          }

          if (etat[i,j-1] == "00") {

            if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(gender)),2]) {

              a01 <- a011[(1+40*(it-1)):(40*it),];
              a02 <- a021[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "21";

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "11";

                } else {

                  etat[i,j] <- "01";

                }

              };

            } else {

              a01 <- a010[(1+40*(it-1)):(40*it),];
              a02 <- a020[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "20";

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "10";

                } else {

                  etat[i,j] <- "00";

                }

              }

            }

          } else {

            if (etat[i,j-1] == "01") {

              a01 <- a011[(1+40*(it-1)):(40*it),];
              a02 <- a021[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "21";

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "11";

                } else {

                  etat[i,j] <- "01";

                }

              }

            } else {

              if (etat[i,j-1] == "10") {

                if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(gender)),2]) {

                  a12 <- a121[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "21";

                  } else {

                    etat[i,j] <- "11";

                  }

                } else {

                  a12 <- a120[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "20";

                  } else {

                    etat[i,j] <- "10";

                  }

                }

              } else {

                if (etat[i,j-1] == "11") {

                  a12 <- a121[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "21";

                  } else {

                    etat[i,j] <- "11";

                  }

                } else {

                  if (etat[i,j-1] == "20") {

                    etat[i,j] <- "20";

                  } else {

                    etat[i,j] <- "21";

                  }

                }

              }

            }

          }

        }

      };

      #####################################################################################
      #####################################################################################
      #####################################################################################
      #####################################################################################
      #####################################################################################

      ### Computation of number of exposed peoples at least one time

      if (age <= 105) {

        s0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("10") | etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

        s1 <- sum(etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

        if (s0 != 0) {

          pr_conso_benzo[age-64,2] <- s1/s0;

        } else {

          pr_conso_benzo[age-64,2] <- pr_conso_benzo[age-64-1,2];

        };

      };

    }

    ### Computation of transition intensities

    ### Incidence of disease on non exposed peoples

    data_a01 <- data_a01[which(data_a01[,3]%in%(gender)),]

    new_data_a01 <- data_a01[(1+40*(it-1)):(40*it),];

    new_data_theta01 <- data_theta01[which(data_theta01[,3]%in%(gender)),]

    new_data_theta01 <- new_data_theta01[(1+41*(it-1)):(41*it),];

    a010[(1+40*(it-1)):(40*it),2] <- as.numeric(new_data_a01[which(new_data_a01[,1] != 65 & new_data_a01[,3]%in%(gender)),2]) / (new_data_theta01[which(new_data_theta01[,1] != 65 & new_data_theta01[,3]%in%(gender)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

    ### Incidence of disease on exposed peoples

    a011[(1+40*(it-1)):(40*it),2] <- a010[(1+40*(it-1)):(40*it),2]*new_data_theta01[which(new_data_theta01[,1] != 65 & new_data_theta01[,3]%in%(gender)),2];

    ### Global incidence of disease

    a01_global[(1+40*(it-1)):(40*it),2] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a010[(1+40*(it-1)):(40*it),2]*new_data_theta01[which(new_data_theta01[,1] != 65 & new_data_theta01[,3]%in%(gender)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a010[(1+40*(it-1)):(40*it),2];

    ### Mortality of healthy subjects on non exposed peoples

    data_a02 <- data_a02[which(data_a02[,2]%in%(gender)),]

    new_data_a02 <- data_a02[(1+40*(it-1)):(40*it),];

    new_data_theta02 <- data_theta02[which(data_theta02[,3]%in%(gender)),]

    for (a in 2:ncol(a020)){ # pour chaque année

      a020[(1+40*(it-1)):(40*it),a] <- as.numeric(new_data_a02[which(new_data_a02[,1] != 65 & new_data_a02[,2]%in%(gender)),a+1]) / (new_data_theta02[which(new_data_theta02[,1] != 65 & new_data_theta02[,3]%in%(gender)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

    }

    ### Mortality of healthy subjects on exposed peoples

    a021[(1+40*(it-1)):(40*it),-1] <- a020[(1+40*(it-1)):(40*it),-1]*new_data_theta02[which(new_data_theta02[,1] != 65 & new_data_theta02[,3]%in%(gender)),2];

    ### Global mortality of healthy subjects

    a02_global[(1+40*(it-1)):(40*it),-1] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a020[(1+40*(it-1)):(40*it),-1]*new_data_theta02[which(new_data_theta02[,1] != 65 & new_data_theta02[,3]%in%(gender)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a020[(1+40*(it-1)):(40*it),-1];

    ### Mortality of diseased subjects on non exposed peoples

    new_data_theta12 <- data_theta12[which(data_theta12[,3]%in%(gender)),]

    for (a in 2:ncol(a020)){ # pour chaque année

      a120[(1+40*(it-1)):(40*it),a] <- as.numeric(RR[(1+40*(it-1)):(40*it),2])*a020[(1+40*(it-1)):(40*it),a] / (new_data_theta12[which(new_data_theta12[,1] != 65 & new_data_theta12[,3]%in%(gender)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

    }

    ### Mortality of diseased subjects on exposed peoples

    a121[(1+40*(it-1)):(40*it),-1] <- a120[(1+40*(it-1)):(40*it),-1]*new_data_theta12[which(new_data_theta12[,1] != 65 & new_data_theta12[,3]%in%(gender)),2];

    ### Global mortality of diseased subjects

    a12_global[(1+40*(it-1)):(40*it),-1] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a120[(1+40*(it-1)):(40*it),-1]*new_data_theta12[which(new_data_theta12[,1] != 65 & new_data_theta12[,3]%in%(gender)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a120[(1+40*(it-1)):(40*it),-1];

    ###############################
    ###          STEP 2         ###
    ###############################

    for (age in 65:105) {

      an_naiss <- t-age;

      an0 <- an_naiss + 65;

      annee <- an0;

      etat <- matrix(c(0),
                     nrow=nb_people,
                     ncol=105-65+1,
                     byrow=T);

      colnames(etat) <- c(65:105);

      for (i in 1:nrow(etat)) {

        alea0 <- runif(1, 0, 1);

        if (intervention == 0) {

          donnees_conso <- pr_conso_benzo

        } else {

          if (intervention == 1) {

            if (annee < year_intervention) {

              donnees_conso <- pr_conso_benzo

            } else {

              donnees_conso <- pr_conso_benzo
              donnees_conso[,2] <- pr_conso_benzo[,2] / 2

            }

          } else {

            if (intervention == 2) {

              if (annee < year_intervention) {

                donnees_conso <- pr_conso_benzo

              } else {

                donnees_conso <- pr_conso_benzo
                donnees_conso[,2] <- 0

              }

            }

          }

        }

        if (alea0 <= donnees_conso[which(donnees_conso[,1]%in%(65)),2]) {

          etat[i,1] <- "01"

        } else {

          etat[i,1] <- "00"

        }

      };

      for (i in 1:nrow(etat)) {

        for (j in 2:ncol(etat)) {

          annee <- an0 + (j-1)

          alea <- runif(1, 0, 1);

          alea0 <- runif(1, 0, 1);

          if (intervention == 0) {

            donnees_incid <- incid_benzo

          } else {

            if (intervention == 1) {

              if (annee < year_intervention) {

                donnees_incid <- incid_benzo

              } else {

                donnees_incid <- incid_benzo
                donnees_incid[,2] <- incid_benzo[,2] / 2

              }

            } else {

              if (intervention == 2) {

                if (annee < year_intervention) {

                  donnees_incid <- incid_benzo

                } else {

                  donnees_incid <- incid_benzo
                  donnees_incid[,2] <- 0

                }

              }

            }

          }

          if (etat[i,j-1] == "00") {

            if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(gender)),2]) {

              a01 <- a011[(1+40*(it-1)):(40*it),];
              a02 <- a021[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "21";

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "11";

                } else {

                  etat[i,j] <- "01";

                }

              };

            } else {

              a01 <- a010[(1+40*(it-1)):(40*it),];
              a02 <- a020[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "20";

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "10";

                } else {

                  etat[i,j] <- "00";

                }

              }

            }

          } else {

            if (etat[i,j-1] == "01") {

              a01 <- a011[(1+40*(it-1)):(40*it),];
              a02 <- a021[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "21";

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "11";

                } else {

                  etat[i,j] <- "01";

                }

              }

            } else {

              if (etat[i,j-1] == "10") {

                if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(gender)),2]) {

                  a12 <- a121[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "21";

                  } else {

                    etat[i,j] <- "11";

                  }

                } else {

                  a12 <- a120[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "20";

                  } else {

                    etat[i,j] <- "10";

                  }

                }

              } else {

                if (etat[i,j-1] == "11") {

                  a12 <- a121[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "21";

                  } else {

                    etat[i,j] <- "11";

                  }

                } else {

                  if (etat[i,j-1] == "20") {

                    etat[i,j] <- "20";

                  } else {

                    etat[i,j] <- "21";

                  }

                }

              }

            }

          }

        }

      };

      #####################################################################################
      #####################################################################################
      #####################################################################################
      #####################################################################################
      #####################################################################################

      ### Computation of health indicators :

      ### Overall life-expectancy

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[,j]%in%("00") | etat[,j]%in%("01") | etat[,j]%in%("10") | etat[,j]%in%("11"))

          };

          result$ev_gen[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_gen[age-64] <- NA;

        }

      };

      ###	Overall life-expectancy on exposed peoples

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("01") | etat[,age-64]%in%("11")),age-64]%in%("01") | etat[which(etat[,age-64]%in%("01") | etat[,age-64]%in%("11")),age-64]%in%("11"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("01") | etat[,age-64]%in%("11")),j]%in%("01") | etat[which(etat[,age-64]%in%("01") | etat[,age-64]%in%("11")),j]%in%("11"))

          };

          result$ev_gen_conso[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_gen_conso[age-64] <- NA;

        }

      };

      ###	Overall life-expectancy on non exposed peoples

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("10")),age-64]%in%("00") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("10")),age-64]%in%("10"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("10")),j]%in%("00") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("10")),j]%in%("10") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("10")),j]%in%("01") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("10")),j]%in%("11"))

          };

          result$ev_gen_nonconso[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_gen_nonconso[age-64] <- NA;

        }

      };

      ###	Life-expectancy without disease

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),age-64]%in%("00") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),age-64]%in%("01"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),j]%in%("00") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),j]%in%("01"))

          };

          result$ev_sans_mal[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_sans_mal[age-64] <- NA;

        }

      };

      ###	Life-expectancy without disease on exposed peoples

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("01")),age-64]%in%("01"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("01")),j]%in%("01"))

          };

          result$ev_sans_mal_conso[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_sans_mal_conso[age-64] <- NA;

        }

      };

      ###	Life-expectancy without disease on non exposed peoples

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("00")),age-64]%in%("00"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("00")),j]%in%("00") | etat[which(etat[,age-64]%in%("00")),j]%in%("01"))

          };

          result$ev_sans_mal_nonconso[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_sans_mal_nonconso[age-64] <- NA;

        }

      };

      ### Life-expectancy for diseased subject

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("10") | etat[,age-64]%in%("11")),age-64]%in%("10") | etat[which(etat[,age-64]%in%("10") | etat[,age-64]%in%("11")),age-64]%in%("11"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("10") | etat[,age-64]%in%("11")),j]%in%("10") | etat[which(etat[,age-64]%in%("10") | etat[,age-64]%in%("11")),j]%in%("11"))

          };

          result$ev_mal[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_mal[age-64] <- NA;

        }

      };

      ### Life-expectancy for diseased subject on exposed peoples

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("11")),age-64]%in%("11"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("11")),j]%in%("11"))

          };

          result$ev_mal_conso[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_mal_conso[age-64] <- NA;

        }

      };

      ### Life-expectancy for diseased subject on non exposed peoples

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("10")),age-64]%in%("10"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("10")),j]%in%("10") | etat[which(etat[,age-64]%in%("10")),j]%in%("11"))

          };

          result$ev_mal_nonconso[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_mal_nonconso[age-64] <- NA;

        }

      };

      ### Life-expectancy for non diseased subject

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),age-64]%in%("00") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),age-64]%in%("01"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),j]%in%("00") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),j]%in%("01") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),j]%in%("10") | etat[which(etat[,age-64]%in%("00") | etat[,age-64]%in%("01")),j]%in%("11"))

          };

          result$ev_non_mal[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_non_mal[age-64] <- NA;

        }

      };

      ### Life-expectancy for non diseased subject on exposed peoples

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("01")),age-64]%in%("01"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("01")),j]%in%("01") | etat[which(etat[,age-64]%in%("01")),j]%in%("11"))

          };

          result$ev_non_mal_conso[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_non_mal_conso[age-64] <- NA;

        }

      };

      ### Life-expectancy for non diseased subject on non exposed peoples

      if (age < 101) {

        n0 <- vector(length = ncol(etat));

        s0 <- sum(etat[which(etat[,age-64]%in%("00")),age-64]%in%("00"));

        if (s0 != 0) {

          for (j in (age-63):ncol(etat)) {

            n0[j] <- sum(etat[which(etat[,age-64]%in%("00")),j]%in%("00") | etat[which(etat[,age-64]%in%("00")),j]%in%("01") | etat[which(etat[,age-64]%in%("00")),j]%in%("10") | etat[which(etat[,age-64]%in%("00")),j]%in%("11"))

          };

          result$ev_non_mal_nonconso[age-64] <- 1 + sum(n0) / s0;

        } else {

          n0 <- NA;

          result$ev_non_mal_nonconso[age-64] <- NA;

        }

      };

      ### Prevalence of disease

      if (age > 65 & age < 100) {

        d0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("00") | etat[,age-65]%in%("01"))) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11")));

        s0 <- sum(etat[,1]%in%("00") | etat[,1]%in%("01") | etat[,1]%in%("10") | etat[,1]%in%("11"));

        s1 <- sum((etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11"))) + 0.5 * sum((etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) & (etat[,age-65]%in%("00") | etat[,age-65]%in%("01"))) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11")));

        if (s0 != 0) {

          p01 <- s1/s0;

          result$tp_dem[age-64] <- s1/d0;

          nb <- p01*data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(gender)),2];

          result$np_age[age-64] <- nb;

        };

      };

      ### Survival

      if (age == 65) {

        result$nsurvie[age-64] <- data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(gender)),2]

      };

      if (age > 65 & age < 100) {

        d0 <- sum(etat[,age-65]%in%("00") | etat[,age-65]%in%("01") | etat[,age-65]%in%("10") | etat[,age-65]%in%("11"));

        s0 <- sum(etat[,1]%in%("00") | etat[,1]%in%("01") | etat[,1]%in%("10") | etat[,1]%in%("11"));

        s1 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11"));

        if (s0 != 0) {

          p01 <- s1/s0;

          result$tsurvie[age-64] <- s1/d0;

          nb <- p01*data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(gender)),2];

          result$nsurvie[age-64] <- nb;

        };

      };

      ### Mean number of years spent with disease

      if (age < 101) {

        result$nm_dem[age-64] <- result$esp_vie_non_mal[age-64] - result$esp_vie_sans_mal[age-64]

      } else {

        n0 <- NA

      };

      ### Mean number of years spent with disease on exposed peoples

      if (age < 101) {

        result$nm_dem_conso[age-64] <- result$esp_vie_non_mal_conso[age-64] - result$esp_vie_sans_mal_conso[age-64];

      } else {

        n0 <- NA;

      };

      ###  Mean number of years spent with disease on non exposed peoples

      if (age < 101) {

        result$nm_dem_nonconso[age-64] <- result$esp_vie_non_mal_nonconso[age-64] - result$esp_vie_sans_mal_nonconso[age-64];

      } else {

        n0 <- NA;

      };

      ### Life-long probability of disease

      if (age == 65) {

        for (i in (age-63):nrow(prb_dem)) {

          result$p_dem[i] <- sum(etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("10") | etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("11"))

        }

      };

      ### Average age at disease onset

      if (age == 65) {

        for (i in (age-63):nrow(age_dem)) {

          result$a_dem[i] <- sum(etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("10") | etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("11"))

        }

      };

      ### Mean number of years of exposition

      if (age == 65) {

        for (i in (age-63):nrow(age_dem)) {

          result$m_conso[i] <- sum(etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),i]%in%("01") | etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),i]%in%("11"))

        }

      };

      ### Number of exposed peoples at least one time

      if (age <= 105) {

        s0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("10") | etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

        s1 <- sum(etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

        if (s0 != 0) {

          result$p_conso[age-64] <- s1/s0;

        } else {

          result$p_conso[age-64] <- NA;

        };

      };

      ### Mortality rate

      if (age > 65 & age < 100) {

        s0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("10") | etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

        s1 <- sum(etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) - sum(etat[,age-65]%in%("20") | etat[,age-65]%in%("21"));

        if (s0 != 0) {

          result$q_mortalite[age-64] <- s1/s0;

        } else {

          result$q_mortalite[age-64] <- NA;

        };

      };

    }

    return(result)

  }

  #stop cluster
  stopCluster(cl)

  return(indicateurs)

}
