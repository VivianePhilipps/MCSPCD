#' Computation of the variability of the Health Indicators
#'
#' This function computes many iterations of the \code{estimHI} function for compute the variability health indicators.
#'
#' @param t year for the projections.
#' @param scenario 0 = pas de changement; 1 = reduction par 2 de la conso de benzo; 2 = reduction totale. Default is \code{0}.
#' @param an_scenario année de mise en place de la reduction de la consommation de benzo.
#' @param nbind nombre d'individus dont on va simuler la trajectoire pour chaque generation.
#' @param nb_iter nombre d'iterations de l'algo.
#' @param data_pop population en entree.
#' @param sexe sexe de la population en entree.
#' @param an_proj ???
#' @param data_conso prevalence de la consommation de benzo selon ag.
#' @param data_incid incidence de la consommation de benzo selon age.
#' @param a010 risques de devenir dement selon l'age.
#' @param a011 risques de devenir dement selon l'age.
#' @param a01_global risques de devenir dement selon l'age.
#' @param a020 risque de deces chez les non-dements par age et par annee.
#' @param a021 risque de deces chez les non-dements par age et par annee.
#' @param a02_global risque de deces chez les non-dements par age et par annee.
#' @param a120 risque de deces chez les dements par age et par annee.
#' @param a121 risque de deces chez les dements par age et par annee.
#' @param a12_global risque de deces chez les dements par age et par annee.
#' @param data_a01 risques de devenir dement selon l'age.
#' @param data_theta01 risque relatif de la demence selon la consommation de benzo.
#' @param data_a02 risque de deces chez les non-dements par age et par annee.
#' @param data_theta02 risque relatif de deces chez les non-dements selon la consommation de benzo.
#' @param data_theta12 risque relatif de deces chez les dements selon la consommation de benzo.
#' @param RR RR de deces pour les dements VS non-dements.
#' @param prb_dem ??
#' @param age_dem ??
#'
#' @return a list containing the variability of health indicators
#'
#' @export
#'
#' @examples
#' varHI(t = t,
#' scenario = scenario,
#' an_scenario = an_scenario,
#' nbind = nbind,
#' nb_iter = nb_iter,
#' data_pop = data_pop,
#' sexe = sexe,
#' an_proj = an_proj,
#' data_conso = data_conso,
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
varHI <- function(t, scenario, an_scenario, nbind, nb_iter, data_pop, sexe,
                  an_proj, data_conso, data_incid,
                  a010, a011, a01_global, a020, a021, a02_global, a120, a121, a12_global,
                  data_a01, data_theta01, data_a02, data_theta02, data_theta12,
                  RR, prb_dem, age_dem)

{

  #setup parallel backend to use many processors
  library(doParallel)
  cores <- detectCores()-1
  cl <- makeCluster(cores) #not to overload your computer
  registerDoParallel(cl)

  indicateurs <- foreach(it=1:nb_iter, .combine='cbind', .verbose=T, .export="multiResultClass") %dopar% { # nombre d'itérations de chaque génération

    result <- multiResultClass()

    ###############################
    ###          ETAPE 1        ###
    ###############################

    ### Recalcul de la prévalence de la conso estimée

    pr_conso_benzo <- matrix(c(0), # matrice de calcul du nombre moyen d'années de conso
                             nrow=105-65+1,
                             ncol=2,
                             byrow=T);
    pr_conso_benzo[,1] <- c(65:105);

    conso_benzo <- data_conso[which(data_conso[,3]%in%(sexe)),];

    conso_benzo <- conso_benzo[(1+41*(it-1)):(41*it),];

    incid_benzo <- data_incid[which(data_incid[,3]%in%(sexe)),];

    incid_benzo <- incid_benzo[(1+40*(it-1)):(40*it),];

    for (age in 65:105) { # suivi de chacune des générations

      an_naiss <- an_proj-age; # année de naissance

      an0 <- an_naiss + 65; # première année à risque (année des 65 ans)

      annee <- an0; # annee d'estimation

      etat <- matrix(c(0), # état initial (non dément)
                     nrow=nbind, # nombre d'individus suivis
                     ncol=105-65+1, # nombre d'années à risque
                     byrow=T);

      colnames(etat) <- c(65:105); # âges

      for (i in 1:nrow(etat)) {

        alea0 <- runif(1, 0, 1); # tirage d'un nombre aléatoire

        if (scenario == 0) {

          donnees_conso <- conso_benzo

        } else {

          if (scenario == 1) {

            if (annee < an_scenario) {

              donnees_conso <- conso_benzo

            } else {

              donnees_conso <- conso_benzo
              donnees_conso[,2] <- conso_benzo[,2] / 2

            }

          } else {

            if (scenario == 2) {

              if (annee < an_scenario) {

                donnees_conso <- conso_benzo

              } else {

                donnees_conso <- conso_benzo
                donnees_conso[,2] <- 0

              }

            }

          }

        }

        if (alea0 <= donnees_conso[which(donnees_conso[,1]%in%(65) & donnees_conso[,3]%in%(sexe)),2]) {

          etat[i,1] <- "01"

        } else {

          etat[i,1] <- "00"

        }

      };

      for (i in 1:nrow(etat)) { # pour chaque individu

        for (j in 2:ncol(etat)) { # pour chaque âge

          annee <- an0 + (j-1) # annee d'estimation

          alea <- runif(1, 0, 1); # tirage d'un nombre aléatoire pour savoir si il devient sain/dément/décès

          alea0 <- runif(1, 0, 1); # tirage d'un nombre aléatoire pour savoir si il est nouveau consommateur

          if (scenario == 0) {

            donnees_incid <- incid_benzo

          } else {

            if (scenario == 1) {

              if (annee < an_scenario) {

                donnees_incid <- incid_benzo

              } else {

                donnees_incid <- incid_benzo
                donnees_incid[,2] <- incid_benzo[,2] / 2

              }

            } else {

              if (scenario == 2) {

                if (annee < an_scenario) {

                  donnees_incid <- incid_benzo

                } else {

                  donnees_incid <- incid_benzo
                  donnees_incid[,2] <- 0

                }

              }

            }

          }

          if (etat[i,j-1] == "00") { # état année précédente est non-dément et non-consommateur

            if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(sexe)),2]) {

              a01 <- a011[(1+40*(it-1)):(40*it),];
              a02 <- a021[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "21"; # individu est devenu consommateur au cours de l'année et décède

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "11"; # individu est devenu consommateur au cours de l'année et devient dément

                } else {

                  etat[i,j] <- "01"; # individu est devenu consommateur au cours de l'année et reste non-malade

                }

              };

            } else {

              a01 <- a010[(1+40*(it-1)):(40*it),];
              a02 <- a020[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "20"; # individu est resté non-consommateur au cours de l'année et décède

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "10"; # individu est resté non-consommateur au cours de l'année et devient dément

                } else {

                  etat[i,j] <- "00"; # individu est resté non-consommateur au cours de l'année et reste non-malade

                }

              }

            }

          } else {

            if (etat[i,j-1] == "01") { # état année précédente est non-dément et consommateur

              a01 <- a011[(1+40*(it-1)):(40*it),];
              a02 <- a021[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "21"; # individu est toujours consommateur au cours de l'année et décède

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "11"; # individu est toujours consommateur au cours de l'année et devient dément

                } else {

                  etat[i,j] <- "01"; # individu est toujours consommateur au cours de l'année et reste non-malade

                }

              }

            } else {

              if (etat[i,j-1] == "10") { # état année précédente est dément et non-consommateur

                if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(sexe)),2]) {

                  a12 <- a121[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "21"; # individu est devenu consommateur au cours de l'année et décède

                  } else {

                    etat[i,j] <- "11"; # individu est devenu consommateur au cours de l'année et reste malade

                  }

                } else {

                  a12 <- a120[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "20"; # individu est resté non-consommateur au cours de l'année et décède

                  } else {

                    etat[i,j] <- "10"; # individu est resté non-consommateur au cours de l'année et reste malade

                  }

                }

              } else {

                if (etat[i,j-1] == "11") { # état année précédente est dément et consommateur

                  a12 <- a121[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "21"; # individu est toujours consommateur au cours de l'année et décède

                  } else {

                    etat[i,j] <- "11"; # individu est toujours consommateur au cours de l'année et reste malade

                  }

                } else {

                  if (etat[i,j-1] == "20") { # individu est décédé et était non-consommateur

                    etat[i,j] <- "20"; # individu est décédé et était non-consommateur

                  } else {

                    etat[i,j] <- "21"; # individu est décédé et était consommateur

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

      ### Calcul d'indicateurs :

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

    # Recalcul des intensités de transition

    ### Définition des risques de démence chez les non-consommateurs

    data_a01 <- data_a01[which(data_a01[,3]%in%(sexe)),]

    new_data_a01 <- data_a01[(1+40*(it-1)):(40*it),];

    new_data_theta01 <- data_theta01[which(data_theta01[,3]%in%(sexe)),]

    new_data_theta01 <- new_data_theta01[(1+41*(it-1)):(41*it),];

    a010[(1+40*(it-1)):(40*it),2] <- as.numeric(new_data_a01[which(new_data_a01[,1] != 65 & new_data_a01[,3]%in%(sexe)),2]) / (new_data_theta01[which(new_data_theta01[,1] != 65 & new_data_theta01[,3]%in%(sexe)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

    ### Définition des risques de démence chez les consommateurs

    a011[(1+40*(it-1)):(40*it),2] <- a010[(1+40*(it-1)):(40*it),2]*new_data_theta01[which(new_data_theta01[,1] != 65 & new_data_theta01[,3]%in%(sexe)),2];

    ### Définition des risques de démence global

    a01_global[(1+40*(it-1)):(40*it),2] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a010[(1+40*(it-1)):(40*it),2]*new_data_theta01[which(new_data_theta01[,1] != 65 & new_data_theta01[,3]%in%(sexe)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a010[(1+40*(it-1)):(40*it),2];

    ### Définition des risques de décès chez les non-consommateurs

    data_a02 <- data_a02[which(data_a02[,2]%in%(sexe)),]

    new_data_a02 <- data_a02[(1+40*(it-1)):(40*it),];

    new_data_theta02 <- data_theta02[which(data_theta02[,3]%in%(sexe)),]

    for (a in 2:ncol(a020)){ # pour chaque année

      a020[(1+40*(it-1)):(40*it),a] <- as.numeric(new_data_a02[which(new_data_a02[,1] != 65 & new_data_a02[,2]%in%(sexe)),a+1]) / (new_data_theta02[which(new_data_theta02[,1] != 65 & new_data_theta02[,3]%in%(sexe)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

    }

    ### Définition des risques de décès chez les consommateurs

    a021[(1+40*(it-1)):(40*it),-1] <- a020[(1+40*(it-1)):(40*it),-1]*new_data_theta02[which(new_data_theta02[,1] != 65 & new_data_theta02[,3]%in%(sexe)),2];

    ### Définition des risques de décès global

    a02_global[(1+40*(it-1)):(40*it),-1] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a020[(1+40*(it-1)):(40*it),-1]*new_data_theta02[which(new_data_theta02[,1] != 65 & new_data_theta02[,3]%in%(sexe)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a020[(1+40*(it-1)):(40*it),-1];

    ### Définition des risques de décès pour un dément chez les non-consommateurs

    new_data_theta12 <- data_theta12[which(data_theta12[,3]%in%(sexe)),]

    for (a in 2:ncol(a020)){ # pour chaque année

      a120[(1+40*(it-1)):(40*it),a] <- as.numeric(RR[(1+40*(it-1)):(40*it),2])*a020[(1+40*(it-1)):(40*it),a] / (new_data_theta12[which(new_data_theta12[,1] != 65 & new_data_theta12[,3]%in%(sexe)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

    }

    ### Définition des risques de décès pour un dément chez les consommateurs

    a121[(1+40*(it-1)):(40*it),-1] <- a120[(1+40*(it-1)):(40*it),-1]*new_data_theta12[which(new_data_theta12[,1] != 65 & new_data_theta12[,3]%in%(sexe)),2];

    ### Définition des risques de décès pour un dément global

    a12_global[(1+40*(it-1)):(40*it),-1] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a120[(1+40*(it-1)):(40*it),-1]*new_data_theta12[which(new_data_theta12[,1] != 65 & new_data_theta12[,3]%in%(sexe)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a120[(1+40*(it-1)):(40*it),-1];

    ###############################
    ###          ETAPE 2        ###
    ###############################

    for (age in 65:105) { # suivi de chacune des générations

      an_naiss <- an_proj-age; # année de naissance

      an0 <- an_naiss + 65; # première année à risque (année des 65 ans)

      annee <- an0; # annee d'estimation

      etat <- matrix(c(0), # état initial (non dément)
                     nrow=nbind, # nombre d'individus suivis
                     ncol=105-65+1, # nombre d'années à risque
                     byrow=T);

      colnames(etat) <- c(65:105); # âges

      for (i in 1:nrow(etat)) {

        alea0 <- runif(1, 0, 1); # tirage d'un nombre aléatoire

        if (scenario == 0) {

          donnees_conso <- pr_conso_benzo

        } else {

          if (scenario == 1) {

            if (annee < an_scenario) {

              donnees_conso <- pr_conso_benzo

            } else {

              donnees_conso <- pr_conso_benzo
              donnees_conso[,2] <- pr_conso_benzo[,2] / 2

            }

          } else {

            if (scenario == 2) {

              if (annee < an_scenario) {

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

      for (i in 1:nrow(etat)) { # pour chaque individu

        for (j in 2:ncol(etat)) { # pour chaque âge

          annee <- an0 + (j-1) # annee d'estimation

          alea <- runif(1, 0, 1); # tirage d'un nombre aléatoire pour savoir si il devient sain/dément/décès

          alea0 <- runif(1, 0, 1); # tirage d'un nombre aléatoire pour savoir si il est nouveau consommateur

          if (scenario == 0) {

            donnees_incid <- incid_benzo

          } else {

            if (scenario == 1) {

              if (annee < an_scenario) {

                donnees_incid <- incid_benzo

              } else {

                donnees_incid <- incid_benzo
                donnees_incid[,2] <- incid_benzo[,2] / 2

              }

            } else {

              if (scenario == 2) {

                if (annee < an_scenario) {

                  donnees_incid <- incid_benzo

                } else {

                  donnees_incid <- incid_benzo
                  donnees_incid[,2] <- 0

                }

              }

            }

          }

          if (etat[i,j-1] == "00") { # état année précédente est non-dément et non-consommateur

            if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(sexe)),2]) {

              a01 <- a011[(1+40*(it-1)):(40*it),];
              a02 <- a021[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "21"; # individu est devenu consommateur au cours de l'année et décède

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "11"; # individu est devenu consommateur au cours de l'année et devient dément

                } else {

                  etat[i,j] <- "01"; # individu est devenu consommateur au cours de l'année et reste non-malade

                }

              };

            } else {

              a01 <- a010[(1+40*(it-1)):(40*it),];
              a02 <- a020[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "20"; # individu est resté non-consommateur au cours de l'année et décède

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "10"; # individu est resté non-consommateur au cours de l'année et devient dément

                } else {

                  etat[i,j] <- "00"; # individu est resté non-consommateur au cours de l'année et reste non-malade

                }

              }

            }

          } else {

            if (etat[i,j-1] == "01") { # état année précédente est non-dément et consommateur

              a01 <- a011[(1+40*(it-1)):(40*it),];
              a02 <- a021[(1+40*(it-1)):(40*it),];

              if (alea <= a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "21"; # individu est toujours consommateur au cours de l'année et décède

              } else {

                if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "11"; # individu est toujours consommateur au cours de l'année et devient dément

                } else {

                  etat[i,j] <- "01"; # individu est toujours consommateur au cours de l'année et reste non-malade

                }

              }

            } else {

              if (etat[i,j-1] == "10") { # état année précédente est dément et non-consommateur

                if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(sexe)),2]) {

                  a12 <- a121[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "21"; # individu est devenu consommateur au cours de l'année et décède

                  } else {

                    etat[i,j] <- "11"; # individu est devenu consommateur au cours de l'année et reste malade

                  }

                } else {

                  a12 <- a120[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "20"; # individu est resté non-consommateur au cours de l'année et décède

                  } else {

                    etat[i,j] <- "10"; # individu est resté non-consommateur au cours de l'année et reste malade

                  }

                }

              } else {

                if (etat[i,j-1] == "11") { # état année précédente est dément et consommateur

                  a12 <- a121[(1+40*(it-1)):(40*it),];

                  if (alea <= a12[j-1,j+(an0-1950)+1]) {

                    etat[i,j] <- "21"; # individu est toujours consommateur au cours de l'année et décède

                  } else {

                    etat[i,j] <- "11"; # individu est toujours consommateur au cours de l'année et reste malade

                  }

                } else {

                  if (etat[i,j-1] == "20") { # individu est décédé et était non-consommateur

                    etat[i,j] <- "20"; # individu est décédé et était non-consommateur

                  } else {

                    etat[i,j] <- "21"; # individu est décédé et était consommateur

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

      ### Calcul d'indicateurs :

      ### Espérance de vie générale

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

      ### Espérance de vie générale chez les consommateurs de benzo

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

      ### Espérance de vie générale chez les non-consommateurs de benzo

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

      ### Espérance de vie sans la maladie

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

      ### Espérance de vie sans la maladie chez les consommateurs de benzo

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

      ### Espérance de vie sans la maladie chez les non-consommateurs de benzo

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

      ### Espérance de vie d'un malade

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

      ### Espérance de vie d'un malade chez les consommateurs de benzo

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

      ### Espérance de vie d'un malade chez les non-consommateurs de benzo

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

      ### Espérance de vie d'un non-malade

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

      ### Espérance de vie d'un non-malade chez les consommateurs de benzo

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

      ### Espérance de vie d'un non-malade chez les non-consommateurs de benzo

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

      ### Prévalence de la démence

      if (age > 65 & age < 100) {

        d0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("00") | etat[,age-65]%in%("01"))) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11")));

        s0 <- sum(etat[,1]%in%("00") | etat[,1]%in%("01") | etat[,1]%in%("10") | etat[,1]%in%("11"));

        s1 <- sum((etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11"))) + 0.5 * sum((etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) & (etat[,age-65]%in%("00") | etat[,age-65]%in%("01"))) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11")));

        if (s0 != 0) {

          p01 <- s1/s0;

          result$tp_dem[age-64] <- s1/d0;

          nb <- p01*data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(sexe)),2];

          result$np_age[age-64] <- nb;

        };

      };

      ### Survie

      if (age == 65) {

        result$nsurvie[age-64] <- data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(sexe)),2]

      };

      if (age > 65 & age < 100) {

        d0 <- sum(etat[,age-65]%in%("00") | etat[,age-65]%in%("01") | etat[,age-65]%in%("10") | etat[,age-65]%in%("11"));

        s0 <- sum(etat[,1]%in%("00") | etat[,1]%in%("01") | etat[,1]%in%("10") | etat[,1]%in%("11"));

        s1 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11"));

        if (s0 != 0) {

          p01 <- s1/s0;

          result$tsurvie[age-64] <- s1/d0;

          nb <- p01*data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(sexe)),2];

          result$nsurvie[age-64] <- nb;

        };

      };

      ### Nombre moyen d'années passées avec la maladie

      if (age < 101) {

        result$nm_dem[age-64] <- result$esp_vie_non_mal[age-64] - result$esp_vie_sans_mal[age-64]

      } else {

        n0 <- NA

      };

      ### Nombre moyen d'années passées avec la maladie chez les consommateurs de benzo

      if (age < 101) {

        result$nm_dem_conso[age-64] <- result$esp_vie_non_mal_conso[age-64] - result$esp_vie_sans_mal_conso[age-64];

      } else {

        n0 <- NA;

      };

      ### Nombre moyen d'années passées avec la maladie chez les non-consommateurs de benzo

      if (age < 101) {

        result$nm_dem_nonconso[age-64] <- result$esp_vie_non_mal_nonconso[age-64] - result$esp_vie_sans_mal_nonconso[age-64];

      } else {

        n0 <- NA;

      };

      ### Proba d'apparition de la démence

      if (age == 65) {

        for (i in (age-63):nrow(prb_dem)) {

          result$p_dem[i] <- sum(etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("10") | etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("11"))

        }

      };

      ### Age moyen d'apparition de la démence

      if (age == 65) {

        for (i in (age-63):nrow(age_dem)) {

          result$a_dem[i] <- sum(etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("10") | etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("11"))

        }

      };

      ### Nombre moyen d'années de consommation de benzo

      if (age == 65) {

        for (i in (age-63):nrow(age_dem)) {

          result$m_conso[i] <- sum(etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),i]%in%("01") | etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),i]%in%("11"))

        }

      };

      ### Prévalence de la conso estimée

      if (age <= 105) {

        s0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("10") | etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

        s1 <- sum(etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

        if (s0 != 0) {

          result$p_conso[age-64] <- s1/s0;

        } else {

          result$p_conso[age-64] <- NA;

        };

      };

      ### Quotient de mortalité estimée

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
