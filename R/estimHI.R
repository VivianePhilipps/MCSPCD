#' Computation of Health Indicators
#'
#' This function computes many health indicators under several scenarios of intervention in benzodiazepines chronic consumption.
#'
#' @param t year for the projections.
#' @param scenario 0 = pas de changement; 1 = reduction par 2 de la conso de benzo; 2 = reduction totale. Default is \code{0}.
#' @param an_scenario année de mise en place de la reduction de la consommation de benzo.
#' @param nbind nombre d'individus dont on va simuler la trajectoire pour chaque generation.
#' @param nb_iter nombre d'iterations de l'algo.
#' @param data_pop population en entree.
#' @param sexe sexe de la population en entree.
#' @param data_a01 risques de devenir dement selon l'age.
#' @param data_a02 risque de deces chez les non-dements par age et par annee.
#' @param data_theta01 risque relatif de la demence selon la consommation de benzo.
#' @param data_theta02 risque relatif de deces chez les non-dements selon la consommation de benzo.
#' @param data_theta12 risque relatif de deces chez les dements selon la consommation de benzo.
#' @param data_conso prevalence de la consommation de benzo selon ag.
#' @param data_incid incidence de la consommation de benzo selon age.
#' @param data_rr_DvsND RR de deces pour les dements VS non-dements.
#' @param data_a01_valeurs variability of risques de devenir dement selon l'age.
#' @param data_a02_valeurs variability of risque de deces chez les non-dements par age et par annee.
#' @param data_theta01_valeurs variability of risque relatif de la demence selon la consommation de benzo.
#' @param data_conso_valeurs variability of prevalence de la consommation de benzo selon ag.
#' @param data_incid_valeurs variability of incidence de la consommation de benzo selon age.
#' @param data_rr_DvsND_valeurs variability of RR de deces pour les dements VS non-dements.
#'
#' @return a list containing the health indicators
#'
#' @export
#'
#' @examples
#' estimHI(t = 2040,
#' scenario = 0,
#' an_scenario = 2020,
#' nbind = 10000,
#' nb_iter = 100,
#' data_pop = pop,
#' sexe = "F",
#' data_a01 = a01,
#' data_a02 = a02,
#' data_theta01 = theta01_cas_1_6,
#' data_theta02 = theta02_1,
#' data_theta12 = theta02_1,
#' data_conso = prevconso,
#' data_incid = incidconso,
#' data_rr_DvsND = rr_DvsND,
#' data_a01_valeurs = a01_valeurs,
#' data_a02_valeurs = a02_valeurs,
#' data_theta01_valeurs = theta01_cas_1_6_valeurs,
#' data_conso_valeurs <- prevconso_valeurs,
#' data_incid_valeurs <- incidconso_valeurs,
#' data_rr_DvsND_valeurs = rr_DvsND_valeurs)
estimHI <- function(t, scenario, an_scenario, nbind, nb_iter, data_pop, sexe,
                    data_a01, data_a02, data_theta01, data_theta02, data_theta12,
                    data_conso, data_incid, data_rr_DvsND, data_a01_valeurs, data_a02_valeurs,
                    data_theta01_valeurs, data_conso_valeurs, data_incid_valeurs, data_rr_DvsND_valeurs)

{

  ### Année de projection

  an_proj <- t

  ### Définition des risques de démence chez les non-consommateurs

  a010 <- matrix(c(0), # matrice vide
                 nrow=40*nb_iter, # âges
                 ncol=2, # âges, a010
                 byrow=T);
  colnames(a010) <- c("age_risque","a010")

  a010[,1] <- c(66:105)
  a010[,2] <- as.numeric(data_a01[which(data_a01[,1] != 65 & data_a01[,3]%in%(sexe)),2]) / (data_theta01[which(data_theta01[,1] != 65 & data_theta01[,3]%in%(sexe)),2]*data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2] - data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2] + 1);

  ### Définition des risques de démence chez les consommateurs

  a011 <- matrix(c(0), # matrice vide
                 nrow=40*nb_iter, # âges
                 ncol=2, # âges, a011
                 byrow=T);
  colnames(a011) <- c("age_risque","a011")

  a011[,1] <- c(66:105)
  a011[,2] <- a010[,2]*data_theta01[which(data_theta01[,1] != 65 & data_theta01[,3]%in%(sexe)),2];

  ### Définition des risques de démence global

  a01_global <- matrix(c(0), # matrice vide
                       nrow=40*nb_iter, # âges
                       ncol=2, # âges, a01_global
                       byrow=T);
  colnames(a01_global) <- c("age_risque","a01_global")

  a01_global[,1] <- c(66:105)
  a01_global[,2] <- data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2]*a010[,2]*data_theta01[which(data_theta01[,1] != 65 & data_theta01[,3]%in%(sexe)),2] + (1-data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2])*a010[,2];

  ### Définition des risques de décès chez les non-consommateurs

  a020 <- matrix(c(0), # matrice vide
                 nrow=40*nb_iter, # âges
                 ncol=132, # âges, années
                 byrow=T);
  colnames(a020) <- c("age_risque",1950:2080)

  a020[,1] <- c(66:105)

  for (a in 2:ncol(a020)){ # pour chaque année

    a020[,a] <- as.numeric(data_a02[which(data_a02[,1] != 65 & data_a02[,2]%in%(sexe)),a+1]) / (data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(sexe)),2]*data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2] - data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2] + 1);

  }

  ### Définition des risques de décès chez les consommateurs

  a021 <- matrix(c(0), # matrice vide
                 nrow=40*nb_iter, # âges
                 ncol=132, # âges, années
                 byrow=T);
  colnames(a021) <- c("age_risque",1950:2080)

  a021[,1] <- c(66:105)
  a021[,-1] <- a020[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(sexe)),2];

  ### Définition des risques de décès global

  a02_global <- matrix(c(0), # matrice vide
                       nrow=40*nb_iter, # âges
                       ncol=132, # âges, années
                       byrow=T);
  colnames(a02_global) <- c("age_risque",1950:2080)

  a02_global[,1] <- c(66:105)
  a02_global[,-1] <- data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2]*a020[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(sexe)),2] + (1-data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2])*a020[,-1];

  ### Définition des risques relatifs de décès pour les déments VS non-déments

  RR <- matrix(c(0), # matrice vide
               nrow=40*nb_iter, # âges
               ncol=2, # âges, RR
               byrow=T);
  colnames(RR) <- c("age_risque","rr_DvsND")

  RR[,1] <- c(66:105)
  RR[,2] <- data_rr_DvsND[which(data_rr_DvsND[,1] != 65 & data_rr_DvsND[,3]%in%(sexe)),2];

  ### Définition des risques de décès pour un dément chez les non-consommateurs

  a120 <- matrix(c(0), # matrice vide
                 nrow=40*nb_iter, # âges
                 ncol=132, # âges, années
                 byrow=T);
  colnames(a120) <- c("age_risque",1950:2080)

  a120[,1] <- c(66:105)

  for (a in 2:ncol(a020)){ # pour chaque année

    a120[,a] <- as.numeric(RR[,2])*a020[,a] / (data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(sexe)),2]*data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2] - data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2] + 1);

  }

  ### Définition des risques de décès pour un dément chez les consommateurs

  a121 <- matrix(c(0), # matrice vide
                 nrow=40*nb_iter, # âges
                 ncol=132, # âges, années
                 byrow=T);
  colnames(a121) <- c("age_risque",1950:2080)

  a121[,1] <- c(66:105)
  a121[,-1] <- a120[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(sexe)),2];

  ### Définition des risques de décès pour un dément global

  a12_global <- matrix(c(0), # matrice vide
                       nrow=40*nb_iter, # âges
                       ncol=132, # âges, années
                       byrow=T);
  colnames(a12_global) <- c("age_risque",1950:2080)

  a12_global[,1] <- c(66:105)
  a12_global[,-1] <- data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2]*a120[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(sexe)),2] + (1-data_conso[which(data_conso[,1] != 65 & data_conso[,3]%in%(sexe)),2])*a120[,-1];

  # Définition des vraies valeurs

  # Définition des risques de démence chez les non-consommateurs

  a010_valeurs <- matrix(c(0), # matrice vide
                         nrow=40, # âges
                         ncol=2, # âges, a010
                         byrow=T);
  colnames(a010_valeurs) <- c("age_risque","a010")

  a010_valeurs[,1] <- c(66:105)
  a010_valeurs[,2] <- as.numeric(data_a01_valeurs[which(data_a01_valeurs[,1] != 65 & data_a01_valeurs[,3]%in%(sexe)),2]) / (data_theta01_valeurs[which(data_theta01_valeurs[,1] != 65 & data_theta01_valeurs[,3]%in%(sexe)),2]*data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2] - data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2] + 1);

  ### Définition des risques de démence chez les consommateurs

  a011_valeurs <- matrix(c(0), # matrice vide
                         nrow=40, # âges
                         ncol=2, # âges, a011
                         byrow=T);
  colnames(a011_valeurs) <- c("age_risque","a011")

  a011_valeurs[,1] <- c(66:105)
  a011_valeurs[,2] <- a010_valeurs[,2]*data_theta01_valeurs[which(data_theta01_valeurs[,1] != 65 & data_theta01_valeurs[,3]%in%(sexe)),2];

  ### Définition des risques de démence global

  a01_global_valeurs <- matrix(c(0), # matrice vide
                               nrow=40, # âges
                               ncol=2, # âges, a01_global
                               byrow=T);
  colnames(a01_global_valeurs) <- c("age_risque","a01_global")

  a01_global_valeurs[,1] <- c(66:105)
  a01_global_valeurs[,2] <- data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2]*a010_valeurs[,2]*data_theta01_valeurs[which(data_theta01_valeurs[,1] != 65 & data_theta01_valeurs[,3]%in%(sexe)),2] + (1-data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2])*a010_valeurs[,2];

  ### Définition des risques de décès chez les non-consommateurs

  a020_valeurs <- matrix(c(0), # matrice vide
                         nrow=40, # âges
                         ncol=132, # âges, années
                         byrow=T);
  colnames(a020_valeurs) <- c("age_risque",1950:2080)

  a020_valeurs[,1] <- c(66:105)

  for (a in 2:ncol(a020_valeurs)){ # pour chaque année

    a020_valeurs[,a] <- as.numeric(data_a02_valeurs[which(data_a02_valeurs[,1] != 65 & data_a02_valeurs[,2]%in%(sexe)),a+1]) / (data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(sexe)),2]*data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2] - data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2] + 1);

  }

  ### Définition des risques de décès chez les consommateurs

  a021_valeurs <- matrix(c(0), # matrice vide
                         nrow=40, # âges
                         ncol=132, # âges, années
                         byrow=T);
  colnames(a021_valeurs) <- c("age_risque",1950:2080)

  a021_valeurs[,1] <- c(66:105)
  a021_valeurs[,-1] <- a020_valeurs[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(sexe)),2];

  ### Définition des risques de décès global

  a02_global_valeurs <- matrix(c(0), # matrice vide
                               nrow=40, # âges
                               ncol=132, # âges, années
                               byrow=T);
  colnames(a02_global_valeurs) <- c("age_risque",1950:2080)

  a02_global_valeurs[,1] <- c(66:105)
  a02_global_valeurs[,-1] <- data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2]*a020_valeurs[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(sexe)),2] + (1-data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2])*a020_valeurs[,-1];

  ### Définition des risques relatifs de décès pour les déments VS non-déments

  RR_valeurs <- matrix(c(0), # matrice vide
                       nrow=40, # âges
                       ncol=2, # âges, RR
                       byrow=T);
  colnames(RR_valeurs) <- c("age_risque","rr_DvsND")

  RR_valeurs[,1] <- c(66:105)
  RR_valeurs[,2] <- data_rr_DvsND_valeurs[which(data_rr_DvsND_valeurs[,1] != 65 & data_rr_DvsND_valeurs[,3]%in%(sexe)),2];

  ### Définition des risques de décès pour un dément chez les non-consommateurs

  a120_valeurs <- matrix(c(0), # matrice vide
                         nrow=40, # âges
                         ncol=132, # âges, années
                         byrow=T);
  colnames(a120_valeurs) <- c("age_risque",1950:2080)

  a120_valeurs[,1] <- c(66:105)

  for (a in 2:ncol(a020_valeurs)){ # pour chaque année

    a120_valeurs[,a] <- as.numeric(RR_valeurs[,2])*a020_valeurs[,a] / (data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(sexe)),2]*data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2] - data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2] + 1);

  }

  ### Définition des risques de décès pour un dément chez les consommateurs

  a121_valeurs <- matrix(c(0), # matrice vide
                         nrow=40, # âges
                         ncol=132, # âges, années
                         byrow=T);
  colnames(a121_valeurs) <- c("age_risque",1950:2080)

  a121_valeurs[,1] <- c(66:105)
  a121_valeurs[,-1] <- a120_valeurs[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(sexe)),2];

  ### Définition des risques de décès pour un dément global

  a12_global_valeurs <- matrix(c(0), # matrice vide
                               nrow=40, # âges
                               ncol=132, # âges, années
                               byrow=T);
  colnames(a12_global_valeurs) <- c("age_risque",1950:2080)

  a12_global_valeurs[,1] <- c(66:105)
  a12_global_valeurs[,-1] <- data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2]*a120_valeurs[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(sexe)),2] + (1-data_conso_valeurs[which(data_conso_valeurs[,1] != 65 & data_conso_valeurs[,3]%in%(sexe)),2])*a120_valeurs[,-1];

  ### Définition des matrices pour les indicateurs de santé

  ### Espérance de vie générale à l'âge a

  esp_vie_gen <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                        nrow=100-65+1,
                        ncol=2+nb_iter,
                        byrow = T);
  esp_vie_gen[,1] <- c(65:100);

  ### Espérance de vie générale chez les consommateurs de benzo à l'âge a

  esp_vie_gen_conso <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                              nrow=100-65+1,
                              ncol=2+nb_iter,
                              byrow = T);
  esp_vie_gen_conso[,1] <- c(65:100);

  ### Espérance de vie générale chez les non-consommateurs de benzo à l'âge a

  esp_vie_gen_nonconso <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                                 nrow=100-65+1,
                                 ncol=2+nb_iter,
                                 byrow = T);
  esp_vie_gen_nonconso[,1] <- c(65:100);

  ### Espérance de vie sans la maladie à l'âge a

  esp_vie_sans_mal <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                             nrow=100-65+1,
                             ncol=2+nb_iter,
                             byrow = T);
  esp_vie_sans_mal[,1] <- c(65:100);

  ### Espérance de vie sans la maladie chez les consommateurs de benzo à l'âge a

  esp_vie_sans_mal_conso <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                                   nrow=100-65+1,
                                   ncol=2+nb_iter,
                                   byrow = T);
  esp_vie_sans_mal_conso[,1] <- c(65:100);

  ### Espérance de vie sans la maladie chez les non-consommateurs de benzo à l'âge a

  esp_vie_sans_mal_nonconso <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                                      nrow=100-65+1,
                                      ncol=2+nb_iter,
                                      byrow = T);
  esp_vie_sans_mal_nonconso[,1] <- c(65:100);

  ### Espérance de vie d'un malade à l'âge a

  esp_vie_mal <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                        nrow=100-65+1,
                        ncol=2+nb_iter,
                        byrow = T);
  esp_vie_mal[,1] <- c(65:100);

  ### Espérance de vie d'un malade chez les consommateurs de benzo à l'âge a

  esp_vie_mal_conso <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                              nrow=100-65+1,
                              ncol=2+nb_iter,
                              byrow = T);
  esp_vie_mal_conso[,1] <- c(65:100);

  ### Espérance de vie d'un malade chez les non-consommateurs de benzo à l'âge a

  esp_vie_mal_nonconso <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                                 nrow=100-65+1,
                                 ncol=2+nb_iter,
                                 byrow = T);
  esp_vie_mal_nonconso[,1] <- c(65:100);

  ### Espérance de vie d'un non-malade à l'âge a

  esp_vie_non_mal <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                            nrow=100-65+1,
                            ncol=2+nb_iter,
                            byrow = T);
  esp_vie_non_mal[,1] <- c(65:100);

  ### Espérance de vie d'un non-malade chez les consommateurs de benzo à l'âge a

  esp_vie_non_mal_conso <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                                  nrow=100-65+1,
                                  ncol=2+nb_iter,
                                  byrow = T);
  esp_vie_non_mal_conso[,1] <- c(65:100);

  ### Espérance de vie d'un non-malade chez les non-consommateurs de benzo à l'âge a

  esp_vie_non_mal_nonconso <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                                     nrow=100-65+1,
                                     ncol=2+nb_iter,
                                     byrow = T);
  esp_vie_non_mal_nonconso[,1] <- c(65:100);

  ### Prévalence pour une année t

  prevalence <- matrix(c(0), # matrice de calcul de la prévalence
                       nrow=99-65+1,
                       ncol=2+nb_iter,
                       byrow=T);
  prevalence[,1] <- c(65:99);

  ### Taux de prévalence pour une année t

  taux_prevalence <- matrix(c(0), # matrice de calcul du taux de prévalence
                            nrow=99-65+1,
                            ncol=2+nb_iter,
                            byrow=T);
  taux_prevalence[,1] <- c(65:99);

  ### Survie pour une année t

  survie <- matrix(c(0), # matrice de calcul de la survie
                   nrow=99-65+1,
                   ncol=2+nb_iter,
                   byrow=T);
  survie[,1] <- c(65:99);

  ### Taux de survivants pour une année t

  taux_survivants <- matrix(c(0), # matrice de calcul du taux de survivants
                            nrow=99-65+1,
                            ncol=2+nb_iter,
                            byrow=T);
  taux_survivants[,1] <- c(65:99);

  ### Nombre moyen d'années passées avec la maladie à l'âge a

  nb_moy_dem <- matrix(c(0), # matrice de calcul du nombre moyen d'années en démence
                       nrow=100-65+1,
                       ncol=2+nb_iter,
                       byrow = T);
  nb_moy_dem[,1] <- c(65:100);

  ### Nombre moyen d'années passées avec la maladie chez les consommateurs de benzo à l'âge a

  nb_moy_dem_conso <- matrix(c(0), # matrice de calcul du nombre moyen d'années en démence
                             nrow=100-65+1,
                             ncol=2+nb_iter,
                             byrow = T);
  nb_moy_dem_conso[,1] <- c(65:100);

  ### Nombre moyen d'années passées avec la maladie chez les non-consommateurs de benzo à l'âge a

  nb_moy_dem_nonconso <- matrix(c(0), # matrice de calcul du nombre moyen d'années en démence
                                nrow=100-65+1,
                                ncol=2+nb_iter,
                                byrow = T);
  nb_moy_dem_nonconso[,1] <- c(65:100);

  ### Proba vie entière d'apparition de la démence

  prb_dem <- matrix(c(0), # matrice de calcul de la proba vie entière d'apparition de la démence
                    nrow=99-65+1,
                    ncol=2+nb_iter,
                    byrow=T);
  prb_dem[,1] <- c(65:99);

  ### Age moyen d'apparition de la démence pour une année t

  age_dem <- matrix(c(0), # matrice de calcul de l'âge moyen d'apparition de la démence
                    nrow=99-65+1,
                    ncol=2+nb_iter,
                    byrow=T);
  age_dem[,1] <- c(65:99);

  ### Nombre moyen d'années de consommation de benzo

  moy_conso <- matrix(c(0), # matrice de calcul du nombre moyen d'années de conso
                      nrow=99-65+1,
                      ncol=2+nb_iter,
                      byrow=T);
  moy_conso[,1] <- c(65:99);

  # Prévalence de la conso estimée

  prevalence_conso <- matrix(c(0), # matrice de calcul du nombre moyen d'années de conso
                             nrow=105-65+1,
                             ncol=2+nb_iter,
                             byrow=T);
  prevalence_conso[,1] <- c(65:105);

  # Quotient de mortalité estimée

  quotient_mortalite <- matrix(c(0), # matrice de calcul du quotient de mortalité
                               nrow=99-65+1,
                               ncol=2+nb_iter,
                               byrow=T);
  quotient_mortalite[,1] <- c(65:99);

  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################

  ### Boucle de calcul pour valeur avec les vrais paramètres

  ###############################
  ###          ETAPE 1        ###
  ###############################

  ### Recalcul de la prévalence de la conso estimée

  pr_conso_benzo <- matrix(c(0), # matrice de calcul du nombre moyen d'années de conso
                           nrow=105-65+1,
                           ncol=2,
                           byrow=T);
  pr_conso_benzo[,1] <- c(65:105);

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

      if (alea0 <= data_conso_valeurs[which(data_conso_valeurs[,1]%in%(65) & data_conso_valeurs[,3]%in%(sexe)),2]) {

        etat[i,1] <- "01"

      } else {

        etat[i,1] <- "00"

      }

    };

    for (i in 1:nrow(etat)) { # pour chaque individu

      for (j in 2:ncol(etat)) { # pour chaque âge

        alea <- runif(1, 0, 1); # tirage d'un nombre aléatoire pour savoir si il devient sain/dément/décès

        alea0 <- runif(1, 0, 1); # tirage d'un nombre aléatoire pour savoir si il est nouveau consommateur

        if (etat[i,j-1] == "00") { # état année précédente est non-dément et non-consommateur

          if (alea0 <= data_incid_valeurs[which(data_incid_valeurs[,1]%in%(j-1+65) & data_incid_valeurs[,3]%in%(sexe)),2]) {

            a01 <- a011_valeurs;
            a02 <- a021_valeurs;

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

            a01 <- a010_valeurs;
            a02 <- a020_valeurs;

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

            a01 <- a011_valeurs;
            a02 <- a021_valeurs;

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

              if (alea0 <= data_incid_valeurs[which(data_incid_valeurs[,1]%in%(j-1+65) & data_incid_valeurs[,3]%in%(sexe)),2]) {

                a12 <- a121_valeurs;

                if (alea <= a12[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "21"; # individu est devenu consommateur au cours de l'année et décède

                } else {

                  etat[i,j] <- "11"; # individu est devenu consommateur au cours de l'année et reste malade

                }

              } else {

                a12 <- a120_valeurs;

                if (alea <= a12[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "20"; # individu est resté non-consommateur au cours de l'année et décède

                } else {

                  etat[i,j] <- "10"; # individu est resté non-consommateur au cours de l'année et reste malade

                }

              }

            } else {

              if (etat[i,j-1] == "11") { # état année précédente est dément et consommateur

                a12 <- a121_valeurs;

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

  a010_valeurs[,2] <- as.numeric(data_a01_valeurs[which(data_a01_valeurs[,1] != 65 & data_a01_valeurs[,3]%in%(sexe)),2]) / (data_theta01_valeurs[which(data_theta01_valeurs[,1] != 65 & data_theta01_valeurs[,3]%in%(sexe)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

  ### Définition des risques de démence chez les consommateurs

  a011_valeurs[,2] <- a010_valeurs[,2]*data_theta01_valeurs[which(data_theta01_valeurs[,1] != 65 & data_theta01_valeurs[,3]%in%(sexe)),2];

  ### Définition des risques de démence global

  a01_global_valeurs[,2] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a010_valeurs[,2]*data_theta01_valeurs[which(data_theta01_valeurs[,1] != 65 & data_theta01_valeurs[,3]%in%(sexe)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a010_valeurs[,2];

  ### Définition des risques de décès chez les non-consommateurs

  for (a in 2:ncol(a020_valeurs)){ # pour chaque année

    a020_valeurs[,a] <- as.numeric(data_a02_valeurs[which(data_a02_valeurs[,1] != 65 & data_a02_valeurs[,2]%in%(sexe)),a+1]) / (data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(sexe)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

  }

  ### Définition des risques de décès chez les consommateurs

  a021_valeurs[,-1] <- a020_valeurs[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(sexe)),2];

  ### Définition des risques de décès global

  a02_global_valeurs[,-1] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a020_valeurs[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(sexe)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a020_valeurs[,-1];

  ### Définition des risques de décès pour un dément chez les non-consommateurs

  for (a in 2:ncol(a020_valeurs)){ # pour chaque année

    a120_valeurs[,a] <- as.numeric(RR_valeurs[,2])*a020_valeurs[,a] / (data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(sexe)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

  }

  ### Définition des risques de décès pour un dément chez les consommateurs

  a121_valeurs[,-1] <- a120_valeurs[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(sexe)),2];

  ### Définition des risques de décès pour un dément global

  a12_global_valeurs[,-1] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a120_valeurs[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(sexe)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a120_valeurs[,-1];

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

          donnees_incid <- data_incid_valeurs

        } else {

          if (scenario == 1) {

            if (annee < an_scenario) {

              donnees_incid <- data_incid_valeurs

            } else {

              donnees_incid <- data_incid_valeurs
              donnees_incid[,2] <- data_incid_valeurs[,2] / 2

            }

          } else {

            if (scenario == 2) {

              if (annee < an_scenario) {

                donnees_incid <- data_incid_valeurs

              } else {

                donnees_incid <- data_incid_valeurs
                donnees_incid[,2] <- 0

              }

            }

          }

        }

        if (etat[i,j-1] == "00") { # état année précédente est non-dément et non-consommateur

          if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(sexe)),2]) {

            a01 <- a011_valeurs;
            a02 <- a021_valeurs;

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

            a01 <- a010_valeurs;
            a02 <- a020_valeurs;

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

            a01 <- a011_valeurs;
            a02 <- a021_valeurs;

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

                a12 <- a121_valeurs;

                if (alea <= a12[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "21"; # individu est devenu consommateur au cours de l'année et décède

                } else {

                  etat[i,j] <- "11"; # individu est devenu consommateur au cours de l'année et reste malade

                }

              } else {

                a12 <- a120_valeurs;

                if (alea <= a12[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "20"; # individu est resté non-consommateur au cours de l'année et décède

                } else {

                  etat[i,j] <- "10"; # individu est resté non-consommateur au cours de l'année et reste malade

                }

              }

            } else {

              if (etat[i,j-1] == "11") { # état année précédente est dément et consommateur

                a12 <- a121_valeurs;

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

    }

    ### Calcul d'indicateurs :

    ### Espérance de vie générale

    if (age < 101) {

      n0 <- vector(length = ncol(etat));

      s0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11"));

      if (s0 != 0) {

        for (j in (age-63):ncol(etat)) {

          n0[j] <- sum(etat[,j]%in%("00") | etat[,j]%in%("01") | etat[,j]%in%("10") | etat[,j]%in%("11"))

        };

        esp_vie_gen[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_gen_conso[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_gen_nonconso[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_sans_mal[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_sans_mal_conso[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_sans_mal_nonconso[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_mal[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_mal_conso[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_mal_nonconso[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_non_mal[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_non_mal_conso[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

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

        esp_vie_non_mal_nonconso[age-64,2] <- 1 + sum(n0) / s0;

      } else {

        n0 <- NA;

      }

    };

    ### Prévalence de la démence

    if (age > 65 & age < 100) {

      d0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("00") | etat[,age-65]%in%("01"))) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11")));

      s0 <- sum(etat[,1]%in%("00") | etat[,1]%in%("01") | etat[,1]%in%("10") | etat[,1]%in%("11"));

      s1 <- sum((etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11"))) + 0.5 * sum((etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) & (etat[,age-65]%in%("00") | etat[,age-65]%in%("01"))) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11")));

      if (s0 != 0) {

        p01 <- s1/s0;

        taux_prevalence[age-64,2] <- s1/d0;

        nb <- p01*data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(sexe)),2];

        prevalence[age-64,2] <- nb;

      };

    };

    ### Survie

    if (age == 65) {

      survie[age-64,2] <- data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(sexe)),2]

    };

    if (age > 65 & age < 100) {

      d0 <- sum(etat[,age-65]%in%("00") | etat[,age-65]%in%("01") | etat[,age-65]%in%("10") | etat[,age-65]%in%("11"));

      s0 <- sum(etat[,1]%in%("00") | etat[,1]%in%("01") | etat[,1]%in%("10") | etat[,1]%in%("11"));

      s1 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11"));

      if (s0 != 0) {

        p01 <- s1/s0;

        taux_survivants[age-64,2] <- s1/d0;

        nb <- p01*data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(sexe)),2];

        survie[age-64,2] <- nb;

      };

    };

    ### Nombre moyen d'années passées avec la maladie

    if (age < 101) {

      nb_moy_dem[age-64,2] <- esp_vie_non_mal[age-64,2] - esp_vie_sans_mal[age-64,2]

    } else {

      n0 <- NA

    };

    ### Nombre moyen d'années passées avec la maladie chez les consommateurs de benzo

    if (age < 101) {

      nb_moy_dem_conso[age-64,2] <- esp_vie_non_mal_conso[age-64,2] - esp_vie_sans_mal_conso[age-64,2];

    } else {

      n0 <- NA;

    };

    ### Nombre moyen d'années passées avec la maladie chez les non-consommateurs de benzo

    if (age < 101) {

      nb_moy_dem_nonconso[age-64,2] <- esp_vie_non_mal_nonconso[age-64,2] - esp_vie_sans_mal_nonconso[age-64,2];

    } else {

      n0 <- NA;

    };

    ### Proba d'apparition de la démence

    if (age == 65) {

      for (i in (age-63):nrow(prb_dem)) {

        prb_dem[i,2] <- sum(etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("10") | etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("11"))

      }

    };

    ### Age moyen d'apparition de la démence

    if (age == 65) {

      for (i in (age-63):nrow(age_dem)) {

        age_dem[i,2] <- sum(etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("10") | etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("11"))

      }

    };

    ### Nombre moyen d'années de consommation de benzo

    if (age == 65) {

      for (i in (age-63):nrow(age_dem)) {

        moy_conso[i,2] <- sum(etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),i]%in%("01") | etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),i]%in%("11"))

      }

    };

    ### Prévalence de la conso estimée

    if (age <= 105) {

      s0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("10") | etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

      s1 <- sum(etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

      if (s0 != 0) {

        prevalence_conso[age-64,2] <- s1/s0;

      } else {

        prevalence_conso[age-64,2] <- NA;

      };

    };

    ### Quotient de mortalité estimée

    if (age > 65 & age < 100) {

      s0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("10") | etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

      s1 <- sum(etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) - sum(etat[,age-65]%in%("20") | etat[,age-65]%in%("21"));

      if (s0 != 0) {

        quotient_mortalite[age-64,2] <- s1/s0;

      } else {

        quotient_mortalite[age-64,2] <- NA;

      };

    };

  }

  ### Boucle de calcul pour variabilité

  if (nb_iter != 0) {

    indicateurs <- varHI(t = t,
                         scenario = scenario,
                         an_scenario = an_scenario,
                         nbind = nbind,
                         nb_iter = nb_iter,
                         data_pop = data_pop,
                         sexe = sexe,
                         an_proj = an_proj,
                         data_conso = data_conso,
                         data_incid = data_incid,
                         a010 = a010,
                         a011 = a011,
                         a01_global = a01_global,
                         a020 = a020,
                         a021 = a021,
                         a02_global = a02_global,
                         a120 = a120,
                         a121 = a121,
                         a12_global = a12_global,
                         data_a01 = data_a01,
                         data_theta01 = data_theta01,
                         data_a02 = data_a02,
                         data_theta02 = data_theta02,
                         data_theta12 = data_theta12,
                         RR = RR,
                         prb_dem = prb_dem,
                         age_dem = age_dem)

    for (i in 1:nb_iter) {
      esp_vie_gen[,i+2] <- indicateurs[,i]$ev_gen
      esp_vie_gen_conso[,i+2] <- indicateurs[,i]$ev_gen_conso
      esp_vie_gen_nonconso[,i+2] <- indicateurs[,i]$ev_gen_nonconso
      esp_vie_sans_mal[,i+2] <- indicateurs[,i]$ev_sans_mal
      esp_vie_sans_mal_conso[,i+2] <- indicateurs[,i]$ev_sans_mal_conso
      esp_vie_sans_mal_nonconso[,i+2] <- indicateurs[,i]$ev_sans_mal_nonconso
      esp_vie_mal[,i+2] <- indicateurs[,i]$ev_mal
      esp_vie_mal_conso[,i+2] <- indicateurs[,i]$ev_mal_conso
      esp_vie_mal_nonconso[,i+2] <- indicateurs[,i]$ev_mal_nonconso
      esp_vie_non_mal[,i+2] <- indicateurs[,i]$ev_non_mal
      esp_vie_non_mal_conso[,i+2] <- indicateurs[,i]$ev_non_mal_conso
      esp_vie_non_mal_nonconso[,i+2] <- indicateurs[,i]$ev_non_mal_nonconso
      prevalence[,i+2] <- indicateurs[,i]$np_age
      taux_prevalence[,i+2] <- indicateurs[,i]$tp_dem
      survie[,i+2] <- indicateurs[,i]$nsurvie
      taux_survivants[,i+2] <- indicateurs[,i]$tsurvie
      nb_moy_dem[,i+2] <- esp_vie_non_mal[,i+2] - esp_vie_sans_mal[,i+2]
      nb_moy_dem_conso [,i+2] <- esp_vie_non_mal_conso[,i+2] - esp_vie_sans_mal_conso[,i+2]
      nb_moy_dem_nonconso [,i+2] <- esp_vie_non_mal_nonconso[,i+2] - esp_vie_sans_mal_nonconso[,i+2]
      prb_dem[,i+2] <- indicateurs[,i]$p_dem
      age_dem[,i+2] <- indicateurs[,i]$a_dem
      moy_conso[,i+2] <- indicateurs[,i]$m_conso
      prevalence_conso[,i+2] <- indicateurs[,i]$p_conso
      quotient_mortalite[,i+2] <- indicateurs[,i]$q_mortalite
    }

  }

    ### Calcul des indicateurs

  ### Espérance de vie générale à l'âge a

  esperance_vie_gen <- matrix(c(0), # matrice des espérances de vie à chaque âge
                              nrow=100-65+1,
                              ncol=5,
                              byrow = T);
  esperance_vie_gen[,1] <- esp_vie_gen[,1]
  esperance_vie_gen[,2] <- esp_vie_gen[,2]
  esperance_vie_gen[,3] <- apply(esp_vie_gen[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_gen[,4] <- esperance_vie_gen[,2] - 1.96*esperance_vie_gen[,3];
  esperance_vie_gen[,5] <- esperance_vie_gen[,2] + 1.96*esperance_vie_gen[,3];

  for (i in 1:nrow(esperance_vie_gen)) {
    if (esperance_vie_gen[i,4]<0  & is.na(esperance_vie_gen[i,4])==F) {
      esperance_vie_gen[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_gen)) {
    if (esperance_vie_gen[i,2]==0 & esperance_vie_gen[i,4]==0 & esperance_vie_gen[i,5]==0 &
        is.na(esperance_vie_gen[i,2])==F & is.na(esperance_vie_gen[i,4])==F & is.na(esperance_vie_gen[i,5])==F) {
      esperance_vie_gen[i,2] <- NA;
      esperance_vie_gen[i,3] <- NA;
      esperance_vie_gen[i,4] <- NA;
      esperance_vie_gen[i,5] <- NA
    }
  }

  colnames(esperance_vie_gen) <- c("Age a","Espérance de vie générale à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie générale chez les consommateurs de benzo à l'âge a

  esperance_vie_gen_conso <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                    nrow=100-65+1,
                                    ncol=5,
                                    byrow = T);
  esperance_vie_gen_conso[,1] <- esp_vie_gen_conso[,1]
  esperance_vie_gen_conso[,2] <- esp_vie_gen_conso[,2]
  esperance_vie_gen_conso[,3] <- apply(esp_vie_gen_conso[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_gen_conso[,4] <- esperance_vie_gen_conso[,2] - 1.96*esperance_vie_gen_conso[,3];
  esperance_vie_gen_conso[,5] <- esperance_vie_gen_conso[,2] + 1.96*esperance_vie_gen_conso[,3];

  for (i in 1:nrow(esperance_vie_gen_conso)) {
    if (esperance_vie_gen_conso[i,4]<0  & is.na(esperance_vie_gen_conso[i,4])==F) {
      esperance_vie_gen_conso[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_gen_conso)) {
    if (esperance_vie_gen_conso[i,2]==0 & esperance_vie_gen_conso[i,4]==0 & esperance_vie_gen_conso[i,5]==0 &
        is.na(esperance_vie_gen_conso[i,2])==F & is.na(esperance_vie_gen_conso[i,4])==F & is.na(esperance_vie_gen_conso[i,5])==F) {
      esperance_vie_gen_conso[i,2] <- NA;
      esperance_vie_gen_conso[i,3] <- NA;
      esperance_vie_gen_conso[i,4] <- NA;
      esperance_vie_gen_conso[i,5] <- NA
    }
  }

  colnames(esperance_vie_gen_conso) <- c("Age a","Espérance de vie générale à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie générale chez les non-consommateurs de benzo à l'âge a

  esperance_vie_gen_nonconso <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                       nrow=100-65+1,
                                       ncol=5,
                                       byrow = T);
  esperance_vie_gen_nonconso[,1] <- esp_vie_gen_nonconso[,1]
  esperance_vie_gen_nonconso[,2] <- esp_vie_gen_nonconso[,2]
  esperance_vie_gen_nonconso[,3] <- apply(esp_vie_gen_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_gen_nonconso[,4] <- esperance_vie_gen_nonconso[,2] - 1.96*esperance_vie_gen_nonconso[,3];
  esperance_vie_gen_nonconso[,5] <- esperance_vie_gen_nonconso[,2] + 1.96*esperance_vie_gen_nonconso[,3];

  for (i in 1:nrow(esperance_vie_gen_nonconso)) {
    if (esperance_vie_gen_nonconso[i,4]<0 & is.na(esperance_vie_gen_nonconso[i,4])==F) {
      esperance_vie_gen_nonconso[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_gen_nonconso)) {
    if (esperance_vie_gen_nonconso[i,2]==0 & esperance_vie_gen_nonconso[i,4]==0 & esperance_vie_gen_nonconso[i,5]==0 &
        is.na(esperance_vie_gen_nonconso[i,2])==F & is.na(esperance_vie_gen_nonconso[i,4])==F & is.na(esperance_vie_gen_nonconso[i,5])==F) {
      esperance_vie_gen_nonconso[i,2] <- NA;
      esperance_vie_gen_nonconso[i,3] <- NA;
      esperance_vie_gen_nonconso[i,4] <- NA;
      esperance_vie_gen_nonconso[i,5] <- NA
    }
  }

  colnames(esperance_vie_gen_nonconso) <- c("Age a","Espérance de vie générale à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie sans la maladie à l'âge a

  esperance_vie_sans_mal <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                   nrow=100-65+1,
                                   ncol=5,
                                   byrow = T);
  esperance_vie_sans_mal[,1] <- esp_vie_sans_mal[,1]
  esperance_vie_sans_mal[,2] <- esp_vie_sans_mal[,2]
  esperance_vie_sans_mal[,3] <- apply(esp_vie_sans_mal[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_sans_mal[,4] <- esperance_vie_sans_mal[,2] - 1.96*esperance_vie_sans_mal[,3];
  esperance_vie_sans_mal[,5] <- esperance_vie_sans_mal[,2] + 1.96*esperance_vie_sans_mal[,3];

  for (i in 1:nrow(esperance_vie_sans_mal)) {
    if (esperance_vie_sans_mal[i,4]<0  & is.na(esperance_vie_sans_mal[i,4])==F) {
      esperance_vie_sans_mal[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_sans_mal)) {
    if (esperance_vie_sans_mal[i,2]==0 & esperance_vie_sans_mal[i,4]==0 & esperance_vie_sans_mal[i,5]==0 &
        is.na(esperance_vie_sans_mal[i,2])==F & is.na(esperance_vie_sans_mal[i,4])==F & is.na(esperance_vie_sans_mal[i,5])==F) {
      esperance_vie_sans_mal[i,2] <- NA;
      esperance_vie_sans_mal[i,3] <- NA;
      esperance_vie_sans_mal[i,4] <- NA;
      esperance_vie_sans_mal[i,5] <- NA
    }
  }

  colnames(esperance_vie_sans_mal) <- c("Age a","Espérance de vie sans la maladie à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie sans la maladie chez les consommateurs à l'âge a

  esperance_vie_sans_mal_conso <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                         nrow=100-65+1,
                                         ncol=5,
                                         byrow = T);
  esperance_vie_sans_mal_conso[,1] <- esp_vie_sans_mal_conso[,1]
  esperance_vie_sans_mal_conso[,2] <- esp_vie_sans_mal_conso[,2]
  esperance_vie_sans_mal_conso[,3] <- apply(esp_vie_sans_mal_conso[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_sans_mal_conso[,4] <- esperance_vie_sans_mal_conso[,2] - 1.96*esperance_vie_sans_mal_conso[,3];
  esperance_vie_sans_mal_conso[,5] <- esperance_vie_sans_mal_conso[,2] + 1.96*esperance_vie_sans_mal_conso[,3];

  for (i in 1:nrow(esperance_vie_sans_mal_conso)) {
    if (esperance_vie_sans_mal_conso[i,4]<0  & is.na(esperance_vie_sans_mal_conso[i,4])==F) {
      esperance_vie_sans_mal_conso[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_sans_mal_conso)) {
    if (esperance_vie_sans_mal_conso[i,2]==0 & esperance_vie_sans_mal_conso[i,4]==0 & esperance_vie_sans_mal_conso[i,5]==0 &
        is.na(esperance_vie_sans_mal_conso[i,2])==F & is.na(esperance_vie_sans_mal_conso[i,4])==F & is.na(esperance_vie_sans_mal_conso[i,5])==F) {
      esperance_vie_sans_mal_conso[i,2] <- NA;
      esperance_vie_sans_mal_conso[i,3] <- NA;
      esperance_vie_sans_mal_conso[i,4] <- NA;
      esperance_vie_sans_mal_conso[i,5] <- NA
    }
  }

  colnames(esperance_vie_sans_mal_conso) <- c("Age a","Espérance de vie sans la maladie à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie sans la maladie chez les non-consommateurs à l'âge a

  esperance_vie_sans_mal_nonconso <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                            nrow=100-65+1,
                                            ncol=5,
                                            byrow = T);
  esperance_vie_sans_mal_nonconso[,1] <- esp_vie_sans_mal_nonconso[,1]
  esperance_vie_sans_mal_nonconso[,2] <- esp_vie_sans_mal_nonconso[,2]
  esperance_vie_sans_mal_nonconso[,3] <- apply(esp_vie_sans_mal_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_sans_mal_nonconso[,4] <- esperance_vie_sans_mal_nonconso[,2] - 1.96*esperance_vie_sans_mal_nonconso[,3];
  esperance_vie_sans_mal_nonconso[,5] <- esperance_vie_sans_mal_nonconso[,2] + 1.96*esperance_vie_sans_mal_nonconso[,3];

  for (i in 1:nrow(esperance_vie_sans_mal_nonconso)) {
    if (esperance_vie_sans_mal_nonconso[i,4]<0  & is.na(esperance_vie_sans_mal_nonconso[i,4])==F) {
      esperance_vie_sans_mal_nonconso[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_sans_mal_nonconso)) {
    if (esperance_vie_sans_mal_nonconso[i,2]==0 & esperance_vie_sans_mal_nonconso[i,4]==0 & esperance_vie_sans_mal_nonconso[i,5]==0 &
        is.na(esperance_vie_sans_mal_nonconso[i,2])==F & is.na(esperance_vie_sans_mal_nonconso[i,4])==F & is.na(esperance_vie_sans_mal_nonconso[i,5])==F) {
      esperance_vie_sans_mal_nonconso[i,2] <- NA;
      esperance_vie_sans_mal_nonconso[i,3] <- NA;
      esperance_vie_sans_mal_nonconso[i,4] <- NA;
      esperance_vie_sans_mal_nonconso[i,5] <- NA
    }
  }

  colnames(esperance_vie_sans_mal_nonconso) <- c("Age a","Espérance de vie sans la maladie à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie d'un malade à l'âge a

  esperance_vie_mal <- matrix(c(0), # matrice des espérances de vie à chaque âge
                              nrow=100-65+1,
                              ncol=5,
                              byrow = T);
  esperance_vie_mal[,1] <- esp_vie_mal[,1]
  esperance_vie_mal[,2] <- esp_vie_mal[,2]
  esperance_vie_mal[,3] <- apply(esp_vie_mal[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_mal[,4] <- esperance_vie_mal[,2] - 1.96*esperance_vie_mal[,3];
  esperance_vie_mal[,5] <- esperance_vie_mal[,2] + 1.96*esperance_vie_mal[,3];

  for (i in 1:nrow(esperance_vie_mal)) {
    if (esperance_vie_mal[i,4]<0  & is.na(esperance_vie_mal[i,4])==F) {
      esperance_vie_mal[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_mal)) {
    if (esperance_vie_mal[i,2]==0 & esperance_vie_mal[i,4]==0 & esperance_vie_mal[i,5]==0 &
        is.na(esperance_vie_mal[i,2])==F & is.na(esperance_vie_mal[i,4])==F & is.na(esperance_vie_mal[i,5])==F) {
      esperance_vie_mal[i,2] <- NA;
      esperance_vie_mal[i,3] <- NA;
      esperance_vie_mal[i,4] <- NA;
      esperance_vie_mal[i,5] <- NA
    }
  }

  colnames(esperance_vie_mal) <- c("Age a","Espérance de vie d'un malade à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie d'un malade chez les consommateurs à l'âge a

  esperance_vie_mal_conso <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                    nrow=100-65+1,
                                    ncol=5,
                                    byrow = T);
  esperance_vie_mal_conso[,1] <- esp_vie_mal_conso[,1]
  esperance_vie_mal_conso[,2] <- esp_vie_mal_conso[,2]
  esperance_vie_mal_conso[,3] <- apply(esp_vie_mal_conso[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_mal_conso[,4] <- esperance_vie_mal_conso[,2] - 1.96*esperance_vie_mal_conso[,3];
  esperance_vie_mal_conso[,5] <- esperance_vie_mal_conso[,2] + 1.96*esperance_vie_mal_conso[,3];

  for (i in 1:nrow(esperance_vie_mal_conso)) {
    if (esperance_vie_mal_conso[i,4]<0 & is.na(esperance_vie_mal_conso[i,4])==F) {
      esperance_vie_mal_conso[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_mal_conso)) {
    if (esperance_vie_mal_conso[i,2]==0 & esperance_vie_mal_conso[i,4]==0 & esperance_vie_mal_conso[i,5]==0 &
        is.na(esperance_vie_mal_conso[i,2])==F & is.na(esperance_vie_mal_conso[i,4])==F & is.na(esperance_vie_mal_conso[i,5])==F) {
      esperance_vie_mal_conso[i,2] <- NA;
      esperance_vie_mal_conso[i,3] <- NA;
      esperance_vie_mal_conso[i,4] <- NA;
      esperance_vie_mal_conso[i,5] <- NA
    }
  }

  colnames(esperance_vie_mal_conso) <- c("Age a","Espérance de vie d'un malade à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie d'un malade chez les non-consommateurs à l'âge a

  esperance_vie_mal_nonconso <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                       nrow=100-65+1,
                                       ncol=5,
                                       byrow = T);
  esperance_vie_mal_nonconso[,1] <- esp_vie_mal_nonconso[,1]
  esperance_vie_mal_nonconso[,2] <- esp_vie_mal_nonconso[,2]
  esperance_vie_mal_nonconso[,3] <- apply(esp_vie_mal_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_mal_nonconso[,4] <- esperance_vie_mal_nonconso[,2] - 1.96*esperance_vie_mal_nonconso[,3];
  esperance_vie_mal_nonconso[,5] <- esperance_vie_mal_nonconso[,2] + 1.96*esperance_vie_mal_nonconso[,3];

  for (i in 1:nrow(esperance_vie_mal_nonconso)) {
    if (esperance_vie_mal_nonconso[i,4]<0  & is.na(esperance_vie_mal_nonconso[i,4])==F) {
      esperance_vie_mal_nonconso[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_mal_nonconso)) {
    if (esperance_vie_mal_nonconso[i,2]==0 & esperance_vie_mal_nonconso[i,4]==0 & esperance_vie_mal_nonconso[i,5]==0 &
        is.na(esperance_vie_mal_nonconso[i,2])==F & is.na(esperance_vie_mal_nonconso[i,4])==F & is.na(esperance_vie_mal_nonconso[i,5])==F) {
      esperance_vie_mal_nonconso[i,2] <- NA;
      esperance_vie_mal_nonconso[i,3] <- NA;
      esperance_vie_mal_nonconso[i,4] <- NA;
      esperance_vie_mal_nonconso[i,5] <- NA
    }
  }

  colnames(esperance_vie_mal_nonconso) <- c("Age a","Espérance de vie d'un malade à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie d'un non-malade à l'âge a

  esperance_vie_non_mal <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                  nrow=100-65+1,
                                  ncol=5,
                                  byrow = T);
  esperance_vie_non_mal[,1] <- esp_vie_non_mal[,1]
  esperance_vie_non_mal[,2] <- esp_vie_non_mal[,2]
  esperance_vie_non_mal[,3] <- apply(esp_vie_non_mal[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_non_mal[,4] <- esperance_vie_non_mal[,2] - 1.96*esperance_vie_non_mal[,3];
  esperance_vie_non_mal[,5] <- esperance_vie_non_mal[,2] + 1.96*esperance_vie_non_mal[,3];

  for (i in 1:nrow(esperance_vie_non_mal)) {
    if (esperance_vie_non_mal[i,4]<0  & is.na(esperance_vie_non_mal[i,4])==F) {
      esperance_vie_non_mal[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_non_mal)) {
    if (esperance_vie_non_mal[i,2]==0 & esperance_vie_non_mal[i,4]==0 & esperance_vie_non_mal[i,5]==0 &
        is.na(esperance_vie_non_mal[i,2])==F & is.na(esperance_vie_non_mal[i,4])==F & is.na(esperance_vie_non_mal[i,5])==F) {
      esperance_vie_non_mal[i,2] <- NA;
      esperance_vie_non_mal[i,3] <- NA;
      esperance_vie_non_mal[i,4] <- NA;
      esperance_vie_non_mal[i,5] <- NA
    }
  }

  colnames(esperance_vie_non_mal) <- c("Age a","Espérance de vie d'un non-malade à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie d'un non-malade chez les consommateurs de benzo à l'âge a

  esperance_vie_non_mal_conso <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                        nrow=100-65+1,
                                        ncol=5,
                                        byrow = T);
  esperance_vie_non_mal_conso[,1] <- esp_vie_non_mal_conso[,1]
  esperance_vie_non_mal_conso[,2] <- esp_vie_non_mal_conso[,2]
  esperance_vie_non_mal_conso[,3] <- apply(esp_vie_non_mal_conso[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_non_mal_conso[,4] <- esperance_vie_non_mal_conso[,2] - 1.96*esperance_vie_non_mal_conso[,3];
  esperance_vie_non_mal_conso[,5] <- esperance_vie_non_mal_conso[,2] + 1.96*esperance_vie_non_mal_conso[,3];

  for (i in 1:nrow(esperance_vie_non_mal_conso)) {
    if (esperance_vie_non_mal_conso[i,4]<0  & is.na(esperance_vie_non_mal_conso[i,4])==F) {
      esperance_vie_non_mal_conso[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_non_mal_conso)) {
    if (esperance_vie_non_mal_conso[i,2]==0 & esperance_vie_non_mal_conso[i,4]==0 & esperance_vie_non_mal_conso[i,5]==0 &
        is.na(esperance_vie_non_mal_conso[i,2])==F & is.na(esperance_vie_non_mal_conso[i,4])==F & is.na(esperance_vie_non_mal_conso[i,5])==F) {
      esperance_vie_non_mal_conso[i,2] <- NA;
      esperance_vie_non_mal_conso[i,3] <- NA;
      esperance_vie_non_mal_conso[i,4] <- NA;
      esperance_vie_non_mal_conso[i,5] <- NA
    }
  }

  colnames(esperance_vie_non_mal_conso) <- c("Age a","Espérance de vie d'un non-malade à l'âge a","std","Borne inf","Borne sup");

  ### Espérance de vie d'un non-malade chez les non-consommateurs de benzo à l'âge a

  esperance_vie_non_mal_nonconso <- matrix(c(0), # matrice des espérances de vie à chaque âge
                                           nrow=100-65+1,
                                           ncol=5,
                                           byrow = T);
  esperance_vie_non_mal_nonconso[,1] <- esp_vie_non_mal_nonconso[,1]
  esperance_vie_non_mal_nonconso[,2] <- esp_vie_non_mal_nonconso[,2]
  esperance_vie_non_mal_nonconso[,3] <- apply(esp_vie_non_mal_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  esperance_vie_non_mal_nonconso[,4] <- esperance_vie_non_mal_nonconso[,2] - 1.96*esperance_vie_non_mal_nonconso[,3];
  esperance_vie_non_mal_nonconso[,5] <- esperance_vie_non_mal_nonconso[,2] + 1.96*esperance_vie_non_mal_nonconso[,3];

  for (i in 1:nrow(esperance_vie_non_mal_nonconso)) {
    if (esperance_vie_non_mal_nonconso[i,4]<0  & is.na(esperance_vie_non_mal_nonconso[i,4])==F) {
      esperance_vie_non_mal_nonconso[i,4] <- 0
    }
  }

  for (i in 1:nrow(esperance_vie_non_mal_nonconso)) {
    if (esperance_vie_non_mal_nonconso[i,2]==0 & esperance_vie_non_mal_nonconso[i,4]==0 & esperance_vie_non_mal_nonconso[i,5]==0 &
        is.na(esperance_vie_non_mal_nonconso[i,2])==F & is.na(esperance_vie_non_mal_nonconso[i,4])==F & is.na(esperance_vie_non_mal_nonconso[i,5])==F) {
      esperance_vie_non_mal_nonconso[i,2] <- NA;
      esperance_vie_non_mal_nonconso[i,3] <- NA;
      esperance_vie_non_mal_nonconso[i,4] <- NA;
      esperance_vie_non_mal_nonconso[i,5] <- NA
    }
  }

  colnames(esperance_vie_non_mal_nonconso) <- c("Age a","Espérance de vie d'un non-malade à l'âge a","std","Borne inf","Borne sup");

  ### Prévalence de la démence

  prevalence <- prevalence[-1,]

  prev <- colSums(prevalence[,-1]); # nombre de prévalents estimés à chaque itération

  nombre_prevalence <- matrix(c(0), # matrice de prévalence
                              nrow=1,
                              ncol=5,
                              byrow = T);
  nombre_prevalence[1] <- as.integer(t);
  nombre_prevalence[2] <- prev[1];
  nombre_prevalence[3] <- sd(prev[-1], na.rm=T);
  nombre_prevalence[4] <- nombre_prevalence[2] - 1.96*nombre_prevalence[3];
  nombre_prevalence[5] <- nombre_prevalence[2] + 1.96*nombre_prevalence[3];

  if (nombre_prevalence[4]<0 & is.na(nombre_prevalence[4])==F) {
    nombre_prevalence[4] <- 0
  };

  if (nombre_prevalence[2]==0 & nombre_prevalence[4]==0 & nombre_prevalence[5]==0 &
      is.na(nombre_prevalence[2])==F & is.na(nombre_prevalence[4])==F & is.na(nombre_prevalence[5])==F) {
    nombre_prevalence[2] <- NA;
    nombre_prevalence[3] <- NA;
    nombre_prevalence[4] <- NA;
    nombre_prevalence[5] <- NA
  };

  colnames(nombre_prevalence) <- c("Année de projection","Nombre de prévalents","std","Borne inf","Borne sup");
  nombre_prevalence <- nombre_prevalence[,-1];

  ### Prévalence de la démence par âge

  nombre_prev_age <- matrix(c(0), # matrice du nombre de prévalents par âge
                            nrow=99-66+1,
                            ncol=5,
                            byrow = T);
  nombre_prev_age[,1] <- prevalence[,1]
  nombre_prev_age[,2] <- prevalence[,2]
  nombre_prev_age[,3] <- apply(prevalence[,-c(1:2)], 1, sd, na.rm=T);
  nombre_prev_age[,4] <- nombre_prev_age[,2] - 1.96*nombre_prev_age[,3];
  nombre_prev_age[,5] <- nombre_prev_age[,2] + 1.96*nombre_prev_age[,3];

  for (i in 1:nrow(nombre_prev_age)) {
    if (nombre_prev_age[i,4]<0 & is.na(nombre_prev_age[i,4])==F) {
      nombre_prev_age[i,4] <- 0
    }
  }

  for (i in 1:nrow(nombre_prev_age)) {
    if (nombre_prev_age[i,2]==0 & nombre_prev_age[i,4]==0 & nombre_prev_age[i,5]==0 &
        is.na(nombre_prev_age[i,2])==F & is.na(nombre_prev_age[i,4])==F & is.na(nombre_prev_age[i,5])==F) {
      nombre_prev_age[i,2] <- NA;
      nombre_prev_age[i,3] <- NA;
      nombre_prev_age[i,4] <- NA;
      nombre_prev_age[i,5] <- NA
    }
  }

  colnames(nombre_prev_age) <- c("Age a","Nombre de prévalents de l'âge a","std","Borne inf","Borne sup");

  ### Taux de prévalence de la démence estimé

  taux_prevalence <- taux_prevalence[-1,]

  taux_prevalence_dem <- matrix(c(0), # matrice des taux de prévalence estimé
                                nrow=99-66+1,
                                ncol=5,
                                byrow = T);
  taux_prevalence_dem[,1] <- taux_prevalence[,1]
  taux_prevalence_dem[,2] <- taux_prevalence[,2]
  taux_prevalence_dem[,3] <- apply(taux_prevalence[,-c(1:2)], 1, sd, na.rm=T);
  taux_prevalence_dem[,4] <- taux_prevalence_dem[,2] - 1.96*taux_prevalence_dem[,3];
  taux_prevalence_dem[,5] <- taux_prevalence_dem[,2] + 1.96*taux_prevalence_dem[,3];

  for (i in 1:nrow(taux_prevalence_dem)) {
    if (taux_prevalence_dem[i,4]<0 & is.na(taux_prevalence_dem[i,4])==F) {
      taux_prevalence_dem[i,4] <- 0
    }
  }

  for (i in 1:nrow(taux_prevalence_dem)) {
    if (taux_prevalence_dem[i,5]>1 & is.na(taux_prevalence_dem[i,5])==F) {
      taux_prevalence_dem[i,5] <- 1
    }
  }

  for (i in 1:nrow(taux_prevalence_dem)) {
    if (taux_prevalence_dem[i,2]==0 & taux_prevalence_dem[i,4]==0 & taux_prevalence_dem[i,5]==0 &
        is.na(taux_prevalence_dem[i,2])==F & is.na(taux_prevalence_dem[i,4])==F & is.na(taux_prevalence_dem[i,5])==F) {
      taux_prevalence_dem[i,2] <- NA;
      taux_prevalence_dem[i,3] <- NA;
      taux_prevalence_dem[i,4] <- NA;
      taux_prevalence_dem[i,5] <- NA
    }
  }

  colnames(taux_prevalence_dem) <- c("Age a","Taux de prévalence de la démence à l'âge a","std","Borne inf","Borne sup");

  ### Survie

  surv <- colSums(survie[,-1]); # nombre de survivants estimés à chaque itération

  nombre_survie <- matrix(c(0), # matrice de survie
                          nrow=1,
                          ncol=5,
                          byrow = T);
  nombre_survie[1] <- as.integer(t);
  nombre_survie[2] <- surv[1]
  nombre_survie[3] <- sd(surv[-1], na.rm=T);
  nombre_survie[4] <- nombre_survie[2] - 1.96*nombre_survie[3];
  nombre_survie[5] <- nombre_survie[2] + 1.96*nombre_survie[3];

  if (nombre_survie[4]<0 & is.na(nombre_survie[4])==F) {
    nombre_survie[4] <- 0
  };

  if (nombre_survie[2]==0 & nombre_survie[4]==0 & nombre_survie[5]==0 &
      is.na(nombre_survie[2])==F & is.na(nombre_survie[4])==F & is.na(nombre_survie[5])==F) {
    nombre_survie[2] <- NA;
    nombre_survie[3] <- NA;
    nombre_survie[4] <- NA;
    nombre_survie[5] <- NA
  };

  colnames(nombre_survie) <- c("Année de projection","Nombre de survivants","std","Borne inf","Borne sup");
  nombre_survie <- nombre_survie[,-1];

  ### Survie par âge

  nombre_surv_age <- matrix(c(0), # matrice des survivants à chaque âge
                            nrow=99-65+1,
                            ncol=5,
                            byrow = T);
  nombre_surv_age[,1] <- survie[,1]
  nombre_surv_age[,2] <- survie[,2]
  nombre_surv_age[,3] <- apply(survie[,-c(1:2)], 1, sd, na.rm=T);
  nombre_surv_age[,4] <- nombre_surv_age[,2] - 1.96*nombre_surv_age[,3];
  nombre_surv_age[,5] <- nombre_surv_age[,2] + 1.96*nombre_surv_age[,3];

  for (i in 1:nrow(nombre_surv_age)) {
    if (nombre_surv_age[i,4]<0 & is.na(nombre_surv_age[i,4])==F) {
      nombre_surv_age[i,4] <- 0
    }
  }

  for (i in 1:nrow(nombre_surv_age)) {
    if (nombre_surv_age[i,2]==0 & nombre_surv_age[i,4]==0 & nombre_surv_age[i,5]==0 &
        is.na(nombre_surv_age[i,2])==F & is.na(nombre_surv_age[i,4])==F & is.na(nombre_surv_age[i,5])==F) {
      nombre_surv_age[i,2] <- NA;
      nombre_surv_age[i,3] <- NA;
      nombre_surv_age[i,4] <- NA;
      nombre_surv_age[i,5] <- NA
    }
  }

  colnames(nombre_surv_age) <- c("Age a","Nombre de survivants de l'âge a","std","Borne inf","Borne sup");

  ### Taux de survie estimé

  taux_survivants <- taux_survivants[-1,]

  taux_survie <- matrix(c(0), # matrice des taux de survie estimé
                        nrow=99-66+1,
                        ncol=5,
                        byrow = T);
  taux_survie[,1] <- taux_survivants[,1]
  taux_survie[,2] <- taux_survivants[,2]
  taux_survie[,3] <- apply(taux_survivants[,-c(1:2)], 1, sd, na.rm=T);
  taux_survie[,4] <- taux_survie[,2] - 1.96*taux_survie[,3];
  taux_survie[,5] <- taux_survie[,2] + 1.96*taux_survie[,3];

  for (i in 1:nrow(taux_survie)) {
    if (taux_survie[i,4]<0 & is.na(taux_survie[i,4])==F) {
      taux_survie[i,4] <- 0
    }
  }

  for (i in 1:nrow(taux_survie)) {
    if (taux_survie[i,5]>1 & is.na(taux_survie[i,5])==F) {
      taux_survie[i,5] <- 1
    }
  }

  for (i in 1:nrow(taux_survie)) {
    if (taux_survie[i,2]==0 & taux_survie[i,4]==0 & taux_survie[i,5]==0 &
        is.na(taux_survie[i,2])==F & is.na(taux_survie[i,4])==F & is.na(taux_survie[i,5])==F) {
      taux_survie[i,2] <- NA;
      taux_survie[i,3] <- NA;
      taux_survie[i,4] <- NA;
      taux_survie[i,5] <- NA
    }
  }

  colnames(taux_survie) <- c("Age a","Taux de survie à l'âge a","std","Borne inf","Borne sup");

  ### Taux de prévalence pour une année t

  taux_prev_demence <- matrix(c(0), # matrice des taux de prévalence de la démence
                              nrow=1,
                              ncol=length(prev)+1,
                              byrow = T);

  taux_prev_demence[1] <- t;

  for (i in 2:ncol(taux_prev_demence)) {

    taux_prev_demence[1,i] <- prev[i-1] / surv[i-1]; # proba vie entière d'apparition de la démence

  };

  taux_moyen_demence <- matrix(c(0), # matrice de moyenne des taux de prévalence de la démence
                               nrow=1,
                               ncol=5,
                               byrow = T);
  taux_moyen_demence[1] <- as.integer(t);
  taux_moyen_demence[2] <- taux_prev_demence[,2];
  taux_moyen_demence[3] <- sd(taux_prev_demence[,-c(1:2)], na.rm=T);
  taux_moyen_demence[4] <- taux_moyen_demence[2] - 1.96*taux_moyen_demence[3];
  taux_moyen_demence[5] <- taux_moyen_demence[2] + 1.96*taux_moyen_demence[3];

  if (taux_moyen_demence[4]<0 & is.na(taux_moyen_demence[4])==F) {
    taux_moyen_demence[4] <- 0
  };

  if (taux_moyen_demence[2]==0 & taux_moyen_demence[4]==0 & taux_moyen_demence[5]==0 &
      is.na(taux_moyen_demence[2])==F & is.na(taux_moyen_demence[4])==F & is.na(taux_moyen_demence[5])==F) {
    taux_moyen_demence[2] <- NA;
    taux_moyen_demence[3] <- NA;
    taux_moyen_demence[4] <- NA;
    taux_moyen_demence[5] <- NA
  };

  colnames(taux_moyen_demence) <- c("Année de projection","Taux prévalence de la démence","std","Borne inf","Borne sup");
  taux_moyen_demence <- taux_moyen_demence[,-1];

  ### Nombre moyen d'années passées avec la maladie à l'âge a

  nombre_moy_dem <- matrix(c(0), # matrice du nombre moyen d'années passées en démence
                           nrow=100-65+1,
                           ncol=5,
                           byrow = T);
  nombre_moy_dem[,1] <- nb_moy_dem[,1]
  nombre_moy_dem[,2] <- nb_moy_dem[,2]
  nombre_moy_dem[,3] <- apply(nb_moy_dem[,-c(1:2)], 1, sd, na.rm=T);
  nombre_moy_dem[,4] <- nombre_moy_dem[,2] - 1.96*nombre_moy_dem[,3];
  nombre_moy_dem[,5] <- nombre_moy_dem[,2] + 1.96*nombre_moy_dem[,3];

  for (i in 1:nrow(nombre_moy_dem)) {
    if (nombre_moy_dem[i,4]<0 & is.na(nombre_moy_dem[i,4])==F) {
      nombre_moy_dem[i,4] <- 0
    }
  }

  for (i in 1:nrow(nombre_moy_dem)) {
    if (nombre_moy_dem[i,2]==0 & nombre_moy_dem[i,4]==0 & nombre_moy_dem[i,5]==0 &
        is.na(nombre_moy_dem[i,2])==F  & is.na(nombre_moy_dem[i,4])==F  & is.na(nombre_moy_dem[i,5])==F) {
      nombre_moy_dem[i,2] <- NA;
      nombre_moy_dem[i,3] <- NA;
      nombre_moy_dem[i,4] <- NA;
      nombre_moy_dem[i,5] <- NA
    }
  }

  colnames(nombre_moy_dem) <- c("Age a","Nombre moyen d'années en démence à l'âge a","std","Borne inf","Borne sup");

  ### Nombre moyen d'années passées avec la maladie chez les consommateurs de benzo à l'âge a

  nombre_moy_dem_conso <- matrix(c(0), # matrice du nombre moyen d'années passées en démence
                                 nrow=100-65+1,
                                 ncol=5,
                                 byrow = T);
  nombre_moy_dem_conso[,1] <- nb_moy_dem_conso[,1]
  nombre_moy_dem_conso[,2] <- nb_moy_dem_conso[,2]
  nombre_moy_dem_conso[,3] <- apply(nb_moy_dem_conso[,-c(1:2)], 1, sd, na.rm=T);
  nombre_moy_dem_conso[,4] <- nombre_moy_dem_conso[,2] - 1.96*nombre_moy_dem_conso[,3];
  nombre_moy_dem_conso[,5] <- nombre_moy_dem_conso[,2] + 1.96*nombre_moy_dem_conso[,3];

  for (i in 1:nrow(nombre_moy_dem_conso)) {
    if (nombre_moy_dem_conso[i,4]<0 & is.na(nombre_moy_dem_conso[i,4])==F) {
      nombre_moy_dem_conso[i,4] <- 0
    }
  }

  for (i in 1:nrow(nombre_moy_dem_conso)) {
    if (nombre_moy_dem_conso[i,2]==0 & nombre_moy_dem_conso[i,4]==0 & nombre_moy_dem_conso[i,5]==0 &
        is.na(nombre_moy_dem_conso[i,2])==F & is.na(nombre_moy_dem_conso[i,4])==F & is.na(nombre_moy_dem_conso[i,5])==F) {
      nombre_moy_dem_conso[i,2] <- NA;
      nombre_moy_dem_conso[i,3] <- NA;
      nombre_moy_dem_conso[i,4] <- NA;
      nombre_moy_dem_conso[i,5] <- NA
    }
  }

  colnames(nombre_moy_dem_conso) <- c("Age a","Nombre moyen d'années en démence à l'âge a","std","Borne inf","Borne sup");

  ### Nombre moyen d'années passées avec la maladie chez les non-consommateurs de benzo à l'âge a

  nombre_moy_dem_nonconso <- matrix(c(0), # matrice du nombre moyen d'années passées en démence
                                    nrow=100-65+1,
                                    ncol=5,
                                    byrow = T);
  nombre_moy_dem_nonconso[,1] <- nb_moy_dem_nonconso[,1]
  nombre_moy_dem_nonconso[,2] <- nb_moy_dem_nonconso[,2]
  nombre_moy_dem_nonconso[,3] <- apply(nb_moy_dem_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  nombre_moy_dem_nonconso[,4] <- nombre_moy_dem_nonconso[,2] - 1.96*nombre_moy_dem_nonconso[,3];
  nombre_moy_dem_nonconso[,5] <- nombre_moy_dem_nonconso[,2] + 1.96*nombre_moy_dem_nonconso[,3];

  for (i in 1:nrow(nombre_moy_dem_nonconso)) {
    if (nombre_moy_dem_nonconso[i,4]<0 & is.na(nombre_moy_dem_nonconso[i,4])==F) {
      nombre_moy_dem_nonconso[i,4] <- 0
    }
  }

  for (i in 1:nrow(nombre_moy_dem_nonconso)) {
    if (nombre_moy_dem_nonconso[i,2]==0 & nombre_moy_dem_nonconso[i,4]==0 & nombre_moy_dem_nonconso[i,5]==0 &
        is.na(nombre_moy_dem_nonconso[i,2])==F & is.na(nombre_moy_dem_nonconso[i,4])==F & is.na(nombre_moy_dem_nonconso[i,5])==F) {
      nombre_moy_dem_nonconso[i,2] <- NA;
      nombre_moy_dem_nonconso[i,3] <- NA;
      nombre_moy_dem_nonconso[i,4] <- NA;
      nombre_moy_dem_nonconso[i,5] <- NA
    }
  }

  colnames(nombre_moy_dem_nonconso) <- c("Age a","Nombre moyen d'années en démence à l'âge a","std","Borne inf","Borne sup");

  ### Proba vie entière d'apparition de la démence

  prb_dem <- prb_dem[-1,]

  prb_demence <- matrix(c(0), # matrice d'âge moyen d'apparition de la démence
                        nrow=1,
                        ncol=ncol(prb_dem),
                        byrow = T);

  prb_demence[1] <- t;

  for (i in 2:ncol(prb_dem)) {

    prb_demence[i] <- sum(prb_dem[,i]) / nrow(etat); # proba vie entière d'apparition de la démence

  };

  prb_moyen_demence <- matrix(c(0), # matrice proba vie entière d'âge moyen d'apparition de la démence
                              nrow=1,
                              ncol=5,
                              byrow = T);
  prb_moyen_demence[1] <- as.integer(t);
  prb_moyen_demence[2] <- prb_demence[,2];
  prb_moyen_demence[3] <- sd(prb_demence[,-c(1:2)], na.rm=T);
  prb_moyen_demence[4] <- prb_moyen_demence[2] - 1.96*prb_moyen_demence[3];
  prb_moyen_demence[5] <- prb_moyen_demence[2] + 1.96*prb_moyen_demence[3];

  if (prb_moyen_demence[4]<0 & is.na(prb_moyen_demence[4])==F) {
    prb_moyen_demence[4] <- 0
  };

  if (prb_moyen_demence[2]==0 & prb_moyen_demence[4]==0 & prb_moyen_demence[5]==0 &
      is.na(prb_moyen_demence[2])==F & is.na(prb_moyen_demence[4])==F & is.na(prb_moyen_demence[5])==F) {
    prb_moyen_demence[2] <- NA;
    prb_moyen_demence[3] <- NA;
    prb_moyen_demence[4] <- NA;
    prb_moyen_demence[5] <- NA
  };

  colnames(prb_moyen_demence) <- c("Année de projection","Proba d'apparition de la démence","std","Borne inf","Borne sup");
  prb_moyen_demence <- prb_moyen_demence[,-1];

  ### Age moyen d'apparition de la démence

  age_dem <- age_dem[-1,]

  age_demence <- matrix(c(0), # matrice d'âge moyen d'apparition de la démence
                        nrow=1,
                        ncol=ncol(age_dem),
                        byrow = T);

  age_demence[1] <- t;

  n0 <- vector(length = nrow(age_dem));

  for (i in 2:ncol(age_dem)) {

    for (j in 1:nrow(age_dem)) {

      n0[j] <- age_dem[j,1]*age_dem[j,i]

    };

    age_demence[i] <- sum(n0)/sum(age_dem[,i]); # âge moyen d'apparition de la démence

  }

  age_moyen_demence <- matrix(c(0), # matrice d'âge moyen d'apparition de la démence
                              nrow=1,
                              ncol=5,
                              byrow = T);
  age_moyen_demence[1] <- as.integer(t);
  age_moyen_demence[2] <- age_demence[,2]
  age_moyen_demence[3] <- sd(age_demence[,-c(1:2)], na.rm=T);
  age_moyen_demence[4] <- age_moyen_demence[2] - 1.96*age_moyen_demence[3];
  age_moyen_demence[5] <- age_moyen_demence[2] + 1.96*age_moyen_demence[3];

  if (age_moyen_demence[4]<0 & is.na(age_moyen_demence[4])==F) {
    age_moyen_demence[4] <- 0
  };

  if (age_moyen_demence[2]==0 & age_moyen_demence[4]==0 & age_moyen_demence[5]==0 &
      is.na(age_moyen_demence[2])==F & is.na(age_moyen_demence[4])==F & is.na(age_moyen_demence[5])==F) {
    age_moyen_demence[2] <- NA;
    age_moyen_demence[3] <- NA;
    age_moyen_demence[4] <- NA;
    age_moyen_demence[5] <- NA
  };

  colnames(age_moyen_demence) <- c("Année de projection","Age moyen d'apparition de la démence","std","Borne inf","Borne sup");
  age_moyen_demence <- age_moyen_demence[,-1];

  ### Nombre moyen d'années de consommation de benzo

  moy_conso <- moy_conso[-1,]

  moyenne_conso <- matrix(c(0), # matrice du nombre d'années de consommation de benzo
                          nrow=1,
                          ncol=ncol(moy_conso),
                          byrow = T);

  moyenne_conso[1] <- t;

  for (i in 2:ncol(moy_conso)) {

    moyenne_conso[i] <- sum(moy_conso[,i]) / nrow(etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),]); # nombre moyen d'années de consommation de benzo

  }

  moyenne_conso_benzo <- matrix(c(0), # matrice du nombre moyen d'années de consommation de benzo
                                nrow=1,
                                ncol=5,
                                byrow = T);
  moyenne_conso_benzo[1] <- as.integer(t);
  moyenne_conso_benzo[2] <- moyenne_conso[,2]
  moyenne_conso_benzo[3] <- sd(moyenne_conso[,-c(1:2)], na.rm=T);
  moyenne_conso_benzo[4] <- moyenne_conso_benzo[2] - 1.96*moyenne_conso_benzo[3];
  moyenne_conso_benzo[5] <- moyenne_conso_benzo[2] + 1.96*moyenne_conso_benzo[3];

  if (moyenne_conso_benzo[4]<0 & is.na(moyenne_conso_benzo[4])==F) {
    moyenne_conso_benzo[4] <- 0
  };

  if (moyenne_conso_benzo[2]==0 & moyenne_conso_benzo[4]==0 & moyenne_conso_benzo[5]==0 &
      is.na(moyenne_conso_benzo[2])==F & is.na(moyenne_conso_benzo[4])==F & is.na(moyenne_conso_benzo[5])==F) {
    moyenne_conso_benzo[2] <- NA;
    moyenne_conso_benzo[3] <- NA;
    moyenne_conso_benzo[4] <- NA;
    moyenne_conso_benzo[5] <- NA
  };

  colnames(moyenne_conso_benzo) <- c("Année de projection","Nombre moyen d'années de consommation de benzo","std","Borne inf","Borne sup");
  moyenne_conso_benzo <- moyenne_conso_benzo[,-1];

  ### Prévalence de la conso estimée

  prevalence_conso_benzo <- matrix(c(0), # matrice des prévalence de la conso estimée
                                   nrow=105-65+1,
                                   ncol=5,
                                   byrow = T);
  prevalence_conso_benzo[,1] <- prevalence_conso[,1]
  prevalence_conso_benzo[,2] <- prevalence_conso[,2]
  prevalence_conso_benzo[,3] <- apply(prevalence_conso[,-c(1:2)], 1, sd, na.rm=T);
  prevalence_conso_benzo[,4] <- prevalence_conso_benzo[,2] - 1.96*prevalence_conso_benzo[,3];
  prevalence_conso_benzo[,5] <- prevalence_conso_benzo[,2] + 1.96*prevalence_conso_benzo[,3];

  for (i in 1:nrow(prevalence_conso_benzo)) {
    if (prevalence_conso_benzo[i,4]<0 & is.na(prevalence_conso_benzo[i,4])==F) {
      prevalence_conso_benzo[i,4] <- 0
    } else {
      if (prevalence_conso_benzo[i,5]>1 & is.na(prevalence_conso_benzo[i,5])==F) {
        prevalence_conso_benzo[i,5] <- 1
      }
    }
  }

  for (i in 1:nrow(prevalence_conso_benzo)) {
    if (prevalence_conso_benzo[i,2]==0 & prevalence_conso_benzo[i,4]==0 & prevalence_conso_benzo[i,5]==0 &
        is.na(prevalence_conso_benzo[i,2])==F & is.na(prevalence_conso_benzo[i,4])==F & is.na(prevalence_conso_benzo[i,5])==F) {
      prevalence_conso_benzo[i,2] <- NA;
      prevalence_conso_benzo[i,3] <- NA;
      prevalence_conso_benzo[i,4] <- NA;
      prevalence_conso_benzo[i,5] <- NA
    }
  }

  colnames(prevalence_conso_benzo) <- c("Age a","Prévalence de la consommation de benzo à l'âge a","std","Borne inf","Borne sup");

  ### Quotient de mortalité estimé

  quotient_mortalite <- quotient_mortalite[-1,]

  quotient_de_mortalite <- matrix(c(0), # matrice des quotients de mortalité estimé
                                  nrow=99-66+1,
                                  ncol=5,
                                  byrow = T);
  quotient_de_mortalite[,1] <- quotient_mortalite[,1]
  quotient_de_mortalite[,2] <- quotient_mortalite[,2]
  quotient_de_mortalite[,3] <- apply(quotient_mortalite[,-c(1:2)], 1, sd, na.rm=T);
  quotient_de_mortalite[,4] <- quotient_de_mortalite[,2] - 1.96*quotient_de_mortalite[,3];
  quotient_de_mortalite[,5] <- quotient_de_mortalite[,2] + 1.96*quotient_de_mortalite[,3];

  for (i in 1:nrow(quotient_de_mortalite)) {
    if (quotient_de_mortalite[i,4]<0 & is.na(quotient_de_mortalite[i,4])==F) {
      quotient_de_mortalite[i,4] <- 0
    }
  }

  for (i in 1:nrow(quotient_de_mortalite)) {
    if (quotient_de_mortalite[i,5]>1 & is.na(quotient_de_mortalite[i,5])==F) {
      quotient_de_mortalite[i,5] <- 1
    }
  }

  for (i in 1:nrow(quotient_de_mortalite)) {
    if (quotient_de_mortalite[i,2]==0 & quotient_de_mortalite[i,4]==0 & quotient_de_mortalite[i,5]==0 &
        is.na(quotient_de_mortalite[i,2])==F & is.na(quotient_de_mortalite[i,4])==F & is.na(quotient_de_mortalite[i,5])==F) {
      quotient_de_mortalite[i,2] <- NA;
      quotient_de_mortalite[i,3] <- NA;
      quotient_de_mortalite[i,4] <- NA;
      quotient_de_mortalite[i,5] <- NA
    }
  }

  colnames(quotient_de_mortalite) <- c("Age a","Quotient de mortalité à l'âge a","std","Borne inf","Borne sup");

  ### Listes de sorties de l'algo

  ### Risques de devenir déments et de dédéder chez les non-consommateurs de benzo

  list_risques <- list(a010, a020, a120)
  names(list_risques) <- c("a010", "a020", "a120")

  ### Espérance de vie générale

  list_esp_vie <- list(esperance_vie_gen, esperance_vie_gen_conso, esperance_vie_gen_nonconso)
  names(list_esp_vie) <- c("esperance_vie_gen", "esperance_vie_gen_conso", "esperance_vie_gen_nonconso")

  ### Espérance de vie sans la maladie

  list_esp_vie_sans_mal <- list(esperance_vie_sans_mal, esperance_vie_sans_mal_conso, esperance_vie_sans_mal_nonconso)
  names(list_esp_vie_sans_mal) <- c("esperance_vie_sans_mal", "esperance_vie_sans_mal_conso", "esperance_vie_sans_mal_nonconso")

  ### Espérance de vie d'un malade

  list_esp_vie_mal <- list(esperance_vie_mal, esperance_vie_mal_conso, esperance_vie_mal_nonconso)
  names(list_esp_vie_mal) <- c("esperance_vie_mal", "esperance_vie_mal_conso", "esperance_vie_mal_nonconso")

  ### Espérance de vie d'un non-malade

  list_esp_vie_non_mal <- list(esperance_vie_non_mal, esperance_vie_non_mal_conso, esperance_vie_non_mal_nonconso)
  names(list_esp_vie_non_mal) <- c("esperance_vie_non_mal", "esperance_vie_non_mal_conso", "esperance_vie_non_mal_nonconso")

  ### Prévalence de la démence

  list_prevalence <- list(nombre_prev_age, nombre_prevalence, taux_prevalence_dem, taux_moyen_demence)
  names(list_prevalence) <- c("nombre_prev_age", "nombre_prevalence", "taux_prevalence_dem", "taux_moyen_demence")

  ### Survie

  list_survie <- list(nombre_surv_age, nombre_survie, taux_survie)
  names(list_survie) <- c("nombre_surv_age", "nombre_survie", "taux_survie")

  ### Nombre moyen d'années passées en démence

  list_nombre_moy_dem <- list(nombre_moy_dem, nombre_moy_dem_conso, nombre_moy_dem_nonconso)
  names(list_nombre_moy_dem) <- c("nombre_moy_dem", "nombre_moy_dem_conso", "nombre_moy_dem_nonconso")

  ### Liste finale

  list_sortie <- list(list_risques, list_esp_vie, list_esp_vie_sans_mal, list_esp_vie_mal, list_esp_vie_non_mal, list_prevalence, list_survie, list_nombre_moy_dem, prb_moyen_demence, age_moyen_demence, moyenne_conso_benzo, prevalence_conso_benzo, quotient_de_mortalite);
  names(list_sortie) <- c("list_risques", "list_esp_vie", "list_esp_vie_sans_mal", "list_esp_vie_mal", "list_esp_vie_non_mal", "list_prevalence", "list_survie", "list_nombre_moy_dem", "prb_moyen_demence", "age_moyen_demence", "moyenne_conso_benzo", "prevalence_conso_benzo", "quotient_de_mortalite");

  return(list_sortie)

}
