#' Computation of Health Indicators
#'
#' This function computes many health indicators under several interventions of intervention in risk factor distribution for a given year.
#'
#' @param t year of the projections for health indicators.
#' @param intervention 0 = no change; 1 = reduction by two of risk factor distribution; 2 = risk factor distribution considered as null. Default is \code{0}.
#' @param year_intervention year of the intervention in risk factor distribution takes place. Default is \code{NULL}.
#' @param nb_people number of people whose trajectory will be simulated for each generation. Default is \code{100}.
#' @param nb_iter number of iterations for the algorithm. Default is \code{0}.
#' @param data_pop data source for demographics data.
#' @param gender gender for computation. \code{"W"} for women and \code{"M"} for men. Default is \code{"W"}.
#' @param data_a01 data source for the incidence of disease.
#' @param data_a02 data source for the mortality of healthy subjects.
#' @param data_theta01 data source for the relative risks associated with the exposure for disease.
#' @param data_theta02 data source for the relative risks associated with the exposure for mortality among healthy subjects.
#' @param data_theta12 data source for the relative risks associated with the exposure for mortality among diseased subjects.
#' @param data_prev data source for the prevalence of the exposition.
#' @param data_incid data source for the incidence of the exposition.
#' @param data_rr_DvsND data source for the relative risks associated with the disease for mortality.
#' @param data_a01_values variability of data source for the incidence of disease.
#' @param data_a02_values variability of data source for the mortality of healthy subjects.
#' @param data_theta01_values variability of data source for the relative risks associated with the exposure for disease.
#' @param data_prev_values variability of data source for the prevalence of the exposition.
#' @param data_incid_values variability of data source for the incidence of the exposition.
#' @param data_rr_DvsND_values variability of data source for the relative risks associated with the disease for mortality.
#'
#' @return a list containing the health indicators
#'
#' @export
#'
#' @examples
#' estimHI(t = 2040,
#' intervention = 1,
#' year_intervention = 2020,
#' nb_people = 10000,
#' nb_iter = 100,
#' data_pop = pop,
#' gender = "W",
#' data_a01 = a01,
#' data_a02 = a02,
#' data_theta01 = theta01_cas_1_6,
#' data_theta02 = theta02_1,
#' data_theta12 = theta02_1,
#' data_prev = prevconso,
#' data_incid = incidconso,
#' data_rr_DvsND = rr_DvsND,
#' data_a01_values = a01_values,
#' data_a02_values = a02_values,
#' data_theta01_values = theta01_cas_1_6_values,
#' data_prev_values <- prevconso_values,
#' data_incid_values <- incidconso_values,
#' data_rr_DvsND_values = rr_DvsND_values)
estimHI <- function(t,
                    intervention = 0,
                    year_intervention = NULL,
                    nb_people = 100,
                    nb_iter = 0,
                    data_pop,
                    gender = "W",
                    data_a01,
                    data_a02,
                    data_theta01,
                    data_theta02,
                    data_theta12,
                    data_prev,
                    data_incid,
                    data_rr_DvsND,
                    data_a01_values,
                    data_a02_values,
                    data_theta01_values,
                    data_prev_values,
                    data_incid_values,
                    data_rr_DvsND_values)

{

  ### Year of projection

  year_proj <- t

  ### Incidence of disease on non exposed peoples

  a010 <- matrix(c(0),
                 nrow=40*nb_iter,
                 ncol=2,
                 byrow=T);
  colnames(a010) <- c("age","a010")

  a010[,1] <- c(66:105)
  a010[,2] <- as.numeric(data_a01[which(data_a01[,1] != 65 & data_a01[,3]%in%(gender)),2]) / (data_theta01[which(data_theta01[,1] != 65 & data_theta01[,3]%in%(gender)),2]*data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2] - data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2] + 1);

  ### Incidence of disease on exposed peoples

  a011 <- matrix(c(0),
                 nrow=40*nb_iter,
                 ncol=2,
                 byrow=T);
  colnames(a011) <- c("age","a011")

  a011[,1] <- c(66:105)
  a011[,2] <- a010[,2]*data_theta01[which(data_theta01[,1] != 65 & data_theta01[,3]%in%(gender)),2];

  ### Global incidence of disease

  a01_global <- matrix(c(0),
                       nrow=40*nb_iter,
                       ncol=2,
                       byrow=T);
  colnames(a01_global) <- c("age","a01_global")

  a01_global[,1] <- c(66:105)
  a01_global[,2] <- data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2]*a010[,2]*data_theta01[which(data_theta01[,1] != 65 & data_theta01[,3]%in%(gender)),2] + (1-data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2])*a010[,2];

  ### Mortality of healthy subjects on non exposed peoples

  a020 <- matrix(c(0),
                 nrow=40*nb_iter,
                 ncol=132,
                 byrow=T);
  colnames(a020) <- c("age",1950:2080)

  a020[,1] <- c(66:105)

  for (a in 2:ncol(a020)){

    a020[,a] <- as.numeric(data_a02[which(data_a02[,1] != 65 & data_a02[,2]%in%(gender)),a+1]) / (data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(gender)),2]*data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2] - data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2] + 1);

  }

  ### Mortality of healthy subjects on exposed peoples

  a021 <- matrix(c(0),
                 nrow=40*nb_iter,
                 ncol=132,
                 byrow=T);
  colnames(a021) <- c("age",1950:2080)

  a021[,1] <- c(66:105)
  a021[,-1] <- a020[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(gender)),2];

  ### Global mortality of healthy subjects

  a02_global <- matrix(c(0),
                       nrow=40*nb_iter,
                       ncol=132,
                       byrow=T);
  colnames(a02_global) <- c("age",1950:2080)

  a02_global[,1] <- c(66:105)
  a02_global[,-1] <- data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2]*a020[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(gender)),2] + (1-data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2])*a020[,-1];

  ### Relative risks associated with the disease for mortality

  RR <- matrix(c(0),
               nrow=40*nb_iter,
               ncol=2,
               byrow=T);
  colnames(RR) <- c("age","rr_DvsND")

  RR[,1] <- c(66:105)
  RR[,2] <- data_rr_DvsND[which(data_rr_DvsND[,1] != 65 & data_rr_DvsND[,3]%in%(gender)),2];

  ### Mortality of diseased subjects on non exposed peoples

  a120 <- matrix(c(0),
                 nrow=40*nb_iter,
                 ncol=132,
                 byrow=T);
  colnames(a120) <- c("age",1950:2080)

  a120[,1] <- c(66:105)

  for (a in 2:ncol(a020)){

    a120[,a] <- as.numeric(RR[,2])*a020[,a] / (data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(gender)),2]*data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2] - data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2] + 1);

  }

  ### Mortality of diseased subjects on exposed peoples

  a121 <- matrix(c(0),
                 nrow=40*nb_iter,
                 ncol=132,
                 byrow=T);
  colnames(a121) <- c("age",1950:2080)

  a121[,1] <- c(66:105)
  a121[,-1] <- a120[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(gender)),2];

  ### Global mortality of diseased subjects

  a12_global <- matrix(c(0),
                       nrow=40*nb_iter,
                       ncol=132,
                       byrow=T);
  colnames(a12_global) <- c("age",1950:2080)

  a12_global[,1] <- c(66:105)
  a12_global[,-1] <- data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2]*a120[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(gender)),2] + (1-data_prev[which(data_prev[,1] != 65 & data_prev[,3]%in%(gender)),2])*a120[,-1];

  ### Variability of parameters

  ### Incidence of disease on non exposed peoples

  a010_values <- matrix(c(0),
                         nrow=40,
                         ncol=2,
                         byrow=T);
  colnames(a010_values) <- c("age","a010")

  a010_values[,1] <- c(66:105)
  a010_values[,2] <- as.numeric(data_a01_values[which(data_a01_values[,1] != 65 & data_a01_values[,3]%in%(gender)),2]) / (data_theta01_values[which(data_theta01_values[,1] != 65 & data_theta01_values[,3]%in%(gender)),2]*data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2] - data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2] + 1);

  ### Incidence of disease on exposed peoples

  a011_values <- matrix(c(0),
                         nrow=40,
                         ncol=2,
                         byrow=T);
  colnames(a011_values) <- c("age","a011")

  a011_values[,1] <- c(66:105)
  a011_values[,2] <- a010_values[,2]*data_theta01_values[which(data_theta01_values[,1] != 65 & data_theta01_values[,3]%in%(gender)),2];

  ### Global incidence of disease

  a01_global_values <- matrix(c(0),
                               nrow=40,
                               ncol=2,
                               byrow=T);
  colnames(a01_global_values) <- c("age","a01_global")

  a01_global_values[,1] <- c(66:105)
  a01_global_values[,2] <- data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2]*a010_values[,2]*data_theta01_values[which(data_theta01_values[,1] != 65 & data_theta01_values[,3]%in%(gender)),2] + (1-data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2])*a010_values[,2];

  ### Mortality of healthy subjects on non exposed peoples

  a020_values <- matrix(c(0),
                         nrow=40,
                         ncol=132,
                         byrow=T);
  colnames(a020_values) <- c("age",1950:2080)

  a020_values[,1] <- c(66:105)

  for (a in 2:ncol(a020_values)){ # pour chaque année

    a020_values[,a] <- as.numeric(data_a02_values[which(data_a02_values[,1] != 65 & data_a02_values[,2]%in%(gender)),a+1]) / (data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(gender)),2]*data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2] - data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2] + 1);

  }

  ### Mortality of healthy subjects on exposed peoples

  a021_values <- matrix(c(0),
                         nrow=40,
                         ncol=132,
                         byrow=T);
  colnames(a021_values) <- c("age",1950:2080)

  a021_values[,1] <- c(66:105)
  a021_values[,-1] <- a020_values[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(gender)),2];

  ### Global mortality of healthy subjects

  a02_global_values <- matrix(c(0),
                               nrow=40,
                               ncol=132,
                               byrow=T);
  colnames(a02_global_values) <- c("age",1950:2080)

  a02_global_values[,1] <- c(66:105)
  a02_global_values[,-1] <- data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(gender)),2] + (1-data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2])*a020_values[,-1];

  ### Relative risks associated with the disease for mortality

  RR_values <- matrix(c(0),
                       nrow=40,
                       ncol=2,
                       byrow=T);
  colnames(RR_values) <- c("age","rr_DvsND")

  RR_values[,1] <- c(66:105)
  RR_values[,2] <- data_rr_DvsND_values[which(data_rr_DvsND_values[,1] != 65 & data_rr_DvsND_values[,3]%in%(gender)),2];

  ### Mortality of diseased subjects on non exposed peoples

  a120_values <- matrix(c(0),
                         nrow=40,
                         ncol=132,
                         byrow=T);
  colnames(a120_values) <- c("age",1950:2080)

  a120_values[,1] <- c(66:105)

  for (a in 2:ncol(a020_values)){

    a120_values[,a] <- as.numeric(RR_values[,2])*a020_values[,a] / (data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(gender)),2]*data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2] - data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2] + 1);

  }

  ### Mortality of diseased subjects on exposed peoples

  a121_values <- matrix(c(0),
                         nrow=40,
                         ncol=132,
                         byrow=T);
  colnames(a121_values) <- c("age",1950:2080)

  a121_values[,1] <- c(66:105)
  a121_values[,-1] <- a120_values[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(gender)),2];

  ### Global mortality of diseased subjects

  a12_global_values <- matrix(c(0),
                               nrow=40,
                               ncol=132,
                               byrow=T);
  colnames(a12_global_values) <- c("age",1950:2080)

  a12_global_values[,1] <- c(66:105)
  a12_global_values[,-1] <- data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(gender)),2] + (1-data_prev_values[which(data_prev_values[,1] != 65 & data_prev_values[,3]%in%(gender)),2])*a120_values[,-1];

  ### Matrix for results of health indicators

  ### Overall life-expectancy

  esp_vie_gen <- matrix(c(0),
                        nrow=100-65+1,
                        ncol=2+nb_iter,
                        byrow = T);
  esp_vie_gen[,1] <- c(65:100);

  ### Overall life-expectancy on exposed peoples

  esp_vie_gen_conso <- matrix(c(0),
                              nrow=100-65+1,
                              ncol=2+nb_iter,
                              byrow = T);
  esp_vie_gen_conso[,1] <- c(65:100);

  ### Overall life-expectancy on non exposed peoples

  esp_vie_gen_nonconso <- matrix(c(0),
                                 nrow=100-65+1,
                                 ncol=2+nb_iter,
                                 byrow = T);
  esp_vie_gen_nonconso[,1] <- c(65:100);

  ### Life-expectancy without disease

  esp_vie_sans_mal <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                             nrow=100-65+1,
                             ncol=2+nb_iter,
                             byrow = T);
  esp_vie_sans_mal[,1] <- c(65:100);

  ### Life-expectancy without disease on exposed peoples

  esp_vie_sans_mal_conso <- matrix(c(0),
                                   nrow=100-65+1,
                                   ncol=2+nb_iter,
                                   byrow = T);
  esp_vie_sans_mal_conso[,1] <- c(65:100);

  ### Life-expectancy without disease on non exposed peoples

  esp_vie_sans_mal_nonconso <- matrix(c(0),
                                      nrow=100-65+1,
                                      ncol=2+nb_iter,
                                      byrow = T);
  esp_vie_sans_mal_nonconso[,1] <- c(65:100);

  ### Life-expectancy for diseased subject

  esp_vie_mal <- matrix(c(0),
                        nrow=100-65+1,
                        ncol=2+nb_iter,
                        byrow = T);
  esp_vie_mal[,1] <- c(65:100);

  ### Life-expectancy for diseased subject on exposed peoples

  esp_vie_mal_conso <- matrix(c(0),
                              nrow=100-65+1,
                              ncol=2+nb_iter,
                              byrow = T);
  esp_vie_mal_conso[,1] <- c(65:100);

  ### Life-expectancy for diseased subject on non exposed peoples

  esp_vie_mal_nonconso <- matrix(c(0),
                                 nrow=100-65+1,
                                 ncol=2+nb_iter,
                                 byrow = T);
  esp_vie_mal_nonconso[,1] <- c(65:100);

  ### Life-expectancy for non-diseased subject

  esp_vie_non_mal <- matrix(c(0),
                            nrow=100-65+1,
                            ncol=2+nb_iter,
                            byrow = T);
  esp_vie_non_mal[,1] <- c(65:100);

  ### Life-expectancy for non-diseased subject on exposed peoples

  esp_vie_non_mal_conso <- matrix(c(0),
                                  nrow=100-65+1,
                                  ncol=2+nb_iter,
                                  byrow = T);
  esp_vie_non_mal_conso[,1] <- c(65:100);

  ### Life-expectancy for non-diseased subject on non exposed peoples

  esp_vie_non_mal_nonconso <- matrix(c(0),
                                     nrow=100-65+1,
                                     ncol=2+nb_iter,
                                     byrow = T);
  esp_vie_non_mal_nonconso[,1] <- c(65:100);

  ### Number of exposed peoples at least one time

  prevalence <- matrix(c(0),
                       nrow=99-65+1,
                       ncol=2+nb_iter,
                       byrow=T);
  prevalence[,1] <- c(65:99);

  ### Rate of exposed peoples at least one time

  taux_prevalence <- matrix(c(0),
                            nrow=99-65+1,
                            ncol=2+nb_iter,
                            byrow=T);
  taux_prevalence[,1] <- c(65:99);

  ### Number of living peoples

  survie <- matrix(c(0),
                   nrow=99-65+1,
                   ncol=2+nb_iter,
                   byrow=T);
  survie[,1] <- c(65:99);

  ### Rate of living peoples a given

  taux_survivants <- matrix(c(0),
                            nrow=99-65+1,
                            ncol=2+nb_iter,
                            byrow=T);
  taux_survivants[,1] <- c(65:99);

  ### Mean number of years spent with disease

  nb_moy_dem <- matrix(c(0),
                       nrow=100-65+1,
                       ncol=2+nb_iter,
                       byrow = T);
  nb_moy_dem[,1] <- c(65:100);

  ### Mean number of years spent with disease on exposed peoples

  nb_moy_dem_conso <- matrix(c(0),
                             nrow=100-65+1,
                             ncol=2+nb_iter,
                             byrow = T);
  nb_moy_dem_conso[,1] <- c(65:100);

  ### Mean number of years spent with disease on non exposed peoples

  nb_moy_dem_nonconso <- matrix(c(0),
                                nrow=100-65+1,
                                ncol=2+nb_iter,
                                byrow = T);
  nb_moy_dem_nonconso[,1] <- c(65:100);

  ### Life-long probability of disease

  prb_dem <- matrix(c(0),
                    nrow=99-65+1,
                    ncol=2+nb_iter,
                    byrow=T);
  prb_dem[,1] <- c(65:99);

  ###	Average age at disease onset

  age_dem <- matrix(c(0),
                    nrow=99-65+1,
                    ncol=2+nb_iter,
                    byrow=T);
  age_dem[,1] <- c(65:99);

  ### Mean number of years of exposition

  moy_conso <- matrix(c(0),
                      nrow=99-65+1,
                      ncol=2+nb_iter,
                      byrow=T);
  moy_conso[,1] <- c(65:99);

  ### Number of exposed peoples at least one time

  prevalence_conso <- matrix(c(0),
                             nrow=105-65+1,
                             ncol=2+nb_iter,
                             byrow=T);
  prevalence_conso[,1] <- c(65:105);

  ### Mortality rate

  quotient_mortalite <- matrix(c(0),
                               nrow=99-65+1,
                               ncol=2+nb_iter,
                               byrow=T);
  quotient_mortalite[,1] <- c(65:99);

  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################

  ### Computation with estimated parameters

  ###############################
  ###          STEP 1         ###
  ###############################

  ### Computation of number of exposed peoples at least one time

  pr_conso_benzo <- matrix(c(0),
                           nrow=105-65+1,
                           ncol=2,
                           byrow=T);
  pr_conso_benzo[,1] <- c(65:105);

  for (age in 65:105) { # for each generation

    an_naiss <- year_proj-age; # year of birth

    an0 <- an_naiss + 65; # first year at risk

    annee <- an0; # year of estimation

    etat <- matrix(c(0), # initial state (non diseased)
                   nrow=nb_people,
                   ncol=105-65+1,
                   byrow=T);

    colnames(etat) <- c(65:105);

    for (i in 1:nrow(etat)) {

      alea0 <- runif(1, 0, 1);

      if (alea0 <= data_prev_values[which(data_prev_values[,1]%in%(65) & data_prev_values[,3]%in%(gender)),2]) {

        etat[i,1] <- "01" # non diseased and exposed

      } else {

        etat[i,1] <- "00" # non diseased and non exposed

      }

    };

    for (i in 1:nrow(etat)) { # for each people

      for (j in 2:ncol(etat)) { # for each age

        alea <- runif(1, 0, 1);

        alea0 <- runif(1, 0, 1);

        if (etat[i,j-1] == "00") {

          if (alea0 <= data_incid_values[which(data_incid_values[,1]%in%(j-1+65) & data_incid_values[,3]%in%(gender)),2]) {

            a01 <- a011_values;
            a02 <- a021_values;

            if (alea <= a02[j-1,j+(an0-1950)+1]) {

              etat[i,j] <- "21"; # dead (with exposed state)

            } else {

              if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "11"; # diseased and exposed

              } else {

                etat[i,j] <- "01"; # non diseased and exposed

              }

            };

          } else {

            a01 <- a010_values;
            a02 <- a020_values;

            if (alea <= a02[j-1,j+(an0-1950)+1]) {

              etat[i,j] <- "20"; # dead (with non exposed state)

            } else {

              if (alea <= a01[j-1,2] + a02[j-1,j+(an0-1950)+1]) {

                etat[i,j] <- "10"; # diseased and non exposed

              } else {

                etat[i,j] <- "00"; # non diseased and non exposed

              }

            }

          }

        } else {

          if (etat[i,j-1] == "01") {

            a01 <- a011_values;
            a02 <- a021_values;

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

              if (alea0 <= data_incid_values[which(data_incid_values[,1]%in%(j-1+65) & data_incid_values[,3]%in%(gender)),2]) {

                a12 <- a121_values;

                if (alea <= a12[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "21";

                } else {

                  etat[i,j] <- "11";

                }

              } else {

                a12 <- a120_values;

                if (alea <= a12[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "20";

                } else {

                  etat[i,j] <- "10";

                }

              }

            } else {

              if (etat[i,j-1] == "11") {

                a12 <- a121_values;

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

  a010_values[,2] <- as.numeric(data_a01_values[which(data_a01_values[,1] != 65 & data_a01_values[,3]%in%(gender)),2]) / (data_theta01_values[which(data_theta01_values[,1] != 65 & data_theta01_values[,3]%in%(gender)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

  ### Incidence of disease on exposed peoples

  a011_values[,2] <- a010_values[,2]*data_theta01_values[which(data_theta01_values[,1] != 65 & data_theta01_values[,3]%in%(gender)),2];

  ### Global incidence of disease

  a01_global_values[,2] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a010_values[,2]*data_theta01_values[which(data_theta01_values[,1] != 65 & data_theta01_values[,3]%in%(gender)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a010_values[,2];

  ### Mortality of healthy subjects on non exposed peoples

  for (a in 2:ncol(a020_values)){

    a020_values[,a] <- as.numeric(data_a02_values[which(data_a02_values[,1] != 65 & data_a02_values[,2]%in%(gender)),a+1]) / (data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(gender)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

  }

  ### Mortality of healthy subjects on exposed peoples

  a021_values[,-1] <- a020_values[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(gender)),2];

  ### Global mortality of healthy subjects

  a02_global_values[,-1] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a020_values[,-1]*data_theta02[which(data_theta02[,1] != 65 & data_theta02[,3]%in%(gender)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a020_values[,-1];

  ### Mortality of diseased subjects on non exposed peoples

  for (a in 2:ncol(a020_values)){

    a120_values[,a] <- as.numeric(RR_values[,2])*a020_values[,a] / (data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(gender)),2]*pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] - pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2] + 1);

  }

  ### Mortality of diseased subjects on exposed peoples

  a121_values[,-1] <- a120_values[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(gender)),2];

  ### Global mortality of diseased subjects

  a12_global_values[,-1] <- pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2]*a120_values[,-1]*data_theta12[which(data_theta12[,1] != 65 & data_theta12[,3]%in%(gender)),2] + (1-pr_conso_benzo[which(pr_conso_benzo[,1] != 65),2])*a120_values[,-1];

  ###############################
  ###          STEP 2         ###
  ###############################

  for (age in 65:105) {

    an_naiss <- year_proj-age;

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

          donnees_incid <- data_incid_values

        } else {

          if (intervention == 1) {

            if (annee < year_intervention) {

              donnees_incid <- data_incid_values

            } else {

              donnees_incid <- data_incid_values
              donnees_incid[,2] <- data_incid_values[,2] / 2

            }

          } else {

            if (intervention == 2) {

              if (annee < year_intervention) {

                donnees_incid <- data_incid_values

              } else {

                donnees_incid <- data_incid_values
                donnees_incid[,2] <- 0

              }

            }

          }

        }

        if (etat[i,j-1] == "00") {

          if (alea0 <= donnees_incid[which(donnees_incid[,1]%in%(j-1+65) & donnees_incid[,3]%in%(gender)),2]) {

            a01 <- a011_values;
            a02 <- a021_values;

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

            a01 <- a010_values;
            a02 <- a020_values;

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

            a01 <- a011_values;
            a02 <- a021_values;

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

                a12 <- a121_values;

                if (alea <= a12[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "21";

                } else {

                  etat[i,j] <- "11";

                }

              } else {

                a12 <- a120_values;

                if (alea <= a12[j-1,j+(an0-1950)+1]) {

                  etat[i,j] <- "20";

                } else {

                  etat[i,j] <- "10";

                }

              }

            } else {

              if (etat[i,j-1] == "11") {

                a12 <- a121_values;

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

    }

    ### Computation of health indicators :

    ### Overall life-expectancy

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

    ###	Overall life-expectancy on exposed peoples

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

    ### Overall life-expectancy on non exposed peoples

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

    ###	Life-expectancy without disease

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

    ###	Life-expectancy without disease on exposed peoples

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

    ###	Life-expectancy without disease on non exposed peoples

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

    ### Life-expectancy for diseased subject

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

    ### Life-expectancy for diseased subject on exposed peoples

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

    ### Life-expectancy for diseased subject on non exposed peoples

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

    ### Life-expectancy for non diseased subject

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

    ### Life-expectancy for non diseased subject on exposed peoples

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

    ### Life-expectancy for non diseased subject on non exposed peoples

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

    ### Prevalence of disease

    if (age > 65 & age < 100) {

      d0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("00") | etat[,age-65]%in%("01"))) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11")));

      s0 <- sum(etat[,1]%in%("00") | etat[,1]%in%("01") | etat[,1]%in%("10") | etat[,1]%in%("11"));

      s1 <- sum((etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11"))) + 0.5 * sum((etat[,age-64]%in%("10") | etat[,age-64]%in%("11")) & (etat[,age-65]%in%("00") | etat[,age-65]%in%("01"))) + 0.5 * sum((etat[,age-64]%in%("20") | etat[,age-64]%in%("21")) & (etat[,age-65]%in%("10") | etat[,age-65]%in%("11")));

      if (s0 != 0) {

        p01 <- s1/s0;

        taux_prevalence[age-64,2] <- s1/d0;

        nb <- p01*data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(gender)),2];

        prevalence[age-64,2] <- nb;

      };

    };

    ### Survival

    if (age == 65) {

      survie[age-64,2] <- data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(gender)),2]

    };

    if (age > 65 & age < 100) {

      d0 <- sum(etat[,age-65]%in%("00") | etat[,age-65]%in%("01") | etat[,age-65]%in%("10") | etat[,age-65]%in%("11"));

      s0 <- sum(etat[,1]%in%("00") | etat[,1]%in%("01") | etat[,1]%in%("10") | etat[,1]%in%("11"));

      s1 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("01") | etat[,age-64]%in%("10") | etat[,age-64]%in%("11"));

      if (s0 != 0) {

        p01 <- s1/s0;

        taux_survivants[age-64,2] <- s1/d0;

        nb <- p01*data_pop[which(data_pop[,1]%in%(an0) & data_pop[,3]%in%(gender)),2];

        survie[age-64,2] <- nb;

      };

    };

    ### Mean number of years spent with disease

    if (age < 101) {

      nb_moy_dem[age-64,2] <- esp_vie_non_mal[age-64,2] - esp_vie_sans_mal[age-64,2]

    } else {

      n0 <- NA

    };

    ### Mean number of years spent with disease on exposed peoples

    if (age < 101) {

      nb_moy_dem_conso[age-64,2] <- esp_vie_non_mal_conso[age-64,2] - esp_vie_sans_mal_conso[age-64,2];

    } else {

      n0 <- NA;

    };

    ### Mean number of years spent with disease on non exposed peoples

    if (age < 101) {

      nb_moy_dem_nonconso[age-64,2] <- esp_vie_non_mal_nonconso[age-64,2] - esp_vie_sans_mal_nonconso[age-64,2];

    } else {

      n0 <- NA;

    };

    ###	Life-long probability of disease

    if (age == 65) {

      for (i in (age-63):nrow(prb_dem)) {

        prb_dem[i,2] <- sum(etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("10") | etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("11"))

      }

    };

    ### Average age at disease onset

    if (age == 65) {

      for (i in (age-63):nrow(age_dem)) {

        age_dem[i,2] <- sum(etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("10") | etat[which(etat[,i-1]%in%("00") | etat[,i-1]%in%("01")),i]%in%("11"))

      }

    };

    ### Mean number of years of exposition

    if (age == 65) {

      for (i in (age-63):nrow(age_dem)) {

        moy_conso[i,2] <- sum(etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),i]%in%("01") | etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),i]%in%("11"))

      }

    };

    ### Number of exposed peoples at least one time

    if (age <= 105) {

      s0 <- sum(etat[,age-64]%in%("00") | etat[,age-64]%in%("10") | etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

      s1 <- sum(etat[,age-64]%in%("01") | etat[,age-64]%in%("11"));

      if (s0 != 0) {

        prevalence_conso[age-64,2] <- s1/s0;

      } else {

        prevalence_conso[age-64,2] <- NA;

      };

    };

    ### Mortality rate

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

  ### Computation for variability

  if (nb_iter != 0) {

    indicateurs <- varHI(t = t,
                         intervention = intervention,
                         year_intervention = year_intervention,
                         nb_people = nb_people,
                         nb_iter = nb_iter,
                         data_pop = data_pop,
                         gender = gender,
                         data_prev = data_prev,
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

  ### Computation of health indicators :

  ### Overall life-expectancy

  life_expectancy <- matrix(c(0),
                            nrow=100-65+1,
                            ncol=5,
                            byrow = T);
  life_expectancy[,1] <- esp_vie_gen[,1]
  life_expectancy[,2] <- esp_vie_gen[,2]
  life_expectancy[,3] <- apply(esp_vie_gen[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy[,4] <- life_expectancy[,2] - 1.96*life_expectancy[,3];
  life_expectancy[,5] <- life_expectancy[,2] + 1.96*life_expectancy[,3];

  for (i in 1:nrow(life_expectancy)) {
    if (life_expectancy[i,4]<0  & is.na(life_expectancy[i,4])==F) {
      life_expectancy[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy)) {
    if (life_expectancy[i,2]==0 & life_expectancy[i,4]==0 & life_expectancy[i,5]==0 &
        is.na(life_expectancy[i,2])==F & is.na(life_expectancy[i,4])==F & is.na(life_expectancy[i,5])==F) {
      life_expectancy[i,2] <- NA;
      life_expectancy[i,3] <- NA;
      life_expectancy[i,4] <- NA;
      life_expectancy[i,5] <- NA
    }
  }

  colnames(life_expectancy) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Overall life-expectancy on exposed peoples

  life_expectancy_exp <- matrix(c(0),
                                nrow=100-65+1,
                                ncol=5,
                                byrow = T);
  life_expectancy_exp[,1] <- esp_vie_gen_conso[,1]
  life_expectancy_exp[,2] <- esp_vie_gen_conso[,2]
  life_expectancy_exp[,3] <- apply(esp_vie_gen_conso[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_exp[,4] <- life_expectancy_exp[,2] - 1.96*life_expectancy_exp[,3];
  life_expectancy_exp[,5] <- life_expectancy_exp[,2] + 1.96*life_expectancy_exp[,3];

  for (i in 1:nrow(life_expectancy_exp)) {
    if (life_expectancy_exp[i,4]<0  & is.na(life_expectancy_exp[i,4])==F) {
      life_expectancy_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_exp)) {
    if (life_expectancy_exp[i,2]==0 & life_expectancy_exp[i,4]==0 & life_expectancy_exp[i,5]==0 &
        is.na(life_expectancy_exp[i,2])==F & is.na(life_expectancy_exp[i,4])==F & is.na(life_expectancy_exp[i,5])==F) {
      life_expectancy_exp[i,2] <- NA;
      life_expectancy_exp[i,3] <- NA;
      life_expectancy_exp[i,4] <- NA;
      life_expectancy_exp[i,5] <- NA
    }
  }

  colnames(life_expectancy_exp) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Overall life-expectancy on non exposed peoples

  life_expectancy_n_exp <- matrix(c(0),
                                  nrow=100-65+1,
                                  ncol=5,
                                  byrow = T);
  life_expectancy_n_exp[,1] <- esp_vie_gen_nonconso[,1]
  life_expectancy_n_exp[,2] <- esp_vie_gen_nonconso[,2]
  life_expectancy_n_exp[,3] <- apply(esp_vie_gen_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_n_exp[,4] <- life_expectancy_n_exp[,2] - 1.96*life_expectancy_n_exp[,3];
  life_expectancy_n_exp[,5] <- life_expectancy_n_exp[,2] + 1.96*life_expectancy_n_exp[,3];

  for (i in 1:nrow(life_expectancy_n_exp)) {
    if (life_expectancy_n_exp[i,4]<0 & is.na(life_expectancy_n_exp[i,4])==F) {
      life_expectancy_n_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_n_exp)) {
    if (life_expectancy_n_exp[i,2]==0 & life_expectancy_n_exp[i,4]==0 & life_expectancy_n_exp[i,5]==0 &
        is.na(life_expectancy_n_exp[i,2])==F & is.na(life_expectancy_n_exp[i,4])==F & is.na(life_expectancy_n_exp[i,5])==F) {
      life_expectancy_n_exp[i,2] <- NA;
      life_expectancy_n_exp[i,3] <- NA;
      life_expectancy_n_exp[i,4] <- NA;
      life_expectancy_n_exp[i,5] <- NA
    }
  }

  colnames(life_expectancy_n_exp) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Life-expectancy without disease

  life_expectancy_w_dis <- matrix(c(0),
                                  nrow=100-65+1,
                                  ncol=5,
                                  byrow = T);
  life_expectancy_w_dis[,1] <- esp_vie_sans_mal[,1]
  life_expectancy_w_dis[,2] <- esp_vie_sans_mal[,2]
  life_expectancy_w_dis[,3] <- apply(esp_vie_sans_mal[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_w_dis[,4] <- life_expectancy_w_dis[,2] - 1.96*life_expectancy_w_dis[,3];
  life_expectancy_w_dis[,5] <- life_expectancy_w_dis[,2] + 1.96*life_expectancy_w_dis[,3];

  for (i in 1:nrow(life_expectancy_w_dis)) {
    if (life_expectancy_w_dis[i,4]<0  & is.na(life_expectancy_w_dis[i,4])==F) {
      life_expectancy_w_dis[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_w_dis)) {
    if (life_expectancy_w_dis[i,2]==0 & life_expectancy_w_dis[i,4]==0 & life_expectancy_w_dis[i,5]==0 &
        is.na(life_expectancy_w_dis[i,2])==F & is.na(life_expectancy_w_dis[i,4])==F & is.na(life_expectancy_w_dis[i,5])==F) {
      life_expectancy_w_dis[i,2] <- NA;
      life_expectancy_w_dis[i,3] <- NA;
      life_expectancy_w_dis[i,4] <- NA;
      life_expectancy_w_dis[i,5] <- NA
    }
  }

  colnames(life_expectancy_w_dis) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Life-expectancy without disease on exposed peoples

  life_expectancy_w_dis_exp <- matrix(c(0),
                                      nrow=100-65+1,
                                      ncol=5,
                                      byrow = T);
  life_expectancy_w_dis_exp[,1] <- esp_vie_sans_mal_conso[,1]
  life_expectancy_w_dis_exp[,2] <- esp_vie_sans_mal_conso[,2]
  life_expectancy_w_dis_exp[,3] <- apply(esp_vie_sans_mal_conso[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_w_dis_exp[,4] <- life_expectancy_w_dis_exp[,2] - 1.96*life_expectancy_w_dis_exp[,3];
  life_expectancy_w_dis_exp[,5] <- life_expectancy_w_dis_exp[,2] + 1.96*life_expectancy_w_dis_exp[,3];

  for (i in 1:nrow(life_expectancy_w_dis_exp)) {
    if (life_expectancy_w_dis_exp[i,4]<0  & is.na(life_expectancy_w_dis_exp[i,4])==F) {
      life_expectancy_w_dis_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_w_dis_exp)) {
    if (life_expectancy_w_dis_exp[i,2]==0 & life_expectancy_w_dis_exp[i,4]==0 & life_expectancy_w_dis_exp[i,5]==0 &
        is.na(life_expectancy_w_dis_exp[i,2])==F & is.na(life_expectancy_w_dis_exp[i,4])==F & is.na(life_expectancy_w_dis_exp[i,5])==F) {
      life_expectancy_w_dis_exp[i,2] <- NA;
      life_expectancy_w_dis_exp[i,3] <- NA;
      life_expectancy_w_dis_exp[i,4] <- NA;
      life_expectancy_w_dis_exp[i,5] <- NA
    }
  }

  colnames(life_expectancy_w_dis_exp) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Life-expectancy without disease on non exposed peoples

  life_expectancy_w_dis_n_exp <- matrix(c(0),
                                        nrow=100-65+1,
                                        ncol=5,
                                        byrow = T);
  life_expectancy_w_dis_n_exp[,1] <- esp_vie_sans_mal_nonconso[,1]
  life_expectancy_w_dis_n_exp[,2] <- esp_vie_sans_mal_nonconso[,2]
  life_expectancy_w_dis_n_exp[,3] <- apply(esp_vie_sans_mal_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_w_dis_n_exp[,4] <- life_expectancy_w_dis_n_exp[,2] - 1.96*life_expectancy_w_dis_n_exp[,3];
  life_expectancy_w_dis_n_exp[,5] <- life_expectancy_w_dis_n_exp[,2] + 1.96*life_expectancy_w_dis_n_exp[,3];

  for (i in 1:nrow(life_expectancy_w_dis_n_exp)) {
    if (life_expectancy_w_dis_n_exp[i,4]<0  & is.na(life_expectancy_w_dis_n_exp[i,4])==F) {
      life_expectancy_w_dis_n_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_w_dis_n_exp)) {
    if (life_expectancy_w_dis_n_exp[i,2]==0 & life_expectancy_w_dis_n_exp[i,4]==0 & life_expectancy_w_dis_n_exp[i,5]==0 &
        is.na(life_expectancy_w_dis_n_exp[i,2])==F & is.na(life_expectancy_w_dis_n_exp[i,4])==F & is.na(life_expectancy_w_dis_n_exp[i,5])==F) {
      life_expectancy_w_dis_n_exp[i,2] <- NA;
      life_expectancy_w_dis_n_exp[i,3] <- NA;
      life_expectancy_w_dis_n_exp[i,4] <- NA;
      life_expectancy_w_dis_n_exp[i,5] <- NA
    }
  }

  colnames(life_expectancy_w_dis_n_exp) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Life-expectancy for diseased subject

  life_expectancy_dis <- matrix(c(0),
                                nrow=100-65+1,
                                ncol=5,
                                byrow = T);
  life_expectancy_dis[,1] <- esp_vie_mal[,1]
  life_expectancy_dis[,2] <- esp_vie_mal[,2]
  life_expectancy_dis[,3] <- apply(esp_vie_mal[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_dis[,4] <- life_expectancy_dis[,2] - 1.96*life_expectancy_dis[,3];
  life_expectancy_dis[,5] <- life_expectancy_dis[,2] + 1.96*life_expectancy_dis[,3];

  for (i in 1:nrow(life_expectancy_dis)) {
    if (life_expectancy_dis[i,4]<0  & is.na(life_expectancy_dis[i,4])==F) {
      life_expectancy_dis[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_dis)) {
    if (life_expectancy_dis[i,2]==0 & life_expectancy_dis[i,4]==0 & life_expectancy_dis[i,5]==0 &
        is.na(life_expectancy_dis[i,2])==F & is.na(life_expectancy_dis[i,4])==F & is.na(life_expectancy_dis[i,5])==F) {
      life_expectancy_dis[i,2] <- NA;
      life_expectancy_dis[i,3] <- NA;
      life_expectancy_dis[i,4] <- NA;
      life_expectancy_dis[i,5] <- NA
    }
  }

  colnames(life_expectancy_dis) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Life-expectancy for diseased subject on exposed peoples

  life_expectancy_dis_exp <- matrix(c(0),
                                    nrow=100-65+1,
                                    ncol=5,
                                    byrow = T);
  life_expectancy_dis_exp[,1] <- esp_vie_mal_conso[,1]
  life_expectancy_dis_exp[,2] <- esp_vie_mal_conso[,2]
  life_expectancy_dis_exp[,3] <- apply(esp_vie_mal_conso[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_dis_exp[,4] <- life_expectancy_dis_exp[,2] - 1.96*life_expectancy_dis_exp[,3];
  life_expectancy_dis_exp[,5] <- life_expectancy_dis_exp[,2] + 1.96*life_expectancy_dis_exp[,3];

  for (i in 1:nrow(life_expectancy_dis_exp)) {
    if (life_expectancy_dis_exp[i,4]<0 & is.na(life_expectancy_dis_exp[i,4])==F) {
      life_expectancy_dis_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_dis_exp)) {
    if (life_expectancy_dis_exp[i,2]==0 & life_expectancy_dis_exp[i,4]==0 & life_expectancy_dis_exp[i,5]==0 &
        is.na(life_expectancy_dis_exp[i,2])==F & is.na(life_expectancy_dis_exp[i,4])==F & is.na(life_expectancy_dis_exp[i,5])==F) {
      life_expectancy_dis_exp[i,2] <- NA;
      life_expectancy_dis_exp[i,3] <- NA;
      life_expectancy_dis_exp[i,4] <- NA;
      life_expectancy_dis_exp[i,5] <- NA
    }
  }

  colnames(life_expectancy_dis_exp) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Life-expectancy for diseased subject on non exposed peoples

  life_expectancy_dis_n_exp <- matrix(c(0),
                                      nrow=100-65+1,
                                      ncol=5,
                                      byrow = T);
  life_expectancy_dis_n_exp[,1] <- esp_vie_mal_nonconso[,1]
  life_expectancy_dis_n_exp[,2] <- esp_vie_mal_nonconso[,2]
  life_expectancy_dis_n_exp[,3] <- apply(esp_vie_mal_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_dis_n_exp[,4] <- life_expectancy_dis_n_exp[,2] - 1.96*life_expectancy_dis_n_exp[,3];
  life_expectancy_dis_n_exp[,5] <- life_expectancy_dis_n_exp[,2] + 1.96*life_expectancy_dis_n_exp[,3];

  for (i in 1:nrow(life_expectancy_dis_n_exp)) {
    if (life_expectancy_dis_n_exp[i,4]<0  & is.na(life_expectancy_dis_n_exp[i,4])==F) {
      life_expectancy_dis_n_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_dis_n_exp)) {
    if (life_expectancy_dis_n_exp[i,2]==0 & life_expectancy_dis_n_exp[i,4]==0 & life_expectancy_dis_n_exp[i,5]==0 &
        is.na(life_expectancy_dis_n_exp[i,2])==F & is.na(life_expectancy_dis_n_exp[i,4])==F & is.na(life_expectancy_dis_n_exp[i,5])==F) {
      life_expectancy_dis_n_exp[i,2] <- NA;
      life_expectancy_dis_n_exp[i,3] <- NA;
      life_expectancy_dis_n_exp[i,4] <- NA;
      life_expectancy_dis_n_exp[i,5] <- NA
    }
  }

  colnames(life_expectancy_dis_n_exp) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Life-expectancy for non diseased subject

  life_expectancy_n_dis <- matrix(c(0),
                                  nrow=100-65+1,
                                  ncol=5,
                                  byrow = T);
  life_expectancy_n_dis[,1] <- esp_vie_non_mal[,1]
  life_expectancy_n_dis[,2] <- esp_vie_non_mal[,2]
  life_expectancy_n_dis[,3] <- apply(esp_vie_non_mal[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_n_dis[,4] <- life_expectancy_n_dis[,2] - 1.96*life_expectancy_n_dis[,3];
  life_expectancy_n_dis[,5] <- life_expectancy_n_dis[,2] + 1.96*life_expectancy_n_dis[,3];

  for (i in 1:nrow(life_expectancy_n_dis)) {
    if (life_expectancy_n_dis[i,4]<0  & is.na(life_expectancy_n_dis[i,4])==F) {
      life_expectancy_n_dis[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_n_dis)) {
    if (life_expectancy_n_dis[i,2]==0 & life_expectancy_n_dis[i,4]==0 & life_expectancy_n_dis[i,5]==0 &
        is.na(life_expectancy_n_dis[i,2])==F & is.na(life_expectancy_n_dis[i,4])==F & is.na(life_expectancy_n_dis[i,5])==F) {
      life_expectancy_n_dis[i,2] <- NA;
      life_expectancy_n_dis[i,3] <- NA;
      life_expectancy_n_dis[i,4] <- NA;
      life_expectancy_n_dis[i,5] <- NA
    }
  }

  colnames(life_expectancy_n_dis) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Life-expectancy for non diseased subject on exposed peoples

  life_expectancy_n_dis_exp <- matrix(c(0),
                                      nrow=100-65+1,
                                      ncol=5,
                                      byrow = T);
  life_expectancy_n_dis_exp[,1] <- esp_vie_non_mal_conso[,1]
  life_expectancy_n_dis_exp[,2] <- esp_vie_non_mal_conso[,2]
  life_expectancy_n_dis_exp[,3] <- apply(esp_vie_non_mal_conso[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_n_dis_exp[,4] <- life_expectancy_n_dis_exp[,2] - 1.96*life_expectancy_n_dis_exp[,3];
  life_expectancy_n_dis_exp[,5] <- life_expectancy_n_dis_exp[,2] + 1.96*life_expectancy_n_dis_exp[,3];

  for (i in 1:nrow(life_expectancy_n_dis_exp)) {
    if (life_expectancy_n_dis_exp[i,4]<0  & is.na(life_expectancy_n_dis_exp[i,4])==F) {
      life_expectancy_n_dis_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_n_dis_exp)) {
    if (life_expectancy_n_dis_exp[i,2]==0 & life_expectancy_n_dis_exp[i,4]==0 & life_expectancy_n_dis_exp[i,5]==0 &
        is.na(life_expectancy_n_dis_exp[i,2])==F & is.na(life_expectancy_n_dis_exp[i,4])==F & is.na(life_expectancy_n_dis_exp[i,5])==F) {
      life_expectancy_n_dis_exp[i,2] <- NA;
      life_expectancy_n_dis_exp[i,3] <- NA;
      life_expectancy_n_dis_exp[i,4] <- NA;
      life_expectancy_n_dis_exp[i,5] <- NA
    }
  }

  colnames(life_expectancy_n_dis_exp) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Life-expectancy for non diseased subject on non exposed peoples

  life_expectancy_n_dis_n_exp <- matrix(c(0),
                                        nrow=100-65+1,
                                        ncol=5,
                                        byrow = T);
  life_expectancy_n_dis_n_exp[,1] <- esp_vie_non_mal_nonconso[,1]
  life_expectancy_n_dis_n_exp[,2] <- esp_vie_non_mal_nonconso[,2]
  life_expectancy_n_dis_n_exp[,3] <- apply(esp_vie_non_mal_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  life_expectancy_n_dis_n_exp[,4] <- life_expectancy_n_dis_n_exp[,2] - 1.96*life_expectancy_n_dis_n_exp[,3];
  life_expectancy_n_dis_n_exp[,5] <- life_expectancy_n_dis_n_exp[,2] + 1.96*life_expectancy_n_dis_n_exp[,3];

  for (i in 1:nrow(life_expectancy_n_dis_n_exp)) {
    if (life_expectancy_n_dis_n_exp[i,4]<0  & is.na(life_expectancy_n_dis_n_exp[i,4])==F) {
      life_expectancy_n_dis_n_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(life_expectancy_n_dis_n_exp)) {
    if (life_expectancy_n_dis_n_exp[i,2]==0 & life_expectancy_n_dis_n_exp[i,4]==0 & life_expectancy_n_dis_n_exp[i,5]==0 &
        is.na(life_expectancy_n_dis_n_exp[i,2])==F & is.na(life_expectancy_n_dis_n_exp[i,4])==F & is.na(life_expectancy_n_dis_n_exp[i,5])==F) {
      life_expectancy_n_dis_n_exp[i,2] <- NA;
      life_expectancy_n_dis_n_exp[i,3] <- NA;
      life_expectancy_n_dis_n_exp[i,4] <- NA;
      life_expectancy_n_dis_n_exp[i,5] <- NA
    }
  }

  colnames(life_expectancy_n_dis_n_exp) <- c("age","life-expectancy","std","CI95_low","CI95_upp");

  ### Prevalence of disease

  prevalence <- prevalence[-1,]

  if (nb_iter != 0) {

    prev <- colSums(prevalence[,-1])

  } else {

    prev <- sum(prevalence[,-1])

  }

  number_prevalence <- matrix(c(0),
                              nrow=1,
                              ncol=5,
                              byrow = T);
  number_prevalence[1] <- as.integer(t);
  number_prevalence[2] <- prev[1];
  number_prevalence[3] <- sd(prev[-1], na.rm=T);
  number_prevalence[4] <- number_prevalence[2] - 1.96*number_prevalence[3];
  number_prevalence[5] <- number_prevalence[2] + 1.96*number_prevalence[3];

  if (number_prevalence[4]<0 & is.na(number_prevalence[4])==F) {
    number_prevalence[4] <- 0
  };

  if (number_prevalence[2]==0 & number_prevalence[4]==0 & number_prevalence[5]==0 &
      is.na(number_prevalence[2])==F & is.na(number_prevalence[4])==F & is.na(number_prevalence[5])==F) {
    number_prevalence[2] <- NA;
    number_prevalence[3] <- NA;
    number_prevalence[4] <- NA;
    number_prevalence[5] <- NA
  };

  colnames(number_prevalence) <- c("year","prevalence","std","CI95_low","CI95_upp");
  number_prevalence <- number_prevalence[,-1];

  ### Prevalence of disease by age

  number_prev_age <- matrix(c(0),
                            nrow=99-66+1,
                            ncol=5,
                            byrow = T);
  number_prev_age[,1] <- prevalence[,1]
  number_prev_age[,2] <- prevalence[,2]
  number_prev_age[,3] <- apply(prevalence[,-c(1:2)], 1, sd, na.rm=T);
  number_prev_age[,4] <- number_prev_age[,2] - 1.96*number_prev_age[,3];
  number_prev_age[,5] <- number_prev_age[,2] + 1.96*number_prev_age[,3];

  for (i in 1:nrow(number_prev_age)) {
    if (number_prev_age[i,4]<0 & is.na(number_prev_age[i,4])==F) {
      number_prev_age[i,4] <- 0
    }
  }

  for (i in 1:nrow(number_prev_age)) {
    if (number_prev_age[i,2]==0 & number_prev_age[i,4]==0 & number_prev_age[i,5]==0 &
        is.na(number_prev_age[i,2])==F & is.na(number_prev_age[i,4])==F & is.na(number_prev_age[i,5])==F) {
      number_prev_age[i,2] <- NA;
      number_prev_age[i,3] <- NA;
      number_prev_age[i,4] <- NA;
      number_prev_age[i,5] <- NA
    }
  }

  colnames(number_prev_age) <- c("age","prevalence","std","CI95_low","CI95_upp");

  ### Prevalence rate of disease by age

  taux_prevalence <- taux_prevalence[-1,]

  prev_rate_disease <- matrix(c(0),
                              nrow=99-66+1,
                              ncol=5,
                              byrow = T);
  prev_rate_disease[,1] <- taux_prevalence[,1]
  prev_rate_disease[,2] <- taux_prevalence[,2]
  prev_rate_disease[,3] <- apply(taux_prevalence[,-c(1:2)], 1, sd, na.rm=T);
  prev_rate_disease[,4] <- prev_rate_disease[,2] - 1.96*prev_rate_disease[,3];
  prev_rate_disease[,5] <- prev_rate_disease[,2] + 1.96*prev_rate_disease[,3];

  for (i in 1:nrow(prev_rate_disease)) {
    if (prev_rate_disease[i,4]<0 & is.na(prev_rate_disease[i,4])==F) {
      prev_rate_disease[i,4] <- 0
    }
  }

  for (i in 1:nrow(prev_rate_disease)) {
    if (prev_rate_disease[i,5]>1 & is.na(prev_rate_disease[i,5])==F) {
      prev_rate_disease[i,5] <- 1
    }
  }

  for (i in 1:nrow(prev_rate_disease)) {
    if (prev_rate_disease[i,2]==0 & prev_rate_disease[i,4]==0 & prev_rate_disease[i,5]==0 &
        is.na(prev_rate_disease[i,2])==F & is.na(prev_rate_disease[i,4])==F & is.na(prev_rate_disease[i,5])==F) {
      prev_rate_disease[i,2] <- NA;
      prev_rate_disease[i,3] <- NA;
      prev_rate_disease[i,4] <- NA;
      prev_rate_disease[i,5] <- NA
    }
  }

  colnames(prev_rate_disease) <- c("age","prevalence_rate","std","CI95_low","CI95_upp");

  ### Survival

  if (nb_iter != 0) {

    surv <- colSums(survie[,-1])

  } else {

    surv <- sum(survie[,-1])

  }

  number_survival <- matrix(c(0),
                            nrow=1,
                            ncol=5,
                            byrow = T);
  number_survival[1] <- as.integer(t);
  number_survival[2] <- surv[1]
  number_survival[3] <- sd(surv[-1], na.rm=T);
  number_survival[4] <- number_survival[2] - 1.96*number_survival[3];
  number_survival[5] <- number_survival[2] + 1.96*number_survival[3];

  if (number_survival[4]<0 & is.na(number_survival[4])==F) {
    number_survival[4] <- 0
  };

  if (number_survival[2]==0 & number_survival[4]==0 & number_survival[5]==0 &
      is.na(number_survival[2])==F & is.na(number_survival[4])==F & is.na(number_survival[5])==F) {
    number_survival[2] <- NA;
    number_survival[3] <- NA;
    number_survival[4] <- NA;
    number_survival[5] <- NA
  };

  colnames(number_survival) <- c("year","survival","std","CI95_low","CI95_upp");
  number_survival <- number_survival[,-1];

  ### Survival by age

  number_survival_age <- matrix(c(0),
                                nrow=99-65+1,
                                ncol=5,
                                byrow = T);
  number_survival_age[,1] <- survie[,1]
  number_survival_age[,2] <- survie[,2]
  number_survival_age[,3] <- apply(survie[,-c(1:2)], 1, sd, na.rm=T);
  number_survival_age[,4] <- number_survival_age[,2] - 1.96*number_survival_age[,3];
  number_survival_age[,5] <- number_survival_age[,2] + 1.96*number_survival_age[,3];

  for (i in 1:nrow(number_survival_age)) {
    if (number_survival_age[i,4]<0 & is.na(number_survival_age[i,4])==F) {
      number_survival_age[i,4] <- 0
    }
  }

  for (i in 1:nrow(number_survival_age)) {
    if (number_survival_age[i,2]==0 & number_survival_age[i,4]==0 & number_survival_age[i,5]==0 &
        is.na(number_survival_age[i,2])==F & is.na(number_survival_age[i,4])==F & is.na(number_survival_age[i,5])==F) {
      number_survival_age[i,2] <- NA;
      number_survival_age[i,3] <- NA;
      number_survival_age[i,4] <- NA;
      number_survival_age[i,5] <- NA
    }
  }

  colnames(number_survival_age) <- c("age","survival","std","CI95_low","CI95_upp");

  ### Survival rate

  taux_survivants <- taux_survivants[-1,]

  survival_rate <- matrix(c(0),
                          nrow=99-66+1,
                          ncol=5,
                          byrow = T);
  survival_rate[,1] <- taux_survivants[,1]
  survival_rate[,2] <- taux_survivants[,2]
  survival_rate[,3] <- apply(taux_survivants[,-c(1:2)], 1, sd, na.rm=T);
  survival_rate[,4] <- survival_rate[,2] - 1.96*survival_rate[,3];
  survival_rate[,5] <- survival_rate[,2] + 1.96*survival_rate[,3];

  for (i in 1:nrow(survival_rate)) {
    if (survival_rate[i,4]<0 & is.na(survival_rate[i,4])==F) {
      survival_rate[i,4] <- 0
    }
  }

  for (i in 1:nrow(survival_rate)) {
    if (survival_rate[i,5]>1 & is.na(survival_rate[i,5])==F) {
      survival_rate[i,5] <- 1
    }
  }

  for (i in 1:nrow(survival_rate)) {
    if (survival_rate[i,2]==0 & survival_rate[i,4]==0 & survival_rate[i,5]==0 &
        is.na(survival_rate[i,2])==F & is.na(survival_rate[i,4])==F & is.na(survival_rate[i,5])==F) {
      survival_rate[i,2] <- NA;
      survival_rate[i,3] <- NA;
      survival_rate[i,4] <- NA;
      survival_rate[i,5] <- NA
    }
  }

  colnames(survival_rate) <- c("age","survival_rate","std","CI95_low","CI95_upp");

  ### Global prevalence rate of disease

  taux_prev_demence <- matrix(c(0),
                              nrow=1,
                              ncol=length(prev)+1,
                              byrow = T);

  taux_prev_demence[1] <- t;

  for (i in 2:ncol(taux_prev_demence)) {

    taux_prev_demence[1,i] <- prev[i-1] / surv[i-1];

  };

  prev_rate_disease <- matrix(c(0),
                              nrow=1,
                              ncol=5,
                              byrow = T);
  prev_rate_disease[1] <- as.integer(t);
  prev_rate_disease[2] <- taux_prev_demence[,2];
  prev_rate_disease[3] <- sd(taux_prev_demence[,-c(1:2)], na.rm=T);
  prev_rate_disease[4] <- prev_rate_disease[2] - 1.96*prev_rate_disease[3];
  prev_rate_disease[5] <- prev_rate_disease[2] + 1.96*prev_rate_disease[3];

  if (prev_rate_disease[4]<0 & is.na(prev_rate_disease[4])==F) {
    prev_rate_disease[4] <- 0
  };

  if (prev_rate_disease[2]==0 & prev_rate_disease[4]==0 & prev_rate_disease[5]==0 &
      is.na(prev_rate_disease[2])==F & is.na(prev_rate_disease[4])==F & is.na(prev_rate_disease[5])==F) {
    prev_rate_disease[2] <- NA;
    prev_rate_disease[3] <- NA;
    prev_rate_disease[4] <- NA;
    prev_rate_disease[5] <- NA
  };

  colnames(prev_rate_disease) <- c("year","prevalence_rate","std","CI95_low","CI95_upp");
  prev_rate_disease <- prev_rate_disease[,-1];

  ### Mean number of years spent with disease

  number_years_disease <- matrix(c(0),
                                 nrow=100-65+1,
                                 ncol=5,
                                 byrow = T);
  number_years_disease[,1] <- nb_moy_dem[,1]
  number_years_disease[,2] <- nb_moy_dem[,2]
  number_years_disease[,3] <- apply(nb_moy_dem[,-c(1:2)], 1, sd, na.rm=T);
  number_years_disease[,4] <- number_years_disease[,2] - 1.96*number_years_disease[,3];
  number_years_disease[,5] <- number_years_disease[,2] + 1.96*number_years_disease[,3];

  for (i in 1:nrow(number_years_disease)) {
    if (number_years_disease[i,4]<0 & is.na(number_years_disease[i,4])==F) {
      number_years_disease[i,4] <- 0
    }
  }

  for (i in 1:nrow(number_years_disease)) {
    if (number_years_disease[i,2]==0 & number_years_disease[i,4]==0 & number_years_disease[i,5]==0 &
        is.na(number_years_disease[i,2])==F  & is.na(number_years_disease[i,4])==F  & is.na(number_years_disease[i,5])==F) {
      number_years_disease[i,2] <- NA;
      number_years_disease[i,3] <- NA;
      number_years_disease[i,4] <- NA;
      number_years_disease[i,5] <- NA
    }
  }

  colnames(number_years_disease) <- c("age","number_years","std","CI95_low","CI95_upp");

  ### Mean number of years spent with disease on exposed peoples

  number_years_disease_exp <- matrix(c(0),
                                     nrow=100-65+1,
                                     ncol=5,
                                     byrow = T);
  number_years_disease_exp[,1] <- nb_moy_dem_conso[,1]
  number_years_disease_exp[,2] <- nb_moy_dem_conso[,2]
  number_years_disease_exp[,3] <- apply(nb_moy_dem_conso[,-c(1:2)], 1, sd, na.rm=T);
  number_years_disease_exp[,4] <- number_years_disease_exp[,2] - 1.96*number_years_disease_exp[,3];
  number_years_disease_exp[,5] <- number_years_disease_exp[,2] + 1.96*number_years_disease_exp[,3];

  for (i in 1:nrow(number_years_disease_exp)) {
    if (number_years_disease_exp[i,4]<0 & is.na(number_years_disease_exp[i,4])==F) {
      number_years_disease_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(number_years_disease_exp)) {
    if (number_years_disease_exp[i,2]==0 & number_years_disease_exp[i,4]==0 & number_years_disease_exp[i,5]==0 &
        is.na(number_years_disease_exp[i,2])==F & is.na(number_years_disease_exp[i,4])==F & is.na(number_years_disease_exp[i,5])==F) {
      number_years_disease_exp[i,2] <- NA;
      number_years_disease_exp[i,3] <- NA;
      number_years_disease_exp[i,4] <- NA;
      number_years_disease_exp[i,5] <- NA
    }
  }

  colnames(number_years_disease_exp) <- c("age","number_years","std","CI95_low","CI95_upp");

  ### Mean number of years spent with disease on non exposed peoples

  number_years_disease_n_exp <- matrix(c(0),
                                       nrow=100-65+1,
                                       ncol=5,
                                       byrow = T);
  number_years_disease_n_exp[,1] <- nb_moy_dem_nonconso[,1]
  number_years_disease_n_exp[,2] <- nb_moy_dem_nonconso[,2]
  number_years_disease_n_exp[,3] <- apply(nb_moy_dem_nonconso[,-c(1:2)], 1, sd, na.rm=T);
  number_years_disease_n_exp[,4] <- number_years_disease_n_exp[,2] - 1.96*number_years_disease_n_exp[,3];
  number_years_disease_n_exp[,5] <- number_years_disease_n_exp[,2] + 1.96*number_years_disease_n_exp[,3];

  for (i in 1:nrow(number_years_disease_n_exp)) {
    if (number_years_disease_n_exp[i,4]<0 & is.na(number_years_disease_n_exp[i,4])==F) {
      number_years_disease_n_exp[i,4] <- 0
    }
  }

  for (i in 1:nrow(number_years_disease_n_exp)) {
    if (number_years_disease_n_exp[i,2]==0 & number_years_disease_n_exp[i,4]==0 & number_years_disease_n_exp[i,5]==0 &
        is.na(number_years_disease_n_exp[i,2])==F & is.na(number_years_disease_n_exp[i,4])==F & is.na(number_years_disease_n_exp[i,5])==F) {
      number_years_disease_n_exp[i,2] <- NA;
      number_years_disease_n_exp[i,3] <- NA;
      number_years_disease_n_exp[i,4] <- NA;
      number_years_disease_n_exp[i,5] <- NA
    }
  }

  colnames(number_years_disease_n_exp) <- c("age","number_years","std","CI95_low","CI95_upp");

  ### Life-long probability of disease

  prb_dem <- prb_dem[-1,]

  prb_demence <- matrix(c(0),
                        nrow=1,
                        ncol=ncol(prb_dem),
                        byrow = T);

  prb_demence[1] <- t;

  for (i in 2:ncol(prb_dem)) {

    prb_demence[i] <- sum(prb_dem[,i]) / nrow(etat);

  };

  ll_prob_disease <- matrix(c(0),
                            nrow=1,
                            ncol=5,
                            byrow = T);
  ll_prob_disease[1] <- as.integer(t);
  ll_prob_disease[2] <- prb_demence[,2];
  ll_prob_disease[3] <- sd(prb_demence[,-c(1:2)], na.rm=T);
  ll_prob_disease[4] <- ll_prob_disease[2] - 1.96*ll_prob_disease[3];
  ll_prob_disease[5] <- ll_prob_disease[2] + 1.96*ll_prob_disease[3];

  if (ll_prob_disease[4]<0 & is.na(ll_prob_disease[4])==F) {
    ll_prob_disease[4] <- 0
  };

  if (ll_prob_disease[2]==0 & ll_prob_disease[4]==0 & ll_prob_disease[5]==0 &
      is.na(ll_prob_disease[2])==F & is.na(ll_prob_disease[4])==F & is.na(ll_prob_disease[5])==F) {
    ll_prob_disease[2] <- NA;
    ll_prob_disease[3] <- NA;
    ll_prob_disease[4] <- NA;
    ll_prob_disease[5] <- NA
  };

  colnames(ll_prob_disease) <- c("year","probability","std","CI95_low","CI95_upp");
  ll_prob_disease <- ll_prob_disease[,-1];

  ### Average age at disease onset

  age_dem <- age_dem[-1,]

  age_demence <- matrix(c(0),
                        nrow=1,
                        ncol=ncol(age_dem),
                        byrow = T);

  age_demence[1] <- t;

  n0 <- vector(length = nrow(age_dem));

  for (i in 2:ncol(age_dem)) {

    for (j in 1:nrow(age_dem)) {

      n0[j] <- age_dem[j,1]*age_dem[j,i]

    };

    age_demence[i] <- sum(n0)/sum(age_dem[,i]);

  }

  average_age_disease <- matrix(c(0),
                                nrow=1,
                                ncol=5,
                                byrow = T);
  average_age_disease[1] <- as.integer(t);
  average_age_disease[2] <- age_demence[,2]
  average_age_disease[3] <- sd(age_demence[,-c(1:2)], na.rm=T);
  average_age_disease[4] <- average_age_disease[2] - 1.96*average_age_disease[3];
  average_age_disease[5] <- average_age_disease[2] + 1.96*average_age_disease[3];

  if (average_age_disease[4]<0 & is.na(average_age_disease[4])==F) {
    average_age_disease[4] <- 0
  };

  if (average_age_disease[2]==0 & average_age_disease[4]==0 & average_age_disease[5]==0 &
      is.na(average_age_disease[2])==F & is.na(average_age_disease[4])==F & is.na(average_age_disease[5])==F) {
    average_age_disease[2] <- NA;
    average_age_disease[3] <- NA;
    average_age_disease[4] <- NA;
    average_age_disease[5] <- NA
  };

  colnames(average_age_disease) <- c("year","mean_age","std","CI95_low","CI95_upp");
  average_age_disease <- average_age_disease[,-1];

  ### Mean number of years of exposition

  moy_conso <- moy_conso[-1,]

  moyenne_conso <- matrix(c(0),
                          nrow=1,
                          ncol=ncol(moy_conso),
                          byrow = T);

  moyenne_conso[1] <- t;

  for (i in 2:ncol(moy_conso)) {

    moyenne_conso[i] <- sum(moy_conso[,i]) / nrow(etat[which(etat[,1]%in%("00") | etat[,1]%in%("10")),]); # nombre moyen d'années de consommation de benzo

  }

  number_years_exp <- matrix(c(0),
                             nrow=1,
                             ncol=5,
                             byrow = T);
  number_years_exp[1] <- as.integer(t);
  number_years_exp[2] <- moyenne_conso[,2]
  number_years_exp[3] <- sd(moyenne_conso[,-c(1:2)], na.rm=T);
  number_years_exp[4] <- number_years_exp[2] - 1.96*number_years_exp[3];
  number_years_exp[5] <- number_years_exp[2] + 1.96*number_years_exp[3];

  if (number_years_exp[4]<0 & is.na(number_years_exp[4])==F) {
    number_years_exp[4] <- 0
  };

  if (number_years_exp[2]==0 & number_years_exp[4]==0 & number_years_exp[5]==0 &
      is.na(number_years_exp[2])==F & is.na(number_years_exp[4])==F & is.na(number_years_exp[5])==F) {
    number_years_exp[2] <- NA;
    number_years_exp[3] <- NA;
    number_years_exp[4] <- NA;
    number_years_exp[5] <- NA
  };

  colnames(number_years_exp) <- c("year","mean_number","std","CI95_low","CI95_upp");
  number_years_exp <- number_years_exp[,-1];

  ### Number of exposed peoples at least one time

  exposition <- matrix(c(0),
                       nrow=105-65+1,
                       ncol=5,
                       byrow = T);
  exposition[,1] <- prevalence_conso[,1]
  exposition[,2] <- prevalence_conso[,2]
  exposition[,3] <- apply(prevalence_conso[,-c(1:2)], 1, sd, na.rm=T);
  exposition[,4] <- exposition[,2] - 1.96*exposition[,3];
  exposition[,5] <- exposition[,2] + 1.96*exposition[,3];

  for (i in 1:nrow(exposition)) {
    if (exposition[i,4]<0 & is.na(exposition[i,4])==F) {
      exposition[i,4] <- 0
    }
  }

  for (i in 1:nrow(exposition)) {
    if (exposition[i,5]>1 & is.na(exposition[i,5])==F) {
      exposition[i,5] <- 1
    }
  }

  for (i in 1:nrow(exposition)) {
    if (exposition[i,2]==0 & exposition[i,4]==0 & exposition[i,5]==0 &
        is.na(exposition[i,2])==F & is.na(exposition[i,4])==F & is.na(exposition[i,5])==F) {
      exposition[i,2] <- NA;
      exposition[i,3] <- NA;
      exposition[i,4] <- NA;
      exposition[i,5] <- NA
    }
  }

  colnames(exposition) <- c("age","exposition","std","CI95_low","CI95_upp");

  ### Mortality rate

  quotient_mortalite <- quotient_mortalite[-1,]

  mortality_rate <- matrix(c(0),
                           nrow=99-66+1,
                           ncol=5,
                           byrow = T);
  mortality_rate[,1] <- quotient_mortalite[,1]
  mortality_rate[,2] <- quotient_mortalite[,2]
  mortality_rate[,3] <- apply(quotient_mortalite[,-c(1:2)], 1, sd, na.rm=T);
  mortality_rate[,4] <- mortality_rate[,2] - 1.96*mortality_rate[,3];
  mortality_rate[,5] <- mortality_rate[,2] + 1.96*mortality_rate[,3];

  for (i in 1:nrow(mortality_rate)) {
    if (mortality_rate[i,4]<0 & is.na(mortality_rate[i,4])==F) {
      mortality_rate[i,4] <- 0
    }
  }

  for (i in 1:nrow(mortality_rate)) {
    if (mortality_rate[i,5]>1 & is.na(mortality_rate[i,5])==F) {
      mortality_rate[i,5] <- 1
    }
  }

  for (i in 1:nrow(mortality_rate)) {
    if (mortality_rate[i,2]==0 & mortality_rate[i,4]==0 & mortality_rate[i,5]==0 &
        is.na(mortality_rate[i,2])==F & is.na(mortality_rate[i,4])==F & is.na(mortality_rate[i,5])==F) {
      mortality_rate[i,2] <- NA;
      mortality_rate[i,3] <- NA;
      mortality_rate[i,4] <- NA;
      mortality_rate[i,5] <- NA
    }
  }

  colnames(mortality_rate) <- c("age","mortality_rate","std","CI95_low","CI95_upp");

  ### Output of the algorithm

  ### Overall life-expectancy

  list_overall_LE <- list(life_expectancy, life_expectancy_exp, life_expectancy_n_exp)
  names(list_overall_LE) <- c("life_expectancy", "life_expectancy_exp", "life_expectancy_n_exp")

  ###	Life-expectancy without disease

  list_LE_without_disease <- list(life_expectancy_w_dis, life_expectancy_w_dis_exp, life_expectancy_w_dis_n_exp)
  names(list_LE_without_disease) <- c("life_expectancy_w_dis", "life_expectancy_w_dis_exp", "life_expectancy_w_dis_n_exp")

  ### Life-expectancy for diseased subject

  list_LE_diseased <- list(life_expectancy_dis, life_expectancy_dis_exp, life_expectancy_dis_n_exp)
  names(list_LE_diseased) <- c("life_expectancy_dis", "life_expectancy_dis_exp", "life_expectancy_dis_n_exp")

  ###	Life-expectancy for non-diseased subject

  list_LE_non_diseased <- list(life_expectancy_n_dis, life_expectancy_n_dis_exp, life_expectancy_n_dis_n_exp)
  names(list_LE_non_diseased) <- c("life_expectancy_n_dis", "life_expectancy_n_dis_exp", "life_expectancy_n_dis_n_exp")

  ### Prevalence of disease

  list_prevalence_disease <- list(number_prev_age, number_prevalence, prev_rate_disease, prev_rate_disease)
  names(list_prevalence_disease) <- c("number_prev_age", "number_prevalence", "prev_rate_disease", "prev_rate_disease")

  ### Survival

  list_survival <- list(number_survival_age, number_survival, survival_rate)
  names(list_survival) <- c("number_survival_age", "number_survival", "survival_rate")

  ### Mean number of years spent with disease

  list_number_years_disease <- list(number_years_disease, number_years_disease_exp, number_years_disease_n_exp)
  names(list_number_years_disease) <- c("number_years_disease", "number_years_disease_exp", "number_years_disease_n_exp")

  ### Output list

  # With exposition

  #HI_output <- list(list_overall_LE, list_LE_without_disease, list_LE_diseased, list_LE_non_diseased, list_prevalence_disease, list_survival, list_number_years_disease, ll_prob_disease, average_age_disease, number_years_exp, exposition, mortality_rate);
  #names(HI_output) <- c(list_overall_LE", "list_LE_without_disease", "list_LE_diseased", "list_LE_non_diseased", "list_prevalence_disease", "list_survival", "list_number_years_disease", "ll_prob_disease", "average_age_disease", "number_years_exp", "exposition", "mortality_rate");

  # Without exposition

  HI_output <- list(list_overall_LE, list_LE_without_disease, list_LE_diseased, list_LE_non_diseased, list_prevalence_disease, list_survival, list_number_years_disease, ll_prob_disease, average_age_disease, number_years_exp, mortality_rate);
  names(HI_output) <- c("list_overall_LE", "list_LE_without_disease", "list_LE_diseased", "list_LE_non_diseased", "list_prevalence_disease", "list_survival", "list_number_years_disease", "ll_prob_disease", "average_age_disease", "number_years_exp", "mortality_rate");

  return(HI_output)

}
