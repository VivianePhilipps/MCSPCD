#' Computation of Health Indicators
#'
#' This function computes many health indicators under several scenarios of intervention in risk factor distribution for a given year.
#'
#' @param t year of the projections for health indicators.
#' @param intervention 0 = no change; 1 = risk factor prevalence and incidence considered as null. Default is \code{0}.
#' @param year_intervention year of the intervention in risk factor distribution takes place. Default is \code{NULL}.
#' @param nb_people number of people whose trajectory will be simulated for each generation. Default is \code{100}.
#' @param nb_iter number of iterations for the algorithm. Default is \code{0}.
#' @param data_pop data source for demographics data.
#' @param gender gender for computation. \code{"W"} for women and \code{"M"} for men. Default is \code{"W"}.
#' @param data_a01_values data source for the incidence of disease.
#' @param data_a02_values data source for the mortality of healthy subjects.
#' @param data_theta01_1_values data source for the relative risks associated with the exposure 1 for disease.
#' @param data_theta01_2_values data source for the relative risks associated with the exposure 2 for disease.
#' @param data_theta01_3_values data source for the relative risks associated with the exposure 3 for disease.
#' @param data_theta02_1_values data source for the relative risks associated with the exposure 1 for mortality among healthy subjects.
#' @param data_theta02_2_values data source for the relative risks associated with the exposure 2 for mortality among healthy subjects.
#' @param data_theta02_3_values data source for the relative risks associated with the exposure 3 for mortality among healthy subjects.
#' @param data_theta12_1_values data source for the relative risks associated with the exposure 1 for mortality among diseased subjects.
#' @param data_theta12_2_values data source for the relative risks associated with the exposure 2 for mortality among diseased subjects.
#' @param data_theta12_3_values data source for the relative risks associated with the exposure 3 for mortality among diseased subjects.
#' @param data_prev_0_values data source for the prevalence of the exposition 0.
#' @param data_prev_1_values data source for the prevalence of the exposition 1.
#' @param data_prev_2_values data source for the prevalence of the exposition 2.
#' @param data_prev_3_values data source for the prevalence of the exposition 3.
#' @param data_prev_4_values data source for the prevalence of the exposition 4.
#' @param data_prev_5_values data source for the prevalence of the exposition 5.
#' @param data_prev_6_values data source for the prevalence of the exposition 6.
#' @param data_prev_7_values data source for the prevalence of the exposition 7.
#' @param data_incid_0_values data source for the incidence of the exposition 0.
#' @param data_incid_1_values data source for the incidence of the exposition 1.
#' @param data_incid_3_values data source for the incidence of the exposition 3.
#' @param data_incid_5_values data source for the incidence of the exposition 5.
#' @param incid_global data source for the global diabete incidence.
#' @param theta1 HR for HTA in the diabete risk
#' @param theta3 HR for inact in the diabete risk 
#' @param data_rr_DvsND_values data source for the relative risks associated with the disease for mortality.
#'
#' @return a list containing the health indicators
#'
#' @export
#'
estimHI_3FDR <- function(t,
                         intervention = 0,
                         year_intervention = NULL,
                         nb_people = 100,
                         data_pop,
                         gender = "W",
                         data_a01_values,
                         data_a02_values,
                         data_theta01_1_values,
                         data_theta01_2_values,
                         data_theta01_3_values,
                         data_theta02_1_values,
                         data_theta02_2_values,
                         data_theta02_3_values,
                         data_theta12_1_values,
                         data_theta12_2_values,
                         data_theta12_3_values,
                         data_prev_0_values,
                         data_prev_1_values,
                         data_prev_2_values,
                         data_prev_3_values,
                         data_prev_4_values,
                         data_prev_5_values,
                         data_prev_6_values,
                         data_prev_7_values,
                         data_incid_0_values,
                         data_incid_1_values,
                         data_incid_3_values,
                         data_incid_5_values,
                         incid_global,theta1,theta3,
                         data_rr_DvsND_values)
{

### Year of projection

    set.seed(0)

    year_proj <- t


### Incidence of disease on non exposed peoples

    a010_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a010_values) <- c("age",1950:2080)

    a010_values[,1] <- c(66:105)

    for (a in 2:ncol(a010_values)){

        a010_values[,a] <- as.numeric(data_a01_values[which(data_a01_values[,1] != 65 & data_a01_values[,2]%in%(gender)),a+1]) / (data_prev_0_values[which(data_prev_0_values[,1] != 65 & data_prev_0_values[,3]%in%(gender)),2] + data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_prev_1_values[which(data_prev_1_values[,1] != 65 & data_prev_1_values[,3]%in%(gender)),2] + data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_prev_2_values[which(data_prev_2_values[,1] != 65 & data_prev_2_values[,3]%in%(gender)),2] + data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2]*data_prev_3_values[which(data_prev_3_values[,1] != 65 & data_prev_3_values[,3]%in%(gender)),2] + data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_prev_4_values[which(data_prev_4_values[,1] != 65 & data_prev_4_values[,3]%in%(gender)),2] + data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2]*data_prev_5_values[which(data_prev_5_values[,1] != 65 & data_prev_5_values[,3]%in%(gender)),2] + data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2]*data_prev_6_values[which(data_prev_6_values[,1] != 65 & data_prev_6_values[,3]%in%(gender)),2] + data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2]*data_prev_7_values[which(data_prev_7_values[,1] != 65 & data_prev_7_values[,3]%in%(gender)),2]);

    }

### Incidence of disease on exposed peoples

    a011_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a011_values) <- c("age",1950:2080)

    a011_values[,1] <- c(66:105)
    a011_values[,-1] <- a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2];

    a012_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a012_values) <- c("age",1950:2080)

    a012_values[,1] <- c(66:105)
    a012_values[,-1] <- a010_values[,-1]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2];

    a013_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a013_values) <- c("age",1950:2080)

    a013_values[,1] <- c(66:105)
    a013_values[,-1] <- a010_values[,-1]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

    a014_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a014_values) <- c("age",1950:2080)

    a014_values[,1] <- c(66:105)
    a014_values[,-1] <- a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2];

    a015_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a015_values) <- c("age",1950:2080)

    a015_values[,1] <- c(66:105)
    a015_values[,-1] <- a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

    a016_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a016_values) <- c("age",1950:2080)

    a016_values[,1] <- c(66:105)
    a016_values[,-1] <- a010_values[,-1]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

    a017_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a017_values) <- c("age",1950:2080)

    a017_values[,1] <- c(66:105)
    a017_values[,-1] <- a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

### Global incidence of disease

    a01_global_values <- matrix(c(0),
                                nrow=40,
                                ncol=132,
                                byrow=T);
    colnames(a01_global_values) <- c("age",1950:2080)

    a01_global_values[,1] <- c(66:105)
    a01_global_values[,-1] <- data_prev_0_values[which(data_prev_0_values[,1] != 65 & data_prev_0_values[,3]%in%(gender)),2]*a010_values[,-1] + data_prev_1_values[which(data_prev_1_values[,1] != 65 & data_prev_1_values[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2] + data_prev_2_values[which(data_prev_2_values[,1] != 65 & data_prev_2_values[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2] + data_prev_3_values[which(data_prev_3_values[,1] != 65 & data_prev_3_values[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2] + data_prev_4_values[which(data_prev_4_values[,1] != 65 & data_prev_4_values[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2] + data_prev_5_values[which(data_prev_5_values[,1] != 65 & data_prev_5_values[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2] + data_prev_6_values[which(data_prev_6_values[,1] != 65 & data_prev_6_values[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2] + data_prev_7_values[which(data_prev_7_values[,1] != 65 & data_prev_7_values[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

### Mortality of healthy subjects on non exposed peoples

    a020_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a020_values) <- c("age",1950:2080)

    a020_values[,1] <- c(66:105)

    for (a in 2:ncol(a020_values)){

        a020_values[,a] <- as.numeric(data_a02_values[which(data_a02_values[,1] != 65 & data_a02_values[,2]%in%(gender)),a+1]) / (data_prev_0_values[which(data_prev_0_values[,1] != 65 & data_prev_0_values[,3]%in%(gender)),2] + data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_prev_1_values[which(data_prev_1_values[,1] != 65 & data_prev_1_values[,3]%in%(gender)),2] + data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_prev_2_values[which(data_prev_2_values[,1] != 65 & data_prev_2_values[,3]%in%(gender)),2] + data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2]*data_prev_3_values[which(data_prev_3_values[,1] != 65 & data_prev_3_values[,3]%in%(gender)),2] + data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_prev_4_values[which(data_prev_4_values[,1] != 65 & data_prev_4_values[,3]%in%(gender)),2] + data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2]*data_prev_5_values[which(data_prev_5_values[,1] != 65 & data_prev_5_values[,3]%in%(gender)),2] + data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2]*data_prev_6_values[which(data_prev_6_values[,1] != 65 & data_prev_6_values[,3]%in%(gender)),2] + data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2]*data_prev_7_values[which(data_prev_7_values[,1] != 65 & data_prev_7_values[,3]%in%(gender)),2]);

    }

### Mortality of healthy subjects on exposed peoples

    a021_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a021_values) <- c("age",1950:2080)

    a021_values[,1] <- c(66:105)
    a021_values[,-1] <- a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2];

    a022_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a022_values) <- c("age",1950:2080)

    a022_values[,1] <- c(66:105)
    a022_values[,-1] <- a020_values[,-1]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2];

    a023_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a023_values) <- c("age",1950:2080)

    a023_values[,1] <- c(66:105)
    a023_values[,-1] <- a020_values[,-1]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

    a024_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a024_values) <- c("age",1950:2080)

    a024_values[,1] <- c(66:105)
    a024_values[,-1] <- a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2];

    a025_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a025_values) <- c("age",1950:2080)

    a025_values[,1] <- c(66:105)
    a025_values[,-1] <- a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

    a026_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a026_values) <- c("age",1950:2080)

    a026_values[,1] <- c(66:105)
    a026_values[,-1] <- a020_values[,-1]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

    a027_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a027_values) <- c("age",1950:2080)

    a027_values[,1] <- c(66:105)
    a027_values[,-1] <- a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

### Global mortality of healthy subjects

    a02_global_values <- matrix(c(0),
                                nrow=40,
                                ncol=132,
                                byrow=T);
    colnames(a02_global_values) <- c("age",1950:2080)

    a02_global_values[,1] <- c(66:105)
    a02_global_values[,-1] <- data_prev_0_values[which(data_prev_0_values[,1] != 65 & data_prev_0_values[,3]%in%(gender)),2]*a020_values[,-1] + data_prev_1_values[which(data_prev_1_values[,1] != 65 & data_prev_1_values[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2] + data_prev_2_values[which(data_prev_2_values[,1] != 65 & data_prev_2_values[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2] + data_prev_3_values[which(data_prev_3_values[,1] != 65 & data_prev_3_values[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2] + data_prev_4_values[which(data_prev_4_values[,1] != 65 & data_prev_4_values[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2] + data_prev_5_values[which(data_prev_5_values[,1] != 65 & data_prev_5_values[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2] + data_prev_6_values[which(data_prev_6_values[,1] != 65 & data_prev_6_values[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2] + data_prev_7_values[which(data_prev_7_values[,1] != 65 & data_prev_7_values[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

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

    for (a in 2:ncol(a120_values)){

        a120_values[,a] <- as.numeric(RR_values[,2])*as.numeric(data_a02_values[which(data_a02_values[,1] != 65 & data_a02_values[,2]%in%(gender)),a+1]) / (data_prev_0_values[which(data_prev_0_values[,1] != 65 & data_prev_0_values[,3]%in%(gender)),2] + data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_prev_1_values[which(data_prev_1_values[,1] != 65 & data_prev_1_values[,3]%in%(gender)),2] + data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_prev_2_values[which(data_prev_2_values[,1] != 65 & data_prev_2_values[,3]%in%(gender)),2] + data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2]*data_prev_3_values[which(data_prev_3_values[,1] != 65 & data_prev_3_values[,3]%in%(gender)),2] + data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_prev_4_values[which(data_prev_4_values[,1] != 65 & data_prev_4_values[,3]%in%(gender)),2] + data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2]*data_prev_5_values[which(data_prev_5_values[,1] != 65 & data_prev_5_values[,3]%in%(gender)),2] + data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2]*data_prev_6_values[which(data_prev_6_values[,1] != 65 & data_prev_6_values[,3]%in%(gender)),2] + data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2]*data_prev_7_values[which(data_prev_7_values[,1] != 65 & data_prev_7_values[,3]%in%(gender)),2]);

    }

### Mortality of diseased subjects on exposed peoples

    a121_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a121_values) <- c("age",1950:2080)

    a121_values[,1] <- c(66:105)
    a121_values[,-1] <- a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2];

    a122_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a122_values) <- c("age",1950:2080)

    a122_values[,1] <- c(66:105)
    a122_values[,-1] <- a120_values[,-1]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2];

    a123_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a123_values) <- c("age",1950:2080)

    a123_values[,1] <- c(66:105)
    a123_values[,-1] <- a120_values[,-1]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

    a124_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a124_values) <- c("age",1950:2080)

    a124_values[,1] <- c(66:105)
    a124_values[,-1] <- a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2];

    a125_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a125_values) <- c("age",1950:2080)

    a125_values[,1] <- c(66:105)
    a125_values[,-1] <- a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

    a126_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a126_values) <- c("age",1950:2080)

    a126_values[,1] <- c(66:105)
    a126_values[,-1] <- a120_values[,-1]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

    a127_values <- matrix(c(0),
                          nrow=40,
                          ncol=132,
                          byrow=T);
    colnames(a127_values) <- c("age",1950:2080)

    a127_values[,1] <- c(66:105)
    a127_values[,-1] <- a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

### Global mortality of diseased subjects

    a12_global_values <- matrix(c(0),
                                nrow=40,
                                ncol=132,
                                byrow=T);
    colnames(a12_global_values) <- c("age",1950:2080)

    a12_global_values[,1] <- c(66:105)
    a12_global_values[,-1] <- data_prev_0_values[which(data_prev_0_values[,1] != 65 & data_prev_0_values[,3]%in%(gender)),2]*a120_values[,-1] + data_prev_1_values[which(data_prev_1_values[,1] != 65 & data_prev_1_values[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2] + data_prev_2_values[which(data_prev_2_values[,1] != 65 & data_prev_2_values[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2] + data_prev_3_values[which(data_prev_3_values[,1] != 65 & data_prev_3_values[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2] + data_prev_4_values[which(data_prev_4_values[,1] != 65 & data_prev_4_values[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2] + data_prev_5_values[which(data_prev_5_values[,1] != 65 & data_prev_5_values[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2] + data_prev_6_values[which(data_prev_6_values[,1] != 65 & data_prev_6_values[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2] + data_prev_7_values[which(data_prev_7_values[,1] != 65 & data_prev_7_values[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

    ## diabete incidence
    if(!missing(incid_global))
    {
        data_incid_0_values <- incid_global
        data_incid_1_values <- incid_global
        data_incid_3_values <- incid_global
        data_incid_5_values <- incid_global
    }
    
### Matrix for results of health indicators

### Overall life-expectancy

    esp_vie_gen <- matrix(c(0),
                          nrow=100-65+1,
                          ncol=2,
                          byrow = T);
    esp_vie_gen[,1] <- c(65:100);

### Overall life-expectancy on exposed peoples

    esp_vie_gen_conso <- matrix(c(0),
                                nrow=100-65+1,
                                ncol=2,
                                byrow = T);
    esp_vie_gen_conso[,1] <- c(65:100);

### Overall life-expectancy on non exposed peoples

    esp_vie_gen_nonconso <- matrix(c(0),
                                   nrow=100-65+1,
                                   ncol=2,
                                   byrow = T);
    esp_vie_gen_nonconso[,1] <- c(65:100);

### Life-expectancy without disease

    esp_vie_sans_mal <- matrix(c(0), # matrice de calcul des espérances de vie à chaque âge
                               nrow=100-65+1,
                               ncol=2,
                               byrow = T);
    esp_vie_sans_mal[,1] <- c(65:100);

### Life-expectancy without disease on exposed peoples

    esp_vie_sans_mal_conso <- matrix(c(0),
                                     nrow=100-65+1,
                                     ncol=2,
                                     byrow = T);
    esp_vie_sans_mal_conso[,1] <- c(65:100);

### Life-expectancy without disease on non exposed peoples

    esp_vie_sans_mal_nonconso <- matrix(c(0),
                                        nrow=100-65+1,
                                        ncol=2,
                                        byrow = T);
    esp_vie_sans_mal_nonconso[,1] <- c(65:100);

### Life-expectancy for diseased subject

    esp_vie_mal <- matrix(c(0),
                          nrow=100-65+1,
                          ncol=2,
                          byrow = T);
    esp_vie_mal[,1] <- c(65:100);

### Life-expectancy for diseased subject on exposed peoples

    esp_vie_mal_conso <- matrix(c(0),
                                nrow=100-65+1,
                                ncol=2,
                                byrow = T);
    esp_vie_mal_conso[,1] <- c(65:100);

### Life-expectancy for diseased subject on non exposed peoples

    esp_vie_mal_nonconso <- matrix(c(0),
                                   nrow=100-65+1,
                                   ncol=2,
                                   byrow = T);
    esp_vie_mal_nonconso[,1] <- c(65:100);

### Life-expectancy for non-diseased subject

    esp_vie_non_mal <- matrix(c(0),
                              nrow=100-65+1,
                              ncol=2,
                              byrow = T);
    esp_vie_non_mal[,1] <- c(65:100);

### Life-expectancy for non-diseased subject on exposed peoples

    esp_vie_non_mal_conso <- matrix(c(0),
                                    nrow=100-65+1,
                                    ncol=2,
                                    byrow = T);
    esp_vie_non_mal_conso[,1] <- c(65:100);

### Life-expectancy for non-diseased subject on non exposed peoples

    esp_vie_non_mal_nonconso <- matrix(c(0),
                                       nrow=100-65+1,
                                       ncol=2,
                                       byrow = T);
    esp_vie_non_mal_nonconso[,1] <- c(65:100);

### Number of exposed peoples at least one time

    prevalence <- matrix(c(0),
                         nrow=99-65+1,
                         ncol=2,
                         byrow=T);
    prevalence[,1] <- c(65:99);

### Rate of exposed peoples at least one time

    taux_prevalence <- matrix(c(0),
                              nrow=99-65+1,
                              ncol=2,
                              byrow=T);
    taux_prevalence[,1] <- c(65:99);

### Number of living peoples

    survie <- matrix(c(0),
                     nrow=99-65+1,
                     ncol=2,
                     byrow=T);
    survie[,1] <- c(65:99);

### Rate of living peoples a given

    taux_survivants <- matrix(c(0),
                              nrow=99-65+1,
                              ncol=2,
                              byrow=T);
    taux_survivants[,1] <- c(65:99);

### Mean number of years spent with disease

    nb_moy_dem <- matrix(c(0),
                         nrow=100-65+1,
                         ncol=2,
                         byrow = T);
    nb_moy_dem[,1] <- c(65:100);

### Mean number of years spent with disease on exposed peoples

    nb_moy_dem_conso <- matrix(c(0),
                               nrow=100-65+1,
                               ncol=2,
                               byrow = T);
    nb_moy_dem_conso[,1] <- c(65:100);

### Mean number of years spent with disease on non exposed peoples

    nb_moy_dem_nonconso <- matrix(c(0),
                                  nrow=100-65+1,
                                  ncol=2,
                                  byrow = T);
    nb_moy_dem_nonconso[,1] <- c(65:100);

### Life-long probability of disease

    prb_dem <- matrix(c(0),
                      nrow=99-65+1,
                      ncol=2,
                      byrow=T);
    prb_dem[,1] <- c(65:99);

###	Average age at disease onset

    age_dem <- matrix(c(0),
                      nrow=99-65+1,
                      ncol=2,
                      byrow=T);
    age_dem[,1] <- c(65:99);

    
###	Average age of incident cases 
    age_cas_incid <- matrix(c(0),
                            nrow=99-65+1,
                            ncol=2,
                            byrow=T);
    age_cas_incid[,1] <- c(65:99);

    
### Mean number of years of exposition

    moy_conso <- matrix(c(0),
                        nrow=99-65+1,
                        ncol=2,
                        byrow=T);
    moy_conso[,1] <- c(65:99);

### Number of exposed peoples at least one time

    prevalence_conso <- matrix(c(0),
                               nrow=105-65+1,
                               ncol=2,
                               byrow=T);
    prevalence_conso[,1] <- c(65:105);

### Mortality rate

    quotient_mortalite <- matrix(c(0),
                                 nrow=99-65+1,
                                 ncol=2,
                                 byrow=T);
    quotient_mortalite[,1] <- c(65:99);


##############################
########## ESSAIS ############
##############################


### Dementia incidence

    incid_demence <- matrix(c(0),
                            nrow=105-66+1,
                            ncol=105-65+2,
                            byrow=T);
    incid_demence[,1] <- c(66:105);

### Mortality incidence

    incid_mort_Dem <- matrix(c(0),
                             nrow=105-66+1,
                             ncol=105-65+2,
                             byrow=T);
    incid_mort_Dem[,1] <- c(66:105);

    incid_mort_NoD <- matrix(c(0),
                             nrow=105-66+1,
                             ncol=105-65+2,
                             byrow=T);
    incid_mort_NoD[,1] <- c(66:105);

### Non exposed prevalence

    prevalence_noe <- matrix(c(0),
                             nrow=105-65+1,
                             ncol=105-65+2,
                             byrow=T);
    prevalence_noe[,1] <- c(65:105);

### Physical inactivity prevalence

    prevalence_ina <- matrix(c(0),
                             nrow=105-65+1,
                             ncol=105-65+2,
                             byrow=T);
    prevalence_ina[,1] <- c(65:105);

### Hypertension prevalence

    prevalence_hta <- matrix(c(0),
                             nrow=105-65+1,
                             ncol=105-65+2,
                             byrow=T);
    prevalence_hta[,1] <- c(65:105);

### Diabete prevalence

    prevalence_dia <- matrix(c(0),
                             nrow=105-65+1,
                             ncol=105-65+2,
                             byrow=T);
    prevalence_dia[,1] <- c(65:105);

### Expositions prevalence

    prevalence_exp0 <- matrix(c(0),
                              nrow=105-65+1,
                              ncol=105-65+2,
                              byrow=T);
    prevalence_exp0[,1] <- c(65:105);

    prevalence_exp1 <- matrix(c(0),
                              nrow=105-65+1,
                              ncol=105-65+2,
                              byrow=T);
    prevalence_exp1[,1] <- c(65:105);

    prevalence_exp2 <- matrix(c(0),
                              nrow=105-65+1,
                              ncol=105-65+2,
                              byrow=T);
    prevalence_exp2[,1] <- c(65:105);

    prevalence_exp3 <- matrix(c(0),
                              nrow=105-65+1,
                              ncol=105-65+2,
                              byrow=T);
    prevalence_exp3[,1] <- c(65:105);

    prevalence_exp4 <- matrix(c(0),
                              nrow=105-65+1,
                              ncol=105-65+2,
                              byrow=T);
    prevalence_exp4[,1] <- c(65:105);

    prevalence_exp5 <- matrix(c(0),
                              nrow=105-65+1,
                              ncol=105-65+2,
                              byrow=T);
    prevalence_exp5[,1] <- c(65:105);

    prevalence_exp6 <- matrix(c(0),
                              nrow=105-65+1,
                              ncol=105-65+2,
                              byrow=T);
    prevalence_exp6[,1] <- c(65:105);

    prevalence_exp7 <- matrix(c(0),
                              nrow=105-65+1,
                              ncol=105-65+2,
                              byrow=T);
    prevalence_exp7[,1] <- c(65:105);

    prop_dem_diabet <- matrix(c(0),
                              nrow=105-66+1,
                              ncol=105-65+2,
                              byrow=T);
    prop_dem_diabet[,1] <- c(66:105);

    prop_dem_hypert <- matrix(c(0),
                              nrow=105-66+1,
                              ncol=105-65+2,
                              byrow=T);
    prop_dem_hypert[,1] <- c(66:105);

    prop_dem_inacti <- matrix(c(0),
                              nrow=105-66+1,
                              ncol=105-65+2,
                              byrow=T);
    prop_dem_inacti[,1] <- c(66:105);

    prop_dem_global <- matrix(c(0),
                              nrow=105-66+1,
                              ncol=105-65+2,
                              byrow=T);
    prop_dem_global[,1] <- c(66:105);

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

### Computation with estimated parameters

    for (age in 65:105) { # for each generation

        an_naiss <- year_proj-age; # year of birth

        
        an0 <- an_naiss + 65; # first year at risk

        annee <- an0; # year of estimation

        data_prev_0_values_ND <- data_prev_0_values
        data_prev_0_values_D <- data_prev_0_values
        data_prev_1_values_ND <- data_prev_1_values
        data_prev_1_values_D <- data_prev_1_values
        data_prev_2_values_ND <- data_prev_2_values
        data_prev_2_values_D <- data_prev_2_values
        data_prev_3_values_ND <- data_prev_3_values
        data_prev_3_values_D <- data_prev_3_values
        data_prev_4_values_ND <- data_prev_4_values
        data_prev_4_values_D <- data_prev_4_values
        data_prev_5_values_ND <- data_prev_5_values
        data_prev_5_values_D <- data_prev_5_values
        data_prev_6_values_ND <- data_prev_6_values
        data_prev_6_values_D <- data_prev_6_values
        data_prev_7_values_ND <- data_prev_7_values
        data_prev_7_values_D <- data_prev_7_values

        a010_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a010_values) <- c("age",1950:2080)

        a010_values[,1] <- c(66:105)

        a011_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a011_values) <- c("age",1950:2080)

        a011_values[,1] <- c(66:105)

        a012_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a012_values) <- c("age",1950:2080)

        a012_values[,1] <- c(66:105)

        a013_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a013_values) <- c("age",1950:2080)

        a013_values[,1] <- c(66:105)

        a014_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a014_values) <- c("age",1950:2080)

        a014_values[,1] <- c(66:105)

        a015_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a015_values) <- c("age",1950:2080)

        a015_values[,1] <- c(66:105)

        a016_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a016_values) <- c("age",1950:2080)

        a016_values[,1] <- c(66:105)

        a017_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a017_values) <- c("age",1950:2080)

        a017_values[,1] <- c(66:105)

        a01_global_values <- matrix(c(0),
                                    nrow=40,
                                    ncol=132,
                                    byrow=T);
        colnames(a01_global_values) <- c("age",1950:2080)

        a01_global_values[,1] <- c(66:105)

        a020_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a020_values) <- c("age",1950:2080)

        a020_values[,1] <- c(66:105)

        a021_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a021_values) <- c("age",1950:2080)

        a021_values[,1] <- c(66:105)

        a022_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a022_values) <- c("age",1950:2080)

        a022_values[,1] <- c(66:105)

        a023_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a023_values) <- c("age",1950:2080)

        a023_values[,1] <- c(66:105)

        a024_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a024_values) <- c("age",1950:2080)

        a024_values[,1] <- c(66:105)

        a025_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a025_values) <- c("age",1950:2080)

        a025_values[,1] <- c(66:105)

        a026_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a026_values) <- c("age",1950:2080)

        a026_values[,1] <- c(66:105)

        a027_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a027_values) <- c("age",1950:2080)

        a027_values[,1] <- c(66:105)

        a02_global_values <- matrix(c(0),
                                    nrow=40,
                                    ncol=132,
                                    byrow=T);
        colnames(a02_global_values) <- c("age",1950:2080)

        a02_global_values[,1] <- c(66:105)

        a120_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a120_values) <- c("age",1950:2080)

        a120_values[,1] <- c(66:105)

        a121_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a121_values) <- c("age",1950:2080)

        a121_values[,1] <- c(66:105)

        a122_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a122_values) <- c("age",1950:2080)

        a122_values[,1] <- c(66:105)

        a123_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a123_values) <- c("age",1950:2080)

        a123_values[,1] <- c(66:105)

        a124_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a124_values) <- c("age",1950:2080)

        a124_values[,1] <- c(66:105)

        a125_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a125_values) <- c("age",1950:2080)

        a125_values[,1] <- c(66:105)

        a126_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a126_values) <- c("age",1950:2080)

        a126_values[,1] <- c(66:105)

        a127_values <- matrix(c(0),
                              nrow=40,
                              ncol=132,
                              byrow=T);
        colnames(a127_values) <- c("age",1950:2080)

        a127_values[,1] <- c(66:105)

        a12_global_values <- matrix(c(0),
                                    nrow=40,
                                    ncol=132,
                                    byrow=T);
        colnames(a12_global_values) <- c("age",1950:2080)

        a12_global_values[,1] <- c(66:105)

        RR_values <- matrix(c(0),
                            nrow=40,
                            ncol=2,
                            byrow=T);
        colnames(RR_values) <- c("age","rr_DvsND")

        RR_values[,1] <- c(66:105)

        etat <- matrix(c(0), # initial state (non diseased)
                       nrow=nb_people,
                       ncol=105-65+1,
                       byrow=T);

        colnames(etat) <- c(65:105);

                                        # RUN 1

        for (i in 1:nrow(etat)) {

            alea0 <- runif(1, 0, 1);

            if (alea0 <= data_prev_0_values[which(data_prev_0_values[,1]%in%(65) & data_prev_0_values[,3]%in%(gender)),2]) {

                etat[i,1] <- "00" # non diseased and non exposed

            } else {

                if (alea0 <= data_prev_0_values[which(data_prev_0_values[,1]%in%(65) & data_prev_0_values[,3]%in%(gender)),2]+data_prev_1_values[which(data_prev_1_values[,1]%in%(65) & data_prev_1_values[,3]%in%(gender)),2]) {

                    etat[i,1] <- "01" # non diseased and exposed 1

                } else {

                    if (alea0 <= data_prev_0_values[which(data_prev_0_values[,1]%in%(65) & data_prev_0_values[,3]%in%(gender)),2]+data_prev_1_values[which(data_prev_1_values[,1]%in%(65) & data_prev_1_values[,3]%in%(gender)),2]+data_prev_2_values[which(data_prev_2_values[,1]%in%(65) & data_prev_2_values[,3]%in%(gender)),2]) {

                        etat[i,1] <- "02" # non diseased and exposed 2

                    } else {

                        if (alea0 <= data_prev_0_values[which(data_prev_0_values[,1]%in%(65) & data_prev_0_values[,3]%in%(gender)),2]+data_prev_1_values[which(data_prev_1_values[,1]%in%(65) & data_prev_1_values[,3]%in%(gender)),2]+data_prev_2_values[which(data_prev_2_values[,1]%in%(65) & data_prev_2_values[,3]%in%(gender)),2]+data_prev_3_values[which(data_prev_3_values[,1]%in%(65) & data_prev_3_values[,3]%in%(gender)),2]) {

                            etat[i,1] <- "03" # non diseased and exposed 3

                        } else {

                            if (alea0 <= data_prev_0_values[which(data_prev_0_values[,1]%in%(65) & data_prev_0_values[,3]%in%(gender)),2]+data_prev_1_values[which(data_prev_1_values[,1]%in%(65) & data_prev_1_values[,3]%in%(gender)),2]+data_prev_2_values[which(data_prev_2_values[,1]%in%(65) & data_prev_2_values[,3]%in%(gender)),2]+data_prev_3_values[which(data_prev_3_values[,1]%in%(65) & data_prev_3_values[,3]%in%(gender)),2]+data_prev_4_values[which(data_prev_4_values[,1]%in%(65) & data_prev_4_values[,3]%in%(gender)),2]) {

                                etat[i,1] <- "04" # non diseased and exposed 4

                            } else {

                                if (alea0 <= data_prev_0_values[which(data_prev_0_values[,1]%in%(65) & data_prev_0_values[,3]%in%(gender)),2]+data_prev_1_values[which(data_prev_1_values[,1]%in%(65) & data_prev_1_values[,3]%in%(gender)),2]+data_prev_2_values[which(data_prev_2_values[,1]%in%(65) & data_prev_2_values[,3]%in%(gender)),2]+data_prev_3_values[which(data_prev_3_values[,1]%in%(65) & data_prev_3_values[,3]%in%(gender)),2]+data_prev_4_values[which(data_prev_4_values[,1]%in%(65) & data_prev_4_values[,3]%in%(gender)),2]+data_prev_5_values[which(data_prev_5_values[,1]%in%(65) & data_prev_5_values[,3]%in%(gender)),2]) {

                                    etat[i,1] <- "05" # non diseased and exposed 5

                                } else {

                                    if (alea0 <= data_prev_0_values[which(data_prev_0_values[,1]%in%(65) & data_prev_0_values[,3]%in%(gender)),2]+data_prev_1_values[which(data_prev_1_values[,1]%in%(65) & data_prev_1_values[,3]%in%(gender)),2]+data_prev_2_values[which(data_prev_2_values[,1]%in%(65) & data_prev_2_values[,3]%in%(gender)),2]+data_prev_3_values[which(data_prev_3_values[,1]%in%(65) & data_prev_3_values[,3]%in%(gender)),2]+data_prev_4_values[which(data_prev_4_values[,1]%in%(65) & data_prev_4_values[,3]%in%(gender)),2]+data_prev_5_values[which(data_prev_5_values[,1]%in%(65) & data_prev_5_values[,3]%in%(gender)),2]+data_prev_6_values[which(data_prev_6_values[,1]%in%(65) & data_prev_6_values[,3]%in%(gender)),2]) {

                                        etat[i,1] <- "06" # non diseased and exposed 6

                                    } else {

                                        etat[i,1] <- "07" # non diseased and exposed 7

                                    }

                                }

                            }

                        }

                    }

                }

            }

        };

        for (j in 2:ncol(etat)) { # for each age
            
            if(!missing(incid_global)){
                ## proportion of status 0 (non exposed), 1 (HTA), 3 (inact) and 5 (HTA+inact) among non diabetic subjects :
                nondiab <- sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))
                
                if(nondiab != 0){
                    prop0 <- sum(etat[,j-1] %in% c("00","10")) / sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))
                    prop1 <- sum(etat[,j-1] %in% c("01","11")) / sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))
                    prop3 <- sum(etat[,j-1] %in% c("03","13")) / sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))
                    prop5 <- sum(etat[,j-1] %in% c("05","15")) / sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))
                    
                    ## diabete incidence for each status :
                    jj <- which(data_incid_0_values[,1]%in%(j-1+65) & data_incid_0_values[,3]%in%(gender))
                    data_incid_0_values[jj,2] <- incid_global[jj,2]/(prop0 + prop1*theta1 + prop3*theta3 + prop5*theta1*theta3)
                    data_incid_1_values[jj,2] <- data_incid_0_values[jj,2]*theta1
                    data_incid_3_values[jj,2] <- data_incid_0_values[jj,2]*theta3
                    data_incid_5_values[jj,2] <- data_incid_0_values[jj,2]*theta1*theta3
                } else {
                    ## plus aucun non diabetique vivant donc data_incid ne sera plus utilise
                    cat("Nondiab=0 pour j=",j," annee=",annee," etats en j-1 : \n");print(table(etat[,j-1]))
                    data_incid_0_values[jj,2] <- -Inf
                    data_incid_1_values[jj,2] <- -Inf
                    data_incid_3_values[jj,2] <- -Inf
                    data_incid_5_values[jj,2] <- -Inf
                    
                }
            }
            
            for (i in 1:nrow(etat)) { # for each subject

                alea <- runif(1, 0, 1);

                alea0 <- runif(1, 0, 1);

                ## ## verif
                ## if(i==1){
                ## cat("annee=",annee," j=",j,"\n")
                ## print(which(data_incid_1_values[,1]%in%(j-1+65) & data_incid_1_values[,3]%in%(gender)))
                ## print(data_incid_0_values[which(data_incid_0_values[,1]%in%(j-1+65) & data_incid_0_values[,3]%in%(gender)),2])
                ## print(data_incid_1_values[which(data_incid_0_values[,1]%in%(j-1+65) & data_incid_0_values[,3]%in%(gender)),2])
                ## print(data_incid_3_values[which(data_incid_0_values[,1]%in%(j-1+65) & data_incid_0_values[,3]%in%(gender)),2])
                ## print(data_incid_5_values[which(data_incid_0_values[,1]%in%(j-1+65) & data_incid_0_values[,3]%in%(gender)),2])
                ## }


                
                if (etat[i,j-1] == "00") {

                    if (alea0 <= data_incid_0_values[which(data_incid_0_values[,1]%in%(j-1+65) & data_incid_0_values[,3]%in%(gender)),2]) {

                        etat[i,j] <- "02"; # non diseased and exposed 2

                    } else {

                        etat[i,j] <- "00"; # non diseased and non exposed

                    }

                } else {

                    if (etat[i,j-1] == "01") {
                        
                        if (alea0 <= data_incid_1_values[which(data_incid_1_values[,1]%in%(j-1+65) & data_incid_1_values[,3]%in%(gender)),2]) {

                            etat[i,j] <- "04"; # non diseased and exposed 4

                        } else {

                            etat[i,j] <- "01"; # non diseased and exposed 1

                        }

                    } else {

                        if (etat[i,j-1] == "03") {

                            if (alea0 <= data_incid_3_values[which(data_incid_3_values[,1]%in%(j-1+65) & data_incid_3_values[,3]%in%(gender)),2]) {

                                etat[i,j] <- "06"; # non diseased and exposed 6

                            } else {

                                etat[i,j] <- "03"; # non diseased and exposed 3

                            }

                        } else {

                            if (etat[i,j-1] == "05") {

                                if (alea0 <= data_incid_5_values[which(data_incid_5_values[,1]%in%(j-1+65) & data_incid_5_values[,3]%in%(gender)),2]) {

                                    etat[i,j] <- "07"; # non diseased and exposed 7

                                } else {

                                    etat[i,j] <- "05"; # non diseased and exposed 5

                                }

                            } else {

                                if (etat[i,j-1] == "02") {

                                    etat[i,j] <- "02"; # non diseased and exposed 2

                                } else {

                                    if (etat[i,j-1] == "04") {

                                        etat[i,j] <- "04"; # non diseased and exposed 4

                                    } else {

                                        if (etat[i,j-1] == "06") {

                                            etat[i,j] <- "06"; # non diseased and exposed 6

                                        } else {

                                            if (etat[i,j-1] == "07") {

                                                etat[i,j] <- "07"; # non diseased and exposed 7

                                            } else {

                                                if (etat[i,j-1] == "10") {

                                                    if (alea0 <= data_incid_0_values[which(data_incid_0_values[,1]%in%(j-1+65) & data_incid_0_values[,3]%in%(gender)),2]) {

                                                        etat[i,j] <- "12"; # diseased and exposed 2

                                                    } else {

                                                        etat[i,j] <- "10"; # diseased and non exposed

                                                    }

                                                } else {

                                                    if (etat[i,j-1] == "11") {
                                                        if(is.na(alea0 <= data_incid_1_values[which(data_incid_1_values[,1]%in%(j-1+65) & data_incid_1_values[,3]%in%(gender)),2])){
                                                            jj <- which(data_incid_0_values[,1]%in%(j-1+65) & data_incid_0_values[,3]%in%(gender))
                                                            cat("Pb NA pour i=",i," which=",which(data_incid_1_values[,1]%in%(j-1+65) & data_incid_1_values[,3]%in%(gender)),
                                                                " incid=",data_incid_1_values[which(data_incid_1_values[,1]%in%(j-1+65) & data_incid_1_values[,3]%in%(gender)),2],
                                                                " alea0=",alea0,
                                                                " nondiab=",nondiab,
                                                                " prop0=",prop0,
                                                                " prop1=",prop1,
                                                                " prop3=",prop3,
                                                                " prop5=",prop5,
                                                                " theta1=",theta1,
                                                                " theta3=",theta3,
                                                                " jj=",jj,
                                                                " incid=",incid_global[jj,2],
                                                                " incid0=",data_incid_0_values[jj,2],
                                                                "\n")
                                                            print(table(etat[,j-1]))
                                                        }
                                                        if (alea0 <= data_incid_1_values[which(data_incid_1_values[,1]%in%(j-1+65) & data_incid_1_values[,3]%in%(gender)),2]) {

                                                            etat[i,j] <- "14"; # diseased and exposed 4

                                                        } else {

                                                            etat[i,j] <- "11"; # diseased and exposed 1

                                                        }

                                                    } else {

                                                        if (etat[i,j-1] == "13") {

                                                            if (alea0 <= data_incid_3_values[which(data_incid_3_values[,1]%in%(j-1+65) & data_incid_3_values[,3]%in%(gender)),2]) {

                                                                etat[i,j] <- "16"; # diseased and exposed 6

                                                            } else {

                                                                etat[i,j] <- "13"; # diseased and exposed 3

                                                            }

                                                        } else {

                                                            if (etat[i,j-1] == "15") {

                                                                if (alea0 <= data_incid_5_values[which(data_incid_5_values[,1]%in%(j-1+65) & data_incid_5_values[,3]%in%(gender)),2]) {

                                                                    etat[i,j] <- "17"; # diseased and exposed 7

                                                                } else {

                                                                    etat[i,j] <- "15"; # diseased and exposed 5

                                                                }

                                                            } else {

                                                                if (etat[i,j-1] == "12") {

                                                                    etat[i,j] <- "12"; # diseased and exposed 2

                                                                } else {

                                                                    if (etat[i,j-1] == "14") {

                                                                        etat[i,j] <- "14"; # diseased and exposed 4

                                                                    } else {

                                                                        if (etat[i,j-1] == "16") {

                                                                            etat[i,j] <- "16"; # diseased and exposed 6

                                                                        } else {

                                                                            if (etat[i,j-1] == "17") {

                                                                                etat[i,j] <- "17"; # diseased and exposed 7

                                                                            } else {

                                                                                etat[i,j] <- etat[i,j-1]

                                                                            }

                                                                        }

                                                                    }

                                                                }

                                                            }

                                                        }

                                                    }

                                                }

                                            }

                                        }

                                    }

                                }

                            }

                        }

                    }

                }

            }

            if (sum(etat[,j]%in%("00"))!=0) {

                data_prev_0_values_ND[which(data_prev_0_values_ND[,1]%in%(j+64) & data_prev_0_values_ND[,3]%in%(gender)),2] <- sum(etat[,j]%in%("00"))/sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")))

            }

            if (sum(etat[,j]%in%("10"))!=0) {

                data_prev_0_values_D[which(data_prev_0_values_D[,1]%in%(j+64) & data_prev_0_values_D[,3]%in%(gender)),2] <- sum(etat[,j]%in%("10"))/sum(etat[,j]%in%(c("10","11","12","13","14","15","16","17")))

            }

            if (sum(etat[,j]%in%("01"))!=0) {

                data_prev_1_values_ND[which(data_prev_1_values_ND[,1]%in%(j+64) & data_prev_1_values_ND[,3]%in%(gender)),2] <- sum(etat[,j]%in%("01"))/sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")))

            }

            if (sum(etat[,j]%in%("11"))!=0) {

                data_prev_1_values_D[which(data_prev_1_values_D[,1]%in%(j+64) & data_prev_1_values_D[,3]%in%(gender)),2] <- sum(etat[,j]%in%("11"))/sum(etat[,j]%in%(c("10","11","12","13","14","15","16","17")))

            }

            if (sum(etat[,j]%in%("02"))!=0) {

                data_prev_2_values_ND[which(data_prev_2_values_ND[,1]%in%(j+64) & data_prev_2_values_ND[,3]%in%(gender)),2] <- sum(etat[,j]%in%("02"))/sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")))

            }

            if (sum(etat[,j]%in%("12"))!=0) {

                data_prev_2_values_D[which(data_prev_2_values_D[,1]%in%(j+64) & data_prev_2_values_D[,3]%in%(gender)),2] <- sum(etat[,j]%in%("12"))/sum(etat[,j]%in%(c("10","11","12","13","14","15","16","17")))

            }

            if (sum(etat[,j]%in%("03"))!=0) {

                data_prev_3_values_ND[which(data_prev_3_values_ND[,1]%in%(j+64) & data_prev_3_values_ND[,3]%in%(gender)),2] <- sum(etat[,j]%in%("03"))/sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")))

            }

            if (sum(etat[,j]%in%("13"))!=0) {

                data_prev_3_values_D[which(data_prev_3_values_D[,1]%in%(j+64) & data_prev_3_values_D[,3]%in%(gender)),2] <- sum(etat[,j]%in%("13"))/sum(etat[,j]%in%(c("10","11","12","13","14","15","16","17")))

            }

            if (sum(etat[,j]%in%("04"))!=0) {

                data_prev_4_values_ND[which(data_prev_4_values_ND[,1]%in%(j+64) & data_prev_4_values_ND[,3]%in%(gender)),2] <- sum(etat[,j]%in%("04"))/sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")))

            }

            if (sum(etat[,j]%in%("14"))!=0) {

                data_prev_4_values_D[which(data_prev_4_values_D[,1]%in%(j+64) & data_prev_4_values_D[,3]%in%(gender)),2] <- sum(etat[,j]%in%("14"))/sum(etat[,j]%in%(c("10","11","12","13","14","15","16","17")))

            }

            if (sum(etat[,j]%in%("05"))!=0) {

                data_prev_5_values_ND[which(data_prev_5_values_ND[,1]%in%(j+64) & data_prev_5_values_ND[,3]%in%(gender)),2] <- sum(etat[,j]%in%("05"))/sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")))

            }

            if (sum(etat[,j]%in%("15"))!=0) {

                data_prev_5_values_D[which(data_prev_5_values_D[,1]%in%(j+64) & data_prev_5_values_D[,3]%in%(gender)),2] <- sum(etat[,j]%in%("15"))/sum(etat[,j]%in%(c("10","11","12","13","14","15","16","17")))

            }

            if (sum(etat[,j]%in%("06"))!=0) {

                data_prev_6_values_ND[which(data_prev_6_values_ND[,1]%in%(j+64) & data_prev_6_values_ND[,3]%in%(gender)),2] <- sum(etat[,j]%in%("06"))/sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")))

            }

            if (sum(etat[,j]%in%("16"))!=0) {

                data_prev_6_values_D[which(data_prev_6_values_D[,1]%in%(j+64) & data_prev_6_values_D[,3]%in%(gender)),2] <- sum(etat[,j]%in%("16"))/sum(etat[,j]%in%(c("10","11","12","13","14","15","16","17")))

            }

            if (sum(etat[,j]%in%("07"))!=0) {

                data_prev_7_values_ND[which(data_prev_7_values_ND[,1]%in%(j+64) & data_prev_7_values_ND[,3]%in%(gender)),2] <- sum(etat[,j]%in%("07"))/sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")))

            }

            if (sum(etat[,j]%in%("17"))!=0) {

                data_prev_7_values_D[which(data_prev_7_values_D[,1]%in%(j+64) & data_prev_7_values_D[,3]%in%(gender)),2] <- sum(etat[,j]%in%("17"))/sum(etat[,j]%in%(c("10","11","12","13","14","15","16","17")))

            }

### Incidence of disease on non exposed peoples

            for (a in 2:ncol(a010_values)){

                a010_values[,a] <- as.numeric(data_a01_values[which(data_a01_values[,1] != 65 & data_a01_values[,2]%in%(gender)),a+1]) / (data_prev_0_values_ND[which(data_prev_0_values_ND[,1] != 65 & data_prev_0_values_ND[,3]%in%(gender)),2] + data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_prev_1_values_ND[which(data_prev_1_values_ND[,1] != 65 & data_prev_1_values_ND[,3]%in%(gender)),2] + data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_prev_2_values_ND[which(data_prev_2_values_ND[,1] != 65 & data_prev_2_values_ND[,3]%in%(gender)),2] + data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2]*data_prev_3_values_ND[which(data_prev_3_values_ND[,1] != 65 & data_prev_3_values_ND[,3]%in%(gender)),2] + data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_prev_4_values_ND[which(data_prev_4_values_ND[,1] != 65 & data_prev_4_values_ND[,3]%in%(gender)),2] + data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2]*data_prev_5_values_ND[which(data_prev_5_values_ND[,1] != 65 & data_prev_5_values_ND[,3]%in%(gender)),2] + data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2]*data_prev_6_values_ND[which(data_prev_6_values_ND[,1] != 65 & data_prev_6_values_ND[,3]%in%(gender)),2] + data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2]*data_prev_7_values_ND[which(data_prev_7_values_ND[,1] != 65 & data_prev_7_values_ND[,3]%in%(gender)),2]);

            }

### Incidence of disease on exposed peoples

            a011_values[,-1] <- a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2];

            a012_values[,-1] <- a010_values[,-1]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2];

            a013_values[,-1] <- a010_values[,-1]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

            a014_values[,-1] <- a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2];

            a015_values[,-1] <- a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

            a016_values[,-1] <- a010_values[,-1]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

            a017_values[,-1] <- a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

### Global incidence of disease

            a01_global_values[,-1] <- data_prev_0_values_ND[which(data_prev_0_values_ND[,1] != 65 & data_prev_0_values_ND[,3]%in%(gender)),2]*a010_values[,-1] + data_prev_1_values_ND[which(data_prev_1_values_ND[,1] != 65 & data_prev_1_values_ND[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2] + data_prev_2_values_ND[which(data_prev_2_values_ND[,1] != 65 & data_prev_2_values_ND[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2] + data_prev_3_values_ND[which(data_prev_3_values_ND[,1] != 65 & data_prev_3_values_ND[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2] + data_prev_4_values_ND[which(data_prev_4_values_ND[,1] != 65 & data_prev_4_values_ND[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2] + data_prev_5_values_ND[which(data_prev_5_values_ND[,1] != 65 & data_prev_5_values_ND[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2] + data_prev_6_values_ND[which(data_prev_6_values_ND[,1] != 65 & data_prev_6_values_ND[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2] + data_prev_7_values_ND[which(data_prev_7_values_ND[,1] != 65 & data_prev_7_values_ND[,3]%in%(gender)),2]*a010_values[,-1]*data_theta01_1_values[which(data_theta01_1_values[,1] != 65 & data_theta01_1_values[,3]%in%(gender)),2]*data_theta01_2_values[which(data_theta01_2_values[,1] != 65 & data_theta01_2_values[,3]%in%(gender)),2]*data_theta01_3_values[which(data_theta01_3_values[,1] != 65 & data_theta01_3_values[,3]%in%(gender)),2];

### Mortality of healthy subjects on non exposed peoples

            for (a in 2:ncol(a020_values)){

                a020_values[,a] <- as.numeric(data_a02_values[which(data_a02_values[,1] != 65 & data_a02_values[,2]%in%(gender)),a+1]) / (data_prev_0_values_ND[which(data_prev_0_values_ND[,1] != 65 & data_prev_0_values_ND[,3]%in%(gender)),2] + data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_prev_1_values_ND[which(data_prev_1_values_ND[,1] != 65 & data_prev_1_values_ND[,3]%in%(gender)),2] + data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_prev_2_values_ND[which(data_prev_2_values_ND[,1] != 65 & data_prev_2_values_ND[,3]%in%(gender)),2] + data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2]*data_prev_3_values_ND[which(data_prev_3_values_ND[,1] != 65 & data_prev_3_values_ND[,3]%in%(gender)),2] + data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_prev_4_values_ND[which(data_prev_4_values_ND[,1] != 65 & data_prev_4_values_ND[,3]%in%(gender)),2] + data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2]*data_prev_5_values_ND[which(data_prev_5_values_ND[,1] != 65 & data_prev_5_values_ND[,3]%in%(gender)),2] + data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2]*data_prev_6_values_ND[which(data_prev_6_values_ND[,1] != 65 & data_prev_6_values_ND[,3]%in%(gender)),2] + data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2]*data_prev_7_values_ND[which(data_prev_7_values_ND[,1] != 65 & data_prev_7_values_ND[,3]%in%(gender)),2]);

            }

### Mortality of healthy subjects on exposed peoples

            a021_values[,-1] <- a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2];

            a022_values[,-1] <- a020_values[,-1]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2];

            a023_values[,-1] <- a020_values[,-1]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

            a024_values[,-1] <- a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2];

            a025_values[,-1] <- a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

            a026_values[,-1] <- a020_values[,-1]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

            a027_values[,-1] <- a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

### Global mortality of healthy subjects

            a02_global_values[,-1] <- data_prev_0_values_ND[which(data_prev_0_values_ND[,1] != 65 & data_prev_0_values_ND[,3]%in%(gender)),2]*a020_values[,-1] + data_prev_1_values_ND[which(data_prev_1_values_ND[,1] != 65 & data_prev_1_values_ND[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2] + data_prev_2_values_ND[which(data_prev_2_values_ND[,1] != 65 & data_prev_2_values_ND[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2] + data_prev_3_values_ND[which(data_prev_3_values_ND[,1] != 65 & data_prev_3_values_ND[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2] + data_prev_4_values_ND[which(data_prev_4_values_ND[,1] != 65 & data_prev_4_values_ND[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2] + data_prev_5_values_ND[which(data_prev_5_values_ND[,1] != 65 & data_prev_5_values_ND[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2] + data_prev_6_values_ND[which(data_prev_6_values_ND[,1] != 65 & data_prev_6_values_ND[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2] + data_prev_7_values_ND[which(data_prev_7_values_ND[,1] != 65 & data_prev_7_values_ND[,3]%in%(gender)),2]*a020_values[,-1]*data_theta02_1_values[which(data_theta02_1_values[,1] != 65 & data_theta02_1_values[,3]%in%(gender)),2]*data_theta02_2_values[which(data_theta02_2_values[,1] != 65 & data_theta02_2_values[,3]%in%(gender)),2]*data_theta02_3_values[which(data_theta02_3_values[,1] != 65 & data_theta02_3_values[,3]%in%(gender)),2];

### Relative risks associated with the disease for mortality

            RR_values[,2] <- data_rr_DvsND_values[which(data_rr_DvsND_values[,1] != 65 & data_rr_DvsND_values[,3]%in%(gender)),2];

### Mortality of diseased subjects on non exposed peoples

            for (a in 2:ncol(a120_values)){

                a120_values[,a] <- as.numeric(RR_values[,2])*as.numeric(data_a02_values[which(data_a02_values[,1] != 65 & data_a02_values[,2]%in%(gender)),a+1]) / (data_prev_0_values_D[which(data_prev_0_values_D[,1] != 65 & data_prev_0_values_D[,3]%in%(gender)),2] + data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_prev_1_values_D[which(data_prev_1_values_D[,1] != 65 & data_prev_1_values_D[,3]%in%(gender)),2] + data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_prev_2_values_D[which(data_prev_2_values_D[,1] != 65 & data_prev_2_values_D[,3]%in%(gender)),2] + data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2]*data_prev_3_values_D[which(data_prev_3_values_D[,1] != 65 & data_prev_3_values_D[,3]%in%(gender)),2] + data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_prev_4_values_D[which(data_prev_4_values_D[,1] != 65 & data_prev_4_values_D[,3]%in%(gender)),2] + data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2]*data_prev_5_values_D[which(data_prev_5_values_D[,1] != 65 & data_prev_5_values_D[,3]%in%(gender)),2] + data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2]*data_prev_6_values_D[which(data_prev_6_values_D[,1] != 65 & data_prev_6_values_D[,3]%in%(gender)),2] + data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2]*data_prev_7_values_D[which(data_prev_7_values_D[,1] != 65 & data_prev_7_values_D[,3]%in%(gender)),2]);

            }

### Mortality of diseased subjects on exposed peoples

            a121_values[,-1] <- a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2];

            a122_values[,-1] <- a120_values[,-1]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2];

            a123_values[,-1] <- a120_values[,-1]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

            a124_values[,-1] <- a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2];

            a125_values[,-1] <- a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

            a126_values[,-1] <- a120_values[,-1]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

            a127_values[,-1] <- a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

### Global mortality of diseased subjects

            a12_global_values[,-1] <- data_prev_0_values_D[which(data_prev_0_values_D[,1] != 65 & data_prev_0_values_D[,3]%in%(gender)),2]*a120_values[,-1] + data_prev_1_values_D[which(data_prev_1_values_D[,1] != 65 & data_prev_1_values_D[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2] + data_prev_2_values_D[which(data_prev_2_values_D[,1] != 65 & data_prev_2_values_D[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2] + data_prev_3_values_D[which(data_prev_3_values_D[,1] != 65 & data_prev_3_values_D[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2] + data_prev_4_values_D[which(data_prev_4_values_D[,1] != 65 & data_prev_4_values_D[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2] + data_prev_5_values_D[which(data_prev_5_values_D[,1] != 65 & data_prev_5_values_D[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2] + data_prev_6_values_D[which(data_prev_6_values_D[,1] != 65 & data_prev_6_values_D[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2] + data_prev_7_values_D[which(data_prev_7_values_D[,1] != 65 & data_prev_7_values_D[,3]%in%(gender)),2]*a120_values[,-1]*data_theta12_1_values[which(data_theta12_1_values[,1] != 65 & data_theta12_1_values[,3]%in%(gender)),2]*data_theta12_2_values[which(data_theta12_2_values[,1] != 65 & data_theta12_2_values[,3]%in%(gender)),2]*data_theta12_3_values[which(data_theta12_3_values[,1] != 65 & data_theta12_3_values[,3]%in%(gender)),2];

            for (i in 1:nrow(etat)) { # for each age

                alea <- runif(1, 0, 1);

                alea0 <- runif(1, 0, 1);

                if (etat[i,j] == "00") {

                    a01 <- a010_values;
                    a02 <- a020_values;

                    if (alea <= a02[j-1,(an0-1950)+1]) {

                        etat[i,j] <- "20"; # dead (with non exposed state)

                    } else {

                        if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                            etat[i,j] <- "10"; # diseased and non exposed

                        } else {

                            etat[i,j] <- "00"; # non diseased and non exposed

                        }

                    }

                } else {

                    if (etat[i,j] == "01") {

                        a01 <- a011_values;
                        a02 <- a021_values;

                        if (alea <= a02[j-1,(an0-1950)+1]) {

                            etat[i,j] <- "21"; # dead (with exposed state 1)

                        } else {

                            if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                etat[i,j] <- "11"; # diseased and exposed 1

                            } else {

                                etat[i,j] <- "01"; # non diseased and exposed 1

                            }

                        }

                    } else {

                        if (etat[i,j] == "02") {

                            a01 <- a012_values;
                            a02 <- a022_values;

                            if (alea <= a02[j-1,(an0-1950)+1]) {

                                etat[i,j] <- "22"; # dead (with exposed state 2)

                            } else {

                                if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                    etat[i,j] <- "12"; # diseased and exposed 2

                                } else {

                                    etat[i,j] <- "02"; # non diseased and exposed 2

                                }

                            }

                        } else {

                            if (etat[i,j] == "03") {

                                a01 <- a013_values;
                                a02 <- a023_values;

                                if (alea <= a02[j-1,(an0-1950)+1]) {

                                    etat[i,j] <- "23"; # dead (with exposed state 3)

                                } else {

                                    if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                        etat[i,j] <- "13"; # diseased and exposed 3

                                    } else {

                                        etat[i,j] <- "03"; # non diseased and exposed 3

                                    }

                                }

                            } else {

                                if (etat[i,j] == "04") {

                                    a01 <- a014_values;
                                    a02 <- a024_values;

                                    if (alea <= a02[j-1,(an0-1950)+1]) {

                                        etat[i,j] <- "24"; # dead (with exposed state 4)

                                    } else {

                                        if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                            etat[i,j] <- "14"; # diseased and exposed 4

                                        } else {

                                            etat[i,j] <- "04"; # non diseased and exposed 4

                                        }

                                    };

                                } else {

                                    if (etat[i,j] == "05") {

                                        a01 <- a015_values;
                                        a02 <- a025_values;

                                        if (alea <= a02[j-1,(an0-1950)+1]) {

                                            etat[i,j] <- "25"; # dead (with exposed state 5)

                                        } else {

                                            if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                                etat[i,j] <- "15"; # diseased and exposed 5

                                            } else {

                                                etat[i,j] <- "05"; # non diseased and exposed 5

                                            }

                                        }

                                    } else {

                                        if (etat[i,j] == "06") {

                                            a01 <- a016_values;
                                            a02 <- a026_values;

                                            if (alea <= a02[j-1,(an0-1950)+1]) {

                                                etat[i,j] <- "26"; # dead (with exposed state 6)

                                            } else {

                                                if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                                    etat[i,j] <- "16"; # diseased and exposed 6

                                                } else {

                                                    etat[i,j] <- "06"; # non diseased and exposed 6

                                                }

                                            };

                                        } else {

                                            if (etat[i,j] == "07") {

                                                a01 <- a017_values;
                                                a02 <- a027_values;

                                                if (alea <= a02[j-1,(an0-1950)+1]) {

                                                    etat[i,j] <- "27"; # dead (with exposed state 7)

                                                } else {

                                                    if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                                        etat[i,j] <- "17"; # diseased and exposed 7

                                                    } else {

                                                        etat[i,j] <- "07"; # non diseased and exposed 7

                                                    }

                                                };

                                            } else {

                                                if (etat[i,j] == "10") {

                                                    a12 <- a120_values;

                                                    if (alea <= a12[j-1,(an0-1950)+1]) {

                                                        etat[i,j] <- "20";

                                                    } else {

                                                        etat[i,j] <- "10";

                                                    }

                                                } else {

                                                    if (etat[i,j] == "11") {

                                                        a12 <- a121_values;

                                                        if (alea <= a12[j-1,(an0-1950)+1]) {

                                                            etat[i,j] <- "21";

                                                        } else {

                                                            etat[i,j] <- "11";

                                                        }

                                                    } else {

                                                        if (etat[i,j] == "12") {

                                                            a12 <- a122_values;

                                                            if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                etat[i,j] <- "22";

                                                            } else {

                                                                etat[i,j] <- "12";

                                                            }

                                                        } else {

                                                            if (etat[i,j] == "13") {

                                                                a12 <- a123_values;

                                                                if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                    etat[i,j] <- "23";

                                                                } else {

                                                                    etat[i,j] <- "13";

                                                                }

                                                            } else {

                                                                if (etat[i,j] == "14") {

                                                                    a12 <- a124_values;

                                                                    if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                        etat[i,j] <- "24";

                                                                    } else {

                                                                        etat[i,j] <- "14";

                                                                    }

                                                                } else {

                                                                    if (etat[i,j] == "15") {

                                                                        a12 <- a125_values;

                                                                        if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                            etat[i,j] <- "25";

                                                                        } else {

                                                                            etat[i,j] <- "15";

                                                                        }

                                                                    } else {

                                                                        if (etat[i,j] == "16") {

                                                                            a12 <- a126_values;

                                                                            if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                                etat[i,j] <- "26";

                                                                            } else {

                                                                                etat[i,j] <- "16";

                                                                            }

                                                                        } else {

                                                                            if (etat[i,j] == "17") {

                                                                                a12 <- a127_values;

                                                                                if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                                    etat[i,j] <- "27";

                                                                                } else {

                                                                                    etat[i,j] <- "17";

                                                                                }

                                                                            } else {

                                                                                etat[i,j] <- etat[i,j]

                                                                            }

                                                                        }

                                                                    }

                                                                }

                                                            }

                                                        }

                                                    }

                                                }

                                            }

                                        }

                                    }

                                }

                            }

                        }

                    }

                }

            }

        }


############ RUN 2 ############
        
        if(intervention!=0){
            for (i in 1:nrow(etat)) {
                intervention_prev_0_values <- data_prev_0_values
                intervention_prev_1_values <- data_prev_1_values
                intervention_prev_2_values <- data_prev_2_values
                intervention_prev_3_values <- data_prev_3_values
                intervention_prev_4_values <- data_prev_4_values
                intervention_prev_5_values <- data_prev_5_values
                intervention_prev_6_values <- data_prev_6_values
                intervention_prev_7_values <- data_prev_7_values


                if (intervention %in% c(1,3,4)) {

                    if (annee >= year_intervention) {
                        intervention_prev_0_values[,2] <- 1
                        intervention_prev_1_values[,2] <- 0
                        intervention_prev_2_values[,2] <- 0
                        intervention_prev_3_values[,2] <- 0
                        intervention_prev_4_values[,2] <- 0
                        intervention_prev_5_values[,2] <- 0
                        intervention_prev_6_values[,2] <- 0
                        intervention_prev_7_values[,2] <- 0
                    }

                }

                if(intervention==2) {
                    intervention_prev_0_values[,2] <- 1
                    intervention_prev_1_values[,2] <- 0
                    intervention_prev_2_values[,2] <- 0
                    intervention_prev_3_values[,2] <- 0
                    intervention_prev_4_values[,2] <- 0
                    intervention_prev_5_values[,2] <- 0
                    intervention_prev_6_values[,2] <- 0
                    intervention_prev_7_values[,2] <- 0            
                }

                if(intervention=="HTA"){
                    if (annee >= year_intervention) {
                        intervention_prev_0_values[,2] <- intervention_prev_0_values[,2] + intervention_prev_1_values[,2]
                        intervention_prev_1_values[,2] <- 0
                        intervention_prev_2_values[,2] <- intervention_prev_2_values[,2] + intervention_prev_4_values[,2]
                        intervention_prev_3_values[,2] <- intervention_prev_3_values[,2] + intervention_prev_5_values[,2]
                        intervention_prev_4_values[,2] <- 0
                        intervention_prev_5_values[,2] <- 0
                        intervention_prev_6_values[,2] <- intervention_prev_6_values[,2] + intervention_prev_7_values[,2]
                        intervention_prev_7_values[,2] <- 0
                    }
                }

                if(intervention=="diab"){
                    if (annee >= year_intervention) {
                        intervention_prev_0_values[,2] <- intervention_prev_0_values[,2] + intervention_prev_2_values[,2]
                        intervention_prev_1_values[,2] <- intervention_prev_1_values[,2] + intervention_prev_4_values[,2]
                        intervention_prev_2_values[,2] <- 0
                        intervention_prev_3_values[,2] <- intervention_prev_3_values[,2] + intervention_prev_6_values[,2]
                        intervention_prev_4_values[,2] <- 0
                        intervention_prev_5_values[,2] <- intervention_prev_5_values[,2] + intervention_prev_7_values[,2]
                        intervention_prev_6_values[,2] <- 0
                        intervention_prev_7_values[,2] <- 0
                    }
                }

                if(intervention=="inact"){
                    if (annee >= year_intervention) {
                        intervention_prev_0_values[,2] <- intervention_prev_0_values[,2] + intervention_prev_3_values[,2]
                        intervention_prev_1_values[,2] <- intervention_prev_1_values[,2] + intervention_prev_5_values[,2]
                        intervention_prev_2_values[,2] <- intervention_prev_2_values[,2] + intervention_prev_6_values[,2]
                        intervention_prev_3_values[,2] <- 0
                        intervention_prev_4_values[,2] <- intervention_prev_4_values[,2] + intervention_prev_7_values[,2]
                        intervention_prev_5_values[,2] <- 0
                        intervention_prev_6_values[,2] <- 0
                        intervention_prev_7_values[,2] <- 0
                    }
                }


                alea0 <- runif(1, 0, 1);

                if (alea0 <= intervention_prev_0_values[which(intervention_prev_0_values[,1]%in%(65) & intervention_prev_0_values[,3]%in%(gender)),2]) {

                    etat[i,1] <- "00" # non diseased and non exposed

                } else {

                    if (alea0 <= intervention_prev_0_values[which(intervention_prev_0_values[,1]%in%(65) & intervention_prev_0_values[,3]%in%(gender)),2]+intervention_prev_1_values[which(intervention_prev_1_values[,1]%in%(65) & intervention_prev_1_values[,3]%in%(gender)),2]) {

                        etat[i,1] <- "01" # non diseased and exposed 1

                    } else {

                        if (alea0 <= intervention_prev_0_values[which(intervention_prev_0_values[,1]%in%(65) & intervention_prev_0_values[,3]%in%(gender)),2]+intervention_prev_1_values[which(intervention_prev_1_values[,1]%in%(65) & intervention_prev_1_values[,3]%in%(gender)),2]+intervention_prev_2_values[which(intervention_prev_2_values[,1]%in%(65) & intervention_prev_2_values[,3]%in%(gender)),2]) {

                            etat[i,1] <- "02" # non diseased and exposed 2

                        } else {

                            if (alea0 <= intervention_prev_0_values[which(intervention_prev_0_values[,1]%in%(65) & intervention_prev_0_values[,3]%in%(gender)),2]+intervention_prev_1_values[which(intervention_prev_1_values[,1]%in%(65) & intervention_prev_1_values[,3]%in%(gender)),2]+intervention_prev_2_values[which(intervention_prev_2_values[,1]%in%(65) & intervention_prev_2_values[,3]%in%(gender)),2]+intervention_prev_3_values[which(intervention_prev_3_values[,1]%in%(65) & intervention_prev_3_values[,3]%in%(gender)),2]) {

                                etat[i,1] <- "03" # non diseased and exposed 3

                            } else {

                                if (alea0 <= intervention_prev_0_values[which(intervention_prev_0_values[,1]%in%(65) & intervention_prev_0_values[,3]%in%(gender)),2]+intervention_prev_1_values[which(intervention_prev_1_values[,1]%in%(65) & intervention_prev_1_values[,3]%in%(gender)),2]+intervention_prev_2_values[which(intervention_prev_2_values[,1]%in%(65) & intervention_prev_2_values[,3]%in%(gender)),2]+intervention_prev_3_values[which(intervention_prev_3_values[,1]%in%(65) & intervention_prev_3_values[,3]%in%(gender)),2]+intervention_prev_4_values[which(intervention_prev_4_values[,1]%in%(65) & intervention_prev_4_values[,3]%in%(gender)),2]) {

                                    etat[i,1] <- "04" # non diseased and exposed 4

                                } else {

                                    if (alea0 <= intervention_prev_0_values[which(intervention_prev_0_values[,1]%in%(65) & intervention_prev_0_values[,3]%in%(gender)),2]+intervention_prev_1_values[which(intervention_prev_1_values[,1]%in%(65) & intervention_prev_1_values[,3]%in%(gender)),2]+intervention_prev_2_values[which(intervention_prev_2_values[,1]%in%(65) & intervention_prev_2_values[,3]%in%(gender)),2]+intervention_prev_3_values[which(intervention_prev_3_values[,1]%in%(65) & intervention_prev_3_values[,3]%in%(gender)),2]+intervention_prev_4_values[which(intervention_prev_4_values[,1]%in%(65) & intervention_prev_4_values[,3]%in%(gender)),2]+intervention_prev_5_values[which(intervention_prev_5_values[,1]%in%(65) & intervention_prev_5_values[,3]%in%(gender)),2]) {

                                        etat[i,1] <- "05" # non diseased and exposed 5

                                    } else {

                                        if (alea0 <= intervention_prev_0_values[which(intervention_prev_0_values[,1]%in%(65) & intervention_prev_0_values[,3]%in%(gender)),2]+intervention_prev_1_values[which(intervention_prev_1_values[,1]%in%(65) & intervention_prev_1_values[,3]%in%(gender)),2]+intervention_prev_2_values[which(intervention_prev_2_values[,1]%in%(65) & intervention_prev_2_values[,3]%in%(gender)),2]+intervention_prev_3_values[which(intervention_prev_3_values[,1]%in%(65) & intervention_prev_3_values[,3]%in%(gender)),2]+intervention_prev_4_values[which(intervention_prev_4_values[,1]%in%(65) & intervention_prev_4_values[,3]%in%(gender)),2]+intervention_prev_5_values[which(intervention_prev_5_values[,1]%in%(65) & intervention_prev_5_values[,3]%in%(gender)),2]+intervention_prev_6_values[which(intervention_prev_6_values[,1]%in%(65) & intervention_prev_6_values[,3]%in%(gender)),2]) {

                                            etat[i,1] <- "06" # non diseased and exposed 6

                                        } else {

                                            etat[i,1] <- "07" # non diseased and exposed 7

                                        }

                                    }

                                }

                            }

                        }

                    }

                }

            };

            

            ## verif si intervention monofacteur
            nb1 <- sum(etat[,1]=="01")
            nb2 <- sum(etat[,1]=="02")
            nb3 <- sum(etat[,1]=="03")
            nb4 <- sum(etat[,1]=="04")
            nb5 <- sum(etat[,1]=="05")
            nb6 <- sum(etat[,1]=="06")
            nb7 <- sum(etat[,1]=="07")

            if(annee>=year_intervention){
                if(intervention=="HTA"){
                    if(nb1 != 0) stop(paste("No subject exposed 1 is expected and there are",nb1))
                    if(nb4 != 0) stop(paste("No subject exposed 4 is expected and there are",nb4))
                    if(nb5 != 0) stop(paste("No subject exposed 5 is expected and there are",nb5))
                    if(nb7 != 0) stop(paste("No subject exposed 7 is expected and there are",nb7))
                }
                if(intervention=="diab"){
                    if(nb2 != 0) stop(paste("No subject exposed 2 is expected and there are",nb2))
                    if(nb4 != 0) stop(paste("No subject exposed 4 is expected and there are",nb4))
                    if(nb6 != 0) stop(paste("No subject exposed 6 is expected and there are",nb6))
                    if(nb7 != 0) stop(paste("No subject exposed 7 is expected and there are",nb7))
                }
                if(intervention=="inact"){
                    if(nb3 != 0) stop(paste("No subject exposed 3 is expected and there are",nb3))
                    if(nb5 != 0) stop(paste("No subject exposed 5 is expected and there are",nb5))
                    if(nb6 != 0) stop(paste("No subject exposed 6 is expected and there are",nb6))
                    if(nb7 != 0) stop(paste("No subject exposed 7 is expected and there are",nb7))
                }
            }

            for (j in 2:ncol(etat)) {
                
                if(!missing(incid_global)){
                    ## proportion of status 0 (non exposed), 1 (HTA), 3 (inact) and 5 (HTA+inact) among non diabetic subjects :
                    
                    nondiab <- sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))
                    
                    if(nondiab != 0){
                        prop0 <- sum(etat[,j-1] %in% c("00","10")) / sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))
                        prop1 <- sum(etat[,j-1] %in% c("01","11")) / sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))
                        prop3 <- sum(etat[,j-1] %in% c("03","13")) / sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))
                        prop5 <- sum(etat[,j-1] %in% c("05","15")) / sum(etat[,j-1] %in% c("00","10","01","11","03","13","05","15"))

                        if(annee+j-1 == year_intervention){
                                        #cat("Changement prop, annee=",annee," j=",j,"\n")
                            if(intervention=="HTA"){
                                prop0 <- prop0 + prop1
                                prop1 <- 0
                                prop3 <- prop3 + prop5
                                prop5 <- 0
                                        #theta1 <- 0
                            }
                            if(intervention=="inact"){
                                prop0 <- prop0 + prop3
                                prop1 <- prop1 + prop5
                                prop3 <- 0
                                prop5 <- 0
                                        #theta3 <- 0
                            }
                        } 
                        
                        ## diabete incidence for each status :
                        jj <- which(data_incid_0_values[,1]%in%(j-1+65) & data_incid_0_values[,3]%in%(gender))
                        data_incid_0_values[jj,2] <- incid_global[jj,2]/(prop0 + prop1*theta1 + prop3*theta3 + prop5*theta1*theta3)
                        data_incid_1_values[jj,2] <- data_incid_0_values[jj,2]*theta1
                        data_incid_3_values[jj,2] <- data_incid_0_values[jj,2]*theta3
                        data_incid_5_values[jj,2] <- data_incid_0_values[jj,2]*theta1*theta3
                    } else {
                        ## plus aucun non diabetique vivant donc data_incid ne sera plus utilise
                        cat("Nondiab=0 pour j=",j," annee=",annee," etats en j-1 : \n");print(table(etat[,j-1]))
                        data_incid_0_values[jj,2] <- -Inf
                        data_incid_1_values[jj,2] <- -Inf
                        data_incid_3_values[jj,2] <- -Inf
                        data_incid_5_values[jj,2] <- -Inf
                    }
                    
                    ## verif si intervention monofacteur
                    if(annee+j-1 > year_intervention){
                        if(intervention=="HTA"){
                            if(prop1 != 0)  warning(paste("There are",prop1,"% subjects exposed 1 at time",j-1))
                            if(prop5 != 0) stop(paste("There are",prop5,"% subjects exposed 5 at time",j-1))
                        }
                        if(intervention=="inact"){
                            if(prop3 != 0) stop(paste("There are",prop3,"% subjects exposed 3 at time",j-1))
                            if(prop5 != 0) stop(paste("There are",prop5,"% subjects exposed 5 at time",j-1))
                        }
                    }
                }
                
                
                intervention_incid_0_values <- data_incid_0_values
                intervention_incid_1_values <- data_incid_1_values
                intervention_incid_3_values <- data_incid_3_values
                intervention_incid_5_values <- data_incid_5_values


                if(intervention=="diab"){
                    intervention_incid_0_values[,2] <- 0
                    intervention_incid_1_values[,2] <- 0 
                    intervention_incid_3_values[,2] <- 0
                    intervention_incid_5_values[,2] <- 0
                }
                
                if (intervention==1) {
                    ## pas de nouveaux diabetiques a partir de l'annee d'intervention
                    if (annee+j-1 >= year_intervention) {
                        intervention_incid_0_values <- data_incid_0_values
                        intervention_incid_0_values[,2] <- 0

                        intervention_incid_1_values <- data_incid_1_values
                        intervention_incid_1_values[,2] <- 0

                        intervention_incid_3_values <- data_incid_3_values
                        intervention_incid_3_values[,2] <- 0

                        intervention_incid_5_values <- data_incid_5_values
                        intervention_incid_5_values[,2] <- 0
                    }
                }

                if(intervention==2) {
                    intervention_incid_0_values[,2] <- 0
                    intervention_incid_1_values[,2] <- 0
                    intervention_incid_3_values[,2] <- 0
                    intervention_incid_5_values[,2] <- 0            
                }

                for (i in 1:nrow(etat)) { 
                    ## simuler l'exposition
                    alea0 <- runif(1, 0, 1);

                    ## etat au temps d'avant
                    etatavt <- etat[i,j-1]
                    exposavt <- strsplit(etat[i,j-1],split="")[[1]][2]
                    demdcavt <- strsplit(etat[i,j-1],split="")[[1]][1]

                    ## Rappel : 0=aucun, 1=HTA, 2=diab, 3=inact, 4=HTA+diab, 5=HTA+inact, 6=diab+inact, 7=HTA+diab+inact

                    ## l'annee de l'intervention, exposition change :
                    if(intervention == "HTA") {
                        exposavt <- switch(exposavt,"0"="0", "1"="0", "2"="2", "3"="3", "4"="2", "5"="3", "6"="6", "7"="6")
                    }
                    if(intervention == "diab") {
                        exposavt <- switch(exposavt,"0"="0", "1"="1", "2"="0", "3"="3" , "4"="1", "5"="5", "6"="3", "7"="5")
                    }
                    if(intervention == "inact") {
                        exposavt <- switch(exposavt,"0"="0", "1"="1", "2"="2", "3"="0", "4"="4", "5"="1", "6"="2", "7"="4")
                    }

                    
                    if((intervention==3) & (annee+j-1>=year_intervention)) {
                        ## plus aucun expose a partir de l'annee d'intervention
                        ## donc statut en j est dem/dc en j-1 et 0 pour expo

                        etat[i,j] <- paste(strsplit(etat[i,j-1],split="")[[1]][1],"0",sep="")

                    } else {
                        if((intervention==4) & (annee+j-1>=year_intervention)) {
                            ## plus aucun HTA ni diab a partir de l'annee d'intervention,  mais inact pas supprime
                            ## donc statut en j est dem/dc en j-1 et 0 ou 3 pour expo

                            exposavt <- strsplit(etat[i,j-1],split="")[[1]][2]
                            expos <- "0"
                            if(exposavt %in% c(3,5,6,7)) expos <- "3"
                            etat[i,j] <- paste(strsplit(etat[i,j-1],split="")[[1]][1],expos,sep="")

                        } else {
                            
                            if (exposavt == "0") {

                                if (alea0 <= intervention_incid_0_values[which(intervention_incid_0_values[,1]%in%(j-1+65) & intervention_incid_0_values[,3]%in%(gender)),2]) {

                                    etat[i,j] <- paste(demdcavt,"2",sep=""); # preceding disease status and exposed 2

                                } else {

                                    etat[i,j] <- paste(demdcavt,"0",sep=""); # preceding disease status and non exposed

                                }

                            } else {

                                if (exposavt == "1") {

                                    if (alea0 <= intervention_incid_1_values[which(intervention_incid_1_values[,1]%in%(j-1+65) & intervention_incid_1_values[,3]%in%(gender)),2]) {

                                        etat[i,j] <- paste(demdcavt,"4",sep=""); # newly exposed 4

                                    } else {

                                        etat[i,j] <- paste(demdcavt,"1",sep=""); # still exposed 1

                                    }

                                } else {

                                    if (exposavt == "3") {

                                        if (alea0 <= intervention_incid_3_values[which(intervention_incid_3_values[,1]%in%(j-1+65) & intervention_incid_3_values[,3]%in%(gender)),2]) {

                                            etat[i,j] <- paste(demdcavt,"6",sep=""); #newly exposed 6

                                        } else {

                                            etat[i,j] <- paste(demdcavt,"3",sep=""); # still exposed 3

                                        }

                                    } else {

                                        if (exposavt == "5") {

                                            if (alea0 <= intervention_incid_5_values[which(intervention_incid_5_values[,1]%in%(j-1+65) & intervention_incid_5_values[,3]%in%(gender)),2]) {

                                                etat[i,j] <- paste(demdcavt,"7",sep=""); # newly exposed 7

                                            } else {

                                                etat[i,j] <- paste(demdcavt,"5",sep=""); # still exposed 5

                                            }

                                        } else {

                                            etat[i,j] <- paste(demdcavt,exposavt,sep="");

                                        }

                                    }
                                    
                                }
                                
                            }
                            
                        }
                    }
                }
                
                for (i in 1:nrow(etat)) {
                    ## simuler deces ou demence

                    alea <- runif(1, 0, 1);

                    if (etat[i,j] == "00") {

                        a01 <- a010_values;
                        a02 <- a020_values;

                        if (alea <= a02[j-1,(an0-1950)+1]) {

                            etat[i,j] <- "20"; # dead (with non exposed state)

                        } else {

                            if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                etat[i,j] <- "10"; # diseased and non exposed

                            } else {

                                etat[i,j] <- "00"; # non diseased and non exposed

                            }

                        }

                    } else {

                        if (etat[i,j] == "01") {

                            a01 <- a011_values;
                            a02 <- a021_values;

                            if (alea <= a02[j-1,(an0-1950)+1]) {

                                etat[i,j] <- "21"; # dead (with exposed state 1)

                            } else {

                                if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                    etat[i,j] <- "11"; # diseased and exposed 1

                                } else {

                                    etat[i,j] <- "01"; # non diseased and exposed 1

                                }

                            }

                        } else {

                            if (etat[i,j] == "02") {

                                a01 <- a012_values;
                                a02 <- a022_values;

                                if (alea <= a02[j-1,(an0-1950)+1]) {

                                    etat[i,j] <- "22"; # dead (with exposed state 2)

                                } else {

                                    if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                        etat[i,j] <- "12"; # diseased and exposed 2

                                    } else {

                                        etat[i,j] <- "02"; # non diseased and exposed 2

                                    }

                                }

                            } else {

                                if (etat[i,j] == "03") {

                                    a01 <- a013_values;
                                    a02 <- a023_values;

                                    if (alea <= a02[j-1,(an0-1950)+1]) {

                                        etat[i,j] <- "23"; # dead (with exposed state 3)

                                    } else {

                                        if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                            etat[i,j] <- "13"; # diseased and exposed 3

                                        } else {

                                            etat[i,j] <- "03"; # non diseased and exposed 3

                                        }

                                    }

                                } else {

                                    if (etat[i,j] == "04") {

                                        a01 <- a014_values;
                                        a02 <- a024_values;

                                        if (alea <= a02[j-1,(an0-1950)+1]) {

                                            etat[i,j] <- "24"; # dead (with exposed state 4)

                                        } else {

                                            if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                                etat[i,j] <- "14"; # diseased and exposed 4

                                            } else {

                                                etat[i,j] <- "04"; # non diseased and exposed 4

                                            }

                                        };

                                    } else {

                                        if (etat[i,j] == "05") {

                                            a01 <- a015_values;
                                            a02 <- a025_values;

                                            if (alea <= a02[j-1,(an0-1950)+1]) {

                                                etat[i,j] <- "25"; # dead (with exposed state 5)

                                            } else {

                                                if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                                    etat[i,j] <- "15"; # diseased and exposed 5

                                                } else {

                                                    etat[i,j] <- "05"; # non diseased and exposed 5

                                                }

                                            }

                                        } else {

                                            if (etat[i,j] == "06") {

                                                a01 <- a016_values;
                                                a02 <- a026_values;

                                                if (alea <= a02[j-1,(an0-1950)+1]) {

                                                    etat[i,j] <- "26"; # dead (with exposed state 6)

                                                } else {

                                                    if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                                        etat[i,j] <- "16"; # diseased and exposed 6

                                                    } else {

                                                        etat[i,j] <- "06"; # non diseased and exposed 6

                                                    }

                                                };

                                            } else {

                                                if (etat[i,j] == "07") {

                                                    a01 <- a017_values;
                                                    a02 <- a027_values;

                                                    if (alea <= a02[j-1,(an0-1950)+1]) {

                                                        etat[i,j] <- "27"; # dead (with exposed state 7)

                                                    } else {

                                                        if (alea <= a01[j-1,(an0-1950)+1] + a02[j-1,(an0-1950)+1]) {

                                                            etat[i,j] <- "17"; # diseased and exposed 7

                                                        } else {

                                                            etat[i,j] <- "07"; # non diseased and exposed 7

                                                        }

                                                    };

                                                } else {

                                                    if (etat[i,j] == "10") {

                                                        a12 <- a120_values;

                                                        if (alea <= a12[j-1,(an0-1950)+1]) {

                                                            etat[i,j] <- "20";

                                                        } else {

                                                            etat[i,j] <- "10";

                                                        }

                                                    } else {

                                                        if (etat[i,j] == "11") {

                                                            a12 <- a121_values;

                                                            if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                etat[i,j] <- "21";

                                                            } else {

                                                                etat[i,j] <- "11";

                                                            }

                                                        } else {

                                                            if (etat[i,j] == "12") {

                                                                a12 <- a122_values;

                                                                if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                    etat[i,j] <- "22";

                                                                } else {

                                                                    etat[i,j] <- "12";

                                                                }

                                                            } else {

                                                                if (etat[i,j] == "13") {

                                                                    a12 <- a123_values;

                                                                    if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                        etat[i,j] <- "23";

                                                                    } else {

                                                                        etat[i,j] <- "13";

                                                                    }

                                                                } else {

                                                                    if (etat[i,j] == "14") {

                                                                        a12 <- a124_values;

                                                                        if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                            etat[i,j] <- "24";

                                                                        } else {

                                                                            etat[i,j] <- "14";

                                                                        }

                                                                    } else {

                                                                        if (etat[i,j] == "15") {

                                                                            a12 <- a125_values;

                                                                            if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                                etat[i,j] <- "25";

                                                                            } else {

                                                                                etat[i,j] <- "15";

                                                                            }

                                                                        } else {

                                                                            if (etat[i,j] == "16") {

                                                                                a12 <- a126_values;

                                                                                if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                                    etat[i,j] <- "26";

                                                                                } else {

                                                                                    etat[i,j] <- "16";

                                                                                }

                                                                            } else {

                                                                                if (etat[i,j] == "17") {

                                                                                    a12 <- a127_values;

                                                                                    if (alea <= a12[j-1,(an0-1950)+1]) {

                                                                                        etat[i,j] <- "27";

                                                                                    } else {

                                                                                        etat[i,j] <- "17";

                                                                                    }

                                                                                } else {

                                                                                    etat[i,j] <- etat[i,j]

                                                                                }

                                                                            }

                                                                        }

                                                                    }

                                                                }

                                                            }

                                                        }

                                                    }

                                                }

                                            }

                                        }

                                    }

                                }

                            }

                        }

                    }

                }
            }
        }
        
### Computation of health indicators :

        etat_vivant <- c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")
        
        etat_vivant_exp <- c("01","02","03","04","05","06","07","11","12","13","14","15","16","17")
        
        etat_vivant_n_exp <- c("00","10")
        
### Overall life-expectancy

        if (age < 101) {

            n0 <- vector(length = ncol(etat));

            s0 <- sum(etat[,age-64]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

            if (s0 != 0) {

                for (j in (age-63):ncol(etat)) {

                    n0[j] <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")))

                };

                esp_vie_gen[age-64,2] <- 1 + sum(n0) / s0;

            } else {

                n0 <- NA;

            }

        };

###	Overall life-expectancy on exposed peoples

        if (age < 101) {

            n0 <- vector(length = ncol(etat));

            s0 <- sum(etat[,age-64] %in% etat_vivant_exp)

            if (s0 != 0) {

                for (j in (age-63):ncol(etat)) {

                    n0[j] <- sum(etat[which(etat[,age-64]%in%etat_vivant_exp),j]%in%etat_vivant_exp)

                };

                esp_vie_gen_conso[age-64,2] <- 1 + sum(n0) / s0;

            } else {

                n0 <- NA;

            }

        };

### Overall life-expectancy on non exposed peoples

        if (age < 101) {

            n0 <- vector(length = ncol(etat));

            s0 <- sum(etat[,age-64]%in%etat_vivant_n_exp)

            if (s0 != 0) {

                for (j in (age-63):ncol(etat)) {

                    n0[j] <- sum(etat[which(etat[,age-64]%in%etat_vivant_n_exp),j]%in%etat_vivant_n_exp)

                };

                esp_vie_gen_nonconso[age-64,2] <- 1 + sum(n0) / s0;

            } else {

                n0 <- NA;

            }

        };

###	Life-expectancy without disease

        if (age < 101) {

            n0 <- vector(length = ncol(etat));

            s0 <- sum(etat[which(etat[,age-64]%in%(c("00","01","02","03","04","05","06","07"))),age-64]%in%(c("00","01","02","03","04","05","06","07")))

            if (s0 != 0) {

                for (j in (age-63):ncol(etat)) {

                    n0[j] <- sum(etat[which(etat[,age-64]%in%(c("00","01","02","03","04","05","06","07"))),j]%in%(c("00","01","02","03","04","05","06","07")))

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

            s0 <- sum(etat[which(etat[,age-64]%in%(c("00","01","02","03","04","05","06","07"))),age-64]%in%(c("00","01","02","03","04","05","06","07")));

            if (s0 != 0) {

                for (j in (age-63):ncol(etat)) {

                    n0[j] <- sum(etat[which(etat[,age-64]%in%(c("00","01","02","03","04","05","06","07"))),j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")))

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

            d0 <- sum(etat[,age-64]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17"))) + 0.5 * sum((etat[,age-64]%in%(c("20","21","22","23","24","25","26","27"))) & (etat[,age-65]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17"))));

            s0 <- sum(etat[,1]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

            s1 <- sum((etat[,age-64]%in%(c("10","11","12","13","14","15","16","17"))) & (etat[,age-65]%in%(c("10","11","12","13","14","15","16","17")))) + 0.5 * sum((etat[,age-64]%in%(c("10","11","12","13","14","15","16","17"))) & (etat[,age-65]%in%(c("00","01","02","03","04","05","06","07")))) + 0.5 * sum((etat[,age-64]%in%(c("20","21","22","23","24","25","26","27"))) & (etat[,age-65]%in%(c("10","11","12","13","14","15","16","17"))));

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

            d0 <- sum(etat[,age-65]%in%etat_vivant)

            s0 <- sum(etat[,1]%in%etat_vivant)
            
            s1 <- sum(etat[,age-64]%in%etat_vivant)

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

                prb_dem[i,2] <- sum(etat[which(etat[,i-1]%in%(c("00","01","02","03","04","05","06","07"))),i]%in%(c("10","11","12","13","14","15","16","17")))

            }

        };

### Average age at disease onset

        if (age == 65) {

            for (i in (age-63):nrow(age_dem)) {

                age_dem[i,2] <- sum(etat[which(etat[,i-1]%in%(c("00","01","02","03","04","05","06","07"))),i]%in%(c("10","11","12","13","14","15","16","17")))

            }

        };
        
### cas incidents l'annee de projection
        if((age>65) & (age<100)) age_cas_incid[age-65+1,2] <- sum(etat[which(etat[,age-65]%in%(c("00","01","02","03","04","05","06","07"))),age-65+1]%in%(c("10","11","12","13","14","15","16","17"))) # on ajoute le nb de sujets dements incidents (en 2040) de cette generation  
        
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


##############################
########## ESSAIS ############
##############################


### Dementia incidence

        if (age >= 65) {

            for (j in 1:nrow(incid_demence)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")));

                s1 <- sum((etat[,j+1]%in%(c("10","11","12","13","14","15","16","17"))) & (etat[,j]%in%(c("00","01","02","03","04","05","06","07"))));

                if (s0 != 0) {

                    incid_demence[j,age-63] <- s1/s0;

                };

            }

        };

### Mortality incidence

        if (age >= 65) {

            for (j in 1:nrow(incid_mort_NoD)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07")));

                s1 <- sum((etat[,j+1]%in%(c("20","21","22","23","24","25","26","27"))) & (etat[,j]%in%(c("00","01","02","03","04","05","06","07"))));

                if (s0 != 0) {

                    incid_mort_NoD[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(incid_mort_Dem)) {

                s0 <- sum(etat[,j]%in%(c("10","11","12","13","14","15","16","17")));

                s1 <- sum((etat[,j+1]%in%(c("20","21","22","23","24","25","26","27"))) & (etat[,j]%in%(c("10","11","12","13","14","15","16","17"))));

                if (s0 != 0) {

                    incid_mort_Dem[j,age-63] <- s1/s0;

                };

            }

        };

### Non exposed prevalence

        if (age >= 65) {

            for (j in 1:nrow(prevalence_noe)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("00","10")));

                if (s0 != 0) {

                    prevalence_noe[j,age-63] <- s1/s0;

                };

            }

        };

### Hypertension prevalence

        if (age >= 65) {

            for (j in 1:nrow(prevalence_hta)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("01","04","05","07","11","14","15","17")));

                if (s0 != 0) {

                    prevalence_hta[j,age-63] <- s1/s0;

                };

            }

        };

### Physical inactivity prevalence

        if (age >= 65) {

            for (j in 1:nrow(prevalence_ina)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("03","05","06","07","13","15","16","17")));

                if (s0 != 0) {

                    prevalence_ina[j,age-63] <- s1/s0;

                };

            }

        };

### Diabete prevalence

        if (age >= 65) {

            for (j in 1:nrow(prevalence_dia)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("02","04","06","07","12","14","16","17")));

                if (s0 != 0) {

                    prevalence_dia[j,age-63] <- s1/s0;

                };

            }

        };

### Expositions prevalence

        if (age >= 65) {

            for (j in 1:nrow(prevalence_exp0)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("00","10")));

                if (s0 != 0) {

                    prevalence_exp0[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prevalence_exp1)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("01","11")));

                if (s0 != 0) {

                    prevalence_exp1[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prevalence_exp2)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("02","12")));

                if (s0 != 0) {

                    prevalence_exp2[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prevalence_exp3)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("03","13")));

                if (s0 != 0) {

                    prevalence_exp3[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prevalence_exp4)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("04","14")));

                if (s0 != 0) {

                    prevalence_exp4[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prevalence_exp5)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("05","15")));

                if (s0 != 0) {

                    prevalence_exp5[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prevalence_exp6)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("06","16")));

                if (s0 != 0) {

                    prevalence_exp6[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prevalence_exp7)) {

                s0 <- sum(etat[,j]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j]%in%(c("07","17")));

                if (s0 != 0) {

                    prevalence_exp7[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prop_dem_global)) {

                s0 <- sum(etat[,j+1]%in%(c("00","01","02","03","04","05","06","07","10","11","12","13","14","15","16","17")));

                s1 <- sum(etat[,j+1]%in%(c("10","11","12","13","14","15","16","17")));

                if (s0 != 0) {

                    prop_dem_global[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prop_dem_diabet)) {

                s0 <- sum(etat[,j+1]%in%(c("02","04","06","07","12","14","16","17")));

                s1 <- sum(etat[,j+1]%in%(c("12","14","16","17")));

                if (s0 != 0) {

                    prop_dem_diabet[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prop_dem_hypert)) {

                s0 <- sum(etat[,j+1]%in%(c("01","04","05","07","11","14","15","17")));

                s1 <- sum(etat[,j+1]%in%(c("11","14","15","17")));

                if (s0 != 0) {

                    prop_dem_hypert[j,age-63] <- s1/s0;

                };

            }

        };

        if (age >= 65) {

            for (j in 1:nrow(prop_dem_inacti)) {

                s0 <- sum(etat[,j+1]%in%(c("03","05","06","07","13","15","16","17")));

                s1 <- sum(etat[,j+1]%in%(c("13","15","16","17")));

                if (s0 != 0) {

                    prop_dem_inacti[j,age-63] <- s1/s0;

                };

            }

        };

    }


### Computation of health indicators :

### Overall life-expectancy
    life_expectancy <- esp_vie_gen

### Overall life-expectancy on exposed peoples
    life_expectancy_exp <- esp_vie_gen_conso
    
### Overall life-expectancy on non exposed peoples
    life_expectancy_n_exp <- esp_vie_gen_nonconso

### Life-expectancy without disease
    life_expectancy_w_dis <- esp_vie_sans_mal
    
### Life-expectancy without disease on exposed peoples
    life_expectancy_w_dis_exp <-esp_vie_sans_mal_conso

### Life-expectancy without disease on non exposed peoples
    life_expectancy_w_dis_n_exp <- esp_vie_sans_mal_nonconso

### Life-expectancy for diseased subject
    life_expectancy_dis <- esp_vie_mal

### Life-expectancy for diseased subject on exposed peoples
    life_expectancy_dis_exp <- esp_vie_mal_conso

### Life-expectancy for diseased subject on non exposed peoples
    life_expectancy_dis_n_exp <- esp_vie_mal_nonconso

### Life-expectancy for non diseased subject
    life_expectancy_n_dis <- esp_vie_non_mal
    
### Life-expectancy for non diseased subject on exposed peoples
    life_expectancy_n_dis_exp <- esp_vie_non_mal_conso

### Life-expectancy for non diseased subject on non exposed peoples
    life_expectancy_n_dis_n_exp <- esp_vie_non_mal_nonconso

### Prevalence of disease
    prevalence <- prevalence[-1,]
    prev <- sum(prevalence[,-1])
    number_prevalence <- prev[1]

### Prevalence of disease by age
    number_prev_age <- prevalence
    
### Prevalence rate of disease by age
    taux_prevalence <- taux_prevalence[-1,]
    prev_rate_disease_age <- taux_prevalence

### Survival
    surv <- sum(survie[,-1])
    number_survival <- surv[1]
    
### Survival by age
    number_survival_age <- survie
    
### Survival rate
    taux_survivants <- taux_survivants[-1,]
    survival_rate <- taux_survivants

### Global prevalence rate of disease
    taux_prev_demence <- matrix(c(0),
                                nrow=1,
                                ncol=length(prev)+1,
                                byrow = T);

    taux_prev_demence[1,1] <- t;

    for (i in 2:ncol(taux_prev_demence)) {

        taux_prev_demence[1,i] <- prev[i-1] / surv[i-1];

    };

    prev_rate_disease <- taux_prev_demence[,2]

### Mean number of years spent with disease
    number_years_disease <- nb_moy_dem

### Mean number of years spent with disease on exposed peoples
    number_years_disease_exp <- nb_moy_dem_conso
    
### Mean number of years spent with disease on non exposed peoples
    number_years_disease_n_exp <- nb_moy_dem_nonconso

### Life-long probability of disease
    prb_dem <- prb_dem[-1,]

    prb_demence <- matrix(c(0),
                          nrow=1,
                          ncol=ncol(prb_dem),
                          byrow = T);

    prb_demence[1,1] <- t;

    for (i in 2:ncol(prb_dem)) {

        prb_demence[i] <- sum(prb_dem[,i]) / nrow(etat);

    };

    ll_prob_disease <- prb_demence[1,2]

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

    average_age_disease <-  age_demence[1,2]
    
### average age of disease for incident cases     
    average_age_incident <- sum(age_cas_incid[,1]*age_cas_incid[,2])/sum(age_cas_incid[,2])

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

    number_years_exp <- matrix(c(as.integer(t), moyenne_conso[,2]),nrow=1)
    
### Number of exposed peoples at least one time
    exposition <- prevalence_conso
    
### Mortality rate

    quotient_mortalite <- quotient_mortalite[-1,]

    mortality_rate <- quotient_mortalite
    
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

    list_prevalence_disease <- list(number_prev_age, number_prevalence, prev_rate_disease_age, prev_rate_disease)
    names(list_prevalence_disease) <- c("number_prev_age", "number_prevalence", "prev_rate_disease_age", "prev_rate_disease")

### Survival

    list_survival <- list(number_survival_age, number_survival, survival_rate)
    names(list_survival) <- c("number_survival_age", "number_survival", "survival_rate")

### Mean number of years spent with disease

    list_number_years_disease <- list(number_years_disease, number_years_disease_exp, number_years_disease_n_exp)
    names(list_number_years_disease) <- c("number_years_disease", "number_years_disease_exp", "number_years_disease_n_exp")

### Summary of all iterations

    list_summary_iterations <- list(esp_vie_sans_mal, prevalence, taux_prevalence, prb_dem, age_dem, nb_moy_dem, prev, surv)
    names(list_summary_iterations) <- c("esp_vie_sans_mal", "prevalence", "taux_prevalence", "prb_dem", "age_dem", "nb_moy_dem", "prev", "surv")

### Output list

                                        # With exposition

                                        #HI_output <- list(list_overall_LE, list_LE_without_disease, list_LE_diseased, list_LE_non_diseased, list_prevalence_disease, list_survival, list_number_years_disease, ll_prob_disease, average_age_disease, number_years_exp, exposition, mortality_rate);
                                        #names(HI_output) <- c(list_overall_LE", "list_LE_without_disease", "list_LE_diseased", "list_LE_non_diseased", "list_prevalence_disease", "list_survival", "list_number_years_disease", "ll_prob_disease", "average_age_disease", "number_years_exp", "exposition", "mortality_rate");

                                        # Without exposition

    HI_output <- list(list_overall_LE, list_LE_without_disease, list_LE_diseased, list_LE_non_diseased, list_prevalence_disease, list_survival, list_number_years_disease, list_summary_iterations, ll_prob_disease, average_age_disease,average_age_incident, number_years_exp, mortality_rate);
    names(HI_output) <- c("list_overall_LE", "list_LE_without_disease", "list_LE_diseased", "list_LE_non_diseased", "list_prevalence_disease", "list_survival", "list_number_years_disease", "list_summary_iterations", "ll_prob_disease", "average_age_disease","average_age_incident,", "number_years_exp", "mortality_rate");

    return(HI_output)

}
