\docType{data}
\name{data_MCSPCD_3FDR}
\title{Input data for estimHI_3FDR}
\format{
  Several matrices.
  Matrices of 80 lines and 3 columns contain an 'age' column ranging
  from 66 to 105, a 'gender' column being either "M" for men or "W" for
  women, and another column containing the requested value.
  Matrices of 80 lines and 133 columns contain the same 'age' and
  'gender' mentionned above, and one column for each generation from 1950
  to  2080.
  In the following, RR stands for relative risk.
\describe{
  \item{a01_constant_values}{matrix of 80 lines and 133 columns
    containing the (initial) dementia indicence for each age from 66 to 105 years and for
    each generation, by gender}
  \item{a01_75pct_values}{matrix of 80 lines and 133 columns,
   same format as a01_constant_values. Alternative incidence
   representing a reduction of 25 \% of the initial scenario.}
 \item{a02_constant_values}{matrix of 80 lines and 133 columns
   containing the mortality of people without dementia, and
   corresponding to the initial dementia incidence.}
 \item{a02_75pct_values}{matrix of 80 lines and 133 columns
   containing the mortality of people without dementia, and
   corresponding to the alternative dementia incidence.}
 \item{rr_DvsND_values}{matrix of 80 lines and 3 colums containing the
   relative risk of death for people with dementia compared to people
   without dementia, at each age from 66 to 105 years, by gender.}
 \item{theta01_HTA}{RR of dementia associated with high blood pressure
   (80 lines and 3 columns)}
 \item{theta01_DIAB}{RR of dementia associated with diabetes (80 lines and 3 columns)}
 \item{theta01_INACT}{RR of dementia associated with physical inactivity (80 lines and 3 columns)}
 \item{theta02_HTA}{RR of death for people without dementia associated
   with high blood pressure (80 lines and 3 columns)}
 \item{theta02_DIAB}{RR of death for people without dementia associated
   with diabetes (80 lines and 3 columns)}
 \item{theta02_INACT}{RR of death for people without dementia associated
   with physical inactivity (80 lines and 3 columns)}
 \item{theta12_HTA}{RR of death for people with dementia associated
   with high blood pressure (80 lines and 3 columns)}
 \item{theta12_DIAB}{RR of death for people with dementia associated
   with diabetes (80 lines and 3 columns)}
 \item{theta12_INACT}{RR of death for people with dementia associated
   with physical inactivity (80 lines and 3 columns)}
 \item{prev_0_values}{prevalence at each age from 65 to 105 of the status
   HTA-/Diab-/INACT- (ie, without high blood pressure, without diabetes
   and without physical inactivity.). Matrix with 82 lines and 3 columns.}
 \item{prev_1_values}{prevalence at 65 years old of the status
   HTA+/Diab-/INACT- (82 lines and 3 columns)}
 \item{prev_2_values}{prevalence at 65 years old of the status
   HTA-/Diab+/INACT- (82 lines and 3 columns)}
 \item{prev_3_values}{prevalence at 65 years old of the status
   HTA-/Diab-/INACT+ (82 lines and 3 columns)}
 \item{prev_4_values}{prevalence at 65 years old of the status
   HTA+/Diab+/INACT- (82 lines and 3 columns)}
 \item{prev_5_values}{prevalence at 65 years old of the status
   HTA+/Diab-/INACT+ (82 lines and 3 columns)}
 \item{prev_6_values}{prevalence at 65 years old of the status
   HTA-/Diab+/INACT+ (82 lines and 3 columns)}
 \item{prev_7_values}{prevalence at 65 years old of the status
   HTA+/Diab+/INACT+ (82 lines and 3 columns)}
 \item{incid_glob}{incidence of diabetes (80 lines and 3 columns)}
 \item{theta1_F}{RR of diabetes for women associated with high blood
   pressure (single value)}
 \item{theta1_H}{RR of diabetes for men associated with high blood
   pressure (single value)}
 \item{theta3_F}{RR of diabetes for women associated with physical
   inactivitiy (single value)}
 \item{theta3_H}{RR of diabetes for men associated with physical
   inactivitiy (single value)}
 \item{pop}{number of people alive at each year from 1950 to 2070, by
   gender (242 lines and 3 columns).}
}
}
\description{
The data were obtained from different sources and were designed for the
franch population.
}
\keyword{datasets}
