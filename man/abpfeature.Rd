% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abpfeature.R
\name{abpfeature}
\alias{abpfeature}
\title{ABPFEATURE}
\usage{
abpfeature(abp, OnsetTimes)
}
\arguments{
\item{abp}{Arterial blood pressure waveform}

\item{OnsetTimes}{Onset times (in indices)}
}
\value{
List of ABP features:
1. Time of systole (samples)
2. Systolic BP (mmHg)
3. Time of diastole (samples)
4. Diastolic BP (mmHg)
5. Pulse Pressure (mmHg)
6. Mean ABP (mmHg)
7. Beat Period (samples)
8. mean_dyneg TODO
9. End of systole time, 0.3*sqrt(RR) method
10. Area under systole, 0.3*sqrt(RR) method TODO
11. End of systole time, 1st min-slope method TODO
12. Area under systole, 1st min-slope method TODO
}
\description{
ABP waveform feature extractor.
Original MATLAB code by James Sun on Nov 19, 2005.
Ported to R by Ran Liu on June 24, 2018.
}
