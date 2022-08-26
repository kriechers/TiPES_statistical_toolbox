
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Quantification of NGRIP dating uncertainties

This directory is comprised of this readme.md (and its readme.Rmd
counterpart), 12 R files and a subdirectory which includes in total 9
data files of various formats. These files can be used to reproduce the
results from Myrvoll-Nilsen et al. (2022). Moreover, this directory
include code to sample chronologies that have been synchronized using
tie-points using a method on which we are currently writing a paper.

## R files

-   *main.R*: Primary file for end-users. Sets options and runs other
    functions.
-   *prepare.R*: Formats input and data
-   *modelfitter.R*: Function that fits the model to layer increments
    using INLA
-   *chronology_simulation.R*: Simulates an ensemble of chronologies.
-   *biased_chronologies.R*: Simulates an ensemble of chronologies where
    a possible stochastic bias has been added.
-   *linrampfitter.R*: Fits the linear ramp model of Erhardt et
    al. (2019) using INLA.
-   *event_depth_to_age.R*: Performs Monte Carlo sampling for estimating
    dating uncertainty of transition onsets given a linear ramp fit on
    onset depth.
-   *rgeneric_model.R*: Specifies the rgeneric model used in
    linrampfitter.R
-   *summary_results.R*: Prints a summary of our fit and sampling to
    screen
-   *plot_results.R*: Plots key figures to screen
-   *mainfunctions.R*: Main wrapper function for fitting model and
    simulating chronologies, also includes the function which produces
    synchronized chronology simulations.
-   *helpfulfunctions.R*: Other functions that are used by other
    functions

## Data sets

The data sets are included in subdirectory *datasets_used*

-   *Adolphi18_pdf_tie_1-3.txt*: Tie-point distribution for the 3 first
    tie-points of Adolphi et al. (2018).
-   *Adolphi18_pdf_tie_4\_update_Muscheler2020.txt*: Tie-point
    distribution for the fourth tie-point of Adolphi et al. (2018),
    updated in Muscheler et al. (2020).
-   *Adolphi18_pdf_tie_5\_update_Muscheler2020.txt*: Tie-point
    distribution for the fifth tie-point of Adolphi et al. (2018),
    updated in Muscheler et al. (2020).
-   *Buizert_onsets.xlsx*: Excel file for onset dates given by Buizert
    et al. (2015).
-   *Capron_onsets.xls*: Excel file for onset dates given by Capron et
    al. (2021).
-   *event_intervals.txt*: Contains transition windows suggested by
    Myrvoll-Nilsen et al. (2022).
-   *GISevents.ods*: Contains transition onsets taken from table 2 in
    Rasmussen et al. (2014).
-   *NGRIP_d18O_and_dust_5cm.xls*: The 5cm depth resolution (0-60ka bp)
    downloaded from <https://www.iceandclimate.nbi.ku.dk/data/>.
    References: NGRIP members (2004), Gkinis et al. (2014) and Ruth et
    al. (2003).
-   *Rasmussen_et_al_2014_QSR_Table_2.xlsx*: Also contains transition
    onsets from table 2 in Rasmussen et al. (2014). Formatted slightly
    different from *GISevents.ods*.

## The model

The base model is detailed in Myrvoll-Nilsen et al. (2022), and is a
Bayesian regression model for the layer increments. The layer increase
per unit of depth is expressed as the sum of one physical component
depending on temperature (through proxy) and the depth of the core, and
one stochastic component which we assume to be an AR(1) process. There
is also limited support for independent and AR(2) processes.

The model is fitted using R-INLA (Rue et al. 2009 and Rue et al. 2017)
which takes advantage of the computational properties of AR processes
and gives full Bayesian inference of all model variables and parameters.

Using Monte Carlo sampling we can produce an ensemble of plausible
chronologies from the fitted layer-increment model. Both the fitting and
sampling procedures are performed in linear time, which is particularly
useful for large datasets such as NGRIP.

We can also generate an ensemble of chronologies sampled from a
synchronized time scale. The theory will be detailed in an upcoming
paper, but we have included this feature in this repository.

It is possible to apply our methodology to estimate the dating
uncertainty of abrupt warming transitions in the ice core. This is done
in two steps. First, we apply the linear ramp model designed by Erhardt
et al. (2019) to the d18O NGRIP data for a window enclosing the
transition, giving posterior distributions for the NGRIP depth
associated with the onset of the transition. Using INLA instead of
Markov chain Monte Carlo sampling will substantially speed up this
process. Then, for each sample drawn from the onset depth posterior
distribution we sample from the associated age-depth uncertainty.

## How to use the model

For the end-user only the file *main.R* is of relevance. Of course, the
user is free and encouraged to experiment and make modifications of all
code in this directory. Within *main.R* there is a number of variables
which determine the settings and what should be computed. Change these
to suit your needs and then run *source(main.R)* to compute everything.

The unsynchronized chronology ensemble can be obtained by
*object$simulation$age* with summary statistics available from
*object$simulation$summary*. The synchronized chronology ensemble will
instead be stored in *object_sync$simulation$age*, with summary
statistics available in *object_sync$simulation$summary*. Summaries of
can be printed to screen using *summary_results(object)* or
*summary_results(object_sync)*, and key figures can be plotted using
*plot_results(object)* or *plot_results(object_sync)*.

Summary statistics of the onset depth and onset age of the abrupt
warming transitions are stored in *depthresults* and *ageresults*,
respectively.

## Requirements

For this code to work you must have installed:

-   *R* (version 4.2.0)
-   *INLA* (version 22.06.03, and its dependencies) This package is not
    available on CRAN, but can be downloaded by running
    *install.packages(“INLA”,repos=c(getOption(“repos”),INLA=“<https://inla.r-inla-download.org/R/testing>”),
    dep=TRUE)*
-   *ggplot2* (version 3.3.6)
-   *zoo* (version 1.8-10)
-   *matrixStats* (version 0.62.0)
-   *numDeriv* (version 2016.8-1.1)
-   *readODS* (version 1.7.0)
-   *readxl* (version 1.4.0)
-   *stringr* (version 1.4.0)

## References

Adolphi, F., Bronk Ramsey, C., Erhardt, T., Edwards, R. L., Cheng, H.,
Turney, C. S. M., Cooper, A., Svensson, A., Rasmussen, S. O., Fischer,
H., and Muscheler, R.: Connecting the Greenland ice-core and U∕Th
timescales via cosmogenic radionuclides: testing the synchroneity of
Dansgaard–Oeschger events, Clim. Past, 14, 1755–1781,
<https://doi.org/10.5194/cp-14-1755-2018>, 2018.

Buizert, C., Cuffey, K. M., Severinghaus, J. P., Baggenstos, D., Fudge,
T. J., Steig, E. J., Markle, B. R., Winstrup, M., Rhodes, R. H., Brook,
E. J., Sowers, T. A., Clow, G. D., Cheng, H., Edwards, R. L., Sigl, M.,
McConnell, J. R., and Taylor, K. C.: The WAIS Divide deep ice core
WD2014 chronology – Part 1: Methane synchronization (68-31 ka BP) and
the gas age-ice age difference, Climate of the Past, 11,520 153–173,
<https://doi.org/10.5194/cp-11-153-2015>, 2015.

Capron, E., Rasmussen, S. O., Popp, T. J., Erhardt, T., Fischer, H.,
Landais, A., Pedro, J. B., Vettoretti, G., Grinsted, A., Gkinis, V.,
Vaughn, B., Svensson, A., Vinther, B. M., and White, J. W.: The anatomy
of past abrupt warmings recorded in Greenland ice, Nature
Communications, 12, <https://doi.org/10.1038/s41467-021-22241-w>, 2021.

Erhardt, T., Capron, E., Rasmussen, S.O., Schüpbach, S., Bigler, M.,
Adolphi, F., and Fischer, H.: Decadal-scale progression of the onset of
Dansgaard-Oeschger warming events, Climate of the Past, 15, 811-825,
<https://doi.org/10.5194/cp-15-811-2019>, 2019.

Gkinis, V., Simonsen, S. B., Buchardt, S. L., White, J. W., and Vinther,
B. M.: Water isotope diffusion rates from the NorthGRIP ice core for the
last 16,000 years - Glaciological and paleoclimatic implications, Earth
and Planetary Science Letters, 405, 132–141,
<https://doi.org/10.1016/j.epsl.2014.08.022>, 2014.

Muscheler, R., Adolphi, F., Heaton, T., Bronk Ramsey, C., Svensson, A.,
Van der Plicht, J., & Reimer, P. (2020). Testing and Improving the
IntCal20 Calibration Curve with Independent Records. Radiocarbon, 62(4),
1079-1094. <doi:10.1017/RDC.2020.54>

Myrvoll-Nilsen, E., Riechers, K., Rypdal, M. & Boers, N.: Comprehensive
uncertainty estimation of the timing of Greenland warmings of the
Greenland Ice core records. Climate of the Past, 18, 1275-1294.
<https://doi.org/10.5194/cp-18-1275-2022>, 2022

North Greenland Ice Core Project members: High-resolution record of
Northern Hemisphere climate extending into the last interglacial period,
Nature, 431, 147–151, <https://doi.org/10.1038/nature02805>, 2004.

Rasmussen, S. O., Bigler, M., Blockley, S. P., Blunier, T., Buchardt, S.
L., Clausen, H. B., Cvijanovic, I., Dahl-Jensen, D., Johnsen, S. J.,
Fischer, H., Gkinis, V., Guillevic, M., Hoek, W. Z., Lowe, J. J., Pedro,
J. B., Popp, T., Seierstad, I. K., Steffensen, J. P., Svensson, A. M.,
Vallelonga, P., Vinther, B. M., Walker, M. J., Wheatley, J. J., and
Winstrup, M.: A stratigraphic framework for abrupt climatic changes
during the Last Glacial period based on three synchronized Greenland
ice-core records: refining and extending the INTIMATE event
stratigraphy, Quaternary Science Reviews, 106, 14–28,
<https://doi.org/https://doi.org/10.1016/j.quascirev.2014.09.007>,
dating, Synthesis, and Interpretation of Palaeoclimatic Records and
Model-data Integration: Advances of the INTIMATE project(INTegration of
Ice core, Marine and TErrestrial records, COST Action ES0907), 2014.

Rue, H., Martino, S., and Chopin, N.: Approximate Bayesian inference for
latent Gaussian models using integrated nested Laplace approximations
(with discussion)., J R Stat Soc Series B, 71, 319–392,
<https://doi.org/10.1111/j.1467-9868.2008.00700.x>, 2009.

Rue, H., Riebler, A., Sørbye, S. H., Illian, J. B., Simpson, D. P., and
Lindgren, F. K.: Bayesian Computing with INLA: A Review, Annu.
Rev. Stat. Appl., 4, 395–421,
<https://doi.org/10.1146/annurev-statistics-060116-054045>, 2017.

Ruth, U., Wagenbach, D., Steffensen, J. P., and Bigler, M.: Continuous
record of microparticle concentration and size distribution in the
central Greenland NGRIP ice core during the last glacial period, J.
Geophys. Res., 108, 4098, <a href="https://doi:10.1029/2002JD002376"
class="uri">https://doi:10.1029/2002JD002376</a>, D3, 2003.
