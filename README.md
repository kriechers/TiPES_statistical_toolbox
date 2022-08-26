# Statistical Toolbox for the analysis of (past) abrupt climate change

All code comprised in this repository was developed in the
context of the TiPES ('Tipping points in the Earth System')
project which was funded by the European Union’s Horizon 2020
research and innovation program under grant agreement no. 820970.

This toolbox aims at facilitating the investigation of past
abrupt climate change as recorded in climate proxy data. All
contributing authors were involved in the TiPES project and
accordingly were involved in the assessment of paleoclimate proxy
data with respect to abrupt shifts.

-----------------------------------------------------------------

## Content

### transition_characterization (python) MCMC-based

A baysian ramp fit to characterize previously detected and
isolated transitions. Applicable to abrupt changes of the mean of
a time series. The method returns uncertainty sensitive estimates
of the transition's starting point, duration, amplitude and
initial value.

### TransitionDetection (julia)

A transition detection algorithm based on the comparison of the
local mean calculated over two adjacent sliding running windows.
Applicable to univariate time series. Returns the time points of
potential abrupt shifts in the local mean of the time series.

### TransitionDetectionSisalv2 (julia)

Readily implemented of the above TransitionDetection algorithm to
the Sisalv2 data base that comprises an extensive collection of
speleothem records. The implementation automatically downloads
the Sisalv2 data base and runs the TransitionDetection on a
selection of records.

### KS_detection (Matlab)

Transition detection algorithm based on comparing the KS
statistic of two adjacent running windows. Returns the time
points of detected transitions.


### RQA_detected (Matlab)

Transition detection algorithm based recurrence plot analysis.
Returns the time points of detected transitions.

### NGRIP_datinguncertainty (R)

Algorithm to efficiently sample realizations of the GICC05
chronology for the NGRIP ice core. The method relies on a
comprehensive statistical model for the respective dating
uncertainties. 

-----------------------------------------------------------

## Related Publications

* Comas-Bru, L., Rehfeld, K., Roesch, C., Amirnezhad-Mozhdehi, S.,
Harrison, S. P., Atsawawaranunt, K., Ahmad, S. M., Brahim, Y. A.,
Baker, A., Bosomworth, M., Breitenbach, S. F. M., Burstyn, Y.,
Columbu, A., Deininger, M., Demény, A., Dixon, B., Fohlmeister,
J., Hatvani, I. G., Hu, J., Kaushal, N., Kern, Z., Labuhn, I.,
Lechleitner, F. A., Lorrey, A., Martrat, B., Novello, V. F.,
Oster, J., Pérez-Mejías, C., Scholz, D., Scroxton, N., Sinha, N.,
Ward, B. M., Warken, S., Zhang, H., and SISAL Working Group
members: SISALv2: a comprehensive speleothem isotope database
with multiple age–depth models, Earth Syst. Sci. Data, 12,
2579–2606, https://doi.org/10.5194/essd-12-2579-2020, 2020.

* Riechers, K. & Boers, N. Significance of uncertain phasing
between the onsets of stadial-interstadial transitions in
different Greenland ice core proxies. Clim. Past 17, 1751–1775,
https://doi.org/10.5194/cp-17-1751-2021, 2021.

* Myrvoll-Nilsen, E., Riechers, K., Rypdal, M. W. & Boers, N.
Comprehensive uncertainty estimation of the timing of Greenland
warmings in the Greenland ice core records. Clim. Past 18,
1275–1294, https://doi.org/10.5194/cp-18-1275-2022, 2022.

* Bagniewski, W., Ghil, M. & Rousseau, D. D. Automatic detection of
abrupt transitions in paleoclimate records. Chaos 31,
https://doi.org/10.1063/5.0062543, 2021.


-----------------------------------------------------------------

## Main Authors:

* Witold Bagniewski
* Eirik Myrvoll-Nilsen
* Jens Fohlmeister
* Keno Riechers
