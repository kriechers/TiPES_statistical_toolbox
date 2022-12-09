# KS transition detection

The directory **KS_detection** is comprised of the following files:

    **KS_detection.m**: this is a function for the detection of transitions based on the Kolmogorov-Smirnov statistic within a given time series. For technical guidance, please see the file itself. For the details of the method, please see Bagniewski, W., Ghil, M. & Rousseau, D. D. Automatic detection of abrupt transitions in paleoclimate records. Chaos 31, (2021).

    **KS_example.m**: This file examplary runs the function defined in **KS_detection.m** on the data contained in **NGRIP.csv** and plots the results. Doing so, it uses the functions defined in **vline.m** and **findgigs.m**.

    **findgigs.m**: a function used by **KS_example.m** for plotting.

    **vline.m**: a function used by **KS_example.m** for plotting.

    **NGRIP.csv**: examplary data on which **KS_example.m** runs the **KS_detection.m** function. For details on the data please see: Rasmussen, S. O. et al. A stratigraphic framework for abrupt climatic changes during the Last Glacial period based on three synchronized Greenland ice-core records: Refining and extending the INTIMATE event stratigraphy. Quat. Sci. Rev. 106, 14–28 (2014).
