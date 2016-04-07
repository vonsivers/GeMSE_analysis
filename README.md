GeMSE_analysis
======
Determines the activity for several isotopes by comparing a sample to a background spectrum
## Usage
```
GeMSE_analysis <parameters_activity_calculation.txt>
```
## Output
The results folder contains the following files:
* sample_activities_summary.txt: The header contains some basic information like time and date, sample name, spectrum files, measurement times etc. For ever isotope the activity and the Bayesfactor are given. If the Bayesfactor is smaller than a specified value the activity is the mode of the marginalized posterior distribution. The uncertainties correspond to the central 68% interval. If the Bayesfactor is larger than a certain value the quantile at a specified confidence level is given.
The following files are written for every isotope both for the signal+background and the background-only fit. The latter are indicated by the subscript „_bckonly“.
* sample_isotope_distributions.pdf: Shows the marginalized posterior distribution for every fit parameter and their correlations.
* sample_isotope_log.txt: Logfile of the BAT fitting algorithm. Any errors that occur during fitting are summarized here.
* sample_isotope_output.root: ROOT file that contains all marginalized distributions and correlations.
* sample_isotope_results.txt: Summary of the results of the fit performed by BAT. 
* sample_stack_bck_<peak_i>.pdf: Shows the fit to the background spectrum for peak nr. i. The individual parameters („signal“, „sample_const“, „bck_const“, „bck_gauss“) are plotted as stacked histograms.
* sample_isotope_stack_sample_peak_i.pdf: Shows the fit to the sample spectrum for peak nr. i.

## parameters_activity_calculation.txt
The file must contain the following information (see „example_parameters_global.txt“):
* sample name: name of the sample used to name the output files. Avoid spaces in the name, use underscore („_“) instead.
* accuracy of MCMC (low/medium/high): accuracy of the Markov Chain Monte-Carlo (MCMC), can be „low“, „medium“ or „high“. „low“ takes a few minutes, „medium“ a few hours, „high“ was never tested
* threshold on Bayes Factor for signal detection: upper limit on Bayes Factor that counts as signal detection. Suggested values are 0.33 (positive evidence) or 0.05 (strong evidence).
* CL for activity limit: confidence level for upper limit on the activity. Suggested value is 0.95.
* fractional uncertainty of efficiencies: systematic error of the simulated efficiencies. Expressed as fractional uncertainty
* isotope parameters folder: folder name of the isotope parameters files (see below)
* sample spectrum file name: name of the ROOT file with the sample spectrum
* background spectrum file name: name of the ROOT file with the background spectrum
* efficiency file name: name of the ROOT file with the simulated efficiencies
* resolution file name: name of the ROOT file with the energy resolution
* results folder: name of the folder where the results should be written
* isotopes to analyze: names of all isotopes to analyze, for every isotope a parameters file must exist in the folder specified above

## parameters_isotope.txt
This file contains the necessary information for every isotope:

* number of peaks: number of peaks that should be used to calculate the activity
* peak energies (keV): energies of the peaks in keV
* background efficiency: accounts for a possible reduction of the peak in the sample spectrum compared to the background spectrum. This factor is usually set to 1, however sometimes the background in the sample spectrum is lower than in the background spectrum. This is for example the case when the sample fills out the whole cavity so there is less background from radon. If it is known how much this reduces the background (e.g. by a factor of 2) the background efficiency can be set to the corresponding value (e.g. 0.5)
* lower fit range (keV): the lower limit of the fit range in keV. The value should correspond to mean-5*sigma
* upper fit range (keV): the upper limit of the fit range in keV. Should correspond to mean+5*sigma
* lower limit of fit parameter „signal“: the lower limit of the (integrated) number of signal counts. The value should always be set to 0.
* upper limit of fit parameter „signal“: upper limit of the (integrated) number of signal counts. Can be increased for samples with very high activity.
* lower limit of fit parameter „sample_const“: the lower limit of the number of counts (per bin) of the constant background in the sample spectrum. The value should always be set to 0. 
* upper limit of fit parameter „sample_const“: upper limit of the number of counts (per bin) of the constant background in the sample spectrum. Can be increased when the background is very high, e.g. due to Compton background from a pronounced line.
* lower limit of fit parameter „bck_const“: the lower limit of the number of counts (per bin) of the constant background in the background spectrum. The value should always be set to 0. 
* upper limit of fit parameter „bck_const“: upper limit of the number of counts (per bin) of the constant background in the background spectrum. Can be increased when the background is very high, e.g. due to Compton background from a pronounced line.
* lower limit of fit parameter „bck_gauss“: the lower limit of the (integrated) number of background counts in the peak. The value should always be set to 0. 
* upper limit of fit parameter „bck_gauss“: upper limit of the (integrated) number of background counts in the peak. Can be increased if the line in the background spectrum has very high activity. 
* integration constant: constant that is added in the calculation of the likelihood. Sometimes the integration of the posterior pdf which is needed to calculate the Bayesfactor returns 0. This happens when the likelihood is so small that values are below the smallest double value. To avoid this an arbitrary constant is added in the calculation of the likelihood. The results of the integration can be found in the <sample_name_isotope>(_bckonly)_results.txt file under „Evidence“. If the value shows „0“ or „inf“ the integration constant has to be increased or decreased.
