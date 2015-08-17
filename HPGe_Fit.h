#ifndef __HPGe_Fit__H
#define __HPGe_Fit__H

// ---------------------------------------------------------

#include <TH1D.h>
#include <TString.h>
#include "BCMTF_HPGe.h"

// ---------------------------------------------------------
class HPGe_Fit
{

public:

    /** \name Constructors and destructors */
    /** @{ */
    
    /**
     * The default constructor. */
    HPGe_Fit();
    
    /**
     * The default destructor. */
    ~HPGe_Fit();
    
    /** @} */
    /** \name Member functions (get) */
    /** @{ */
    
    /**
     * @return The Bayes Factor */
    double GetBayesFactor()
    { return fBayesFactor; };
    
    /**
     * @return The activity as string. */
    TString GetActivity()
    { return fActivityStr; };
    
    /**
     * @return The p-value. */
    double GetPValue()
    { return fPvalue; };
    
    
    /** @} */
    
    /** \name Member functions (set) */
    /** @{ */

    /**
     * Set the folder with isotope parameters files
     */
    void SetParametersFolder(TString folder)
    {fparameters_folder = folder;};
    
    /**
     * Set the precision of the MCMC
     */
    void SetPrecision(TString precision)
    {fprecision = precision;};
    
    /**
     * Set the sample and background histograms
     */
    void SetSpectra(TH1D* hist_sample, TH1D* hist_bck)
    {
        fhist_bck = hist_bck;
        fhist_sample = hist_sample;
    };
    
    /**
     * Set the file with the simulated efficiencies
     */
    void SetEfficiencyFile(TString efficiency_name)
    {fefficiency_name = efficiency_name;};
    
    /**
     * Set the file with the energy resolution
     */
    void SetResolutionFile(TString resolution_name)
    {fresolution_name = resolution_name;};
    
    /**
     * Set the results file name
     */
    void SetResultsFile(TString results_name)
    {fresults_name = results_name;};
    
    /**
     * Set the measurement times of the sample and background spectra
     */
    void SetMeasurementTimes(double t_sample, double t_bck)
    {
        ft_sample = t_sample;
        ft_bck = t_bck;
    };
    
    /**
     * Run the fit for an isotope
     */
    int RunFit(TString isotope_name);

        /** @} */

private:
    
    // read parameters from file 'parameters_'isotope'.txt'
    int read_parameters_isotope(TString isotope_name);
    
    // read simulated efficiencies from file
    int read_efficiencies();
    
    // read energy resolution from file
    int read_resolution();
    
    // get histogram in range
    TH1D* GetHistoRange(TString name, TH1D* hist, double range_low, double range_up);
    
    // create template histograms
    TH1D* CreateHistConst(TString name, double range_low, double range_up, int nBins);
    TH1D* CreateHistGauss(TString name, double mean, double sigma,  double range_low, double range_up, int nBins);

    // get first or last empty bin
    int GetEmptyBin(TH1D* hist, TString option);
    
    // eliminate empty bins
    void SetNewLimits(BCMTF_HPGe* m, TString parameter);
    
    // folder for isotope parameter files
    TString fparameters_folder;
    
    // names of in/output files
    TString fresults_name, fresolution_name, fefficiency_name;
    
    // sample and background spectra
    TH1D* fhist_sample;
    TH1D* fhist_bck;
    
    // precision of MCMC
    TString fprecision;
    
    // p-value
    double fPvalue;
    
    // Bayes Factor
    double fBayesFactor;
    
    // calculated activity as a string
    TString fActivityStr;
    
    /**
     * measurement time of the sample and background spectra (in sec.)*/
    double ft_sample, ft_bck;
    
    
    // parameters to read from file 'parameters_'isotope'.txt'
    double fconstant;
    double flimit_sig_low;
    double flimit_sig_up;
    int fNpeaks;
    std::vector<double> fpeak_energy;
    std::vector<double> fefficiency_bck;
    std::vector<double> ffitRange_low;
    std::vector<double> ffitRange_high;
    std::vector<double> flimit_const_sample_low;
    std::vector<double> flimit_const_sample_up;
    std::vector<double> flimit_const_bck_low;
    std::vector<double> flimit_const_bck_up;
    std::vector<double> flimit_gauss_bck_low;
    std::vector<double> flimit_gauss_bck_up;
    
    // parameters to read from efficiency and resolution file
    std::vector<double> fpeak_sigma;
    std::vector<double> fefficiency_sig;

};
// ---------------------------------------------------------

#endif

