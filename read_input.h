 // ***************************************************************
// This file was created using the CreateProject.sh script
// for project DBD.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>
#include <vector>
#include <fstream>
#include <string>
#include <TVectorD.h>

#include <iostream>
using namespace std;

#endif

// parameters to read from file 'parameters_global.txt'
TString fname, faccuracy, fparameters_folder, fsample_name, fbck_name, fresults_name, fefficiency_name, fresolution_name;
std::vector<TString> fisotope_name;
double ft_sample, ft_bck, fBF_limit, fCL, feff_err;

// sample and background spectra
TH1D* fhist_sample;
TH1D* fhist_bck;


// ----------------------------------------------------
// read global parameters from user specified file
// ----------------------------------------------------


int read_parameters_global(TString FileName) {
    
    // clear isotope name vector
    fisotope_name.clear();
    
    ifstream File;
    File.open(FileName);
    
    if (!File.is_open()) {
        std::cout << "##### ERROR: could not open " << FileName << std::endl;
        return 1;
    }
    
    std::string headerline;
    TString isotope_name;
    
    getline(File, headerline);
    File >> fname;
    getline(File, headerline);
    getline(File, headerline);
    File >> faccuracy;
    getline(File, headerline);
    getline(File, headerline);
    File >> fBF_limit;
    getline(File, headerline);
    getline(File, headerline);
    File >> fCL;
    getline(File, headerline);
    getline(File, headerline);
    File >> feff_err;
    getline(File, headerline);
    getline(File, headerline);
    File >> fparameters_folder;
    getline(File, headerline);
    getline(File, headerline);
    File >> fsample_name;
    getline(File, headerline);
    getline(File, headerline);
    File >> fbck_name;
    getline(File, headerline);
    getline(File, headerline);
    File >> fefficiency_name;
    getline(File, headerline);
    getline(File, headerline);
    File >> fresolution_name;
    getline(File, headerline);
    getline(File, headerline);
    File >> fresults_name;
    getline(File, headerline);
    getline(File, headerline);
    while (true)
    {
        File >> isotope_name;
        
        if(isotope_name=="" ) break;

        fisotope_name.push_back(isotope_name);

    }
    
    File.close();
    
    // print information
    std::cout << "######################################" << std::endl;
    std::cout << "#### Reading "<< FileName << " ..." << std::endl;
    std::cout << "sample name: " << fname << std::endl;
    std::cout << "accuracy of MCMC: " << faccuracy << std::endl;
    std::cout << "threshold for signal detection: " << fBF_limit << std::endl;
    std::cout << "CL for limit on activity: " << fCL << std::endl;
    std::cout << "fractional uncertainty of efficiencies: " << feff_err << std::endl;
    std::cout << "isotope parameters folder: " << fsample_name << std::endl;
    std::cout << "sample spectrum filename: " << fsample_name << std::endl;
    std::cout << "bck spectrum filename: " << fbck_name << std::endl;
    std::cout << "efficiency filename: " << fefficiency_name << std::endl;
    std::cout << "resolution filename: " << fresolution_name << std::endl;
    std::cout << "results folder: " << fresults_name << std::endl;
    std::cout << "isotopes to analyze:" << std::endl;
    for (int i=0; i<fisotope_name.size(); ++i)
    {
        std::cout << fisotope_name[i] << std::endl;
    }
    std::cout << "######################################" << std::endl;
    
    return 0;
}



// ----------------------------------------------------
// open spectrum files
// ----------------------------------------------------

int read_spectra() {
    
    TFile* file_bck = TFile::Open(fbck_name);
    
    if (!file_bck) {
        std::cout << "##### ERROR: could not open " << fbck_name << std::endl;
        return 1;
    }
    
    TFile* file_sample = TFile::Open(fsample_name);
    
    if (!file_sample) {
        std::cout << "##### ERROR: could not open " << fsample_name << std::endl;
        return 1;
    }

    if (!file_bck->GetListOfKeys()->Contains("hist")) {
        std::cout << "##### ERROR: no histogram in file " << fbck_name << std::endl;
        return 1;
    }
    
    if (!file_sample->GetListOfKeys()->Contains("hist")) {
        std::cout << "##### ERROR: no histogram in file " << fsample_name << std::endl;
        return 1;
    }
    
    fhist_bck = (TH1D*) file_bck->Get("hist");
    fhist_sample = (TH1D*) file_sample->Get("hist");
    
    
    if (!file_bck->GetListOfKeys()->Contains("t_live")) {
        std::cout << "##### ERROR: no live time in file " << fbck_name << std::endl;
        return 1;
    }

    if (!file_sample->GetListOfKeys()->Contains("t_live")) {
        std::cout << "##### ERROR: no live time in file " << fsample_name << std::endl;
        return 1;
    }

    
    TVectorD* v_live_bck = (TVectorD*) file_bck->Get("t_live");
    ft_bck = ((*v_live_bck)[0]);
    
    TVectorD* v_live_sample = (TVectorD*) file_sample->Get("t_live");
    ft_sample = ((*v_live_sample)[0]);
    
    return 0;


}


    


