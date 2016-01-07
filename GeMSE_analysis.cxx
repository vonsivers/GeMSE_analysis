/*
 *  W180.c
 *  
 *
 *  Created by sivers on 20.10.13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */
#include "HPGe_Fit.h"
#include "BCMTF_HPGe.h"
#include "read_input.h"

#include <TSystem.h>

#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char *argv[]) {
    
    // argc should be 2 for correct execution
    if ( argc != 2 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <filename>" << std::endl;
        return 1;
    }

    TString FileName = argv[1];
    
    // read global input parameters from file
    if (read_parameters_global(FileName)) return 1;

    // open spectrum files
    if (read_spectra()) return 1;

    // try to open results directory
    if (!gSystem->OpenDirectory(fresults_name)) {
        
        // if directory does not exist make one
        if (gSystem->MakeDirectory(fresults_name)==-1) {
            std::cout << "###### ERROR: could not create directory " << fresults_name << std::endl;
        }
    }
    
    // open results file
    ofstream results_file;
    results_file.open(fresults_name+"/"+fname+"_activities_summary.txt");
    
    // get current date and time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    
    strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
    std::string timestring(buffer);
    
    // write header
    results_file << "#################################" << "\n";
    results_file << timestring << "\n";
    results_file << "sample name: " << fname << "\n";
    results_file << "sample spectrum: " << fsample_name << "\n";
    results_file << "background spectrum: " << fbck_name << "\n";
    results_file << "simulated efficiencies: " << fefficiency_name << "\n";
    results_file << "energy resolution: " << fresolution_name << "\n";
    results_file << "measurement time sample: " << ft_sample << " sec." << "\n";
    results_file << "measurement time background: " << ft_bck << " sec." << "\n";
    results_file << "threshold for signal: " << fBF_limit << "\n";
    results_file << "CL for activity limit: " << fCL << "\n";
    results_file << "#################################" << "\n";
    results_file << "\n";
    results_file << "Isotope \t Activity (Bq) \t Bayes Factor" << "\n";
    
    // create a new fitter
    HPGe_Fit* fit = new HPGe_Fit();
    
    // set parameters
    fit->SetParametersFolder(fparameters_folder);
    fit->SetPrecision(faccuracy);
    fit->SetBFLimit(fBF_limit);
    fit->SetCL(fCL);
    fit->SetSpectra(fhist_sample, fhist_bck);
    fit->SetEfficiencyFile(fefficiency_name);
    fit->SetResolutionFile(fresolution_name);
    fit->SetResultsFile(fresults_name+"/"+fname);
    fit->SetMeasurementTimes(ft_sample,ft_bck);

    // loop over all isotopes
    for (int i=0; i<fisotope_name.size(); ++i) {
        
        // run fit for isotope
        if (fit->RunFit(fisotope_name[i])) return 1;
        
        // print results to txt file
        results_file << fisotope_name[i] << "\t" << fit->GetActivity() << "\t" << fit->GetBayesFactor() << "\n";

    }
    
    results_file.close();

    
    return 0;
}






