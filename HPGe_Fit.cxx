/*
 *  W180.c
 *  
 *
 *  Created by sivers on 20.10.13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */
#include "HPGe_Fit.h"

#include <BAT/BCModelManager.h>
#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCModelOutput.h>
//#include <BAT/BCSummaryTool.h>
#include <BAT/BCH1D.h>
#include <BAT/BCMTFChannel.h>
#include <BAT/BCParameter.h>


#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TROOT.h>

#include <iostream>
#include <fstream>
using namespace std;


// ---------------------------------------------------------
HPGe_Fit::HPGe_Fit()
{}

// ---------------------------------------------------------
HPGe_Fit::~HPGe_Fit()
{}






// ----------------------------------------------------
// run fit for one isotope
// ----------------------------------------------------

int HPGe_Fit::RunFit(TString isotope_name) {
    
    // data spectra
    std::vector<TH1D*> hist_sample_peak;
    std::vector<TH1D*> hist_bck_peak;
    
    // template histograms
    std::vector<TH1D*> temp_const_sample;
    std::vector<TH1D*> temp_const_bck;
    std::vector<TH1D*> temp_gauss_sample;
    std::vector<TH1D*> temp_gauss_bck;
    
    
    // check if all parameters have been set
    if (ft_bck==0||ft_sample==0||fresults_name==""||fefficiency_name==""||fresolution_name==""||fprecision==""||fhist_sample==0||fhist_bck==0||fBF_limit==0||fCL==0) {
        std::cout << "##### ERROR: all parameters must be set before running the fit" << std::endl;
        return 1;
    }
    
    // read input parameters from 'parameters_'isotope'.txt'
    if (read_parameters_isotope(isotope_name)) return 1;
    
    // read efficiencies from file
    if (read_efficiencies()) return 1;
    
    // read resolution from file
    if (read_resolution()) return 1;
    
    // names of all histos
    TString name_hist_sample, name_hist_bck, name_temp_const_sample, name_temp_const_bck, name_temp_gauss_sample, name_temp_gauss_bck;
    
    // set results name
    TString results_name = fresults_name+"_"+isotope_name;
    
    // loop over all peaks
    for (int i=0; i<fNpeaks; ++i) {
        
        // names of data histograms
        name_hist_sample = TString::Format("hist_sample_peak%d",i);
        name_hist_bck = TString::Format("hist_bck_peak%d",i);
        
        // make data histograms
        hist_sample_peak.push_back(GetHistoRange(name_hist_sample,fhist_sample,ffitRange_low[i],ffitRange_high[i]));
        hist_bck_peak.push_back(GetHistoRange(name_hist_bck,fhist_bck,ffitRange_low[i],ffitRange_high[i]));
        
        // get number of bins
        int nbins_sample = hist_sample_peak[i]->GetNbinsX();
        int nbins_bck = hist_bck_peak[i]->GetNbinsX();
        
        // get xbins
        double* xbins_sample = hist_sample_peak[i]->GetXaxis()->GetXbins()->fArray;
        double* xbins_bck = hist_bck_peak[i]->GetXaxis()->GetXbins()->fArray;
        
        // names of template histograms
        name_temp_const_sample = TString::Format("temp_const_sample_peak%d",i);
        name_temp_const_bck = TString::Format("temp_const_bck_peak%d",i);
        name_temp_gauss_sample = TString::Format("temp_gauss_sample_peak%d",i);
        name_temp_gauss_bck = TString::Format("temp_gauss_bck_peak%d",i);

        //std::cout << "##### making templates ..." << std::endl;

        // make template histograms
        temp_const_sample.push_back(CreateHistConst(name_temp_const_sample, nbins_sample, xbins_sample));
        temp_const_bck.push_back(CreateHistConst(name_temp_const_bck, nbins_bck, xbins_bck));
        temp_gauss_sample.push_back(CreateHistGauss(name_temp_gauss_sample, fpeak_energy[i], fpeak_sigma[i], nbins_sample, xbins_sample));
        temp_gauss_bck.push_back(CreateHistGauss(name_temp_gauss_bck, fpeak_energy[i], fpeak_sigma[i], nbins_bck, xbins_bck));
        
        //delete edge_sample;
        //delete edge_bck;

    }
    
   	// ---------------------------------------
    // signal + bck fit
    // ---------------------------------------
    std::cout << "####################################" << std::endl;
    std::cout << "starting signal + bck fit" << std::endl;
    std::cout << "####################################" << std::endl;
    
    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();
    
    // open log file
    BCLog::OpenLog(results_name+"_log.txt");
    BCLog::SetLogLevel(BCLog::warning);
    
	// create new fitter object
	BCMTF_HPGe * m = new BCMTF_HPGe();
    
    // set seed for reproducible results
    m->MCMCSetRandomSeed(21340);
    
    // set Metropolis as marginalization method
    m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    
    // set number of bins for marginalized distributions
    //m->SetNbins(500);
    
    // add signal process
    m->AddProcess("signal", flimit_sig_low,flimit_sig_up);
    
    // names of processes
    TString name_process_sample_const, name_process_bck_const, name_process_bck_gauss;
    
    // names of channels
    TString name_channel_sample, name_channel_bck;
    
    // add systematic uncertainty for efficiencies
    m->AddSystematic("efficiency_err", -5., 5.);
    m->SetPriorGauss("efficiency_err", 0., 1.);

    // loop over all peaks
    for (int i=0; i<fNpeaks; ++i) {
        
        // name processes
        name_process_sample_const = TString::Format("sample_const_peak%d",i);
        name_process_bck_const = TString::Format("bck_const_peak%d",i);
        name_process_bck_gauss = TString::Format("bck_gauss_peak%d",i);
        
        // name channels
        name_channel_sample = TString::Format("sample_peak%d",i);
        name_channel_bck = TString::Format("bck_peak%d",i);
        
        // add channels (the order is important!)
        m->AddChannel(name_channel_sample);
        m->AddChannel(name_channel_bck);
        
        // add processes
        m->AddProcess(name_process_sample_const, flimit_const_sample_low[i], flimit_const_sample_up[i]);
        m->AddProcess(name_process_bck_const, flimit_const_bck_low[i], flimit_const_bck_up[i]);
        m->AddProcess(name_process_bck_gauss, flimit_gauss_bck_low[i], flimit_gauss_bck_up[i]);
        
        // set data
        m->SetData(name_channel_sample, *hist_sample_peak[i]);
        m->SetData(name_channel_bck, *hist_bck_peak[i]);
        
        // set template and histograms
        m->SetTemplate(name_channel_sample, "signal", *temp_gauss_sample[i], fefficiency_sig[i]);
        m->SetTemplate(name_channel_sample, name_process_sample_const, *temp_const_sample[i], 1.);
        m->SetTemplate(name_channel_sample, name_process_bck_gauss, *temp_gauss_sample[i], fefficiency_bck[i]*ft_sample/ft_bck);
        m->SetTemplate(name_channel_bck, name_process_bck_const, *temp_const_bck[i], 1.);
        m->SetTemplate(name_channel_bck, name_process_bck_gauss, *temp_gauss_bck[i], 1.);
        
        // set systematics
        m->SetSystematicVariation(name_channel_sample, "signal", "efficiency_err", feff_err, feff_err);

        // set priors
        m->SetPriorConstant("signal");
        m->SetPriorConstant(name_process_sample_const);
        m->SetPriorConstant(name_process_bck_const);
        m->SetPriorConstant(name_process_bck_gauss);
     
    }
	
    // -----------------------------------------------------------
    // run marginalization to find right parameter limits
    // -----------------------------------------------------------
    std::cout << "######################################" << std::endl;
    std::cout << "#### pre-run fit to set parameter limits ..." << std::endl;
    std::cout << "######################################" << std::endl;

    m->MCMCSetPrecision(BCEngineMCMC::kLow);
    
    // do it 2 times
    for (int i=0; i<2; ++i) {
        
        m->MarginalizeAll();
        
        //----------------------------
        // first comes signal parameter
        //----------------------------
        
        SetNewLimits(m, "signal");
        
        //----------------------------
        // then loop over all peaks
        //----------------------------
        
        for (int i=0; i<fNpeaks; ++i) {
            
            // name processes
            name_process_sample_const = TString::Format("sample_const_peak%d",i);
            name_process_bck_const = TString::Format("bck_const_peak%d",i);
            name_process_bck_gauss = TString::Format("bck_gauss_peak%d",i);
            
            SetNewLimits(m, name_process_sample_const);
            SetNewLimits(m, name_process_bck_const);
            SetNewLimits(m, name_process_bck_gauss);
            
        }

    }
    

    // -----------------------------------------------------------
    // all parameter limits are set now continue with fit
    // -----------------------------------------------------------
    std::cout << "######################################" << std::endl;
    std::cout << "#### now running actual fit ..." << std::endl;
    std::cout << "######################################" << std::endl;
    
    // set constant for integration
    m->SetConstant(fconstant);
    
    // set precision
    if (fprecision=="low") {
        m->MCMCSetPrecision(BCEngineMCMC::kLow);
    }
    else if (fprecision=="medium") {
        m->MCMCSetPrecision(BCEngineMCMC::kMedium);
    }
    else if (fprecision=="high") {
        m->MCMCSetPrecision(BCEngineMCMC::kHigh);
    }
    
    // create new output object
    BCModelOutput* mout = new BCModelOutput(m, results_name+"_output.root");
    
    // create a new summary tool object
    //BCSummaryTool * summary = new BCSummaryTool(m);
    
    // marginalize
    m->MarginalizeAll();
    
    // find global mode
    m->FindMode( m->GetBestFitParameters() );
    
    // calculate normalization
    m->Integrate();
    
    // get signal counts and error
    double counts = m->GetMarginalized("signal")->GetMode();
    double counts_err_low = counts - m->GetMarginalized("signal")->GetQuantile(0.16);
    double counts_err_up = m->GetMarginalized("signal")->GetQuantile(0.84) - counts;
    
    // get upper limit for signal counts
    double counts_limit = m->GetMarginalized("signal")->GetQuantile(fCL);
    
    // print marginalized distributions
    m->PrintAllMarginalized(results_name+"_distributions.pdf");
    
    // print results file
    m->PrintResults(results_name+"_results.txt");
    
    // print summary results
    //summary->PrintKnowledgeUpdatePlots(results_name+"_summary_update.pdf");
    
    // print templates and stacks
    for (int i = 0; i < m->GetNChannels(); ++i) {
        BCMTFChannel * channel = m->GetChannel(i);
        //channel->PrintTemplates(Form("%s_templates_%s.pdf",results_name.Data(),channel->GetName().c_str()));
        m->PrintStack(i, m->GetBestFitParameters(), Form("%s_stack_%s.pdf",results_name.Data(),channel->GetName().c_str()));
    }
    
    // write marginalized distributions to output file
    mout->WriteMarginalizedDistributions();
    
    // close log file
    BCLog::CloseLog();
    
    // close output file
    mout->Close();
    
    // free memory
    //delete m;
    delete mout;
    //delete summary;
    
    // ---------------------------------------
    // bck-only fit
    // ---------------------------------------
    std::cout << "####################################" << std::endl;
    std::cout << "starting bck-only fit" << std::endl;
    std::cout << "####################################" << std::endl;
    
    results_name+="_bckonly";
    
    // open log file
    BCLog::OpenLog(results_name+"_log.txt");
    BCLog::SetLogLevel(BCLog::warning);
    
    // create new fitter object
    BCMTF_HPGe* m_bck = new BCMTF_HPGe();
    
    // set seed for reproducible results
    m_bck->MCMCSetRandomSeed(21340);
    
    // set Metropolis as marginalization method
    m_bck->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    
    // set number of bins for marginalized distributions
    //m_bck->SetNbins(500);
    
    // loop over all peaks
    for (int i=0; i<fNpeaks; ++i) {
        
        // name processes
        name_process_sample_const = TString::Format("sample_const_peak%d",i);
        name_process_bck_const = TString::Format("bck_const_peak%d",i);
        name_process_bck_gauss = TString::Format("bck_gauss_peak%d",i);
        
        // name channels
        name_channel_sample = TString::Format("sample_peak%d",i);
        name_channel_bck = TString::Format("bck_peak%d",i);
        
        // add channels
        m_bck->AddChannel(name_channel_sample);
        m_bck->AddChannel(name_channel_bck);
        
        // add processes
        m_bck->AddProcess(name_process_sample_const, flimit_const_sample_low[i], flimit_const_sample_up[i]);
        m_bck->AddProcess(name_process_bck_const, flimit_const_bck_low[i], flimit_const_bck_up[i]);
        m_bck->AddProcess(name_process_bck_gauss, flimit_gauss_bck_low[i], flimit_gauss_bck_up[i]);
        
        // set data
        m_bck->SetData(name_channel_sample, *hist_sample_peak[i]);
        m_bck->SetData(name_channel_bck, *hist_bck_peak[i]);
        
        // set template and histograms
        m_bck->SetTemplate(name_channel_sample, name_process_sample_const, *temp_const_sample[i], 1.);
        m_bck->SetTemplate(name_channel_sample, name_process_bck_gauss, *temp_gauss_sample[i], fefficiency_bck[i]*ft_sample/ft_bck);
        m_bck->SetTemplate(name_channel_bck, name_process_bck_const, *temp_const_bck[i], 1.);
        m_bck->SetTemplate(name_channel_bck, name_process_bck_gauss, *temp_gauss_bck[i], 1.);
        
        // set priors
        m_bck->SetPriorConstant(name_process_sample_const);
        m_bck->SetPriorConstant(name_process_bck_const);
        m_bck->SetPriorConstant(name_process_bck_gauss);
        
    }
    
    // -----------------------------------------------------------
    // run marginalization to find parameter limits
    // -----------------------------------------------------------
    std::cout << "######################################" << std::endl;
    std::cout << "#### pre-run fit to set parameter limits ..." << std::endl;
    std::cout << "######################################" << std::endl;

    m_bck->MCMCSetPrecision(BCEngineMCMC::kLow);
    
    // do it 2 times
    for (int i=0; i<2; ++i) {

        m_bck->MarginalizeAll();
        
        // loop over all peaks
        for (int i=0; i<fNpeaks; ++i) {
            
            // name processes
            name_process_sample_const = TString::Format("sample_const_peak%d",i);
            name_process_bck_const = TString::Format("bck_const_peak%d",i);
            name_process_bck_gauss = TString::Format("bck_gauss_peak%d",i);
            
            SetNewLimits(m_bck, name_process_sample_const);
            SetNewLimits(m_bck, name_process_bck_const);
            SetNewLimits(m_bck, name_process_bck_gauss);
            
        }
    }
    
    // -----------------------------------------------------------
    // all parameter limits are set now continue with fit
    // -----------------------------------------------------------
    std::cout << "######################################" << std::endl;
    std::cout << "#### now running actual fit ..." << std::endl;
    std::cout << "######################################" << std::endl;

    // set constant for integration
    m_bck->SetConstant(fconstant);
    
    // set precision
    if (fprecision=="low") {
        m_bck->MCMCSetPrecision(BCEngineMCMC::kLow);
    }
    else if (fprecision=="medium") {
        m_bck->MCMCSetPrecision(BCEngineMCMC::kMedium);
    }
    else if (fprecision=="high") {
        m_bck->MCMCSetPrecision(BCEngineMCMC::kHigh);
    }
    
    // create new output object
    BCModelOutput* mout_bck = new BCModelOutput(m_bck, results_name+"_output.root");
    
    // create a new summary tool object
    //BCSummaryTool * summary_bck = new BCSummaryTool(m_bck);
    
    // marginalize
    m_bck->MarginalizeAll();
    
    // find global mode
    m_bck->FindMode( m_bck->GetBestFitParameters() );
    
    // calculate normalization
    m_bck->Integrate();
    
    // print marginalized distributions
    m_bck->PrintAllMarginalized(results_name+"_distributions.pdf");
    
    // print results file
    m_bck->PrintResults(results_name+"_results.txt");
    
    // print summary results
    //summary_bck->PrintKnowledgeUpdatePlots(results_name+"_summary_update.pdf");
    
    // print templates and stacks
    for (int i = 0; i < m_bck->GetNChannels(); ++i) {
        BCMTFChannel * channel = m_bck->GetChannel(i);
        //channel->PrintTemplates(Form("%s_templates_%s.pdf",results_name.Data(),channel->GetName().c_str()));
        m_bck->PrintStack(i, m_bck->GetBestFitParameters(), Form("%s_stack_%s.pdf",results_name.Data(),channel->GetName().c_str()));
    }
    
    // write marginalized distributions to output file
    mout_bck->WriteMarginalizedDistributions();
    
    // close log file
    BCLog::CloseLog();
    
    // close output file
    mout_bck->Close();
    
    // free memory
    //delete m_bck;
    delete mout_bck;
    //delete summary_bck;
    
    // ----------------------------------------------------
    // set up model manager
    // ----------------------------------------------------
    
    // create new BCTemplateFitterManager
    BCModelManager * smm = new BCModelManager();
    smm->AddModel(m_bck);
    smm->AddModel(m);
    
    // compare models
    fBayesFactor = smm->BayesFactor(0,1);
    
    // check for positive signal and calculate activity
    if (fBayesFactor<fBF_limit) {
        fActivityStr = TString::Format("%1.2e - %1.2e + %1.2e",counts/ft_sample,counts_err_low/ft_sample,counts_err_up/ft_sample);
    }
    else {
        fActivityStr = TString::Format("< %1.2e",counts_limit/ft_sample);
    }
    
    
    
    // ----------------------------------------------------
    // clean up
    // ----------------------------------------------------
    
    delete m;
    delete m_bck;
    delete smm;
    
    for(int i = 0; i < fNpeaks; ++i)
    {
        name_hist_sample = TString::Format("hist_sample_peak%d",i);
        name_hist_bck = TString::Format("hist_bck_peak%d",i);
        name_temp_const_sample = TString::Format("temp_const_sample_peak%d",i);
        name_temp_const_bck = TString::Format("temp_const_bck_peak%d",i);
        name_temp_gauss_sample = TString::Format("temp_gauss_sample_peak%d",i);
        name_temp_gauss_bck = TString::Format("temp_gauss_bck_peak%d",i);

        while ((hist_sample_peak[i] = (TH1D*) gROOT->FindObject(name_hist_sample))) {
               delete hist_sample_peak[i];
               hist_sample_peak[i] = NULL;
        }

        while ((hist_bck_peak[i] = (TH1D*) gROOT->FindObject(name_hist_bck))) {
            delete hist_bck_peak[i];
            hist_bck_peak[i] = NULL;
        }
        
        while ((temp_const_sample[i] = (TH1D*) gROOT->FindObject(name_temp_const_sample))) {
            delete temp_const_sample[i];
            temp_const_sample[i] = NULL;
        }
        
        while ((temp_const_bck[i] = (TH1D*) gROOT->FindObject(name_temp_const_bck))) {
            delete temp_const_bck[i];
            temp_const_bck[i] = NULL;
        }
        
        while ((temp_gauss_sample[i] = (TH1D*) gROOT->FindObject(name_temp_gauss_sample))) {
            delete temp_gauss_sample[i];
            temp_gauss_sample[i] = NULL;
        }
        
        while ((temp_gauss_bck[i] = (TH1D*) gROOT->FindObject(name_temp_gauss_bck))) {
            delete temp_gauss_bck[i];
            temp_gauss_bck[i] = NULL;
        }
        

    }
    
    
    return 0;
    
}



// ----------------------------------------------------
// Read Parameters from 'parameters_'isotope'.txt' file
// ----------------------------------------------------

int HPGe_Fit::read_parameters_isotope(TString isotope_name) {
    
    // clear all vectors
    fpeak_energy.clear();
    fefficiency_bck.clear();
    ffitRange_low.clear();
    ffitRange_high.clear();
    flimit_const_sample_low.clear();
    flimit_const_sample_up.clear();
    flimit_const_bck_low.clear();
    flimit_const_bck_up.clear();
    flimit_gauss_bck_low.clear();
    flimit_gauss_bck_up.clear();
    
    // open file
    TString FileName = fparameters_folder+"/parameters_"+isotope_name+".txt";
    
    ifstream File;
    File.open(FileName);
    
    if (!File.is_open()) {
        std::cout << "##### ERROR: could not open " << FileName << std::endl;
        return 1;
    }
    
    std::string headerline;
    
    double peak_energy, efficiency_bck, fitRange_low, fitRange_high, limit_const_sample_low, limit_const_sample_up, limit_const_bck_low, limit_const_bck_up, limit_gauss_bck_low, limit_gauss_bck_up;
    
    getline(File, headerline);
    File >> fNpeaks;
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> peak_energy;
        fpeak_energy.push_back(peak_energy);
    }
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> efficiency_bck;
        fefficiency_bck.push_back(efficiency_bck);
    }
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> fitRange_low;
        ffitRange_low.push_back(fitRange_low);
    }
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> fitRange_high;
        ffitRange_high.push_back(fitRange_high);
    }
    getline(File, headerline);
    getline(File, headerline);
    File >> flimit_sig_low;
    getline(File, headerline);
    getline(File, headerline);
    File >> flimit_sig_up;
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> limit_const_sample_low;
        flimit_const_sample_low.push_back(limit_const_sample_low);
    }
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> limit_const_sample_up;
        flimit_const_sample_up.push_back(limit_const_sample_up);
    }
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> limit_const_bck_low;
        flimit_const_bck_low.push_back(limit_const_bck_low);
    }
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> limit_const_bck_up;
        flimit_const_bck_up.push_back(limit_const_bck_up);
    }
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> limit_gauss_bck_low;
        flimit_gauss_bck_low.push_back(limit_gauss_bck_low);
    }
    getline(File, headerline);
    getline(File, headerline);
    for (int i=0; i<fNpeaks; ++i)
    {
        File >> limit_gauss_bck_up;
        flimit_gauss_bck_up.push_back(limit_gauss_bck_up);
        
    }
    getline(File, headerline);
    getline(File, headerline);
    File >> fconstant;
    
    File.close();
    
    // print information
    std::cout << "######################################" << std::endl;
    std::cout << "#### Reading " << FileName << " ..." << std::endl;
    std::cout << "number of peaks: " << fNpeaks << std::endl;
    std::cout << "peak energies:" << std::endl;
    for (int i=0; i<fNpeaks; ++i)
    {
        std::cout << fpeak_energy[i] << std::endl;
    }
    std::cout << "bck efficiencies:" << std::endl;
    for (int i=0; i<fNpeaks; ++i)
    {
        std::cout << fefficiency_bck[i] << std::endl;
    }
    std::cout << "fit ranges:" << std::endl;
    for (int i=0; i<fNpeaks; ++i)
    {
        std::cout << ffitRange_low[i] << " - " << ffitRange_high[i] << std::endl;
    }
    std::cout << "limits of fit parameter 'signal':" << std::endl;
    std::cout << flimit_sig_low << " - " << flimit_sig_up << std::endl;
    std::cout << "limit of fit parameter 'sample_const':" << std::endl;
    for (int i=0; i<fNpeaks; ++i)
    {
        std::cout << flimit_const_sample_low[i] << " - " << flimit_const_sample_up[i] << std::endl;
    }
    std::cout << "limit of fit parameter 'bck_const':" << std::endl;
    for (int i=0; i<fNpeaks; ++i)
    {
        std::cout << flimit_const_bck_low[i] << " - " << flimit_const_bck_up[i] << std::endl;
    }
    std::cout << "limit of fit parameter 'bck_gauss':" << std::endl;
    for (int i=0; i<fNpeaks; ++i)
    {
        std::cout << flimit_gauss_bck_low[i] << " - " << flimit_gauss_bck_up[i] << std::endl;
    }
    std::cout << "integration constant: " << fconstant << std::endl;
    std::cout << "######################################" << std::endl;
    
    return 0;
}

// ----------------------------------------------------
// Read simulated efficiencies from file
// ----------------------------------------------------

int HPGe_Fit::read_efficiencies() {
    
    // clear efficiency vector
    fefficiency_sig.clear();
    
    
    // open file with simulated efficiencies
    TFile* File = TFile::Open(fefficiency_name);
    
    if (!File) {
        std::cout << "##### ERROR: could not open " << fefficiency_name << std::endl;
        return 1;
    }
    
    std::cout << "#### Reading " << fefficiency_name << " ..." << std::endl;
    
    if (!File->GetListOfKeys()->Contains("tree")) {
        std::cout << "##### ERROR: no tree in file " << fefficiency_name << std::endl;
        return 1;
    }
    
    TTree* tree = (TTree*) File->Get("tree");
    
    double energy, efficiency;
    
    tree->SetBranchAddress("energy",&energy);
    tree->SetBranchAddress("eff_BR",&efficiency);
    
    int nEntries = tree->GetEntries();
    
    // find all efficiencies for one isotope
    for (int j=0; j<fpeak_energy.size(); ++j) {
        
        for (int i=0; i<nEntries; ++i) {
           
            tree->GetEntry(i);
            
            if (energy==fpeak_energy[j]) {
                
                //std::cout << "##### found " << fpeak_energy[j] << " keV peak" << std::endl;
                fefficiency_sig.push_back(efficiency);
                break;
                
            }
            // check if efficiency was found
            if (i==nEntries-1) {
                std::cout << "##### ERROR: could not find efficiency for " << fpeak_energy[j] << " keV peak" << std::endl;
                return 1;
            }
        }
        
    }
    
    File->Close();
    
    return 0;
}

// ----------------------------------------------------
// Read energy resolution from file
// ----------------------------------------------------

int HPGe_Fit::read_resolution() {
    
    // clear resolution vector
    fpeak_sigma.clear();
    
    // open file with energy resolution
    TFile* File = TFile::Open(fresolution_name);
    
    if (!File) {
        std::cout << "##### ERROR: could not open " << fresolution_name << std::endl;
        return 1;
    }
    
    std::cout << "#### Reading " << fresolution_name << " ..." << std::endl;
    
    if (!File->GetListOfKeys()->Contains("c2")) {
        std::cout << "##### ERROR: no canvas in file " << fresolution_name << std::endl;
        return 1;
    }
    
    
    TCanvas* c2 = (TCanvas*) File->Get("c2");
    TGraphErrors* graph = (TGraphErrors*) c2->GetPrimitive("Graph");
    TF1* resolution = graph->GetFunction("fitFunction");
    
    // calculate resolutions for one isotope
    for (int j=0; j<fpeak_energy.size(); ++j) {
        fpeak_sigma.push_back(resolution->Eval(fpeak_energy[j]));
    }
    
    File->Close();
    
    return 0;
}

// ----------------------------------------------------
// Get Histograms in Range
// ----------------------------------------------------

TH1D* HPGe_Fit::GetHistoRange(TString name, TH1D* hist, double range_low, double range_up) {
    
    //int nbins=hist->FindBin(range_up)-hist->FindBin(range_low)+1;
    int bin_low = hist->FindBin(range_low);
    int bin_up = hist->FindBin(range_up);
    
    int nbins = bin_up-bin_low+1;
    double *xbins = new double[nbins+1];
    
    for (int i=0; i<=nbins; ++i) {
        xbins[i] = hist->GetBinLowEdge(bin_low+i);
    }
    
    TH1D* hist_range = new TH1D(name,"",nbins,xbins);
    
    for (int i=0; i<nbins; ++i) {
        hist_range->SetBinContent(i+1,hist->GetBinContent(bin_low+i));
    }
    
    return hist_range;
    
}

// ----------------------------------------------------
// Create Template Histograms
// ----------------------------------------------------

TH1D* HPGe_Fit::CreateHistConst(TString name, int nBins, double* xbins) {
    
    TH1D* histogram = new TH1D(name, ";Energy (keV); Counts (a.u.)", nBins, xbins);
    
    for (int i=1; i<=nBins; i++) {
        histogram->SetBinContent(i,1);
    }
    
    histogram->Scale(1.0/histogram->Integral());
    
    return histogram;
}

TH1D* HPGe_Fit::CreateHistGauss(TString name, double mean, double sigma, int nBins, double* xbins) {
    
    TH1D* histogram = new TH1D(name, ";Energy (keV); Counts (a.u.)", nBins, xbins);
    
    double xBin;
    double yBin;
    
    for (int i=1; i<=nBins; i++) {
        
        xBin = histogram->GetBinCenter(i);
        yBin = TMath::Gaus(xBin, mean, sigma);
        
        histogram->SetBinContent(i,yBin);
        
    }
    
    histogram->Scale(1.0/histogram->Integral());
    
    return histogram;
}




// ----------------------------------------------------
// Get first or last empty bin
// ----------------------------------------------------

int HPGe_Fit::GetEmptyBin(TH1D* hist, TString option) {
    
    int nBins = hist->GetNbinsX();
    int bin;
    double counts;
    
    if (hist->GetEntries()==0) {
        std::cout << "###### ERROR in GetEmptyBin(): histogram is empty" << std::endl;
        return 0;
    }
    
    if (option=="first") {
        bin = 1;
        for (int i=1; i<=nBins; ++i) {
            
            counts = hist->GetBinContent(i);
            
            if (counts>0) {
                bin = i;
                break;
            }
        }
    }

   else if (option=="last") {
       bin = nBins;
       for (int i=nBins; i>0; --i) {
           
           counts = hist->GetBinContent(i);
           
           if (counts>0) {
               bin = i;
               break;
           }
       }
   }
   else {
       std::cout << "###### ERROR in GetEmptyBin(): unknown option" << std::endl;
       return 0;
   }
    
    return bin;

}

// ----------------------------------------------------
// Set new parameter limits (eliminate empty bins)
// ----------------------------------------------------

void HPGe_Fit::SetNewLimits(BCMTF_HPGe* m, TString parameter) {
    
    TH1D* hist = m->GetMarginalized(parameter)->GetHistogram();

    int nBins = hist->GetNbinsX();

    // find first non-empty bin
    int firstBin = GetEmptyBin(hist, "first");
    
    // find last non-empty bin
    int lastBin = GetEmptyBin(hist, "last");
    
    // set new parameter limits
    double limit_low, limit_up;
    
    if (firstBin>1) {
        limit_low = hist->GetBinLowEdge(firstBin-1);
    }
    else {
        limit_low = hist->GetBinLowEdge(1);
    }
    if (lastBin<nBins-1) {
        limit_up = hist->GetBinLowEdge(lastBin+2);
    }
    else {
        limit_up = hist->GetBinLowEdge(lastBin)+hist->GetBinWidth(0);
    }
    
    m->GetParameter(parameter)->SetLimits(limit_low,limit_up);
}

