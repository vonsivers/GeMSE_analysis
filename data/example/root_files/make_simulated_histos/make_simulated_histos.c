// initialize random number generator
int seed = 21430;
TRandom3* fRandom = new TRandom3(seed);


void make_simulated_histos()
{

    double t_bck = 4.e6; // measurement time in sec (background)
    double t_sample = 2.e6; // measurement time in sec (sample)
    double activity_bck = 1.e-3; // activity in Bq (bck)
    double activity_sample = 1.e-2; // activity in Bq (sample)
    
    TH1D* hist_bck = make_histo("simulated_histo_bck",activity_bck,t_bck);
    TH1D* hist_sig = make_histo("simulated_histo_sig",activity_sample+activity_bck,t_sample);
    
    hist_bck->Sumw2();
    hist_sig->Sumw2();
    hist_bck->Scale(1./t_bck);
    hist_sig->Scale(1./t_sample);
    
    TCanvas* c2 = new TCanvas("c2");
    gStyle->SetOptStat(0);
    c2->SetLogy();
    hist_bck->SetLineColor(1);
    hist_bck->SetMarkerColor(1);
    hist_bck->GetYaxis()->SetTitle("Rate (Hz)");
    hist_sig->SetLineColor(2);
    hist_sig->SetMarkerColor(2);
    
    hist_sig->Draw();
    hist_bck->Draw("same");
    
    c2->SaveAs("simulated_histos.pdf");
    c2->SaveAs("simulated_histos.root");
    
}
    // ---------------------------------------------------------
    
TH1D* make_histo(TString resultsname, double activity, double t_meas)
{
    // open histo
    TFile* file = TFile::Open("summed_histogram_convoluted.root");
    TCanvas* c = (TCanvas*) file->Get("c1");
    TH1D* hist_conv = (TH1D*) c->GetPrimitive("hist_conv");
    hist_conv->SetDirectory(0);
    file->Close();
    
    int nBins = hist_conv->GetNbinsX();
    double Bin_min = hist_conv->GetBinLowEdge(1);
    double Bin_max = hist_conv->GetBinLowEdge(nBins)+hist_conv->GetBinWidth(1);

    // scale histogram
    hist_conv->Scale(t_meas*activity);
    
    // new histogram
    TH1D* hist = new TH1D("hist","Simulated Spectrum;Energy (keV);Counts",nBins,Bin_min,Bin_max);
    
    // loop over all bins
    for (int i=1; i<=nBins; ++i) {
        
        int counts = fRandom->Poisson(hist_conv->GetBinContent(i));
        hist->SetBinContent(i,counts);
        
            }
    
    TVectorD v_live(1);
    TVectorD v_real(1);
    v_live[0] = t_meas;
    v_real[0] = t_meas;
    
    // draw histo
    TCanvas* c1 = new TCanvas("c1");
    gStyle->SetOptStat(0);
    c1->SetLogy();
    hist->Draw();
    c1->SaveAs(resultsname+".pdf");

    // write histo to file
    TFile* file_result = new TFile(resultsname+".root","recreate");
    file_result->cd();
    hist->Write();
    v_live.Write("t_live");
    v_real.Write("t_real");
    file_result->Close();
    
    return hist;
}
    







