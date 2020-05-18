#include <vector>
#include <iostream>
#include <iomanip>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"

//#define VALOR_EXPT_LOW_LEVEL_MESG_ENABLED

void ApplyLateralSystematic(TH1D * pdf, double scale)
{
  double orig_pdf_integral = pdf->Integral();
  if(orig_pdf_integral < 1E-6) {
    std::cout << "No need to apply lateral parameter to empty PDF " << pdf->GetName() << std::endl;
    return;
  }

  const unsigned int nbins = pdf->GetNbinsX();
  //const unsigned int n_overflow_bins = 2;

  //get original bin edges
  std::vector<double>            * unscaled_bin_edges         = new std::vector<double>(nbins + 1);
  for(unsigned int i = 0; i <= nbins; i++) {
    unscaled_bin_edges->at(i) = pdf->GetXaxis()->GetBinUpEdge(i);
#ifdef VALOR_EXPT_LOW_LEVEL_MESG_ENABLED
    std::cout << "Bin edge " << i << " value " << unscaled_bin_edges->at(i) << std::endl;
#endif
  }

  //get original bin content
  std::vector<double> pdfArr(nbins);
  for(unsigned int i = 0; i < nbins; i++) {
    pdfArr.at(i) = pdf->GetBinContent(i+1);
  }

  // Scale the bin edges for this kine var. Assuming uniform distribution of events within each bin.
  std::vector<double> scaled_bin_edges = *unscaled_bin_edges;
  for(double & edge : scaled_bin_edges) {
    edge *= scale;
  }

#ifdef VALOR_EXPT_LOW_LEVEL_MESG_ENABLED
  // Copy the un-modified PDF contents into a vector, so we can print comparisons later
  std::vector<double> old_events_per_bin (nbins, 0.0);
  for(unsigned int ireco = 0; ireco < nbins; ireco++) {
    old_events_per_bin.at(ireco) = pdfArr.at(ireco);
  }
#endif

  // Calculate the new number of events per bin by seeing how much overlap there is between scaled and unscaled bins.
  // Note: We have underflow and overflow bins here. This is because as we scale the bin edges, events can drop off
  // either end of the scaled binning. These under/overflow bins catch these events, which are then added back on to
  // the first or last bin at the end, respectively.
  std::vector<double> new_events_per_bin (nbins, 0.0);
  double factor = 0.0;
  for (size_t iedge_s = 0 ; iedge_s < scaled_bin_edges.size() - 1 ; iedge_s++) { // Loop over scaled edges
    const double & bin_min_s = scaled_bin_edges[iedge_s]  ;
    const double & bin_max_s = scaled_bin_edges[iedge_s+1];

    int ibin_s = iedge_s;

    for (size_t iedge_u = 0 ; iedge_u < (*unscaled_bin_edges).size() - 1 ; iedge_u++) {  // Loop over unscaled bins
      const double & bin_min_u = (*unscaled_bin_edges)[iedge_u]  ;
      const double & bin_max_u = (*unscaled_bin_edges)[iedge_u+1];

      int ibin_u = iedge_u;

      //std::cout << iedge_u << "\t" << bin_min_u << "\t" << bin_max_u << std::endl;

      if(bin_min_u > bin_max_s)
      {  // Unscaled bin entirely larger than scaled bin, stop checking higher unscaled
        break;
      }
      else if(bin_min_s > bin_max_u)
      {  // Scaled bin entirely larger than unscaled bin, skip to next unscaled bin
        continue;
      }
      else if(bin_min_u >= bin_min_s && bin_max_u <= bin_max_s)
      {  // Unscaled bin subset of scaled bin
        factor = (bin_max_u - bin_min_u) / (bin_max_s - bin_min_s);
      }
      else if(bin_min_s >= bin_min_u && bin_max_s <= bin_max_u)
      {  // Scaled bin subset of unscaled bin
        // Take the contents of this bin unchanged, there could be additional scaled bins contributing
        factor = 1.0;
      }
      else if(bin_max_u > bin_max_s && bin_min_u >= bin_min_s)
      {  // Scaled bin overlaps low part of unscaled bin
        factor = (bin_max_s - bin_min_u) / (bin_max_s - bin_min_s);
      }
      else if(bin_min_u <  bin_min_s && bin_max_u <= bin_max_s)
      {  // Scaled bin overlaps high part of unscaled bin
        factor = (bin_max_u - bin_min_s) / (bin_max_s - bin_min_s);
      }
      else
      {  // Shouldn't get here
        factor = 0.0;
      }

      if(factor <= 0.0) {
	std::cerr << "Got a factor <= 0 for ibin_s = " << ibin_s << ", ibin_u = "
		  << ibin_u << ". This shouldn't happen!" << std::endl;
        continue;
      }

      new_events_per_bin[ibin_u] += pdfArr[ibin_s] * factor;

#ifdef VALOR_EXPT_LOW_LEVEL_MESG_ENABLED
      std::cout
        << "scaled bin "    << ibin_s << " [" << bin_min_s << ", " << bin_max_s << "]; "
        << "un-scaled bin " << ibin_u << " [" << bin_min_u << ", " << bin_max_u << "]; "
        << "factor = " << factor << "; old Nevents = " << pdfArr[ibin_s] << "; new Nevents = " << new_events_per_bin[ibin_u]
	<< std::endl;
#endif

    }//iedge_u
  }//iedge_s

  //add the underflow/overflow
  //don't add an underflow, because the first couple of bins are actually empty
  if(scale > 1) { //overflow happens
    const int iedge_s = nbins - 1;
    const int ibin_s = iedge_s;
    const double & bin_min_s = scaled_bin_edges[iedge_s]  ;
    const double & bin_max_s = scaled_bin_edges[iedge_s+1];

    const int iedge_u = nbins - 1;
    const int ibin_u = iedge_u;
    const double & bin_min_u = (*unscaled_bin_edges)[iedge_u]  ;
    const double & bin_max_u = (*unscaled_bin_edges)[iedge_u+1];

    // "1-" because this is the bit that is above the maximum (e.g. 30 GeV for 1Rmu)
    // We've already counted the other bit
    factor = 1 - (bin_max_u - bin_min_s) / (bin_max_s - bin_min_s);
    new_events_per_bin[ibin_u] += pdfArr[ibin_s] * factor;
  }//overflow

  // Modify the PDF and print debug info
#ifdef VALOR_EXPT_LOW_LEVEL_MESG_ENABLED
  std::cout << "PDF contents before and after applying lateral systematic:" << std::endl;
#endif
  for (unsigned int ibin = 0 ; ibin < nbins ; ibin++) {  // Loop over global bins
    pdfArr[ibin] = new_events_per_bin.at(ibin);
  }

#ifdef VALOR_EXPT_LOW_LEVEL_MESG_ENABLED
  for (unsigned int ibin = 0 ; ibin < nbins ; ibin++) {  // Loop over global bins
    std::cout
      << "bin " << ibin << " [" << unscaled_bin_edges->at(ibin) << ", " << unscaled_bin_edges->at(ibin + 1) << "]"
      << "; old content = " << old_events_per_bin[ibin] << "; new content = " << pdfArr[ibin] << std::endl;
  }//ibin
#endif

  //set the content of the output histogram
  for(unsigned int i = 0; i < nbins; i++) {
    pdf->SetBinContent(i + 1, pdfArr[i]);
  }

  // Lateral systematics should just migrate events between bins, while keeping the total number of events constant.
  // Check that we haven't changed the total number of events
  if ( fabs((pdf->Integral()/orig_pdf_integral)-1) > 1E-6 ) {
    std::cerr << std::fixed << std::setprecision(9)
	      << "Failed normalisation check after applying lateral systematic for PDF"
	      << pdf->GetName() << ". Original integral = " << orig_pdf_integral
	      << ", new integral = " << pdf->Integral() << std::endl;
    exit(-1);
  }
}

void escale(const double scale=0.024,
	    const char * fname_tweaked = "sigma_var_mode_SKDetFSI_Ereco.root")
{
  TFile f_nominal("write_hk_10yr_sb_run19migration.root");
  TH1D * h_nominal;
  f_nominal.GetObject("valor_pdfts_ob_rf_superk__JPARC_FHC__1rmu_Ereco", h_nominal);
  
  TFile f_tweaked(fname_tweaked);
  TH1D * h_tweaked;

  TCanvas c;
  c.Print("escale.pdf[");

  std::cout << h_nominal << std::endl;
  std::cout << "Integral " << h_nominal->Integral() << std::endl;

  const int nruns = 4;
  int runs[nruns] = {+3, +1, -1, -3};
  string runs_s[nruns] = {"p3", "p1", "m1", "m3"};
  for(int ir = 0; ir < nruns; ir++) {
    TH1D * h_to_tweak = new TH1D(*h_nominal);
    ApplyLateralSystematic(h_to_tweak, 1 + scale * runs[ir]);
    //compare with the internally generated tweaked histogram
    h_to_tweak->SetLineColor(kRed);
    h_to_tweak->Draw();
    f_tweaked.GetObject(TString::Format("env_sb_superk__JPARC_FHC__1rmu_sk_e_scale_%s", runs_s[ir].c_str()), h_tweaked);
    h_tweaked->SetLineStyle(kDashed);
    h_tweaked->DrawCopy("SAME");
    c.Print("escale.pdf");

    for(int i = 1; i <= h_nominal->GetNbinsX(); i++) {
      if(abs(h_to_tweak->GetBinContent(i) - h_tweaked->GetBinContent(i)) > 1E-6)
	std::cout << i << "\t" << h_nominal->GetBinContent(i)
		  << "\t" << h_to_tweak->GetBinContent(i) << "\t"
		  << h_tweaked->GetBinContent(i) << std::endl;
  }

  }
  c.Print("escale.pdf]");

}
