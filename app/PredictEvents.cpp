#include <iostream>
#include <fstream>

#include "tools/CardDealer.h"
#include "event/Flux.h"
#include "event/Reco.h"
#include "physics/Oscillator.h"

int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	CardDealer *cd = new CardDealer(cardFile);

	//flux files
	std::string flxMode0_File, flxModeB_File;

	cd->Get("flux_0", flxMode0_File);
	cd->Get("flux_B", flxModeB_File);

	//flux objects
	Flux *flxMode0 = new Flux(flxMode0_File);	//neutrino mode flux
	Flux *flxModeB = new Flux(flxModeB_File);	//antineut mode flux

	//reco files	
	//neutrino modes
	std::string nuE0_nuE0_m0_File, nuM0_nuM0_m0_File, nuM0_nuE0_m0_File;
	std::string nuEB_nuEB_m0_File, nuMB_nuMB_m0_File, nuMB_nuEB_m0_File;
	//antineutrino modes
	std::string nuE0_nuE0_mB_File, nuM0_nuM0_mB_File, nuM0_nuE0_mB_File;
	std::string nuEB_nuEB_mB_File, nuMB_nuMB_mB_File, nuMB_nuEB_mB_File;

	cd->Get("reco_nuE0_nuE0_m0", nuE0_nuE0_m0_File);
	cd->Get("reco_nuM0_nuM0_m0", nuM0_nuM0_m0_File);
	cd->Get("reco_nuM0_nuE0_m0", nuM0_nuE0_m0_File);
	cd->Get("reco_nuEB_nuEB_m0", nuEB_nuEB_m0_File);
	cd->Get("reco_nuMB_nuMB_m0", nuMB_nuMB_m0_File);
	cd->Get("reco_nuMB_nuEB_m0", nuMB_nuEB_m0_File);

	cd->Get("reco_nuE0_nuE0_mB", nuE0_nuE0_mB_File);
	cd->Get("reco_nuM0_nuM0_mB", nuM0_nuM0_mB_File);
	cd->Get("reco_nuM0_nuE0_mB", nuM0_nuE0_mB_File);
	cd->Get("reco_nuEB_nuEB_mB", nuEB_nuEB_mB_File);
	cd->Get("reco_nuMB_nuMB_mB", nuMB_nuMB_mB_File);
	cd->Get("reco_nuMB_nuEB_mB", nuMB_nuEB_mB_File);

	//reco objects
	//neutrino mode
	//for nue -> nue, num -> num, num -> nue
	Reco *nuE0_nuE0_m0 = new Reco(nuE0_nuE0_m0_File);
	Reco *nuM0_nuM0_m0 = new Reco(nuM0_nuM0_m0_File);
	Reco *nuM0_nuE0_m0 = new Reco(nuM0_nuE0_m0_File);
	//for nueb -> nueb, numb -> numb, numb -> nueb
	Reco *nuEB_nuEB_m0 = new Reco(nuEB_nuEB_m0_File);
	Reco *nuMB_nuMB_m0 = new Reco(nuMB_nuMB_m0_File);
	Reco *nuMB_nuEB_m0 = new Reco(nuMB_nuEB_m0_File);
	//antineutrino mode
	//for nue -> nue, num -> num, num -> nue
	Reco *nuE0_nuE0_mB = new Reco(nuE0_nuE0_mB_File);
	Reco *nuM0_nuM0_mB = new Reco(nuM0_nuM0_mB_File);
	Reco *nuM0_nuE0_mB = new Reco(nuM0_nuE0_mB_File);
	//for nueb -> nueb, numb -> numb, numb -> nueb
	Reco *nuEB_nuEB_mB = new Reco(nuEB_nuEB_mB_File);
	Reco *nuMB_nuMB_mB = new Reco(nuMB_nuMB_mB_File);
	Reco *nuMB_nuEB_mB = new Reco(nuMB_nuEB_mB_File);

	//get POT and scaling factor from SK to HK
	double pot0, potB;
	double scale;

	if (!cd->Get("POT_0", pot0))
		pot0 = 1.0;
	if (!cd->Get("POT_B", potB))
		potB = 1.0;
	if (!cd->Get("scale", scale))
		scale = 1.0;

	//scale fluxes to POT and weight
	flxMode0->Scale(pot0 * scale);
	flxModeB->Scale(potB * scale);

	//set up oscillation
	std::string densityFile;
	cd->Get("density_profile", densityFile);

	double M12, M23;
	double S12, S13, S23;
	double dCP;

	cd->Get("M12", M12);
	cd->Get("M23", M23);
	cd->Get("S12", S12);
	cd->Get("S13", S13);
	cd->Get("S23", S23);
	cd->Get("dCP", dCP);

	Oscillator *osc = new Oscillator(densityFile);
	osc->SetMasses<Oscillator::normal>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP * Const::fPi);

	/*
	std::ofstream out("mmprob.dat");
	for (double en = 0; en < 10.0; en += 0.05)
	{
		out << en << "\t";
		out << osc->Probability(Nu::M_, Nu::M_, en) << "\t";		//2
		out << osc->Probability(Nu::Mb, Nu::Mb, en) << "\t";		//3
		out << osc->Probability(Nu::E_, Nu::E_, en) << "\t";		//4
		out << osc->Probability(Nu::Eb, Nu::Eb, en) << "\t";		//5
		out << osc->Probability(Nu::M_, Nu::E_, en) << "\t";		//6
		out << osc->Probability(Nu::Mb, Nu::Eb, en) << std::endl;	//7
	}
	out.close();
	*/


	/////////////////////////////////////////////

	//set to 0 histograms for NC events in appeareance channel
	nuM0_nuE0_m0->Scale("E_NC", 0.0);
	nuM0_nuE0_m0->Scale("M_NC", 0.0);
	nuMB_nuEB_m0->Scale("E_NC", 0.0);
	nuMB_nuEB_m0->Scale("M_NC", 0.0);

	nuM0_nuE0_mB->Scale("E_NC", 0.0);
	nuM0_nuE0_mB->Scale("M_NC", 0.0);
	nuMB_nuEB_mB->Scale("E_NC", 0.0);
	nuMB_nuEB_mB->Scale("M_NC", 0.0);

	//scale fluxes to oscillation
	TH1D* flx_nuE0_nuE0_m0 = flxMode0->Clone(flxMode0->Get(Flux::E0));
	TH1D* flx_nuM0_nuM0_m0 = flxMode0->Clone(flxMode0->Get(Flux::M0));
	TH1D* flx_nuM0_nuE0_m0 = flxMode0->Clone(flxMode0->Get(Flux::M0));
	TH1D* flx_nuEB_nuEB_m0 = flxMode0->Clone(flxMode0->Get(Flux::EB));
	TH1D* flx_nuMB_nuMB_m0 = flxMode0->Clone(flxMode0->Get(Flux::MB));
	TH1D* flx_nuMB_nuEB_m0 = flxMode0->Clone(flxMode0->Get(Flux::MB));

	//NC doesn't feel oscillation, so apply before
	nuE0_nuE0_m0->Scale("E_NC", flx_nuE0_nuE0_m0);
	nuE0_nuE0_m0->Scale("M_NC", flx_nuE0_nuE0_m0);
	nuM0_nuM0_m0->Scale("E_NC", flx_nuM0_nuM0_m0);
	nuM0_nuM0_m0->Scale("M_NC", flx_nuM0_nuM0_m0);
	//nuM0_nuE0_m0->Scale("E_NC", flx_nuM0_nuE0_m0);
	//nuM0_nuE0_m0->Scale("M_NC", flx_nuM0_nuE0_m0);
	nuEB_nuEB_m0->Scale("E_NC", flx_nuEB_nuEB_m0);
	nuEB_nuEB_m0->Scale("M_NC", flx_nuEB_nuEB_m0);
	nuMB_nuMB_m0->Scale("E_NC", flx_nuMB_nuMB_m0);
	nuMB_nuMB_m0->Scale("M_NC", flx_nuMB_nuMB_m0);
	//nuMB_nuEB_m0->Scale("E_NC", flx_nuMB_nuEB_m0);
	//nuMB_nuEB_m0->Scale("M_NC", flx_nuMB_nuEB_m0);

	osc->Oscillate(Nu::E_, Nu::E_, flx_nuE0_nuE0_m0);
	osc->Oscillate(Nu::M_, Nu::M_, flx_nuM0_nuM0_m0);
	osc->Oscillate(Nu::M_, Nu::E_, flx_nuM0_nuE0_m0);
	osc->Oscillate(Nu::Eb, Nu::Eb, flx_nuEB_nuEB_m0);
	osc->Oscillate(Nu::Mb, Nu::Mb, flx_nuMB_nuMB_m0);
	osc->Oscillate(Nu::Mb, Nu::Eb, flx_nuMB_nuEB_m0);

	nuE0_nuE0_m0->Scale("E_CCQE",  flx_nuE0_nuE0_m0);
	nuE0_nuE0_m0->Scale("E_CCnQE", flx_nuE0_nuE0_m0);
	nuE0_nuE0_m0->Scale("M_CCQE",  flx_nuE0_nuE0_m0);
	nuE0_nuE0_m0->Scale("M_CCnQE", flx_nuE0_nuE0_m0);
	nuM0_nuM0_m0->Scale("E_CCQE",  flx_nuM0_nuM0_m0);
	nuM0_nuM0_m0->Scale("E_CCnQE", flx_nuM0_nuM0_m0);
	nuM0_nuM0_m0->Scale("M_CCQE",  flx_nuM0_nuM0_m0);
	nuM0_nuM0_m0->Scale("M_CCnQE", flx_nuM0_nuM0_m0);
	nuM0_nuE0_m0->Scale("E_CCQE",  flx_nuM0_nuE0_m0);
	nuM0_nuE0_m0->Scale("E_CCnQE", flx_nuM0_nuE0_m0);
	nuM0_nuE0_m0->Scale("M_CCQE",  flx_nuM0_nuE0_m0);
	nuM0_nuE0_m0->Scale("M_CCnQE", flx_nuM0_nuE0_m0);
	nuEB_nuEB_m0->Scale("E_CCQE",  flx_nuEB_nuEB_m0);
	nuEB_nuEB_m0->Scale("E_CCnQE", flx_nuEB_nuEB_m0);
	nuEB_nuEB_m0->Scale("M_CCQE",  flx_nuEB_nuEB_m0);
	nuEB_nuEB_m0->Scale("M_CCnQE", flx_nuEB_nuEB_m0);
	nuMB_nuMB_m0->Scale("E_CCQE",  flx_nuMB_nuMB_m0);
	nuMB_nuMB_m0->Scale("E_CCnQE", flx_nuMB_nuMB_m0);
	nuMB_nuMB_m0->Scale("M_CCQE",  flx_nuMB_nuMB_m0);
	nuMB_nuMB_m0->Scale("M_CCnQE", flx_nuMB_nuMB_m0);
	nuMB_nuEB_m0->Scale("E_CCQE",  flx_nuMB_nuEB_m0);
	nuMB_nuEB_m0->Scale("E_CCnQE", flx_nuMB_nuEB_m0);
	nuMB_nuEB_m0->Scale("M_CCQE",  flx_nuMB_nuEB_m0);
	nuMB_nuEB_m0->Scale("M_CCnQE", flx_nuMB_nuEB_m0);
                                            
	TH1D* flx_nuE0_nuE0_mB = flxModeB->Clone(flxModeB->Get(Flux::E0));
	TH1D* flx_nuM0_nuM0_mB = flxModeB->Clone(flxModeB->Get(Flux::M0));
	TH1D* flx_nuM0_nuE0_mB = flxModeB->Clone(flxModeB->Get(Flux::M0));
	TH1D* flx_nuEB_nuEB_mB = flxModeB->Clone(flxModeB->Get(Flux::EB));
	TH1D* flx_nuMB_nuMB_mB = flxModeB->Clone(flxModeB->Get(Flux::MB));
	TH1D* flx_nuMB_nuEB_mB = flxModeB->Clone(flxModeB->Get(Flux::MB));

	//NC doesn't feel oscillation, so apply before
	nuE0_nuE0_mB->Scale("E_NC", flx_nuE0_nuE0_mB);
	nuE0_nuE0_mB->Scale("M_NC", flx_nuE0_nuE0_mB);
	nuM0_nuM0_mB->Scale("E_NC", flx_nuM0_nuM0_mB);
	nuM0_nuM0_mB->Scale("M_NC", flx_nuM0_nuM0_mB);
	//nuM0_nuE0_mB->Scale("E_NC", flx_nuM0_nuE0_mB);
	//nuM0_nuE0_mB->Scale("M_NC", flx_nuM0_nuE0_mB);
	nuEB_nuEB_mB->Scale("E_NC", flx_nuEB_nuEB_mB);
	nuEB_nuEB_mB->Scale("M_NC", flx_nuEB_nuEB_mB);
	nuMB_nuMB_mB->Scale("E_NC", flx_nuMB_nuMB_mB);
	nuMB_nuMB_mB->Scale("M_NC", flx_nuMB_nuMB_mB);
	//nuMB_nuEB_mB->Scale("E_NC", flx_nuMB_nuEB_mB);
	//nuMB_nuEB_mB->Scale("M_NC", flx_nuMB_nuEB_mB);

	osc->Oscillate(Nu::E_, Nu::E_, flx_nuE0_nuE0_mB);
	osc->Oscillate(Nu::M_, Nu::M_, flx_nuM0_nuM0_mB);
	osc->Oscillate(Nu::M_, Nu::E_, flx_nuM0_nuE0_mB);
	osc->Oscillate(Nu::Eb, Nu::Eb, flx_nuEB_nuEB_mB);
	osc->Oscillate(Nu::Mb, Nu::Mb, flx_nuMB_nuMB_mB);
	osc->Oscillate(Nu::Mb, Nu::Eb, flx_nuMB_nuEB_mB);

	nuE0_nuE0_mB->Scale("E_CCQE",  flx_nuE0_nuE0_mB);
	nuE0_nuE0_mB->Scale("E_CCnQE", flx_nuE0_nuE0_mB);
	nuE0_nuE0_mB->Scale("M_CCQE",  flx_nuE0_nuE0_mB);
	nuE0_nuE0_mB->Scale("M_CCnQE", flx_nuE0_nuE0_mB);
	nuM0_nuM0_mB->Scale("E_CCQE",  flx_nuM0_nuM0_mB);
	nuM0_nuM0_mB->Scale("E_CCnQE", flx_nuM0_nuM0_mB);
	nuM0_nuM0_mB->Scale("M_CCQE",  flx_nuM0_nuM0_mB);
	nuM0_nuM0_mB->Scale("M_CCnQE", flx_nuM0_nuM0_mB);
	nuM0_nuE0_mB->Scale("E_CCQE",  flx_nuM0_nuE0_mB);
	nuM0_nuE0_mB->Scale("E_CCnQE", flx_nuM0_nuE0_mB);
	nuM0_nuE0_mB->Scale("M_CCQE",  flx_nuM0_nuE0_mB);
	nuM0_nuE0_mB->Scale("M_CCnQE", flx_nuM0_nuE0_mB);
	nuEB_nuEB_mB->Scale("E_CCQE",  flx_nuEB_nuEB_mB);
	nuEB_nuEB_mB->Scale("E_CCnQE", flx_nuEB_nuEB_mB);
	nuEB_nuEB_mB->Scale("M_CCQE",  flx_nuEB_nuEB_mB);
	nuEB_nuEB_mB->Scale("M_CCnQE", flx_nuEB_nuEB_mB);
	nuMB_nuMB_mB->Scale("E_CCQE",  flx_nuMB_nuMB_mB);
	nuMB_nuMB_mB->Scale("E_CCnQE", flx_nuMB_nuMB_mB);
	nuMB_nuMB_mB->Scale("M_CCQE",  flx_nuMB_nuMB_mB);
	nuMB_nuMB_mB->Scale("M_CCnQE", flx_nuMB_nuMB_mB);
	nuMB_nuEB_mB->Scale("E_CCQE",  flx_nuMB_nuEB_mB);
	nuMB_nuEB_mB->Scale("E_CCnQE", flx_nuMB_nuEB_mB);
	nuMB_nuEB_mB->Scale("M_CCQE",  flx_nuMB_nuEB_mB);
	nuMB_nuEB_mB->Scale("M_CCnQE", flx_nuMB_nuEB_mB);

	std::string outFile0, outFileB;
	cd->Get("out_0", outFile0);
	cd->Get("out_B", outFileB);

	TFile *out0 = new TFile(outFile0.c_str(), "RECREATE");
	out0->cd();
	
	nuE0_nuE0_m0->SaveProjections("nuE0_nuE0", 'y');
	nuM0_nuM0_m0->SaveProjections("nuM0_nuM0", 'y');
	nuM0_nuE0_m0->SaveProjections("nuM0_nuE0", 'y');
	nuEB_nuEB_m0->SaveProjections("nuEB_nuEB", 'y');
	nuMB_nuMB_m0->SaveProjections("nuMB_nuMB", 'y');
	nuMB_nuEB_m0->SaveProjections("nuMB_nuEB", 'y');

	TFile *outB = new TFile(outFileB.c_str(), "RECREATE");
	outB->cd();

	nuE0_nuE0_mB->SaveProjections("nuE0_nuE0", 'y');
	nuM0_nuM0_mB->SaveProjections("nuM0_nuM0", 'y');
	nuM0_nuE0_mB->SaveProjections("nuM0_nuE0", 'y');
	nuEB_nuEB_mB->SaveProjections("nuEB_nuEB", 'y');
	nuMB_nuMB_mB->SaveProjections("nuMB_nuMB", 'y');
	nuMB_nuEB_mB->SaveProjections("nuMB_nuEB", 'y');

	out0->Close();
	outB->Close();

	delete cd;
	delete flxMode0, flxModeB;
	delete nuE0_nuE0_m0, nuM0_nuM0_m0, nuM0_nuE0_m0;
	delete nuEB_nuEB_m0, nuMB_nuMB_m0, nuMB_nuEB_m0;
	delete nuE0_nuE0_mB, nuM0_nuM0_mB, nuM0_nuE0_mB;
	delete nuEB_nuEB_mB, nuMB_nuMB_mB, nuMB_nuEB_mB;

	return 0;
}
