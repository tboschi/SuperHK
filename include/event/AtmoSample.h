/*
 * This class is meant to build atmospheric observables for fitting
 * It generates them from SK MC, building and populating 19 x 2D histograms
 * in log p and cos zenith of the outgoing lepton
 *
 * Since the SK MC files are massive, this class serves two purpose
 *   1) creates 4D tensors in log Enu, cos zNu  ,  log p,  cos z p
 *      which are then flattened out to 2D. This matrix can be used to 
 *      convert from true neutrino information (necessary for osc prob.) to
 *      reconstructed lepton information.
 *      The advantage is that these 4D tensors must be calculated once!
 *   2) using the 4D tensors the class builds the atmospheric observables
 *      given a combination of oscillation parameters
 *
 * The tensors are saved as two T2HD object for easy binning and one Eigen::MatrixXd
 * that stores the content of the tensor, already flattened out.
 * The global index of each of the T2HD object is related to rows and cols of the matrix object
 * 
 */

#ifndef AtmoSample_H
#define AtmoSample_H

#include <iostream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <memory>

#include "tools/CardDealer.h"
#include "physics/Const.h"
#include "physics/Flavours.h"
#include "physics/Oscillator.h"
#include "physics/Atmosphere.h"

#include "event/Sample.h"

#include "TObject.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TAxis.h"
#include "TMatrixD.h"
#include "TMatrixT.h"

#include "Eigen/Dense"

class AtmoSample : public Sample
{
	public:
		AtmoSample(const std::string card, std::string process = "");
		AtmoSample(const CardDealer &cd, std::string process = "");
		AtmoSample(CardDealer *cd, std::string process = "");
		void Init(const CardDealer &cd, std::string process = "");

		void LoadSimulation();

		//void LoadReconstruction(std::string channel);

		void LoadReconstruction(const CardDealer &cd) override;
		void LoadSystematics(const CardDealer &cd) override;

		std::map<std::string, Eigen::VectorXd>
			BuildSamples(std::shared_ptr<Oscillator> osc = nullptr) override;

		Eigen::VectorXd ConstructSamples(std::shared_ptr<Oscillator> osc = nullptr);


	private:
		// atmospheric oscillation
		std::unique_ptr<Atmosphere> atm_path;

		// binning information is stored as root histograms
		std::map<std::string, TH2D*> _reco_hist, _true_hist;
		std::map<std::string, Eigen::MatrixXd> _bin_contents;
		std::map<int, std::string> _type_names;

		// For root stuff
		std::unique_ptr<TChain> dm, nh, ih;
		int ipnu, mode, itype;
		float dirnu[3], dir[3], flxho[3];
		float pnu, amom, weightx, ErmsHax, nEAveHax;
		int point, bins;
		double data[3000];

		std::map<int, int> pre_point_NH, pre_point_IH;

		double _weight, _reduce;
};

void CreateTensor(CardDealer *cd);

#endif

