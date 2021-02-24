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

#ifndef ATMOSAMPLE_H
#define ATMOSAMPLE_H

#include <iostream>
#include <string>
#include <set>

#include "physics/Atmosphere.h"

#include "event/Sample.h"

#include "TChain.h"

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

		std::unordered_map<std::string, Eigen::VectorXd>
			BuildSamples(std::shared_ptr<Oscillator> osc = nullptr) override;
		Eigen::VectorXd ConstructSamples(std::shared_ptr<Oscillator> osc = nullptr);
		virtual std::unordered_map<std::string, Eigen::VectorXd>
			Unfold(const Eigen::VectorXd &En);


	private:
		// atmospheric oscillation
		std::unique_ptr<Atmosphere> _atm_path;

		// binning information is stored as root histograms
		std::unordered_map<std::string, TH2D*> _reco_hist, _true_hist;
		//std::unordered_map<std::string, Eigen::MatrixXd> _bin_contents;
		std::map<int, std::string> _type_names;	 // order important!

		// For root stuff
		std::unique_ptr<TChain> dm; //, nh, ih;
		int ipnu, mode, itype;
		float dirnu[3], dir[3], flxho[3];
		float pnu, amom, weightx; //, ErmsHax, nEAveHax;
		long int _nentries;

		int point, bins;
		double data[3000];
		std::unordered_map<std::string, std::shared_ptr<TChain> > _pre_tree;
		std::unordered_map<std::string, std::map<size_t, size_t> > _pre_point;

		double _weight, _reduce;
};

//void CreateTensor(CardDealer *cd);

#endif
