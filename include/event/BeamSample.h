/*
 * This class is meant to build beam observables for fitting
 * It generates them from VALOR inputs (HK MC),
 * building and populating 4 x 1D histograms in E_reco
 *
 * The histograms are created from TH2D objects that relate
 * E_true information (needed for osc. prob.) and E_reco
 * Such histograms are stored as Eigen::MatrixXd objects
 */

#ifndef BeamSample_H
#define BeamSample_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "event/Sample.h"

class BeamSample : public Sample
{
	public:
		BeamSample(const std::string &card, std::string process = "");
		BeamSample(const CardDealer &cd, std::string process = "");
		BeamSample(CardDealer *cd, std::string process = "");
		void Init(const CardDealer &cd, std::string process = "");

		void LoadReconstruction(std::string reco_file);

		void LoadReconstruction(const CardDealer &cd) override;
		void LoadSystematics(const CardDealer &cd) override;

		std::unordered_map<std::string, Eigen::VectorXd>
			BuildSamples(std::shared_ptr<Oscillator> osc = nullptr) override;
		virtual std::unordered_map<std::string, Eigen::VectorXd>
			Unfold(const Eigen::VectorXd &En);

		Eigen::SparseMatrix<double> ScaleMatrix(Xi factor, const Eigen::VectorXd &epsil);

	private:
		std::unordered_map<std::string, Eigen::MatrixXd> _reco;
		//std::map<std::string, std::vector<double> > _binX;
		//std::map<std::string, std::vector<double> > _binY;

		//std::vector<std::string> _mode;
		//std::vector<std::string> _chan;
		//std::vector<std::string> _horn;
};

#endif
