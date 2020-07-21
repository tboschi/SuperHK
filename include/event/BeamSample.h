/*
 * BeamSample class, container of various components as 2D root histogram
 * 
 * Author: Tommaso Boschi
 */

#ifndef BeamSample_H
#define BeamSample_H

#include <iostream>
#include <string>
#include <map>

#include "tools/CardDealer.h"

#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"

#include "Eigen/Dense"

class BeamSample : Sample
{
	public:
		BeamSample(CardDealer *card);
		void Init();

		void LoadReconstruction();
		void LoadReconstruction(std::string channel);
		void DefineBinning();
		void LoadSystematics();

		Eigen::VectorXd ConstructSpectrum(Oscillator *osc = 0);
		std::map<std::string, Eigen::VectorXd> BuildSamples(Oscillator *osc = 0);


		/* old stuff */
		Eigen::MatrixXd operator()(std::string name);
		std::vector<double> BinsX(std::string name = "");
		std::vector<double> BinsY(std::string name = "");

		void Scale(double X);
		void Scale(std::string name, double X);

		Eigen::MatrixXd Apply(std::string name, TH1D* h, char axis = 'x');
		Eigen::MatrixXd Apply(std::string name, Eigen::VectorXd &vec, char axis = 'x');

	private:
		std::map<std::string, Eigen::MatrixXd> _reco;
		std::map<std::string, std::vector<double> > _binX;
		std::map<std::string, std::vector<double> > _binY;

		//std::map<std::string, TH2D*> mReco;
		//std::map<std::string, TH2D*>::iterator irh;
};

#endif
