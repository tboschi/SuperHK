
/*
 * Reco class, container of various components as 2D root histogram
 * 
 * Author: Tommaso Boschi
 */

#ifndef Reco_H
#define Reco_H

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

class Reco
{
	public:
		Reco(std::string cardReco);
		//Reco(const Reco & r);	//copy ctor
		//~Reco();

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
