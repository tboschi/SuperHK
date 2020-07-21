/*
 * Sample base class
 * Author: Tommaso Boschi
 */

#ifndef Sample_H
#define Sample_H

#include <iostream>
#include <string>
#include <map>

#include "tools/CardDealer.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"

#include "Eigen/Dense"

class Sample
{
	public:
		Eigen::VectorXd ConstructSpectrum(Oscillator *osc = 0);

	private:
};
