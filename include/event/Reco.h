
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

class Reco
{
	public:
		enum Type
		{
			E0_E0,
			M0_M0,
			M0_E0,
			EB_EB,
			MB_MB,
			MB_EB
		};

		enum Ring
		{
			Elike,
			Mlike,
		};

		enum Interaction
		{
			CCQE,
			CCnQE,
			NC,
		};

		Reco(std::string cardReco);
		Reco(const Reco & r);	//copy ctor
		~Reco();

		//Clone functions, so that an object from this class
		//owns valid copies of the histograms
		template <class T>
		TH2D* Clone(T* x)
		{
			TH2D* h;
			if (x)
			{
				h = dynamic_cast<TH2D*> (x->Clone());
				h->SetDirectory(0);
			}
			else
				h = NULL;

			return h;
		}

		TH2D* Get(std::string name);

		void Scale(std::string name, double X);
		void Scale(double X);
		void Scale(std::string name, TH1D* h, char axis = 'x');
		void Scale(TH1D* h, char axis = 'x');

		TH1D* Project(std::string name, char axis = 'x', std::string app = "");
		void SaveProjections(std::string baseName, char axis = 'x');

		int BinsX(const double *&bins);
		int BinsY(const double *&bins);

	private:
		std::map<std::string, TH2D*> mReco;
		std::map<std::string, TH2D*>::iterator irh;
};

#endif
