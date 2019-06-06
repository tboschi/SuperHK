/*
 * Flux class, container of various components as root histogram
 * 
 * Author: Tommaso Boschi
 */

#ifndef Flux_H
#define Flux_H

#include <iostream>
#include <string>
#include <map>

#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"

#include "tools/CardDealer.h"

class Flux
{
	public:
		enum Type
		{
			E0,
			M0,
			EB,
			MB
		};

	public:
		Flux(std::string cardFlux);
		Flux(const Flux & f);	//copy ctor
		~Flux();

		//Clone functions, so that an object from this class
		//owns valid copies of the histograms
		template <class T>
		TH1D* Clone(T* x)
		{
			TH1D* h;
			if (x)
			{
				h = dynamic_cast<TH1D*> (x->Clone());
				h->SetDirectory(0);
			}
			else
				h = NULL;

			return h;
		}

		TH1D* Get(std::string name);
		TH1D* Get(Type t);

		void Scale(std::string name, double X);
		void Scale(Type t, double X);
		void Scale(double X);

		//void Add();
		//void Add(Hist Name);
		//bool Stretch(Hist Name, double Sx, double Ex);

		//double RangeStart();
		//double RangeEnd();
		//double BinNumber();
		//double BinWidth();

	private:
		std::map<std::string, TH1D*> mFlux;
		std::map<std::string, TH1D*>::iterator ifh;
};

#endif
