#include "Flux.h"

//ctor
Flux::Flux(std::string cardFlux)
{
	CardDealer cd(cardFlux);
	std::string fileFlux;
	cd.Get("flux_path", fileFlux);

	TFile* inFile = new TFile(fileFlux.c_str(), "READ");
	if (inFile->IsZombie())
		std::cerr << "Flux: file " << fileFlux << " does not exist" << std::endl;

	std::map<std::string, std::string> mSD;
	std::map<std::string, std::string>::iterator is;
	//look if there is a particular event
	if (cd.Get("nu_", mSD))
		for (is = mSD.begin(); is != mSD.end(); ++is)
			mFlux[is->first] = Clone(inFile->Get(is->second.c_str()));

	inFile->Close();
}

//copy ctor
Flux::Flux(const Flux & f)
{
	std::map<std::string, TH1D*>::const_iterator ich;
	for (ich = (f.mFlux).begin(); ich != (f.mFlux).end(); ++ich)
		mFlux[ich->first] = Clone(ich->second);
}

//detor
Flux::~Flux()
{
	for (ifh = mFlux.begin(); ifh != mFlux.end(); ++ifh)
		delete ifh->second;
}


//Get functions

TH1D* Flux::Get(std::string name)
{
	return mFlux[name];
}

TH1D* Flux::Get(Type t)
{
	std::string name;
	switch (t)
	{
		case Type::E0:
			name = "E0";
			break;
		case Type::M0:
			name = "M0";
			break;
		case Type::EB:
			name = "EB";
			break;
		case Type::MB:
			name = "MB";
			break;
	}

	return Get(name);
}

void Flux::Scale(std::string name, double X)
{
	if (Get(name))
		Get(name)->Scale(X);
}

void Flux::Scale(Type t, double X)
{
	if (Get(t))
		Get(t)->Scale(X);
}

void Flux::Scale(double X)
{
	for (ifh = mFlux.begin(); ifh != mFlux.end(); ++ifh)
		Scale(ifh->first, X);
}

/*
void Flux::Add()
{
	Get(Total)->Reset("ICES");

	Add(Pion);
	Add(PPion);
	Add(Kaon);
	Add(Kaon0);
	Add(Charm);
	Add(Muon);
	Add(TauE);
	Add(TauM);
}

void Flux::Add(Hist Name)
{
	TH1D* hComponent;

	if (hComponent = Get(Name))
		Get(Total)->Add(hComponent);
}
*/

/*
bool Flux::Stretch(TH1D*& hist, double Sx, double Ex)
{
	if (hist && Sx >= RangeStart() && Ex <= RangeEnd())
	{
		TH1D *hTemp = dynamic_cast<TH1D*> (hist->Clone());
		hist->Reset("ICES");

		double A = hTemp->GetXaxis()->GetBinCenter(hTemp->FindFirstBinAbove(0));
		double B = hTemp->GetXaxis()->GetBinCenter(hTemp->FindLastBinAbove(0));

		double m = (Ex - Sx) / (B - A);
		double Start = RangeStart();
		double End   = RangeEnd();
		double EnStep = (End - Start) / BinNumber();
		for (double Energy = Start; Energy < End; Energy += EnStep)
		{
			double Flux = hTemp->GetBinContent(hTemp->FindBin(Energy));
			hist->Fill( Sx + (Energy - A) * (Ex - Sx) / (B - A), Flux * (Ex - Sx) / (B - A) );
		}

		delete hTemp;
		return true;
	}
	else
		return false;
}

double Flux::RangeStart()
{
	return Get(Total)->GetBinCenter(0) + BinWidth()/2.0;
}

double Flux::RangeEnd()
{
	return Get(Total)->GetBinCenter(BinNumber()) + BinWidth()/2.0;
}

double Flux::BinNumber()
{
	return Get(Total)->GetNbinsX();
}

double Flux::BinWidth()
{
	return (Get(Total)->GetBinCenter(BinNumber()) - Get(Total)->GetBinCenter(0)) / BinNumber();
}
*/
