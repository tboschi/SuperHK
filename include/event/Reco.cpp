#include "Reco.h"

//ctor
Reco::Reco(std::string cardReco)
{
	CardDealer cd(cardReco);
	std::string fileReco;
	double weight;
	bool gt = cd.Get("reco_path", fileReco);
	if (!cd.Get("scale", weight))
		weight = 1.0;

	TFile* inFile = new TFile(fileReco.c_str(), "READ");
	if (inFile->IsZombie())
		std::cerr << "Reco: file " << fileReco << " does not exist" << std::endl;

	std::map<std::string, std::string> mSD;
	std::map<std::string, std::string>::iterator is;
	//look if there is a particular event
	if (cd.Get("Ring_", mSD))
	{
		for (is = mSD.begin(); is != mSD.end(); ++is)
		{
			mReco[is->first] = Clone(inFile->Get(is->second.c_str()));
			Scale(is->first, weight);
		}
	}

	inFile->Close();
}

//copy ctor
Reco::Reco(const Reco & r)
{
	std::map<std::string, TH2D*>::const_iterator ich;
	for (ich = (r.mReco).begin(); ich != (r.mReco).end(); ++ich)
		mReco[ich->first] = Clone(ich->second);
}

//detor
Reco::~Reco()
{
	for (irh = mReco.begin(); irh != mReco.end(); ++irh)
		delete irh->second;
}


//Get functions

TH2D* Reco::Get(std::string name)
{
	return mReco[name];
}

void Reco::Scale(std::string name, double X)
{
	if (Get(name))
		Get(name)->Scale(X);
}

void Reco::Scale(double X)
{
	for (irh = mReco.begin(); irh != mReco.end(); ++irh)
		Scale(irh->first, X);
}

void Reco::Scale(std::string name, TH1D* h, char axis)
{
	TAxis *ax0, *ax1;
	switch (axis)
	{
		case 'x':
			ax0 = Get(name)->GetXaxis();
			ax1 = Get(name)->GetYaxis();
			break;
		case 'y':
			ax0 = Get(name)->GetYaxis();
			ax1 = Get(name)->GetXaxis();
			break;
	}

	for (int ix = 1; ix < ax0->GetNbins()+1; ++ix)
	{
		//std::cout << "scale bin " << ix << "\t" << ax0->GetBinCenter(ix) << " GeV" << std::endl;
		int bin = h->FindBin(ax0->GetBinCenter(ix));
		double cn = h->GetBinContent(bin);

		//std::cout << "range " << lb << "\t" << ub << "\tavg : " << avg << std::endl;
		double cc = 0;
		for (int iy = 1; iy < ax1->GetNbins()+1; ++iy)
		{
			Get(name)->SetBinContent(ix, iy, Get(name)->GetBinContent(ix, iy) * cn);
			//cc += Get(name)->GetBinContent(ix, iy);
		}
	}
}

void Reco::Scale(TH1D* h, char axis)
{
	for (irh = mReco.begin(); irh != mReco.end(); ++irh)
		Scale(irh->first, h, axis);
}

TH1D* Reco::Project(std::string name, char axis, std::string app)
{
	if (Get(name))
	{
		std::string pname = name + app;
		switch (axis)
		{
			case 'x':
				return Get(name)->ProjectionX(pname.c_str());
			case 'y':
				return Get(name)->ProjectionY(pname.c_str());
		}
	}
	else
		return 0;
}

void Reco::SaveProjections(std::string baseName, char axis)
{
	for (irh = mReco.begin(); irh != mReco.end(); ++irh)
	{
		TH1D* hP = Project(irh->first, axis);
		std::string name = baseName + "_" + irh->first;
		if (hP)
		{
			hP->SetName(name.c_str());
			hP->Write(name.c_str(), TObject::kWriteDelete);
		}
	}
}

int Reco::BinsX(const double *&bins)
{
	TAxis *ax = mReco.begin()->second->GetXaxis();
	bins = ax->GetXbins()->GetArray();
	return ax->GetNbins();
}

int Reco::BinsY(const double *&bins)
{
	TAxis *ax = mReco.begin()->second->GetYaxis();
	bins = ax->GetXbins()->GetArray();
	return ax->GetNbins();
}
