#include "EventStacker.h"

EventStacker::EventStacker(CardDealer *& binCard) :
	kVerbosity(Utils::verbosity)
{
	InitMap();	//with eventtypes
	CreateStacks(binCard);
}

EventStacker::~EventStacker()
{
	for (ims = mStack.begin(); ims != mStack.end(); ++ims)
		delete ims->second;
}

void EventStacker::InitMap()
{
	mEvent.clear();

	mEvent[Event::undefined]			= "undefined";
	mEvent[Event::SubGeV_elike_0dcy]		= "SubGeV_elike_0dcy";
	mEvent[Event::SubGeV_elike_1dcy]		= "SubGeV_elike_1dcy";
	mEvent[Event::SubGeV_SingleRing_pi0like]	= "SubGeV_SingleRing_pi0like";
	mEvent[Event::SubGeV_mulike_0dcy]		= "SubGeV_mulike_0dcy";
	mEvent[Event::SubGeV_mulike_1dcy]		= "SubGeV_mulike_1dcy";
	mEvent[Event::SubGeV_mulike_2dcy]		= "SubGeV_mulike_2dcy";
	mEvent[Event::SubGeV_pi0like]			= "SubGeV_pi0like";
	mEvent[Event::MultiGeV_elike_nue]		= "MultiGeV_elike_nue";
	mEvent[Event::MultiGeV_elike_nuebar]		= "MultiGeV_elike_nuebar";
	mEvent[Event::MultiGeV_mulike]			= "MultiGeV_mulike";
	mEvent[Event::MultiRing_elike_nue]		= "MultiRing_elike_nue";
	mEvent[Event::MultiRing_elike_nuebar]		= "MultiRing_elike_nuebar";
	mEvent[Event::MultiRing_mulike]			= "MultiRing_mulike";
	mEvent[Event::MultiRingOther_1]			= "MultiRingOther_1";
	mEvent[Event::PCStop]				= "PCStop";
	mEvent[Event::PCThru]				= "PCThru";
	mEvent[Event::UpStop_mu]			= "UpStop_mu";
	mEvent[Event::UpThruNonShower_mu]		= "UpThruNonShower_mu";
	mEvent[Event::UpThruShower_mu]			= "UpThruShower_mu";

	mEvent[Event::T2Knumu]				= "T2Knumu";
	mEvent[Event::T2Knumu1]				= "T2Knumu1";
	mEvent[Event::T2Knumu2]				= "T2Knumu2";
	mEvent[Event::T2Knumu3]				= "T2Knumu3";
	mEvent[Event::T2Knumu4]				= "T2Knumu4";
	mEvent[Event::T2Knumu5]				= "T2Knumu5";
	mEvent[Event::T2Knue]				= "T2Knue";
	mEvent[Event::T2Knue1]				= "T2Knue1";
	mEvent[Event::T2Knue2]				= "T2Knue2";
	mEvent[Event::T2Knue3]				= "T2Knue3";
	mEvent[Event::T2Knue4]				= "T2Knue4";
	mEvent[Event::T2Knue5]				= "T2Knue5";
	mEvent[Event::T2Knumubar]			= "T2Knumubar";
	mEvent[Event::T2Knumubar1]			= "T2Knumubar1";
	mEvent[Event::T2Knumubar2]			= "T2Knumubar2";
	mEvent[Event::T2Knumubar3]			= "T2Knumubar3";
	mEvent[Event::T2Knumubar4]			= "T2Knumubar4";
	mEvent[Event::T2Knumubar5]			= "T2Knumubar5";
	mEvent[Event::T2Knuebar]			= "T2Knuebar";
	mEvent[Event::T2Knuebar1]			= "T2Knuebar1";
	mEvent[Event::T2Knuebar2]			= "T2Knuebar2";
	mEvent[Event::T2Knuebar3]			= "T2Knuebar3";
	mEvent[Event::T2Knuebar4]			= "T2Knuebar4";
	mEvent[Event::T2Knuebar5]			= "T2Knuebar5";

	mEvent[Event::MINOSnumu]			= "MINOSnumu";
	mEvent[Event::MINOSnue]				= "MINOSnue";
}

void EventStacker::CreateStacks(CardDealer * &binCard)
{
	for (ime = mEvent.begin(); ime != mEvent.end(); ++ime)
	{
		std::map<std::string, std::vector<double> > mSD;
		std::map<std::string, std::vector<double> >::iterator is;
		//look if there is a particular event
		if (binCard->Get(ime->second, mSD))
		{
			std::map<char, std::vector<double> > mAx;
			std::map<char, std::vector<double> >::iterator ix;
			//mStringDouble has n entries, string Name_{X,Y} double {bins}
			while (mSD.size())
			{	
				//pick first and look for entries with same Name
				is = mSD.begin();
				std::string name = is->first.substr(0, is->first.find_first_of('_'));

				//char axis = is->first.at(is->first.find_first_of('_') + 1);
				char axis = is->first.at(name.length()+1);
				if (!mAx.count(axis))
					mAx[axis] = is->second;

				//saving axes with same Name	start from next iterator
				for (++is; is != mSD.end(); ++is)
				{
					std::string _name = is->first.substr(0, is->first.find_first_of('_'));

					if (name == _name)	//same name, so same stack
					{
						//axis = is->first.at(is->first.find_first_of('_') + 1);
						axis = is->first.at(name.length()+1);
						if (!mAx.count(axis))
							mAx[axis] = is->second;
						mSD.erase(is);		//remove, to speed up loop later
					}
				}
				CreateStack(name, ime->first, mAx);	//create T1H or T2H or T3H
			}
		}
	}
}

void EventStacker::CreateStack(std::string name, Event::Type et, std::map<char, std::vector<double> > &mAx)
{
	std::string stackName = mEvent[et] + "_" + name;
	if (!mStack.count(et))
	{
		double *xAx, *yAx, *zAx;
		int xBins, yBins, zBins;	//this inlcudes overflow

		int nDim = mAx.size();
		switch (nDim)
		{
			case 1:
				xAx = &(mAx['X'][0]);		//x-axis array
				xBins = mAx['X'].size() - 1;
				mStack[et] = new TH1D(stackName.c_str(), stackName.c_str(),
							xBins, xAx);
				break;
			case 2:
				xAx = &(mAx['X'][0]);		//x-axis array
				xBins = mAx['X'].size() - 1;
				yAx = &(mAx['Y'][0]);		//y-axis array	
				yBins = mAx['Y'].size() - 1;
				mStack[et] = new TH2D(stackName.c_str(), stackName.c_str(),
							xBins, xAx,
							yBins, yAx);
				break;
			case 3:
				xAx = &(mAx['X'][0]);		//x-axis array
				xBins = mAx['X'].size() - 1;
				yAx = &(mAx['Y'][0]);		//y-axis array	
				yBins = mAx['Y'].size() - 1;
				zAx = &(mAx['Z'][0]);		//z-axis array	
				zBins = mAx['Z'].size() - 1;
				mStack[et] = new TH3D(stackName.c_str(), stackName.c_str(),
							xBins, xAx,
							yBins, yAx,
							zBins, zAx);
				break;
			default:
				if (kVerbosity)
				{
					std::cerr << "EventStacker::CreateStack dimension " << nDim << " too high" << std::endl;
					std::cerr << "                          maximum dimension: 3" << std::endl;
				}
				break;
		}
	}
	else if (kVerbosity)
		std::cerr << "EventStacker::CreateStack name \"" << stackName << "\" unknown" << std::endl;
}

//idea: use weight to distinguish between dt and mc
void EventStacker::Stack(Event::Type et, std::vector<double> &vars, double w)
{
	if (vars.size() == mStack[et]->GetDimension())
	{
		double E = vars.front();		//fron vars is energy
		int bin = FillStack(et, vars, w);

		//pushing back a triplet (row, col, val) to correct stack. Increase column value
		trEnergy[et].push_back(Eigen::Triplet<double>(bin, mColumnIdx[et], E));
		trWeight[et].push_back(Eigen::Triplet<double>(bin, mColumnIdx[et], w));

		++mColumnIdx[et];
	}
	else if (kVerbosity)
	{
		std::cerr << "EventStacker::Stack variables passed: " << vars.size() << " do not match" << std::endl;
		std::cerr << "                    Stack dimension: " <<  mStack[et]->GetDimension() << std::endl;
	}
}

int EventStacker::FillStack(Event::Type et, std::vector<double> &vars, double w)
{
	int bin;
	switch (mStack[et]->GetDimension())
	{
		case 1:
			bin = mStack[et]->Fill(vars.at(1),
					       w);
		case 2:
			bin = (dynamic_cast<TH2*>(mStack[et]))->Fill(vars.at(1),
								     vars.at(2),
								     w);
		case 3:
			bin = (dynamic_cast<TH3*>(mStack[et]))->Fill(vars.at(1),
								     vars.at(2),
								     vars.at(3),
								     w);
	}

	return bin;
}

void EventStacker::ChangeStacks(Oscillator *&l)
{
	for (ims = mStack.begin(); ims != mStack.end(); ++ims)
		ChangeWeights(ims->first, l);
}

void EventStacker::ResetStacks()
{
	for (ims = mStack.begin(); ims != mStack.end(); ++ims)
		ResetWeights(ims->first);
}

void EventStacker::ChangeWeights(Event::Type et, Oscillator *&l)
{
	int bins = mStack[et]->GetNbinsX() *
			    mStack[et]->GetNbinsY() *
			    mStack[et]->GetNbinsZ();
	int cols = mColumnIdx[et];

	Eigen::SparseMatrix<double> spWeight(bins, cols);
	Eigen::SparseMatrix<double> spEnergy(bins, cols);

	spWeight.setFromTriplets(trWeight[et].begin(), trWeight[et].end());
	spEnergy.setFromTriplets(trEnergy[et].begin(), trEnergy[et].end());

	//update values
	for (int k = 0; k < spEnergy.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(spEnergy, k); it; ++it)
			it.valueRef() = l->Probability(Nu::E_, Nu::E_, it.value());

	Eigen::VectorXd noOsc = spWeight			* Eigen::VectorXd::Ones(bins);
	Eigen::VectorXd upOsc = spWeight.cwiseProduct(spEnergy) * Eigen::VectorXd::Ones(bins);

	for (int b = 0; b < bins; ++b)
		mStack[et]->SetBinContent(b, upOsc(b)/noOsc(b));
}

void EventStacker::ResetWeights(Event::Type et)
{
	int bins = mStack[et]->GetNbinsX() *
			    mStack[et]->GetNbinsY() *
			    mStack[et]->GetNbinsZ();
	int cols = mColumnIdx[et];

	Eigen::SparseMatrix<double> spWeight(bins, cols);

	spWeight.setFromTriplets(trWeight[et].begin(), trWeight[et].end());
	Eigen::VectorXd noOsc = spWeight * Eigen::VectorXd::Ones(bins);

	for (int b = 0; b < bins; ++b)
		mStack[et]->SetBinContent(b, noOsc(b));
}

int EventStacker::Bins()
{
	int ret = 0;
	for (ims = mStack.begin(); ims != mStack.end(); ++ims)
		ret += Bins(ims->first);
}

int EventStacker::Bins(Event::Type et)
{
	return mStack[et]->GetNbinsX() *
	       mStack[et]->GetNbinsY() *
	       mStack[et]->GetNbinsZ();
}

//this function return bin content at bin b of stack et
//
double EventStacker::GetBin(Event::Type et, int b) 
{
	mStack[et]->GetBinContent(b+1);
}

//create a vector with all bins, ordered by map internal order
//
std::vector<double> EventStacker::Vectorise()
{
	std::vector<double> vRet;
	for (ims = mStack.begin(); ims != mStack.end(); ++ims)
		for (int bin = 0; bin < Bins(ims->first); ++bin)
			vRet.push_back(GetBin(ims->first, bin));

	return vRet;
}
