#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "physics/Oscillator.h"

#include "TTree.h"
#include "TFile.h"

int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	CardDealer *cd = new CardDealer(cardFile);

	//std::string files;
	//cd->Get("files", files);

	std::vector<double> vValue, vSigma;
	std::vector<std::string> vParm;
	std::map<std::string, std::vector<double> > mvals;
	std::map<std::string, std::vector<double> >::iterator im;

	if (cd->Get("p_", mvals))
	{
		for (im = mvals.begin(); im != mvals.end(); ++im)
		{
			vParm.push_back(im->first);
			vValue.push_back(im->second.front());
			vSigma.push_back(im->second.back());
		}
	}

	if (argc < 3)
	{
		std::cerr << "No files to analyse!\n";
		return 1;
	}

	//std::string cmd = "ls " + files + " > .tmp_list";
	//system(cmd.c_str());

	//std::ifstream in(".tmp_list");
	//std::string name;

	double Epsilons[1000];
	double Errors[1000];
	double X2;
	double OriX2;
	double Time;
	double SysX2;
	double OriSysX2;
	int    Point;
	int    TPoint;
	double CP;
	double TCP;
	double M12;
	double TM12;
	double M23;
	double TM23;
	double S12;
	double TS12;
	double S13;
	double TS13;
	double S23;
	double TS23;

	std::string output = std::string(argv[2]) + "/.tmp_penalise";
	std::string cmd = "find " + std::string(argv[2]) +
			  " -name \"SpaghettiSens.T2HK.*.root\" > " + output;
	std::cout << cmd << std::endl;
	system(cmd.c_str());

	std::string file;
	std::ifstream listExclusion(output.c_str());
	while (std::getline(listExclusion, file))
	//for (int f = 2; f < argc; ++f)
	{
		TFile *inf = new TFile(file.c_str(), "READ");
		if (inf->IsZombie())
		{
			std::cout << "file doesn't exist\n";
			continue;
		}

		TTree* axis    = static_cast<TTree*>(inf->Get("axisTree"));
		TTree* error   = static_cast<TTree*>(inf->Get("errorTree"));
		TTree* sX2_old = static_cast<TTree*>(inf->Get("stepX2Tree"));

		std::cout << "Penalising " << file;
		size_t pos = file.find("SpaghettiSens.");
		if (pos != std::string::npos)
			file.insert(file.find_first_of('.', pos+13), "_penalised"); 
		std::cout << " into " << file << std::endl;

		//file = file.insert(file.find(".0"), "_penalised");
		TFile *out = new TFile(file.c_str(), "RECREATE");

		TTree *sX2_new = sX2_old->CloneTree(0);

		//stepX2Tree *sX2_old = new stepX2Tree(stepX2);
		//stepX2Tree *sX2_new = new stepX2Tree();

		if (sX2_old->GetBranch("Epsilons"))
		{
			sX2_old->SetBranchAddress("Epsilons", Epsilons);
			sX2_new->SetBranchAddress("Epsilons", Epsilons);
		}
		if (sX2_old->GetBranch("Errors"))
		{
			sX2_old->SetBranchAddress("Errors", Errors);
			sX2_new->SetBranchAddress("Errors", Errors);
		}
		if (sX2_old->GetBranch("X2"))
		{
			sX2_old->SetBranchAddress("X2", &X2);
			sX2_new->SetBranchAddress("X2", &X2);
		}
		if (sX2_old->GetBranch("OriX2"))
		{
			sX2_old->SetBranchAddress("OriX2", &OriX2);
			sX2_new->SetBranchAddress("OriX2", &OriX2);
		}
		if (sX2_old->GetBranch("Time"))
		{
			sX2_old->SetBranchAddress("Time", &Time);
			sX2_new->SetBranchAddress("Time", &Time);
		}
		if (sX2_old->GetBranch("SysX2"))
		{
			sX2_old->SetBranchAddress("SysX2", &SysX2);
			sX2_new->SetBranchAddress("SysX2", &SysX2);
		}
		if (sX2_old->GetBranch("OriSysX2"))
		{
			sX2_old->SetBranchAddress("OriSysX2", &OriSysX2);
			sX2_new->SetBranchAddress("OriSysX2", &OriSysX2);
		}
		if (sX2_old->GetBranch("Point"))
		{
			sX2_old->SetBranchAddress("Point", &Point);
			sX2_new->SetBranchAddress("Point", &Point);
		}
		if (sX2_old->GetBranch("TPoint"))
		{
			sX2_old->SetBranchAddress("TPoint", &TPoint);
			sX2_new->SetBranchAddress("TPoint", &TPoint);
		}
		if (sX2_old->GetBranch("CP"))
		{
			sX2_old->SetBranchAddress("CP",   &CP);
			sX2_new->SetBranchAddress("CP",   &CP);
		}
		if (sX2_old->GetBranch("TCP"))
		{
			sX2_old->SetBranchAddress("TCP",  &TCP);
			sX2_new->SetBranchAddress("TCP",  &TCP);
		}
		if (sX2_old->GetBranch("M12"))
		{
			sX2_old->SetBranchAddress("M12",  &M12);
			sX2_new->SetBranchAddress("M12",  &M12);
		}
		if (sX2_old->GetBranch("TM12"))
		{
			sX2_old->SetBranchAddress("TM12", &TM12);
			sX2_new->SetBranchAddress("TM12", &TM12);
		}
		if (sX2_old->GetBranch("M23"))
		{
			sX2_old->SetBranchAddress("M23",  &M23);
			sX2_new->SetBranchAddress("M23",  &M23);
		}
		if (sX2_old->GetBranch("TM23"))
		{
			sX2_old->SetBranchAddress("TM23", &TM23);
			sX2_new->SetBranchAddress("TM23", &TM23);
		}
		if (sX2_old->GetBranch("S12"))
		{
			sX2_old->SetBranchAddress("S12",  &S12);
			sX2_new->SetBranchAddress("S12",  &S12);
		}
		if (sX2_old->GetBranch("TS12"))
		{
			sX2_old->SetBranchAddress("TS12", &TS12);
			sX2_new->SetBranchAddress("TS12", &TS12);
		}
		if (sX2_old->GetBranch("S13"))
		{
			sX2_old->SetBranchAddress("S13",  &S13);
			sX2_new->SetBranchAddress("S13",  &S13);
		}
		if (sX2_old->GetBranch("TS13"))
		{
			sX2_old->SetBranchAddress("TS13", &TS13);
			sX2_new->SetBranchAddress("TS13", &TS13);
		}
		if (sX2_old->GetBranch("S23"))
		{
			sX2_old->SetBranchAddress("S23",  &S23);
			sX2_new->SetBranchAddress("S23",  &S23);
		}
		if (sX2_old->GetBranch("TS23"))
		{
			sX2_old->SetBranchAddress("TS23", &TS23);
			sX2_new->SetBranchAddress("TS23", &TS23);
		}


		for (int i = 0; i < sX2_old->GetEntries(); ++i)
		{
			sX2_old->GetEntry(i);
			double penalty = 0;
			for (size_t p = 0; p < vParm.size(); ++p)
			{
				double nominal = 0;
				if (vParm[p] == "CP")
					nominal = CP;
				else if(vParm[p] == "M12")
					nominal = M12;
				else if(vParm[p] == "M23")
					nominal = M23;
				else if(vParm[p] == "S12")
					nominal = S12;
				else if(vParm[p] == "S13")
					nominal = S13;
				else if(vParm[p] == "S23")
					nominal = S23;
				else
					throw std::invalid_argument("Unknown parameter\n");

				double chi2 = (nominal - vValue[p]) / vSigma[p];
				penalty += chi2*chi2;
			}

			X2 += penalty;
			sX2_new->Fill();
		}

		out->cd();
		if (axis)
			axis->CloneTree()->Write();
		if (error)
			error->CloneTree()->Write();
		sX2_new->Write();


		inf->Close();
		out->Close();
		delete inf;
	}

	return 0;
}
