#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "event/Flux.h"
#include "event/Reco.h"
#include "physics/Oscillator.h"

#include "TTree.h"
#include "TFile.h"

int main(int argc, char** argv)
{
	//std::string files;
	//cd->Get("files", files);

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
	double X2, X2_atmo;
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


	std::string path1(argv[1]);
	std::string path2(argv[2]);

	if (path1.back() != '/')
		path1.append("/");

	if (path2.back() != '/')
		path2.append("/");

	// get beam files
	std::string output = path1 + ".tmp_add_beam";
	std::string cmd = "ls " + path1 + "SpaghettiSens_penalised.T2HK.*.root > " + output;
	system(cmd.c_str());

	std::string file;
	std::vector<std::string> beamFiles;
	std::ifstream allfiles(output.c_str());
	while (std::getline(allfiles, file))
		beamFiles.push_back(file);
	allfiles.close();

	// get atmo files
	output = path2 + ".tmp_add_atmo";
	cmd = "ls " + path2 + "SpaghettiSens_penalised.T2HK.*.root > " + output;
	system(cmd.c_str());

	std::vector<std::string> atmoFiles;
	allfiles.open(output.c_str());
	while (std::getline(allfiles, file))
		atmoFiles.push_back(file);

	if (beamFiles.size() == atmoFiles.size())
		std::cout << "GOOD" << std::endl;

	for (int f = 0; f < beamFiles.size(); ++f)
	{
		file = beamFiles[f];
		TFile *beamf = new TFile(file.c_str(), "READ");
		if (beamf->IsZombie())
		{
			std::cout << "file doesn't exist\n";
			continue;
		}

		TTree* axis    = static_cast<TTree*>(beamf->Get("axisTree"));
		TTree* error   = static_cast<TTree*>(beamf->Get("errorTree"));
		TTree* sX2_old = static_cast<TTree*>(beamf->Get("stepX2Tree"));


		TFile *atmof = new TFile(atmoFiles[f].c_str(), "READ");
		if (beamf->IsZombie())
		{
			std::cout << "file doesn't exist\n";
			continue;
		}

		TTree* sX2_atmo = static_cast<TTree*>(atmof->Get("stepX2Tree"));
		sX2_atmo->SetBranchAddress("X2", &X2_atmo);

		if (sX2_old->GetEntries() != sX2_atmo->GetEntries())
			std::cout << "VERY BAD" << std::endl;

		std::cout << "Combining " << file << " and " << atmoFiles[f] << std::endl;
		int pos = file.find("SpaghettiSens_penalised.");
		if (pos != std::string::npos)
			file.replace(pos+14, 9, "atmo"); 
		std::cout << "\tinto " << file << std::endl;

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
			sX2_atmo->GetEntry(i);

			X2 += X2_atmo;
			sX2_new->Fill();
		}

		out->cd();
		if (axis)
			axis->CloneTree()->Write();
		if (error)
			error->CloneTree()->Write();
		sX2_new->Write();


		beamf->Close();
		atmof->Close();

		out->Close();

		delete beamf;
		delete atmof;
	}

	return 0;
}
