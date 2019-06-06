#include "DataManager.h"

DataManager::DataManager(std::string treeName, int inst, int proc) :
	instance(inst),
	nProcess(proc),
	kVerbosity(Utils::verbosity)
{
	mainChain = new TChain(treeName.c_str());
}

DataManager::~DataManager()
{
	delete mainChain;
	for (unsigned int i = 0; i < vfriendChain.size(); ++i)
		delete vfriendChain.at(i);

	std::map<std::string, double*>::iterator imd;
	for (imd = mBranchDouble.begin(); imd != mBranchDouble.end(); ++imd)
		delete imd->second;


	std::map<std::string, int*>::iterator imi;
	for (imi = mBranchInteger.begin(); imi != mBranchInteger.end(); ++imi)
		delete imi->second;
}

bool DataManager::Add(std::string fileName, std::string treeName)
{
	if (kVerbosity)
	{
		std::cout << "DataManager::Add file(s) \"" << fileName << "\"" << std::endl;
		std::cout << "                 run DataManager::SeedBlossom only after you add all files" << std::endl;;
	}

	mainChain->Add(fileName.c_str());
}

bool DataManager::AddFriend(std::string fileName, std::string treeName)
{
	if (kVerbosity)
	{
		std::cout << "DataManager::Add file(s) \"" << fileName << "\"" << std::endl;
		std::cout << "                 run DataManager::SeedBlossom only after you add all files" << std::endl;;
	}

	vfriendChain.push_back(new TChain(treeName.c_str()));
	vfriendChain.back()->Add(fileName.c_str());	//not sure!!!!

	mainChain->AddFriend(fileName.c_str());
}

void DataManager::SetBranches()		//attach all the branches to a variable
{
	TObjArray *listBranch = mainChain->GetListOfBranches();

	for (unsigned int i = 0; i < listBranch->GetEntries(); ++i)	//loop over the branches of the main tree
	{
		TBranch* br = dynamic_cast<TBranch*>(listBranch->At(i));
		SetBranch(br);
	}

	for (unsigned int f = 0; f < vfriendChain.size(); ++f)		//loop over the friend trees
	{								//and for each loop over the branches
		listBranch = vfriendChain.at(f)->GetListOfBranches();
		for (unsigned int i = 0; i < listBranch->GetEntries(); ++i)
		{
			TBranch* br = dynamic_cast<TBranch*>(listBranch->At(i));
			SetBranch(br);
		}
	}
}

void DataManager::SetBranch(TBranch *& branch)
{
	std::string typeBranch(branch->GetTitle());

	std::string size = "";
	unsigned int dimArray = 1;		//minimum size
	size_t arrayPos = typeBranch.find("[");
	if (arrayPos != std::string::npos)
	{
		size = typeBranch.substr(arrayPos+1, typeBranch.find("]"));
		dimArray = std::stoul(size, NULL, 10);
	}

	size_t typePos = typeBranch.find("/");
	char typeChar;
	if (typePos != std::string::npos)
		typeChar = typeBranch.at(typePos);

	std::string nameBranch(branch->GetName());
	switch (typeChar)
	{
		case 'C':	//a character string terminated by the 0 character
		case 'B':	//an 8 bit signed integer (Char_t)
		case 'b':	//an 8 bit unsigned integer (UChar_t)
		case 'S':	//a 16 bit signed integer (Short_t)
		case 's':	//a 16 bit unsigned integer (UShort_t)
		case 'I':	//a 32 bit signed integer (Int_t)
		case 'i':	//a 32 bit unsigned integer (UInt_t)
		case 'L':	//a 64 bit signed integer (Long64_t)
		case 'l':	//a 64 bit unsigned integer (ULong64_t)
		case 'O':	//[the letter o, not a zero] a boolean (Bool_t)
			mBranchInteger[nameBranch] = new int[dimArray];
			mBrSizeInteger[nameBranch] = dimArray;
			mainChain->SetBranchAddress("nameBranch", mBranchInteger[nameBranch]);
			break;
		case 'F':	//a 32 bit floating point (Float_t)
		case 'D':	//a 64 bit floating point (Double_t)
			mBranchDouble[nameBranch] = new double[dimArray];
			mBrSizeDouble[nameBranch] = dimArray;
			mainChain->SetBranchAddress("nameBranch", mBranchDouble[nameBranch]);
			break;
		default:
			if (kVerbosity)
				std::cerr << "DataManager::SetBranch type \"" << char(typeChar) << "\" unknown" << std::endl;
			break;
	}
}

/*	if branches are neglected while getentry, computation can be speeded up,
 *	but according to manual, SetBranchStatus must be called before assigning addresses
 *
bool DataManager::SetBranchStatus(const std::string branch, double *&va, unsigned int &size)
{
	imd = mBranchDouble.find(branch);

	if (imd != mBranchDouble.end())
	{

	}
}
*/

bool DataManager::GetEntry(unsigned int nEntry)
{
	return mainChain->GetEntry(nEntry);
}

bool DataManager::Get(const std::string branch, double *&va, unsigned int &size)
{
	imd = mBranchDouble.find(branch);

	if (imd != mBranchDouble.end())
	{
		va = imd->second;
		size = mBrSizeDouble[imd->first];
		return true;
	}
	else
	{
		if( kVerbosity )
			std::cerr << "DataManager::Get \"" << branch << "\" not found" << std::endl;
		return false;
	}
}

bool DataManager::Get(const std::string branch, int *&va, unsigned int &size)
{
	imi = mBranchInteger.find(branch);

	if (imi != mBranchInteger.end())
	{
		va = imi->second;
		size = mBrSizeInteger[imd->first];
		return true;
	}
	else
	{
		if( kVerbosity )
			std::cerr << "DataManager::Get \"" << branch << "\" not found" << std::endl;
		return false;
	}
}

unsigned int DataManager::GetInstance()
{
	return instance;
}

unsigned int DataManager::GetEntries()
{
	return mainChain->GetEntries();
}

unsigned int DataManager::GetProcesses()
{
	return nProcess;
}
