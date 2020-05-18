/* DataManager
 * handles root files, trees and their friends
 * GetEntry should be called from here
 * and Get pointers to variables
 */

#ifndef DataManager_H
#define DataManager_H

#include <cstdlib>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Utils.h"

#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

class DataManager
{
	public:
		DataManager(std::string treeName, int inst, int proc);
		~DataManager();

		bool Add(std::string fileName, std::string treeName);
		bool AddFriend(std::string fileName, std::string treeName);

		void SetBranches();
		void SetBranch(TBranch *& branch);

		bool GetEntry(unsigned int nEntry);
		bool Get(const std::string branch, double *&va, unsigned int &size);
		bool Get(const std::string branch, int *&va, unsigned int &size);

		unsigned int GetInstance();
		unsigned int GetEntries();
		unsigned int GetProcesses();

	private:
		std::map<std::string, double*> mBranchDouble;
		std::map<std::string, double*>::iterator imd;
		std::map<std::string, unsigned int> mBrSizeDouble;
		
		std::map<std::string, int*> mBranchInteger;
		std::map<std::string, int*>::iterator imi;
		std::map<std::string, unsigned int> mBrSizeInteger;
		
		//chains, main and list of freinds
		TChain *mainChain;
		std::vector<TChain*> vfriendChain;

		//maps for branch adresses
		unsigned int instance, nProcess;
		bool kVerbosity;
};

#endif
