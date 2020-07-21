std::map<int, TH1D*> AtmoSample::CreateHistograms()
{
	std::map<std::string, int> hist_types;
	std::map<std::string, std::vector<double> > hist_axes;
	cd->Get("BinningType_", hist_types);
	cd->Get("BinAxis_", hist_axes);

	// store binning information as TH2D/TH1D
	// the content will be stored in Eigen matrices
	// and there will be a mapping with row-column and global bins
	//
	std::map<std::string, TH1D*> reco_bins;
	std::map<int, std::string> type_names;

	for (const auto &ih : hist_types) {
		std::string ax0 = ih.first + "_0";
		std::string ax1 = ih.first + "_0";
		if (hist_axes[ax1].size()-1 == 1) {	// only one bin in y-axis
			reco_bins[ih.first] = new TH1D (ih.first.c_str(), ih.first.c_str(),
					hist_axes[ax0].size()-1, &hist_axes[ax0][0]);
		}
		else {
			reco_bins[ih.first] = new TH2D (ih.first.c_str(), ih.first.c_str(),
					hist_axes[ax0].size()-1, &hist_axes[ax0][0],
					hist_axes[ax1].size()-1, &hist_axes[ax1][0]);
		}
		// reverse lookup
		type_names[ih.second] = ih.first;
	}

	return reco_bins, type_names;	// not allowed in C++
}

void AtmoSample::FillHistograms()
{
	std::string chain, friends, tree_name;
	cd->Get("tree_input", chain);
	cd->Get("tree_name", tree_name);
	cd->Get("friend_input", friends);
	cd->Get("friend_name", friend_name);

	TChain *dm = new TChain(tree_name.c_str());
	TChain *fr = new TChain(friend_name.c_str());

	dm->Add(chain.c_str());
	fr->Add(friends.c_str());

	dm->AddFriend(friend_name.c_str());

	int ipnu, mode, itype;
	double drinu[3], dir[3], flxho[3];
	double pnu, amom, weightx, ErmsHax, nEAveHax;

	dm->SetBranchAddress("ipnu",    &ipnu);	//pdg
	dm->SetBranchAddress("mode",    &mode);
	dm->SetBranchAddress("itype",   &itype);

	dm->SetBranchAddress("dirnu",   dirnu);
	dm->SetBranchAddress("pnu",     &pnu);
	dm->SetBranchAddress("dir",     dir);
	dm->SetBranchAddress("amom",    &amom);
	dm->SetBranchAddress("flxho",   &flxho);
	dm->SetBranchAddress("weightx", &weightx);

	// only if friends
	dm->SetBranchAddress("ErmsHax",  ErmsHax  );
	dm->SetBranchAddress("nEAveHax", nEAveHax );

	dm->SetBranchStatus("*", 0);
	dm->SetBranchStatus("dirnu", 1);
	dm->SetBranchStatus("pnu", 1);

	std::map<int, double> E_min, E_max, cosZ_min, cosZ_max;
	for (int i = 0; i < dm->GetEntries(); ++i) {
		dm->GetEntry(i);

		pnu = log10(pnu) + 3;
		if ((E_min.count(ipnu) && E_min[ipnu] > pnu)
			|| !E_min.count(ipnu))
			E_min[ipnu] = pnu;
		if ((E_max.count(ipnu) && E_max[ipnu] < pnu)
			|| !E_max.count(ipnu))
			E_max[ipnu] = pnu;

		if ((cosZ_min.count(ipnu) && cosZ_min[ipnu] > dirnu[2])
			|| !cosZ_min.count(ipnu))
			cosZ_min[ipnu] = dirnu[2];
		if ((cosZ_max.count(ipnu) && cosZ_max[ipnu] < dirnu[2])
			|| !cosZ_max.count(ipnu))
			cosZ_max[ipnu] = dirnu[2];
	}

	// create matrices for bin contents
	std::map<int, TH1D*> true_bins;
	std::map<std::string, Eigen::MatrixXd> bin_contents;
	for (const auto &ir : type_names) {
		for (const auto &it : E_min) {
			std::string name = FromPDGtoString(it.first);
			std::string hname = ir.second + "_" + name;
			// 100 bins in true neutrino energy and 50 bins in cosz
			true_bins = new TH2D(name.c_str(), name.c_str(),
					100, it.second, E_max[it.first],
					50, cosZ_min[it.first], cosZ_max[it.first]);
			// granularity of true bins is fixed
			// but eventually should be defined in card file

			int ny = reco_bins[ir.second]->GetNbinsX() * 
				 reco_bins[ir.second]->GetNbinsY();

			int nx = 50 * 100;

			// one bin content per neutrino initial flavour
			bin_contens[ir.second + "_nuE0_" + name] = Eigen::MatrixXd::Zero(ny, nx);
			bin_contens[ir.second + "_nuEB_" + name] = Eigen::MatrixXd::Zero(ny, nx);
			bin_contens[ir.second + "_nuM0_" + name] = Eigen::MatrixXd::Zero(ny, nx);
			bin_contens[ir.second + "_nuMB_" + name] = Eigen::MatrixXd::Zero(ny, nx);
		}
	}

	// loop through all events and fill histograms
	for (int i = 0; i < dm->GetEntries(); ++i) {
		dm->GetEntry(i, 1);	// all branches

		LeptonP = amom;
		CosineTheta = -1. * dir(2);  // + downgoing , - upgoing 
		NuCosineTheta = -1. * dirnu(2); // + downgoing , - upgoing 

		MCweight = weightx(0);
		Mode     = mode(0);

		NuCosineTheta = -1. * dirnu(2); // + downgoing , - upgoing 
		Energy = pnu;
		PDG = ipnu(0);

		CosineTheta = -1.*dir(2);  // + downgoing , - upgoing 

		SkipFlag = false;
		EventType = itype(0) - 1;
		itype -= 1;

		if (itype < _first_event_type || itype > _last_event_type)
			continue;

		// no NC tau's allowed, so skip
		if (std::abs(mode) > 30 && std::abs(ipnu) == 16)
			continue;

		// this is real data, which we do not want
		if (pnu < 1.0e-7 && ipnu == 0 && mode == 0)
			continue;

		// using these four variables to determine bin and EventType
		double lognu = Energy > 0 ? log10(Energy*1000) : 0;	// GeV
		double loglp = LeptonP > 0 ? log10(LeptonP) : 0;	// already in GeV

		if((useTrue & 2) == 2)
			vars.push_back( GetNuCosineZ() );
		else
			vars.push_back( GetCosineZ() );
		vars.push_back( logp         );


		std::string type = type_names[itype];

		// global bin of event in true and reco histograms
		int glo_nr = reco_hist[type]->FindBin(log10(amom), dir[2]);
		int glo_nt = true_hist[type]->FindBin(log10(pnu)+3, dirnu[2]);

		// total number of bins on x direction to find the absolut bin
		int r_nx = reco_bins[type]->GetNbinsX();
		int t_nx = true_bins[type]->GetNbinsX();

		// absolute bins of the matrix
		// rows of matrix
		int loc_nr = glo_nr - (r_nx + glo_nr / (r_nx + 2) + 1);
		// columns of matrix
		int loc_nt = glo_nt - (t_nx + glo_nt / (t_nx + 2) + 1);

		if (loc_nr < 0 || loc_nr > bin_contents[name].rows() - 1)
			continue;
		if (loc_nt < 0 || loc_nt > bin_contents[name].cols() - 1)
			continue;

		double flxratio = flxho[0] / flxho[1];
		if (std::abs(ipnu) == 12)
			flxratio = 1/flxratio;

		std::string name = FromPDGtoString(ipnu);

		double factor_E = 1, factor_M = 1;
		if (name.find("nuE") != std::string::npos)	// it is e type
      			factor_E = flxho[1] / flxho[0];
		else						// it is mu or tau type
      			factor_M = flxho[0] / flxho[1];

		bin_contens[ir.second + "_nuE0_" + name] += weightx * factor_E;
		bin_contens[ir.second + "_nuEB_" + name] += weightx * factor_E;
		bin_contens[ir.second + "_nuM0_" + name] += weightx * factor_M;
		bin_contens[ir.second + "_nuMB_" + name] += weightx * factor_M;
	}
}

std::string AtmoSample::FromPDGtoString(int pdg)
{
	std::string name;
	switch (pdg) {
		case -12:
			name = "nuEB";
			break;
		case -14:
			name = "nuMB";
			break;
		case -16:
			name = "nuTB";
			break;
		case 12:
			name = "nuE0";
			break;
		case 14:
			name = "nuM0";
			break;
		case 16:
			name = "nuT0";
			break;
	}

	return name;
}


// this things are used
// Energy
// LeptonP
// CosineTheta
// NuCosineTheta
void AtmoSample::Parse(int n)
{
	dm->GetEntry(n);

}

// determine bin of event
int AtmoSample::GetBin(std::vector<double> VarList, int type )
{
	// if type is not recognised as a histogram...
	if( ! _mpBP.count( type ) ) return -1;

	// number of dimensions of histogram bins
	int nDimensions = _mpBP[type]->GetnDims();

	// We will allow the input list to be larger 
	// than the current bin type's allowed dimension
	std::assert(VarList.size() < nDimensions)
	{
		std::cout << "Profiler::GetBin Error. Variable list is too small "  
			<< " size: " << VarList.size() << " expected " << nDimensions << endl 
			<< " for dimensions specified by bin type " << type << ". Exiting."
			<< std::endl ; 
		exit(-1);
	} 

	int DimFactors  [nDimensions];
	int LocalDimNum [nDimensions];
	for(int i = 0 ; i < nDimensions-1; i++ )
	{
		DimFactors[i] = 1;
		for(int j = i + 1; j < nDimensions ; j++ )
			DimFactors[i] *=  _mpBP[type]->GetnBinsPerDim( j );
	}
	DimFactors[nDimensions -1 ] = 1;


	int X  = _mpBP[type]->GetOffset();
	int SubBinOffSet;
	int BaseAddress;
	int ThisAddress;
	int NextAddress;

	// Illustrative examples for bin addressing 
	// in the linearalized schemne:
	//
	// If you have two dimensional bins
	//  (Momentum p, and Zenith Angle Z )
	//  with Np bins in p, and Nz bins in Z:
	// 
	//  The p bin edge ordering is:
	//  p ->  X + p*Nz
	//
	//  The z bin edge ordering is:
	//  z ->  X + p*Nz + z 

	//
	//  If three dimensional 
	//  (Momentum p, Zenith Angle Z , Cats c )
	//  p -> X + p*Nz*Nc
	//  z -> X + p*Nz*Nc + z*Nc
	//  c -> X + p*Nz*Nc + z*Nc + c
	//
	//  Below:
	//  DimFactors for dimension1 (p) is Nz*Nc
	//  DimFactors for dimension2 (z) is Nc
	//  DimFactors for dimension3 (c) is 1
	//

	int nBins;

	// loop over all of the bin dimensions
	for(int i = 0 ; i < nDimensions ; i++ )
	{
		LocalDimNum[i] = -1;

		// compute the base offset for this bin types i_th dimension
		// Start at the global offset
		BaseAddress = X;
		if( i > 0 )
			for(int k = 0 ; k < i-1 ; k ++ )
				BaseAddress += LocalDimNum[k]*DimFactors[k];        

		nBins = _mpBP[type]->GetnBinsPerDim( i );
		for(int j = 0 ; j <  nBins ; j++ )
		{ 

			// increment by one for the current bin edge
			ThisAddress = BaseAddress + j * DimFactors[i];        

			// increment by one more for the next bin edge
			NextAddress = ThisAddress + DimFactors[i];        

			// we are between two edges
			if( j != nBins -1 )
				if( VarList[i] >= _vpID[ ThisAddress ]->Edge(i) && VarList[i] < _vpID[ NextAddress ]->Edge(i) )  
				{ 
					LocalDimNum[i] = j ; 
					break;
				}

			if( j == nBins -1 )
				// we are in the last bin 
				if( VarList[i] >= _vpID[ ThisAddress ]->Edge(i) )
				{
					LocalDimNum[i] = j ; 
					break;
				}

		}// end of loop on number of bins in this dimension


		// we reached the end of the bin edges 
		// and we didn't fall in any of the defined bins -->Error!
		if( LocalDimNum[i] < 0 ) 
			return -1;  


	} // end of loop on dimensions


	ThisAddress = X;
	for(int k = 0 ; k < nDimensions ; k ++ )
		ThisAddress += LocalDimNum[k]*DimFactors[k];        

	return ThisAddress;

}
