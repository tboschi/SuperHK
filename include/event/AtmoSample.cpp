#include "AtmoSample.h"

// process is a string specying what to process of the sample input files
// R -> load reconstruction
// B -> define binning
// S -> load systematics
AtmoSample::AtmoSample(std::string card, std::string process) :
	Sample(card)
{
	CardDealer cd(card);
	Init(cd, process);
}

AtmoSample::AtmoSample(const CardDealer &cd, std::string process) :
	Sample(cd)
{
	Init(cd, process);
}

AtmoSample::AtmoSample(CardDealer *cd, std::string process) :
	Sample(cd)
{
	Init(*cd, process);
}

void AtmoSample::Init(const CardDealer &cd, std::string process)
{
	atm_path = std::unique_ptr<Atmosphere>(new Atmosphere(cd));

	_oscf = { {"nuE0_nuE0", {Nu::E_, Nu::E_}},
		  {"nuE0_nuM0", {Nu::E_, Nu::M_}},
		  {"nuE0_nuT0", {Nu::E_, Nu::T_}},
		  {"nuM0_nuE0", {Nu::M_, Nu::E_}},
		  {"nuM0_nuM0", {Nu::M_, Nu::M_}},
		  {"nuM0_nuT0", {Nu::M_, Nu::T_}},
		  //
		  {"nuEB_nuEB", {Nu::Eb, Nu::Eb}},
		  {"nuEB_nuMB", {Nu::Eb, Nu::Mb}},
		  {"nuEB_nuTB", {Nu::Eb, Nu::Tb}},
		  {"nuMB_nuEB", {Nu::Mb, Nu::Eb}},
		  {"nuMB_nuMB", {Nu::Mb, Nu::Mb}},
		  {"nuMB_nuTB", {Nu::Mb, Nu::Tb}} };

	if (!cd.Get("sample", _type))
		// default value too long..
		_type = { "SubGeVElike0Decay", 
			  "SubGeVElike1Decay", 
			  "SingleRingPi0", 
			  "SubGeVMulike0Decay", 
			  "SubGeVMulike1Decay", 
			  "SubGeVMulike2Decay", 
			  "TwoRingPi0", 
			  "MultiGeV_elike_nue",
			  "MultiGeV_elike_barnue", 
			  "MultiGeV_mulike", 
			  "MultiRing_elike_nue", 
			  "MultiRing_elike_barnue", 
			  "MultiRing_mulike", 
			  "MultiRingOther", 
			  "PCStop", 
			  "PCThru", 
			  "UpStopMu", 
			  "UpThruMuNonShowering", 
			  "UpThruMuShowering" };


	if (!cd.Get("MC_scale", _weight))// adjust MC time and exposure
		_weight = 188.4 / 22.5;	// HK to SK ratio
	if (!cd.Get("reduce", _reduce))	// FV reduction
		_reduce = 0.5;



	if (kVerbosity) {
		std::cout << "AtmoSample: types: ";
		for (const auto &it : _type)
			std::cout << "\t" << it;
		std::cout << std::endl;
		std::cout << "AtmoSample: number of scale systematics " << _nScale << std::endl;
		for (const auto &it : _scale)
			std::cout << "\t" << it.first << " -> " << it.second << std::endl;
	}

	_point = -1;


	if (process.empty())
		process = "RBS";
	else
		std::transform(process.begin(), process.end(), process.begin(),
				[](unsigned char c){ return std::toupper(c); });

	if (process.find('R') != std::string::npos)
		LoadReconstruction(cd);
	if (process.find('B') != std::string::npos)
		DefineBinning();
	if (process.find('S') != std::string::npos)
		LoadSystematics(cd);
}


void AtmoSample::LoadSystematics(const CardDealer &cd)
{
	zeroEpsilons = true;
	std::string file_name, tree_name;
	if (!cd.Get("systematic_file", file_name))
		throw std::invalid_argument("AtmoSample: no systematic file for atmo sample,"
			    "cannot determine number of systematics");

	TFile * mf = new TFile(file_name.c_str());
	if (mf->IsZombie())
		throw std::invalid_argument("WARNING - AtmoSample: file " + file_name
				  + " does not exist");

	if (!cd.Get("systematic_tree", tree_name))
		tree_name = "sigmatree";

	std::set<std::string> skip_sys;
	cd.Get("skip", skip_sys);	// errors to skip

	TString * name_ptr = new TString;
	TTree * syst = static_cast<TTree*>(mf->Get(tree_name.c_str()));
	syst->SetBranchAddress("FijName", &name_ptr);
	int nsyst = syst->GetEntries();
	//double error;
	//syst->SetBranchAddress("Sigma", &error);

	_nSys = 0;
	for (int i = 0; i < nsyst; ++i) {
		syst->GetEntry(i);
		std::string name_err = name_ptr->Data();
		if (skip_sys.find(name_err) == skip_sys.end())
			++_nSys;
	}
	std::cout << "AtmoSample: counted number of systematics " << _nSys << std::endl;

	_corr = Eigen::MatrixXd::Identity(_nSys, _nSys);

	_sysMatrix.clear();
	if (kVerbosity)
		std::cout << "AtmoSample: creating matrix "
			  << _nBin << " x " << _nSys << std::endl;
	_sysMatrix[0] = Eigen::ArrayXXd::Zero(_nBin, _nSys);
	for (int sigma = -3; sigma < 4; sigma += 2)
		_sysMatrix[sigma] = Eigen::ArrayXXd::Zero(_nBin, _nSys);


	// if only stats fit, this function is done here
	if (cd.Get("stats_only")) {
		std::cout << "AtmoSample: systematic errors for atmo sample will not be fitted" << std::endl;
		nsyst = 0;	// with this, the following loop won't work
	}


	int k_err = 0;
	for (int i = 0; i < nsyst; ++i) {
		//int off = 0;	// offset for global bin
		syst->GetEntry(i);
		std::string name_err = name_ptr->Data();
		if (skip_sys.find(name_err) != skip_sys.end()) {
			if (kVerbosity)
				std::cout << "AtmoSample: skipping " << name_err
					  << " systematic" << std::endl;
			continue;
		}
		zeroEpsilons = false;

		TH1D* hsys = static_cast<TH1D*>(mf->Get(name_err.c_str()));
		for (const std::string &it : _type) {
			int i = _offset[it];
			for (int n : _binpos[it]) {
				if (std::abs(hsys->GetBinContent(n+1)) > 1)
					std::cout << k_err << ", " << n << " out of scale\n";
				for (double sigma = -3; sigma < 4; sigma += 2)
					_sysMatrix[sigma](i, k_err)
						= sigma * hsys->GetBinContent(n+1);
						//= sigma * std::min(1., std::max(-1., 
						//	   hsys->GetBinContent(n+1)));
				++i;
			}
		}
		++k_err;
	}

	mf->Close();

	if (kVerbosity) {
		std::cout << "AtmoSample: Importing " << _nSys
			  << " entries from " << file_name
			  << " ( " << tree_name << " ) " << std::endl;
		std::cout << "Total number of systematics is "
			  << _corr.cols() << std::endl;
	}
}


// same of BeamSample _-> move into base class
void AtmoSample::LoadReconstruction(const CardDealer &cd) {

	// check first if precomputed files exist
	std::string chain, tree_name; //, friends, friend_name;
	if (cd.Get("pre_input_NH", chain)) {
		if (!cd.Get("pre_tree_name", tree_name))
			tree_name = "atmoTree";

		if (kVerbosity)
			std::cout << "AtmoSample: opening NH input\n";
		nh = std::unique_ptr<TChain>(new TChain(tree_name.c_str()));

		nh->Add(chain.c_str());

		nh->SetBranchAddress("Point", &point);
		nh->SetBranchAddress("Bins",  &bins);
		nh->SetBranchAddress("Data",  data);

		nh->SetBranchStatus("*", 0);
		nh->SetBranchStatus("Point", 1);
		//linking fitting point to position in tree
		//as they may not be ordered
		for (int i = 0; i < nh->GetEntries(); ++i) {
			nh->GetEntry(i);
			pre_point_NH[point] = i;
		}
	}

	if (cd.Get("pre_input_IH", chain)) {
		if (!cd.Get("pre_tree_name", tree_name))
			tree_name = "atmoTree";
		ih = std::unique_ptr<TChain>(new TChain(tree_name.c_str()));

		if (kVerbosity)
			std::cout << "AtmoSample: opening IH input\n";
		
		ih->Add(chain.c_str());

		ih->SetBranchAddress("Point", &point);
		ih->SetBranchAddress("Bins",  &bins);
		ih->SetBranchAddress("Data",  data);

		ih->SetBranchStatus("*", 0);
		ih->SetBranchStatus("Point", 1);
		//linking fitting point to position in tree
		//as they may not be ordered
		for (int i = 0; i < ih->GetEntries(); ++i) {
			ih->GetEntry(i);
			pre_point_IH[point] = i;
		}
	}


	// extract bining from recostruction matrices
	// binning is compressed, such that nonempty bins are not stored

	// type of histogram, they are numbered with integer
	// to identify the event type in the MC files
	std::map<std::string, int> hist_types;
	std::map<std::string, std::vector<double> > hist_axes;
	cd.Get("bintype_", hist_types);
	cd.Get("binaxis_", hist_axes);

	// store binning information as TH2D
	for (const auto &ih : hist_types) {
		if (_type.find(ih.first) == _type.end())
			continue;
		std::string ax0 = ih.first + "_0";
		std::string ax1 = ih.first + "_1";
		std::string hname = "reco_" + ih.first;
		_reco_hist[ih.first] = new TH2D (hname.c_str(), hname.c_str(),
				hist_axes[ax0].size()-1, &hist_axes[ax0][0],
				hist_axes[ax1].size()-1, &hist_axes[ax1][0]);
		// reverse lookup
		_type_names[ih.second] = ih.first;
	}

	if (!cd.Get("MC_input", chain))
		throw std::invalid_argument("AtmoSample: no reconstruction files in card,"
			       		    "very bad!");
	if (!cd.Get("MC_tree_name", tree_name))
		tree_name = "osc_tuple";
	//cd.Get("friend_input", friends);
	//cd.Get("friend_name", friend_name);

	if (kVerbosity)
		std::cout << "AtmoSample: opening " << chain
			  << " for " << tree_name << std::endl;
	dm = std::unique_ptr<TChain>(new TChain(tree_name.c_str()));
	//TChain *fr = new TChain(friend_name.c_str());

	dm->Add(chain.c_str());
	//fr->Add(friends.c_str());

	//dm->AddFriend(friend_name.c_str());

	//int ipnu, mode, itype;
	//float dirnu[3], dir[3], flxho[3];
	//float pnu, amom, weightx, ErmsHax, nEAveHax;

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
	//dm->SetBranchAddress("ErmsHax",  ErmsHax  );
	//dm->SetBranchAddress("nEAveHax", nEAveHax );

	// loop first...
	// zenith angle will be always between -1 and 1
	// find max and min values of true energy information
	// group by pdg value
	//std::map<int, double> E_min, E_max;
	//for (int i = 0; i < dm->GetEntries(); ++i) {
	//	dm->GetEntry(i);

	//	pnu = log10(pnu) + 3;
	//	if (!E_min.count(ipnu) || E_min[ipnu] > pnu)
	//		E_min[ipnu] = pnu;
	//	if (!E_max.count(ipnu) || E_max[ipnu] < pnu)
	//		E_max[ipnu] = pnu;
	//}
	//for (const auto &it : E_min)
	//	std::cout << "energy lim " << it.first
	//		  << " = " << it.second
	//		  << " : " << E_max[it.first] << std::endl;
	// type_names is map type(int) and hist name(str) pairs

	// ready to fill histograms
	// create matrices for bin contents
}



/*
std::map<std::string, Eigen::VectorXd> AtmoSample::BuildSamples(Oscillator *osc)
{
	std::map<std::string, Eigen::VectorXd> samples;

	for (const std::string &it : _type) {
		const auto &tr = _true_bins[it];

		// map contains the right oscillation flavour combination 
		// for this particular type binning
		// if osc = 0 then probs contains ones only 
		std::map<std::string, Eigen::VectorXd> probs =
						atm_path->Oscillate(_global[it], _oscf, osc);

		for (const auto &ip : probs) {
			std::string binname = it + "_" + ip.first;
			if (samples.count(it))
				samples[it] += _bin_contents[binname] * ip.second;
			else
				samples[it]  = _bin_contents[binname] * ip.second;
		}
	}

	return samples;
}
*/

std::map<std::string, Eigen::VectorXd> AtmoSample::BuildSamples(std::shared_ptr<Oscillator> osc)
{
	std::map<std::string, Eigen::VectorXd> samples;

	for (const auto &ih : _reco_hist)
		ih.second->Reset("ICES");

	dm->SetBranchStatus("*", 1);
	int first_event_type = _type_names.begin()->first;
	int last_event_type = _type_names.rbegin()->first;
	// loop through all events and fill histograms

	//std::ostringstream address;
	//address << "osccalc_" << (void const *)osc << ".root";
	//TFile oscout(address.str().c_str(), "RECREATE");
	//double fE, fM, pE, pM, W;
	//TTree *tt = new TTree("info", "Oscillated info");
	//tt->Branch("ipnu",    &ipnu, "ipnu/I");	//pdg
	//tt->Branch("mode",    &mode, "mode/I");
	//tt->Branch("itype",   &itype, "type/I");
	//tt->Branch("dir",     &dir[2], "dir/F");
	//tt->Branch("amom",    &amom, "amom/F");
	//tt->Branch("W",    &W,  "W/D");
	//tt->Branch("fE",   &fE, "fE/D");
	//tt->Branch("fM",   &fM, "fM/D");
	//tt->Branch("pE",   &pE, "pE/D");
	//tt->Branch("pM",   &pM, "pM/D");
	

	std::cout << "AtmoSample: reading " << dm->GetEntries() << " entries" << std::endl;
	for (int i = 0; i < dm->GetEntries(); ++i) {
	//for (int i = 0; i < 1000; ++i) {
		//std::cout << "entry " << i << std::endl;
		dm->GetEntry(i);	// all branches
		itype -= 1;

		// uknown event types
		if (itype < first_event_type || itype > last_event_type)
			continue;

		// no NC tau's allowed, so skip
		if (std::abs(mode) >= 30 && std::abs(ipnu) == 16)
			continue;

		// this is real data, which we do not want
		if (pnu < 1.0e-7 && ipnu == 0 && mode == 0)
			continue;

		std::string type = _type_names[itype];

		// bin is not in range, just skip it
		int bin = _reco_hist[type]->FindBin(log10(amom), -dir[2]);
		if ( _reco_hist[type]->IsBinOverflow(bin)
		  || _reco_hist[type]->IsBinUnderflow(bin))
			continue;


		weightx *= _weight;
		if (itype > 70)
			weightx *= _reduce;
		//W = weightx;

		double factor_E = 1, factor_M = 1;
		if (std::abs(ipnu) == 12)	// it is e type
      			factor_M = flxho[1] / flxho[0];
		else				// it is mu or tau type
      			factor_E = flxho[0] / flxho[1];
		//fE = factor_E;
		//fM = factor_M;


		Nu::Flavour nu_out = Nu::fromPDG(ipnu);
		//Nu::Flavour nue, num;
		if (osc && std::abs(mode) < 30) {	// NCs have mode >= 30
			// LUT is cleared in this step
			// -1 dirnu[2] ?? ???? ?? ? ? ? ? ? ?
			//double hgt = atm_path->GenerateRandomHeight(nu_out, -dirnu[2], pnu);
			//std::cout << "production height " << hgt << std::endl;
			//Oscillator::Profile ld = atm_path->MatterProfile(-dirnu[2], 15);
			//osc->SetMatterProfile(ld);

			// this automatically get production height
			// if honda flux is defined in card then height is generted
			// otherwise is fixed to value defiend in card
			osc->SetMatterProfile(atm_path->MatterProfile(nu_out, -dirnu[2], pnu));

			// but thanks to LUT, probability is computed only once
			if (ipnu > 0)	// neutrino
			{
				//pE = osc->Probability(Nu::E_, nu_out, pnu);
				//pM = osc->Probability(Nu::M_, nu_out, pnu);

				weightx *= factor_E * osc->Probability(Nu::E_, nu_out, pnu)
					+  factor_M * osc->Probability(Nu::M_, nu_out, pnu);
				//nue = Nu::E_; num = Nu::M_;
			}
			else {		//antineutrino
				//pE = osc->Probability(Nu::Eb, nu_out, pnu);
				//pM = osc->Probability(Nu::Mb, nu_out, pnu);

				weightx *= factor_E * osc->Probability(Nu::Eb, nu_out, pnu)
					+  factor_M * osc->Probability(Nu::Mb, nu_out, pnu);
				//nue = Nu::Eb; num = Nu::Mb;
			}
			//std::cout << "layers at " << -dirnu[2] << " is "
			//	  << ld.size() << ": " << osc->Length() << "\n";
			//for (const auto &l : ld)
			//	std::cout << "\t" << std::get<0>(l) << std::endl;
			//if (itype == 72)
				//std::cout << i << " " << itype << " -> " << factor_E << "\t" << factor_M
				  	//<< "\t" << weightx << std::endl; 
		}

		_reco_hist[type]->Fill(log10(amom), -dir[2], weightx);
	}

	//tt->Write();
	//oscout.Close();

	for (const auto &ir : _reco_hist) {
		int xs = ir.second->GetNbinsX();	// energy
		int ys = ir.second->GetNbinsY();	// cosz
		Eigen::MatrixXd rm = Eigen::Map<Eigen::MatrixXd>
			(ir.second->GetArray() + xs + 2, xs + 2, ys);

		// remove first and last columns
		rm.transposeInPlace();
		rm.leftCols(xs) = rm.middleCols(1, xs);
		rm.conservativeResize(ys, xs);

		//std::cout << ir.first << "\n" << rm.rowwise().sum() << "\n" << std::endl;
		samples[ir.first] = Eigen::Map<Eigen::VectorXd>(rm.data(), rm.cols() * rm.rows());

	}

	return samples;
}

Eigen::VectorXd AtmoSample::ConstructSamples(std::shared_ptr<Oscillator> osc)
{
	if (osc && osc->GetHierarchy() == Oscillator::normal &&
	    pre_point_NH.count(_point)) {

		if (kVerbosity > 1)
			std::cout << "AtmoSample: using precomputed NH at point " << _point << std::endl;

		nh->SetBranchStatus("*", 1);
		nh->GetEntry(pre_point_NH[_point]);

		Eigen::VectorXd samples = Eigen::Map<Eigen::VectorXd>(data, bins);
		return _stats * samples;
	}
	if (osc && osc->GetHierarchy() == Oscillator::inverted &&
	    pre_point_IH.count(_point)) {

		if (kVerbosity > 1)
			std::cout << "AtmoSample: using precomputed IH at point " << _point << std::endl;

		ih->SetBranchStatus("*", 1);
		ih->GetEntry(pre_point_IH[_point]);

		Eigen::VectorXd samples = Eigen::Map<Eigen::VectorXd>(data, bins);
		return _stats * samples;
	}
	else	// no precomputation
		return Sample::ConstructSamples(osc);
}


/*
void AtmoSample::LoadReconstruction(std::string channel)
{
	std::string reco_card;
	if (!cd.Get(channel, reco_card)) {
		if (kVerbosity > 1)
			std::cout << "WARNING - AtmoSample: reconstruction file for channel \""
				  << reco_card << "\" is missing" << std::endl;
		return;
	}

	if (kVerbosity > 2)
		std::cout << "AtmoSample: using reconstruction file " << reco_card << std::endl;

	CardDealer cd(reco_card);

	std::string reco_path;
	if (!cd.Get("reco_path", reco_path)) {
		if (kVerbosity > 0)
			std::cerr << "WARNING - AtmoSample: no reco_path specified in \""
				  << reco_card << "\"" << std::endl;
		return;
	}

	//std::cout << "opening " << fileReco << std::endl;
	TFile* inFile = new TFile(reco_path.c_str(), "READ");
	if (inFile->IsZombie()) {
		if (kVerbosity > 0)
			std::cerr << "WARNING - AtmoSample: file " << reco_path
				  << " does not exist" << std::endl;
		return;
	}


	std::string reconame;
	if (!cd.Get("reco_info", reconame)) {
		if (kVerbosity > 0)
			std::cerr << "WARNING - AtmoSample: no reco_info specified in \""
				  << reco_card << "\"" << std::endl;
		return;
	}

	std::string truename;
	if (!cd.Get("true_info", truename)) {
		if (kVerbosity > 0)
			std::cerr << "WARNING - AtmoSample: no true_info specified in \""
				  << reco_card << "\"" << std::endl;
		return;
	}

	std::string type;
	if (!cd.Get("type", type)) {
		if (kVerbosity > 0)
			std::cerr << "WARNING - AtmoSample: no event type specified in \""
				  << reco_card << "\"" << std::endl;
		return;
	}

	if (inFile->Get(reconame.c_str()))
		_reco_hist[type] = static_cast <TH2D*>(inFile->Get(reconame.c_str()));
	if (inFile->Get(truename.c_str()))
		_true_hist[type] = static_cast <TH2D*>(inFile->Get(truename.c_str()));


	TIter next(inFile->GetListOfKeys());
	TKey *k;
	while ((k = static_cast<TKey*>(next())))
	{
		if ( reconame.compare(k->GetName())
		  || truename.compare(k->GetName()) )
			continue;

		// else it is a matrix object
		TMatrixT<double>* mat = static_cast<TMatrixT<double>*> (k->ReadObj());
		// keyname format is type_oscf
		std::string keyname = k->GetName();
		_bin_contents[keyname] = Eigen::Map<Eigen::MatrixXd> (mat->GetMatrixArray(),
							mat->GetNrows(), mat->GetNcols());
	}

	inFile->Close();
}
*/

// Create the histograms using the binning information from card
// this bins refer to the reconstruction variables, log p and cos zenith
// loop through MC files and defines true histograms
// to add -> create card output
// NOT A MEMBER FUNCTION, to be used as one off
/*
void CreateTensor(CardDealer *cd)
{
	// type of histogram, they are numbered with integer
	// to identify the event type in the MC files
	std::map<std::string, int> hist_types;
	std::map<std::string, std::vector<double> > hist_axes;
	cd.Get("bintype_", hist_types);
	cd.Get("binaxis_", hist_axes);

	// store binning information as TH2D/TH1D
	// the content will be stored in Eigen matrices
	// and there will be a mapping with row-column and global bins
	//
	std::map<std::string, Eigen::MatrixXd> bin_contents;
	std::map<std::string, TH2D*> reco_bins, true_bins;	// saving as TH1D pointer!
	std::map<int, std::string> type_names;

	for (const auto &ih : hist_types) {
		std::cout << ih.first << " and " << ih.second << std::endl;
		std::string ax0 = ih.first + "_0";
		std::string ax1 = ih.first + "_1";
		std::cout << "sizes " << hist_axes[ax0].size() << ", "
			  << hist_axes[ax1].size() << std::endl;
		std::string hname = "reco_" + ih.first;
		reco_bins[ih.first] = new TH2D (hname.c_str(), hname.c_str(),
				hist_axes[ax0].size()-1, &hist_axes[ax0][0],
				hist_axes[ax1].size()-1, &hist_axes[ax1][0]);
		// reverse lookup
		type_names[ih.second] = ih.first;
	}

	std::string chain, friends, tree_name, friend_name;
	cd.Get("tree_input", chain);
	cd.Get("tree_name", tree_name);
	cd.Get("friend_input", friends);
	cd.Get("friend_name", friend_name);

	TChain *dm = new TChain(tree_name.c_str());
	TChain *fr = new TChain(friend_name.c_str());

	dm->Add(chain.c_str());
	fr->Add(friends.c_str());

	dm->AddFriend(friend_name.c_str());

	int ipnu, mode, itype;
	float dirnu[3], dir[3], flxho[3];
	float pnu, amom, weightx, ErmsHax, nEAveHax;

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
	//dm->SetBranchAddress("ErmsHax",  ErmsHax  );
	//dm->SetBranchAddress("nEAveHax", nEAveHax );

	dm->SetBranchStatus("*", 0);
	dm->SetBranchStatus("dirnu", 1);
	dm->SetBranchStatus("pnu", 1);
	dm->SetBranchStatus("ipnu", 1);

	// loop first...
	// zenith angle will be always between -1 and 1
	// find max and min values of true energy information
	// group by pdg value
	//std::map<int, double> E_min, E_max;
	//for (int i = 0; i < dm->GetEntries(); ++i) {
	//	dm->GetEntry(i);

	//	pnu = log10(pnu) + 3;
	//	if (!E_min.count(ipnu) || E_min[ipnu] > pnu)
	//		E_min[ipnu] = pnu;
	//	if (!E_max.count(ipnu) || E_max[ipnu] < pnu)
	//		E_max[ipnu] = pnu;
	//}
	//for (const auto &it : E_min)
	//	std::cout << "energy lim " << it.first
	//		  << " = " << it.second
	//		  << " : " << E_max[it.first] << std::endl;
	// type_names is map type(int) and hist name(str) pairs

	// ready to fill histograms
	// create matrices for bin contents

	std::vector<int> nu_types = {-16, -14, -12, 12, 14, 16};
	int trueX, trueY;
	if (!cd.Get("true_bins_energy", trueX))
		trueX = 100;
	if (!cd.Get("true_bins_cosz", trueY))
		trueY = 50;

	// type names is <int, string>
	for (const auto &ir : type_names) {
		// 100 bins in true neutrino energy and 50 bins in cosz
		std::string hname = "true_" + ir.second;
		true_bins[ir.second] = new TH2D(hname.c_str(), hname.c_str(),
					trueX, 1, 8, //it.second, E_max[it.first],
					trueY, -1, 1);	// -1 < cosz < 1
		// granularity of true bins is fixed
		// but eventually should be defined in card file

		int ny = reco_bins[ir.second]->GetNbinsX() * 
			 reco_bins[ir.second]->GetNbinsY();

		int nx = trueX * trueY;

		std::cout << "hist " << ir.second << "\t" << ny << ", " << nx << std::endl;
		for (int it : nu_types) {
			std::string name = Nu::toString(it);
			// one bin content per neutrino initial flavour
			if (it > 0) {
				bin_contents[ir.second + "_nuE0_" + name]
						= Eigen::MatrixXd::Zero(ny, nx);
				bin_contents[ir.second + "_nuM0_" + name]
						= Eigen::MatrixXd::Zero(ny, nx);
			}
			else {
				bin_contents[ir.second + "_nuEB_" + name]
						= Eigen::MatrixXd::Zero(ny, nx);
				bin_contents[ir.second + "_nuMB_" + name]
						= Eigen::MatrixXd::Zero(ny, nx);
			}
		}
	}

	// maps are sorted
	// type names is <int, string>, makes event by event look up easier
	int first_event_type = type_names.begin()->first;
	int last_event_type = type_names.rbegin()->first;

	double scale;
	if (!cd.Get("MC_scale", scale))
		scale = 188.4 / 22.5;	// HK to SK ratio
	dm->SetBranchStatus("*", 1);
	// loop through all events and fill histograms
	for (int i = 0; i < dm->GetEntries(); ++i) {
		dm->GetEntry(i);	// all branches
		itype -= 1;
		weightx *= scale;

		// uknown event types
		if (itype < first_event_type || itype > last_event_type)
			continue;

		// no NC tau's allowed, so skip
		if (std::abs(mode) > 30 && std::abs(ipnu) == 16)
			continue;

		// this is real data, which we do not want
		if (pnu < 1.0e-7 && ipnu == 0 && mode == 0)
			continue;

		// using these four variables to determine bin and EventType
		//pnu  > 0 ? log10(pnu) + 3 : 0;	// GeV
		//amom > 0 ? log10(amom)    : 0;	// already in GeV
		//cosnu = dirnu[2];
		//coslp = dir[2];

		std::string type = type_names[itype];

		// global bin of event in true and reco histograms
		int glo_nr = reco_bins[type]->FindFixBin(log10(amom), dir[2]);
		int glo_nt = true_bins[type]->FindFixBin(log10(pnu)+3, dirnu[2]);

		// total number of bins on x direction to find the absolut bin
		int r_nx = reco_bins[type]->GetNbinsX() + 2;
		int r_ny = reco_bins[type]->GetNbinsY() + 1;
		int t_nx = true_bins[type]->GetNbinsX() + 2;
		int t_ny = true_bins[type]->GetNbinsY() + 1;

		// means global bin is underflow or overflow, so event is discarded
		if (glo_nr < r_nx || glo_nr > r_nx * r_ny - 1
		 || glo_nr % r_nx == 0 || glo_nr % r_nx == r_nx - 1)
			continue;
		if (glo_nr < t_nx || glo_nr > t_nx * t_ny - 1
		 || glo_nr % t_nx == 0 || glo_nr % t_nx == t_nx - 1)
			continue;

		// absolute bins of the matrix
		// rows of matrix
		int loc_nr = glo_nr - (r_nx - 1 + 2 * (glo_nr / r_nx));
		// columns of matrix
		int loc_nt = glo_nt - (t_nx - 1 + 2 * (glo_nt / t_nx));


		std::string name = Nu::toString(ipnu);

		double factor_E = 1, factor_M = 1;
		if (name.find("nuE") != std::string::npos)	// it is e type
      			factor_E = flxho[1] / flxho[0];
		else						// it is mu or tau type
      			factor_M = flxho[0] / flxho[1];

		if (ipnu > 0) {	// neutrino
			bin_contents[type + "_nuE0_" + name](loc_nr, loc_nt)
							+= weightx * factor_E;
			bin_contents[type + "_nuM0_" + name](loc_nr, loc_nt)
							+= weightx * factor_M;
		}
		else {		//antineutrino
			bin_contents[type + "_nuEB_" + name](loc_nr, loc_nt)
							+= weightx * factor_E;
			bin_contents[type + "_nuMB_" + name](loc_nr, loc_nt)
							+= weightx * factor_M;
		}
	}

	// save histograms and matrices
	std::string outpath, cardpath;
	cd.Get("output_path", outpath);
	cd.Get("card_path", cardpath);
	if (outpath.find(".root") == std::string::npos)
		outpath += ".root";
	if (cardpath.find(".card") == std::string::npos)
		cardpath += ".card";
	// type names is <int, string>
	for (const auto &ir : type_names) {

		std::string outf = outpath;
		outf.insert(outf.find(".root"), ir.second); 
		std::string card = cardpath;
		card.insert(card.find(".card"), ir.second); 

		TFile outp(outf.c_str(), "RECREATE");
		std::ofstream cf(card.c_str());

		// name is histogram type
		reco_bins[ir.second]->Write();
		cf << "reco_info\t\"" << reco_bins[ir.second]->GetName() << "\"\n";
		// name is histogram type + "_" + nu flavour
		true_bins[ir.second]->Write();
		cf << "true_info\t\"" << true_bins[ir.second]->GetName() << "\"\n";

		cf << "type\t\"" << ir.second << "\"\n" << std::endl;
		for (const auto &ib : bin_contents) {
				continue;
			TMatrixD *content = new TMatrixD(ib.second.rows(),
							 ib.second.cols(),
							 ib.second.data());

			content->Write(ib.first.c_str());
		}

		cf << "\n\nreco_path\t\"" << outf << "\"" << std::endl;

		outp.Close();
		cf.close();
	}
	// ROOT should automatically delete memory of pointer
	// when closing file, unless something was not written

	// just to be safe, delete everything
	for (auto &ih : reco_bins)
		delete ih.second;
	for (auto &ih : true_bins)
		delete ih.second;

	bin_contents.clear();
}
*/
