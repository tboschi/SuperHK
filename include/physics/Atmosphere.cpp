#include "Atmosphere.h"

Atmosphere::Atmosphere(const std::string &card) :
	gen(std::random_device()())
{
	CardDealer cd(card),
	LoadDensityProfile(cd);
	LoadProductionHeights(cd);

	if (!cd.Get("verbosity", kVerbosity))
		kVerbosity = 0;
}

Atmosphere::Atmosphere(const CardDealer &cd) :
	gen(std::random_device()())
{
	LoadDensityProfile(cd);
	LoadProductionHeights(cd);

	if (!cd.Get("verbosity", kVerbosity))
		kVerbosity = 0;
}


Atmosphere::Atmosphere(CardDealer *cd) :
	gen(std::random_device()())
{
	LoadDensityProfile(*cd);
	LoadProductionHeights(*cd);

	if (!cd->Get("verbosity", kVerbosity))
		kVerbosity = 0;
}

// loads production heights calculated by Honda (HKKM2014) 
// from http://www.icrr.u-tokyo.ac.jp/~mhonda/
// The tables are per neutrino flavour, production site
// They are matrices in zenith angle and neutrino energy
// giving the cumulative proability distribution for 
// production height. See Honda's website for details
//
// In order to generate a production height, simply generate
// a random number [0, 1] and find the equivalent index in
// _problibs vector -> the respective entry in the (zenith, energy)
// vector gives the production height
void Atmosphere::LoadProductionHeights(const CardDealer &cd)
{
	if (!cd.Get("production_height", _atm))
		_atm = 15.0;	// default 15 km altitude

	// vector with the discretised probabilities
	_problibs.clear();
	// vector with the discretised log energy bins !
	_energies.clear();

	// calling system ls such that table_file can contain wildcards
	std::vector<std::string> prod_files;
	if (!cd.Get("honda_production", prod_files)) {
		std::cerr << "Atmosphere: no production heights in card, please specify production height yourself" << std::endl;
		return;
	}


	// vector with the discretised cos zenith angles
	_zenithas.assign(20, 0);
	std::iota(_zenithas.begin(), _zenithas.end(), -9);
	std::transform(_zenithas.begin(), _zenithas.end(), _zenithas.begin(),
			[](double b) { return -0.1 * b; }); 
	// now zenith_bin is 0.9, 0.8, ... -0.9, -1.0	according to Honda binning

	std::map<int, std::vector<double> > *table;	// "reference"
	for (const std::string &prod : prod_files) {
		std::ifstream ip(prod.c_str());
		if (kVerbosity)
			std::cout << "Loading production file " << prod << std::endl;
		std::string line;

		//int iz = 0;	// zenith index
		int ie = 0;	// energy index

		bool fill_energy = false;
		while (std::getline(ip, line)) {
			// remove hash comments
			if (line.find('#') != std::string::npos)
				line.erase(line.find('#'));
			if (line.empty())
				continue;

			std::stringstream iss(line);

			std::vector<double> row;
			double ent;
			while (iss >> ent)
				row.push_back(ent);

			if (row.size() == 25) {	// header
				fill_energy = !_energies.size();
				if (int(row[0]) != 1)
					continue;
				switch (int(row[1])) {	// nu flavor
					case 1: // nu mu
						table = &nuM0_table;
						break;
					case 2: // nu mu bar
						table = &nuMb_table;
						break;
					case 3: // nu e
						table = &nuE0_table;
						break;
					case 4: // nu e bar
						table = &nuEb_table;
						break;
				}

				//iz = int(row[2]) - 1;

				// row[3] is azimuthal bin, but it is averaged over

				if (!_problibs.size())
					_problibs.assign(row.begin() + 4, row.end());
			}
			else if (row.size() == 22) { // body with energies and cumulative distributions

				if (fill_energy)	// store energy bin information
					_energies.push_back(log10(row[0]));

				std::transform(row.begin()+1, row.end(), row.begin() + 1,
						[](double b) { return b / 1.e3; }); // m -> km
				//int global_bin = iz * (row.size() - 1) + ie;
				(*table)[ie].assign(row.begin()+1, row.end());
				++ie;
			}
			else
				std::cerr << "did you tamper with Honda's inputs?\n";
		}
		ip.close();
	}

	//cmd = "rm .production_files";
	//system(cmd.c_str());
}


// Load matter density profile from file
// the profile is saved as an array[3] of at least 2 elements, but it can be 3 too
// columns above 3 are discarded
// the elements are
// 	radius, which can be normalised
// 	density 
// 	electron density (optional, default 0.5)
void Atmosphere::LoadDensityProfile(const CardDealer &cd, std::string table_file)
{
	_profile.clear();

	if (table_file.empty())
		if (!cd.Get("density_profile", table_file))
			throw std::invalid_argument("Atmosphere: you have not defined any density profile. Too bad!");

	std::ifstream it(table_file.c_str());
	if (!it.good()) {
		std::cerr << "Density profile " << table_file << " cannot be read" << std::endl;
		return;
	}

	// read first line and remember last line
	std::string line;
	double first_rad = -1, last_rad = -1;
	while (std::getline(it, line)) {
		// remove hash comments
		if (line.find('#') != std::string::npos)
			line.erase(line.find('#'));
		if (line.empty())
			continue;

		std::stringstream iss(line);
		if (first_rad < 0)
			iss >> first_rad;
		iss >> last_rad;
	}

	// last radius is smaller than first one, so it is reversed!
	bool kReverse = last_rad < first_rad;
	// last radius is less than one, meaning it is absolute scale
	bool kAbsolute = (last_rad > 1);

	// reset!
	it.clear();
	it.seekg(0);
	while (std::getline(it, line)) {
		// remove hash comments
		if (line.find('#') != std::string::npos)
			line.erase(line.find('#'));
		if (line.empty())
			continue;

		std::stringstream iss(line);

		double ent;
		std::vector<double> row;
		row.reserve(3);	// three columns needed

		while (iss >> ent)
			row.push_back(ent);

		// not absolute scale
		if (!kAbsolute)
			row[0] *= Const::EarthR;
		if (row.size() < 3)
			row[2] = 0.5;	// electron density default

		_profile.push_back({row[0], row[1], row[2]});
	}

	if (kReverse)
		std::reverse(_profile.begin(), _profile.end());

	it.close();
}

// generate a random production height given
// 	flv -> initial neutrino flavour
// 	cosine zenith
// 	energy
double Atmosphere::RandomHeight(Nu::Flavour flv, double cosz, double energy)
{
	// no random heights file loaded
	if (!_problibs.size())
		return -1;

	energy = log10(energy);
	if (energy >= _energies.back()) // error energy >~ 1e4
		return -1;
	if (cosz < _zenithas.back())	// error cosz < -1
		return -1;

	// get reference to correct table
	std::map<int, std::vector<double> > *table;	// "reference"
	switch (flv) {	// nu flavor
		case Nu::M_: // nu mu
		case Nu::T_: // nu mu
			table = &nuM0_table;
			break;
		case Nu::Mb: // nu mu bar
		case Nu::Tb: // nu mu bar
			table = &nuMb_table;
			break;
		case Nu::E_: // nu e
			table = &nuE0_table;
			break;
		case Nu::Eb: // nu e bar
			table = &nuEb_table;
			break;
		default:
			throw std::invalid_argument("Undefine neutrino flavour\n");
			break;
	}

	// define a probability interval first
	std::uniform_real_distribution<> uniform(0, 1);//uniform distribution between 0 and 1
	double prob = uniform(gen);
	auto ip = std::lower_bound(_problibs.begin(), _problibs.end(), prob);
	
	// energy is ordered increasing
	auto ie = std::lower_bound(_energies.begin(), _energies.end(), energy);
	// zenith is ordered decreasing
	auto ic = std::lower_bound(_zenithas.begin(), _zenithas.end(), cosz, std::greater<double>());

	// iterators are the least elements that compare to input
	// if iterator points to begin, move to second
	// if iterator points to end, move them to last valid
	if (ip == _problibs.begin())
		++ip;
	else if (ip == _problibs.end())
		--ip;

	if (ie == _energies.begin())
		++ie;
	else if (ie == _energies.end())
		--ie;

	if (ic == _zenithas.begin())
		++ic;
	else if (ic == _zenithas.end())
		--ic;

	//if (ip == _problibs.end() || ip == std::prev(_problibs.end()))
	//	ip = _problibs.end() - 2;
	//if (ie == _energies.end() || ie == std::prev(_energies.end()))
	//	ie = _energies.end() - 2;
	//if (ic == _zenithas.end() || ic == std::prev(_zenithas.end()))
	//	ic = _zenithas.end() - 2;

	// index positions
	int pn = std::distance(_problibs.begin(), ip);
	int en = std::distance(_energies.begin(), ie);
	int cn = std::distance(_zenithas.begin(), ic);

	// Trilinear interpolation
	// https://en.wikipedia.org/wiki/Trilinear_interpolation
	// with respect to following element
	double pd = (prob   - *ip) / (*std::prev(ip) - *ip);
	double ed = (energy - *ie) / (*std::prev(ie) - *ie);
	double cd = (cosz   - *ic) / (*std::prev(ic) - *ic);
	// meaning of values [0,1] is proximity to following element

	// x direction is prob,
	// y direction is energy
	// z direction is cosz
	// Production height is a function of
	// 	probability (for PDF), energy, cosz

	double h000 = (*table)[_energies.size() *  cn    + en  ][pn  ];
	double h001 = (*table)[_energies.size() *  cn    + en  ][pn-1];
	double h010 = (*table)[_energies.size() *  cn    + en-1][pn  ];
	double h011 = (*table)[_energies.size() *  cn    + en-1][pn-1];
	double h100 = (*table)[_energies.size() * (cn-1) + en  ][pn  ];
	double h101 = (*table)[_energies.size() * (cn-1) + en  ][pn-1];
	double h110 = (*table)[_energies.size() * (cn-1) + en-1][pn  ];
	double h111 = (*table)[_energies.size() * (cn-1) + en-1][pn-1];

	double h00 = h000 * (1 - pd) + h100 * pd;
	double h01 = h001 * (1 - pd) + h101 * pd;
	double h10 = h010 * (1 - pd) + h110 * pd;
	double h11 = h011 * (1 - pd) + h111 * pd;

	double h0 = h00 * (1 - ed) + h10 * ed;
	double h1 = h01 * (1 - ed) + h11 * ed;

	// this is the interpolated value
	return h0 * (1 - cd) + h1 * cd;
}

// calculate the matter profile for neutrino passage
// 	cosz = cosine of zenith angle)
// 	atm  = production height in atmosphere
// return a vector of length-density pairs in order of neutrino travel baseline
// the first entry will be therefore atmospheric propagation
Oscillator::Profile Atmosphere::MatterProfile(Nu::Flavour flv, double cosz, double energy)
{
	return MatterProfile(cosz, RandomHeight(flv, cosz, energy));
}

Oscillator::Profile Atmosphere::MatterProfile(double cosz, double atm)
{
	if (atm < 0)
		atm = _atm;

	if (kVerbosity > 4)
		std::cout << "Atmosphere: generating profile for " << cosz << " from " << atm << " km\n";

	double sinz = sqrt(1 - cosz * cosz);

	// path travellend in air
	double x0 = sqrt(pow(Const::EarthR + atm, 2) - pow(Const::EarthR * sinz, 2))
			- Const::EarthR * std::abs(cosz);
			//- Const::EarthR * (cosz);

	Oscillator::Profile lens_dens = {{x0, 0 /*Const::rhoAir*/, 0.5}};
	// neutrino from above, so no Earth crossing
	if (cosz > -1e-9)
		return lens_dens;


	Oscillator::Profile halves;
	// minimum distance of path to centre of Earth
	double dist = Const::EarthR * sinz;
	double x_prev = 0;
	// find deepest radius touched by neutrino
	auto ir = std::lower_bound(_profile.begin(), _profile.end(), dist,
				  [](const Oscillator::LDY &ldy, double v)
				  	{ return ldy[0] < v; });
	for ( ; ir != _profile.end(); ++ir) {
		Oscillator::LDY &ldy = *ir;
		if (std::abs(ldy[0] - dist) < 1e-9)
			continue;
		// find track length inside this shell
		double x_n = sqrt(pow(ldy[0], 2) - pow(dist, 2)) - x_prev;
		halves.push_back({x_n, ldy[1], ldy[2]});
		x_prev += x_n;
	}

	if (halves.size()) { // double the deepest length
		halves.front()[0] *= 2;

		// halves is 2*xn, xn-1, ... x2, x1
		lens_dens.insert(lens_dens.end(), halves.rbegin(), halves.rend());
		lens_dens.insert(lens_dens.end(), halves.begin() + 1, halves.end());
	} //else shell passed is very very small

	// lens_dens is now x0, x1, x2, ..., xn-1, xn, xn, 
	return lens_dens;
}

/*
std::map<std::string, Eigen::VectorXd> Atmosphere::Oscillate(
			const std::vector<std::pair<double, double> > &bins, 
			const std::map<std::string, std::pair<Nu::Flavour, Nu::Flavour> > &oscf,
			Oscillator *osc)
{
	std::map<std::string, Eigen::VectorXd> allprob;
	for (const auto &io : oscf) {
		if (osc)
			allprob[io.first] = Eigen::VectorXd::Zero(bins.size());
		else
			allprob[io.first] = Eigen::VectorXd::Ones(bins.size());
	}

	if (!osc)
		return allprob;

	// bins is pair of cosz and energy
	int n = 0;
	for (const auto &ib : bins) {
		osc->SetMatterProfile(MatterProfile(ib.first));
		for (const auto &io : oscf) {
			const auto &flav = io.second;
			allprob[io.first](n) = osc->Probability(flav.first, flav.second, //flavs
						     		ib.second);		 //energy
		}
		++n;
	}

	return allprob;
}
*/


////                    WEIRD STUFF HERE BELOW

/*
double Atmosphere::OscillationProbability(Oscillator *osc, EventParser* evt)
{
	// NuOscillatedTo is final state flavour, the one found in Event E
	//  and it is absolute value, so don't dinstigunsh between neutrino anti neutrino
	// MCWeight is weight from montecarlo simulation
	// FactorE and FactorMu are taken from hondafluxratio
	// Converts PDG->internal OscType    e:12->1 , mu: 14->2, tau:16->3 //

	Nu nu_type;
	bool kAnti = evt->GetPDG() > 0;	// true is particle, false otherwise
	switch (evt->GetPDG()) {
		case 12:
			nu_type = Nu::E_;
			break;
		case 14:
			nu_type = Nu::M_;
			break;
		case 16:
			nu_type = Nu::T_;
			break;
		case -12:
		        nu_type = Nu::Eb;
		        break;
		case -14:
		        nu_type = Nu::Mb;
		        break;
		case -16:
			nu_type = Nu::Tb;
			break;
	}

	double MC_weight = evt->GetMCWeight();

	// NC events only recieve their solar flux weighting
	// also zero out any erroneou NC Tau events that may be included
	if ( abs( evt->GetMode() ) >= 30) {
		if (nu_type == Nu::T_ || nu_type == Nu::Tb)
			return 0;
		return MC_weight;
	}

	double cosz = evt->GetNuCosineZ();                  
	double energy = evt->GetEnergy();
	double rms = std::min(0.5*energy, evt->GetEAveRMS());

	std::vector<double> energy_avg = { energy, energy - rms, energy - rms *0.5 , 
					   energy + rms *0.5, energy + rms }; 
	if (evt->GetNEAve()) {	// not sure what this is
		double fac = 1.5;
		for (int i = 1; i < energy_avg.size(); ++i) {
			energy_avg[i] *= fact;
			fact -= i * 0.25;
		}
	}


	//Fortran::nebaseline_h3d_(XXPATH, IDNU2, DirNeu, fortranE, SKheight);
	std::vector<double> heights = GenerateHeights(nu_type, cosz, energy);
	std::vector<double> nupaths = GenerateBaselines(heighs, cosz);

	std::vector<double> path_avg(5); // * 1.0e5  in [cm]
	path_avg[0] = std::accumulate(nupaths.begin(), nupaths.end(), 0.0) / 20.;

	// Original Fortran Indexing -- is this still needed???
	std::vector<int> index_order = {10, 9, 11, 12, 8, 7, 13, 14, 6, 5, 15, 16, 4, 3, 17, 18};
	for(int i = 0 ; i < index_order.size(); ++i)
		avePath[i % 4 + 1] += nupaths[index_order[i]-1] / 4.; // *1.e5 cm

	// with cosz, ProdHeight = 15.
	// above fucntion does this -> 
	double production = 15.0e5; // 15km in cm
	double path_length = sqrt(pow(Const::Rearth + production, 2)
			       - pow(Const::Rearth, 2)*(1 - cosz*cosz)) - REarth*cosz;
	SetDensityProfile(cosz, path_length);
	// varaibles used in earth desnity profiling

	// Mass Hierarchy Flag 
	// ( correction for solar term is done internally )
	// double hFactor           = ( kInverted ? -1.0 : 1.0 ) ; 

	// assuming oscillator osc is already set! 
	//
	//bNu->SetMNS( P.Get("S12") , P.Get("S13") , P.Get("S23") , 
	//		P.Get("M12") ,      hFactor * P.Get("M23") , 
	//		P.Get("CP")  , 
	//		Energy       , kSinX2Vars   ,  evt->GetPDG() );  


	// Flux factor adjustment : 
	// if nue      ,  FactorE   = 1 , otherwise reweight from muonflux to electron flux
	// if numu/tau ,  FactorMu  = 1 , otherwise reweight from electron to muon flux 

	double factor_E = 1, factor_M = 1;
	if (nu_type != Nu::E_ && nu_type != Nu::Eb)	// it is mu or tau type
		factor_E = evt->GetHondaFluxRatio(2);
	else						// it is e type
		factor_M = evt->GetHondaFluxRatio(1);

	// skip energy averaging for 
	// upthrough types
	int nPaths = 5; 	// average on 5 lengths
	int event_type = evt->GetType() % global::nEventTypes;
	if( event_type == "UpThruNonShower_mu" ||
	    event_type == "UpThruShower_mu"  )  
		nPaths = 1 ;

	double prob = 0;	
	Nu nue_in = kAnti ? Nu::E_ : Nu::Eb;
	Nu num_in = kAnti ? Nu::M_ : Nu::Mb;
	for (int i = 0; i < nPaths ; ++i) // energy path length aver
	{
		Eigen::Matrix trans = AverageOscillationMantle(cosz, energy_avg[i], path_avg[i])

		prob += tran(nue_in, nu_type) * factor_E;
		prob += tran(num_in, nu_type) * factor_M;
	}

	prob *= MC_weight / nPaths;

	return prob;
}
*/

/*
// function returns number of layers of Earth the neutrino transverses
std::vector<std::pair<double, double> > Atmosphere::SetDensityProfile
					(double cosz, double path_length)	// should be in cm
{
	std::vector<std::pair<double, double> > lens_dens;	// pathlength - TotalEarthLegth
	//_TraverseRhos[0] = 0.0;
	//_TraverseDistance[0] =  path_length - TotalEarthLength ;

	if (cosz >= 0) {  
		lens_dens.push_back(std::make_pair(0, path_length)); // in [cm]
		return lens_dens;
		//_TraverseDistance[0] =  path_length;
		//return 1; 
	}

	int layers = 0;
	for (const auto &id : _profile)
		// first comparison should be (cosz < 0)
		if (cosz < - sqrt(1 - pow(id.first / Const::Rearth, 2)))
			++layers;

	// add first air component
	lens_dens.push_back(std::make_pair(0, path_length + 2.0 * cosz * Const::Rearth)); // 1e5 in [cm]

	// the zeroth layer is the air!
	for (int i = 0 ; i < layers; ++i) 
	{
		double rho = lens_dens[i].second;
		double rad = lens_dens[i].first;
		double rad_next = lens_dens[i+1].first;

		double layer_distance = sqrt(pow(rad, 2) - pow(Const::Rearth, 2) * (1 - cosz*cosz));

		if (i < layers-1)
			layer_distance -= sqrt(pow(rad_next, 2) - pow(Const::Rearth, 2) * (1 - cosz*cosz));
		else
			layer_distance *= 2.;

		//_TraverseDistance[i+1]  =  0.5*( CrossThis-CrossNext )*1.e5;
		//_TraverseDistance[i+1]  =  CrossThis*1.e5;
		lens_dens.push_back(std::make_pair(rho, layer_distance)); // * 1.e5 in [cm]
	}

	// make it symmetrical by aziumth
	for (int i = layers - 1; i > 0; --i)
		lens_dens.push_back(lens_dens[i-1]);

	return lens_dens;	// which is lens_dens.size() - 1
}

Eigen::MatrixXd Atmosphere::AverageOscillationMantle(double cosz, double energy, double path_length)
{
	//first element should be air length and density
	auto layers_dens = SetDensityProfile(cosz, path_length);

	//Layers = Earth->get_LayersTraversed( );

	int phases = 1;
	double phase_offset = 0;

	if (kMantleAveraging && layers_dens.size() > 2) { // Must both Average and be in Mantle
		// m23 contains delta m_23
		const Eigen::VectorXd &mm = osc->Masses();	//make reference
		double m23 = mm(2);	//only three masses
		if (mm23 < 0)
			m23 += mm(1);
		else
			m23 -= mm(1);
		// check penultimate layer (why?)
		if (layers_dens.end()[-2] > 4.96 * energy / std::abs(m23)) { // * 1.e5 for cm
			phasese = mantle_phases; // defined globally
			phase_offset = Const::pi / 2.0;
		}
	}

	osc->SetMatterProfile(lens_dens);
	Eigen::MatrixXd trans_avg = Eigen::MatriXd::Zero(3,3);
	for (int ip = 0; ip < phases; ++ip) {
		// multiply by some phase
		Eigen::MatrixXd trans = osc->TransitionMatrix(energy, ip * phase_offset);
		trans_avg += trans.cwiseAbs2() / phases;
	}

	return trans_avg;
}
*/






/* Return (random) baseline at given costheta and neutrino energy,
 * for a detector which is a distance d in km above the surface
 * of the earth.
 */
/*
// this if from nebaseline_h3d.F
std::vector<double> Atmosphere::GenerateBaselines(std::vector<double> &heights, double cosz)
{
	const double &R = Const::Rearth;
	double x = R + _tank_h;

	// try solve equaton a x^2 + b x + c = 0
	//double a = 1
	double b = 2 * x * cosz * (cosz < 0 ? 1 : -1); 
	double c = x * x - R * R;

	double Q = (cosz < 0 ? 1 : -1); 
	if (b * b < 4 * c)	// only true very near horizon; use tangent value there
		Q = sqrt(2 * _tank_h * x);
	else
		Q *= (-b + sqrt(b * b - 4 * c) * (cosz > 0 ? -1 : 1)) / 2.;

	double cosalpha = std::abs((Q*Q + - c) / (Q * R)) / 2.;

	std::vector<double> fpath;
	fpath.reserve(heights.size());
	for (double ih : heights) {
		double pp;
		ih = pow(R + ih, 2);

		if (std::abs(costheta) * x > sqrt(c))
			pp = sqrt(pow(R * cosalpha, 2) + pow(R + ih, 2) - R * R)
				- R * cosalpha + Q;
		else
			pp = -x * cosz +
		              sqrt(pow(x * costheta, 2) + pow(R + ih, 2) - x * x);
		fpath.push_back(pp)
	}

	return fpath;
}

// this if from neprodhgt_h3d.f
// too much hardcoded stuff, but don't know the meaning
std::vector<double> AtmospherichPath::GenerateHeights(Nu type, double cosz, double energy)
{
	// the following swtich statement is equivalent to this ifs seuqence
	// if ( E->GetPDG() ==  12 ) IDNU2 = 1;
	// if ( E->GetPDG() == -12 ) IDNU2 = 3;
	// if ( E->GetPDG() ==  14 ) IDNU2 = 2;
	// if ( E->GetPDG() == -14 ) IDNU2 = 4;
	// if ( E->GetPDG() ==  16 ) IDNU2 = 2;
	// if ( E->GetPDG() == -16 ) IDNU2 = 4;

	switch (type) {
		case Nu::E_:	//-12
			nutag = 0;
			break;
		case Nu::M_:	//14
		case Nu::T_:	//16
			nutag = 1;
			break;
		case Nu::Eb:	//-12
			nutag = 2;
			break;
		case Nu::Mb:	//-14
		case Nu::Tb:	//-16	//case 4: case -2:
			nutag = 3;
			break;
	}

	// limited between 0 and 19
	int costag = std::min(19, int((cosz + 1.0) * 10.));
	double delta_cos = cosz - (costag - 11.0) / 10.0 - 0.05;
	
	// limited between 0 and 30
	int enetag = std::max(0, std::min(30, int(10.0*log10(energy) + 10.0)));
	double de = pow(10, (enetag+1)/10. - 1.1) - pow(10, (enetag/10. - 1.1);
	double delta_e =  energy - pow(10, enetag/10. - 1.1));

	// hard coded....
	int kene = 1;
	int kcos = 30;
	int knut = 4;

	int jj = kene * (enetag + kcos * (costag + knut * nutag));
	std::vector<double> heights(20);
	for (int i = 0; i < 19; ++i)
	{
		double delta_h;
		if (enetag < 31 && delta_e > 0) {
			dhe = h3d[jj + kene][i] - h3d[jj][i];
			delta_h = delta_e * dhe / de;
		}

         	if ( (costag == 1 && delta_cos <= 0)
		  || (costag == 20 && delta_cos >= 0))
			dhc = 0.0;
		else {
			if (delta_cos >= 0.)
				dhc = h3d[jj + kcos][i] - h3d[jj][i];
			else
				dhc = h3d[jj - kcos][i] - h3d[jj][i];
			delta_h = std::abs(delta_cos) * dhc * 10.;
		}

		heights[i] = dhc = h3d[jj][i] + delta_h;
	}

	heights[20] = heights[19];

	return hgt;
}
*/
