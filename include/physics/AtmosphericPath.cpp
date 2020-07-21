#include "AtmosphericPath.h"

AtmosphericPath::AtmosphericPath(CardDealer *cd)
{
	LoadTable();

	if (!cd->Get("tank_height", _tank_height))
		_tank_height = 0.76;	// whatever number that is

	if (!cd->Get("MantleAveraging", kMantleAveraging))
		kMantleAveraging = 0;
	if (!cd->Get("mantle_phases", mantle_phases)) {
		if (kMantleAveraging)
			mantle_phases = 4;
		else
			mantle_phases = 1;
	}

	// initialise Earth densities
	_density = { std::make_pair(Const::Rearth,  3.3),	//6371
		     std::make_pair(5701.0,  5.0),
		     std::make_pair(3480.0, 11.3),
		     std::make_pair(1220.0, 13.0),
		     std::make_pair(0 ,     13.0) };
}


void AtmosphericPath::LoadTable()
{
	std::string table_file;
	cd->Get("table_file", table_file);
	std::ifstream it(table_file.c_str());

	while (std::getline(it, line)) {
		std::stringstream iss(line);
		double ent;
		iss >> ent;	// discard first
		std::vector<double> h3d_row;
		while (iss >> ent)
			h3d_row.push_back(ent);
		h3d.push_back(h3d_row);
	}
}


double AtmosphericPath::OscillationProbability(Oscillator *osc, EventParser* evt)
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


// function returns number of layers of Earth the neutrino transverses
std::vector<std::pair<double, double> > AtmosphericPath::SetDensityProfile
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
	for (const auto &id : _density)
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
		lens_dens.push_back(std::make_pair(rho, layer_distance)); //* 1.e5 in [cm]
	}

	// make it symmetrical by aziumth
	for (int i = layers - 1; i > 0; --i)
		lens_dens.push_back(lens_dens[i-1]);

	return lens_dens;	// which is lens_dens.size() - 1
}

Eigen::MatrixXd AtmosphericPath::AverageOscillationMantle(double cosz, double energy, double path_length)
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


/* Return (random) baseline at given costheta and neutrino energy,
 * for a detector which is a distance d in km above the surface
 * of the earth.
 */
// this if from nebaseline_h3d.F
std::vector<double> AtmosphericPath::GenerateBaselines(std::vector<double> &heights, double cosz)
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
