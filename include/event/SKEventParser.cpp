#include "SKEventParser.h"

SKEventParser::SKEventParser()
{
}

void SKEventParser::SetTrueFlag(bool flag)
{
	trueFlag = flag;
}

Event::Type SKEventParser::GetVars(std::vector<double> &pv, double &w)
{
	Event::Type et = *itype - 1;		//event type
	pv.push_back(*pnu);			//energy always first element

	if (et < Event::SubGeV_elike_0dcy || et > Event::UpThruShower_mu)	//event not classified
		return Event::undefined;
	else
	{
		pv.clear();

		if (trueFlag)
		{
			pv.push_back(*pnu);		//nu energy
			pv.push_back(-1.0 * dirnu[2]);	//nu cos(theta)
		}
		else
		{
			pv.push_back(*amom);		//lepton momentum
			pv.push_back(-1.0 * dir[2]);	//lepton cos(theta)
		}

		//if nu energy is 0, has no pdg code or mode, then it is a data file
		if (pnu == NULL && ipnu == NULL && mode == NULL )
			w = -1.0;
		else						//else it is MC
		{
			w = *weightx;			//MC weight
			if (useFSIWeight)
				 w *= *fsiweight;
		}

		return et;
	}
}

/*
	SkipFlag = false;
	PDG = 0;

   // Event(Vertex) Paremeters 
	MCweight = weightx;
	if( kUseFSIWeight )
	{
	   MCweight*= fsiweight;
	   //std::cout << weightx(0) << " " << fsiweight(0) << std::endl;
	}
	
	Mode     = mode;
	
	// Parent Neutrino Parameters 
	NuCosineTheta = -1. * dirnu[2]; // + downgoing , - upgoing 
	PDG           = *ipnu;
	
	// Lepton Parameters 
	LeptonP = amom(0);
	CosineTheta = -1.*dir(2);  // + downgoing , - upgoing 
	
	EventType = itype(0) - 1;	//explain this
	
	// hax rvw 
	//
	  EventType = EventType % global::nEventTypes; 
	// now convert to SK-IV
	  EventType += global::SK4;
	//
	// end hax
	
	
	if ( EventType < 0 || EventType > global::EndOfBinTypes )
	   SkipFlag = true;
	
	// no NC tau's allowed!
	if ( abs(Mode) > 30  && abs(PDG) == 16 ) SkipFlag = true; 
	
	  else if (abs(Mode) < 30 && abs(PDG) == 14)	//CCNumu 
	{
	int elikesample[6]={0,1,7,8,10,11};
	bool elike=false;				//Elike sample
	
	for (int i = 0; i < 6; i++)			
	{
	     if ( EventType % global::nEventTypes == elikesample[i])
	     {
	     	elike = true;
	     	break;
	     }
	}
	
	if (elike)
	{
	     MCweight *= 0.5;	//Half Mu in E sample  	
	}
	}
	
	else if ( abs(PDG) == 16)
	{
		MCweight *= 0.5;	//Half tau
	  	MCweight *= 0;	//No tau
	}
	
	
	
	// Hack to convert old SK3 MC types to SK4 MC Types
	//   if ( MonteCarlo )
	//     {
	//       if ( EventType >= 36 && EventType <= 53 )
	//   	 EventType = EventType + 18;
	//     }
*/


/*
double SKEventParser::GetHondaFluxRatio( int NuType )
{

   if ( NuType == 1 )   
      return (double) flxho(1) / flxho(0);

   if ( NuType == 2 || NuType == 3 )
      return (double) flxho(0) / flxho(1);

   std::cerr << "Returning 0 from SKEventParser::GetHondaFluxRatio for type " << NuType << std::endl;
   return 0;

}
*/


void SKEventParser::GetData(DataManager *&dm)
{
	dm->Get("ipnu",		ipnu);
	dm->Get("mode",		mode);
	dm->Get("ip",		ip);
	dm->Get("itype",	itype);
	
	dm->Get("dirnu",	dirnu);
	dm->Get("pnu",		pnu);
	dm->Get("dir",		dir);
	dm->Get("amom",		amom);
	dm->Get("flxho",	flxho);
	dm->Get("weightx",	weightx);
	
	dm->Get("wall",		wall);
	dm->Get("towall",	towall);
	
	// only if friends
	dm->Get("ErmsHax",	ErmsHax);
	dm->Get("nEAveHax",	nEAveHax);
	
	dm->Get("Erms",		Erms);
	dm->Get("nEAve",	nEAve);
	
	dm->Get("NeighborEHax", NeighborEHax);
	dm->Get("NeighborE",	NeighborE);

/*
	kUseFSIWeight = false ;
	if( dm->Existence( "fsiweight" ) )
	{
		dm->Get("fsiweight" , fsiweight );
//		kUseFSIWeight = true ;
	}
*/
}


/*
void SKEventParser::GetEventBinValues( std::vector<double> & vars)
{
   double logp;
   if( (useTrue & 1) == 1){
	   logp = ( GetEnergy() != 0 ? log10( GetEnergy()*1000 ): 0 );// pnu is saved as GeV
   }
   else
   {
   logp = ( GetLeptonP() != 0 ? log10( GetLeptonP() ): 0 );
   }
   vars.push_back( logp         );

   if((useTrue & 2) == 2)
   {
	   vars.push_back( GetNuCosineZ() );
   }
   else
   {
   vars.push_back( GetCosineZ() );
   }
}          
*/
