#include "EventParser.h"

EventParser::EventParser()
{
}

/* function to be called externally
 * pass a datamanager, containing files and tree (picked up by a GetEntry, external)
 * outputs is a vector of variables to be stacked and the event weight
 * output is event type to fill the correct stack
 */

Event::Type EventParser::Parse(DataManager *&dm, std::vector<double> &pVars, double &w)
{
	GetData(dm);				//load data from tree

	return et = GetVars(pVars, w);		//extract bin variables from data
}

void EventParser::GetData(DataManager *&dm)
{
	/* ...
	 * do somemagic assigning
	 * like
	 * dm->Get("varname", var);
	 */
}

Event::Type EventParser::GetVars(std::vector<double> &pv, double &w)
{
	/* ...
	 * do somemagic assigning
	 * like
	 * pv.push_back(var1);
	 * pv.push_back(var2);
	 *
	 * find the event type and return it
	 */
}

void EventParser::LoadMCWeight( std::string filename)	//used only in beam numu EP
{
	vmcWeight.clear();

	std::ifstream inFile(filename.c_str());		//SK (atmo?) flux data in flux/cm^2 at SK
	std::string line;

	while (std::getline(inFile, line))
	{
		std::stringstream ssl(line);
		std::vector<double> vLine;

		double val;
		while (ssl >> val)
			vLine.push_back(val);

		vmcWeight.push_back(vLine);
	}
}
