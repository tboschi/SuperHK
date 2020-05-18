/* Even parser baseclass
 */

#ifndef EventParser_H
#define EventParser_H

#include <vector>
#include <string>

#include "Event.h"
#include "DataManager.h"

class EventParser
{
	public:
		EventParser();
		Event::Type Parse(DataManager *&dm, std::vector<double> &pVars, double &w);
		void LoadMCWeight( std::string filename);	//used only in beam numu EP

	protected:
		virtual void GetData(DataManager *&dm);
		virtual Event::Type GetVars(std::vector<double> &pv, double &w);

		EventType eventType;
		bool skipFlag;

		std::vector<std::vector<double> > mcWeight;

	private:
};

#endif
