/* Even parser baseclass
 */

#ifndef EventParser_H
#define EventParser_H

#include "Event.h"
#include "EventParser.h"

class SKEventParser : protected EventParser
{
	public:
		SKEventParser();
		SetTrueFlag(bool flag);

	private:
		void GetData(DataManager *&dm);
		Event::Type GetVars(std::vector<double> &pv, double &w);

		double *ipnu,
		       *mode,
		       *ip,
		       *itype,
		       *dirnu,
		       *pnu,
		       *dir,
		       *amom,
		       *flxho,
		       *weightx,
		       *wall,
		       *towall,
		       *ErmsHax,
		       *nEAveHax,
		       *Erms,
		       *nEAve,
		       *NeighborEHax,
		       *NeighborE;
};

#endif
