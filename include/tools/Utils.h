#ifndef Utils_H
#define Utils_H

namespace Utils
{
	static const bool verbosity = false;

	static void SeedBlossom(int & Seed, int & Blossom,
			 unsigned int instance,
			 unsigned int entries,
			 unsigned int processes)
	{
		unsigned int ppf = entries / processes;		//points per process
		unsigned int ppr = entries % processes;		//remainder

		if (instance < ppr)
		{
			Seed    =  instance      * (ppf + 1);
			Blossom = (instance + 1) * (ppf + 1);
		}
		else 
		{
			Seed    =  instance      * ppf + ppr;
			Blossom = (instance + 1) * ppf + ppr;
		}
	}
};

#endif
