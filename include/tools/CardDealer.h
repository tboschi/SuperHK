/* CardDealer
 * handles cards
 */

#ifndef CardDealer_H
#define CardDealer_H

#include <cstdlib>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Utils.h"

class CardDealer
{
	public:
		CardDealer(std::string cardFile);

		bool Status();
		bool ReadCard(const std::string cardFile);

		bool Find(const std::string key);

		bool Get(const std::string key, std::vector<double>& vv);
		bool Get(const std::string key, double &dd);
		bool Get(const std::string key, std::string & ss);
		bool Get(const std::string subKey, std::map<std::string, std::vector<double> >& mv);
		bool Get(const std::string subKey, std::map<std::string, double>& mv);
		bool Get(const std::string subKey, std::map<std::string, std::string>& ms);

		std::vector<std::string> ListKeys();

	private:
		std::map<std::string, std::string> mKeyString;
		std::map<std::string, std::string>::iterator ims;

		std::map<std::string, std::vector<double> > mKeyDouble;
		std::map<std::string, std::vector<double> >::iterator imd;

		bool kVerbosity, kStatus;
};

#endif
