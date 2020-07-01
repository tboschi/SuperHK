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
#include <algorithm>

// trim from start (in place)
static inline void ltrim(std::string &s)
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
				return !std::isspace(ch); }));
}

// trim from end (in place)
static inline void rtrim(std::string &s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
				return !std::isspace(ch); }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s)
{
	ltrim(s);
	rtrim(s);
}

class CardDealer
{
	public:
		CardDealer(std::string cardFile, bool verb = false);
		CardDealer(char* filename, bool verb = true);

		std::string CardName();

		bool Status();
		bool ReadCard(const std::string cardFile);

		bool Find(const std::string key);

		bool Get(const std::string key, int &ii);
		bool Get(const std::string key, double &dd);
		bool Get(const std::string key, std::string & ss);
		bool Get(const std::string key, std::vector<int>& vi);
		bool Get(const std::string key, std::vector<double>& vd);
		bool Get(const std::string key, std::vector<std::string>& vs);
		bool Get(const std::string key, std::map<std::string, int>& mi);
		bool Get(const std::string key, std::map<std::string, double>& mv);
		bool Get(const std::string key, std::map<std::string, std::string>& ms);
		bool Get(const std::string key, std::map<std::string, std::vector<int> >& mi);
		bool Get(const std::string key, std::map<std::string, std::vector<double> >& mv);
		bool Get(const std::string key, std::map<std::string, std::vector<std::string> >& ms);

		std::vector<std::string> ListKeys();
		std::vector<std::string> ListDoubleKeys();
		std::vector<std::string> ListStringKeys();

	private:
		std::map<std::string, std::vector<std::string> > key_strings;
		std::map<std::string, std::vector<std::string> >::iterator ims;

		std::map<std::string, std::vector<double> > key_doubles;
		std::map<std::string, std::vector<double> >::iterator imd;

		bool kVerbosity, kStatus;
		std::string _cardFile;
};

#endif
