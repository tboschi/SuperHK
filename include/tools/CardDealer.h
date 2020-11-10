/* CardDealer
 * handles cards
 * Templated version
 */

#ifndef CardDealer_H
#define CardDealer_H

#include <cstdlib>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <type_traits>

#include "tools/ic.h"

//// trim from start (in place)
//static inline void ltrim(std::string &s)
//{
//	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
//				return !std::isspace(ch); }));
//}
//
//// trim from end (in place)
//static inline void rtrim(std::string &s)
//{
//	s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
//				return !std::isspace(ch); }).base(), s.end());
//}
//
//// trim from both ends (in place)
//static inline void trim(std::string &s)
//{
//	ltrim(s);
//	rtrim(s);
//}

class CardDealer
{
	public:
		template <typename T>
		CardDealer(T cardFile, bool verb = false) :
			_cardFile(cardFile),
			kVerbosity(verb) {
			Parse(_cardFile);
		}

		std::string CardName() {
			return _cardFile;
		}

		bool Parse(const std::string &cardFile) {
			_entries.clear();

			std::ifstream inf(cardFile.c_str());

			if(!inf.is_open()) {
				std::string msg = "ERROR CardDealer: " + cardFile + "does not exist";
				throw std::invalid_argument(msg);
			}
			else if (kVerbosity)
				std::cout << "CardDealer: reading input from " << cardFile << std::endl;

			unsigned int nline = 0;
			std::string line;
			while (std::getline(inf, line))
			{
				++nline;
				if (kVerbosity)
					std::cout << "CardDealer: reading line " << nline << ": " << line << std::endl;

				std::string key, token;
				std::vector<std::string> words;

				//removes all comments (anything after #)
				if (line.find_first_of('#') != std::string::npos)
					line.erase(line.find_first_of('#'));

				if (line.empty())
					continue;

				std::stringstream ssl(line);
				ssl >> key;	//key is first word

				if (!ssl && kVerbosity)	//end of line
				{
					std::cout << "CardDealer: key " << key << " has empty argument" << std::endl;
					continue;
				}

				while (!ssl.eof()) {
					// string delimited by "..." or '...'
					if (ssl.peek() == '\'' || ssl.peek() == '\"')
					{
						char end = ssl.peek();
						ssl.ignore();

						std::string word;
						std::getline(ssl, word, end);
						if (word.length())
							words.push_back(word);
					}
					// no delimiters, so expect words separated by commas
					else if (!std::isspace(ssl.peek())) {
						std::string word;
						std::getline(ssl, word, ',');
						//ssl >> word;
						if (word.length())
							words.push_back(word);
					}
					else
						ssl.ignore();
				}
				//make key lower case
				_entries[key] = words;

				if (kVerbosity)
					std::cout << "CardDealer: obtained key " << key << std::endl;
			}
		}

		
		// empty template
		//template <typename T>
		//bool Get(const std::string &key, T &ret) {}

		template <typename T>
		struct is_sequence : std::false_type {};

		template <typename ... Args>
		struct is_sequence <std::vector<Args...> > : std::true_type {};

		template <typename ... Args>
		struct is_sequence <std::set<Args...> > : std::true_type {};

		template <typename ... Args>
		struct is_sequence <std::list<Args...> > : std::true_type {};

		template <typename ... Args>
		struct is_sequence <std::deque<Args...> > : std::true_type {};


		// generic one, uses stringstream to cast string to correct type
		// good for fundamental type
		template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type* = nullptr>
		bool Get(const std::string &key, T &ret) {
			auto id = _entries.find(key);
			if (id != _entries.end()) {
				std::stringstream st(id->second.front());
				return static_cast<bool>(st >> ret);	// successful casting
			}
			else
				return false;
		}

		// good for insert types
		template <typename C, typename std::enable_if<is_sequence<C>::value>::type* = nullptr>
		bool Get(const std::string &key, C &ret) {
			auto id = _entries.find(key);
			if (id != _entries.end()) {
				ret.clear();
				for (const std::string &it : id->second) {
					typename C::value_type ty;
					std::stringstream st(it);
					if (st >> ty)
						ret.insert(ret.end(), ty);
				}
				return bool(ret.size());
			}
			else
				return false;
		}


		// return pieces of maps
		template <typename T>
		bool Get(const std::string &key, std::map<std::string, T> &ret) {
			ret.clear();
			for (const auto &id : _entries) {
				if (id.first.find(key) != std::string::npos) {
					T ty;
 					if (Get(id.first, ty))
						ret[id.first.substr(id.first.find(key)
								  + key.length())] = ty;
				}
			}

			return bool(ret.size());
		}

		// for strings, there is no need to cast

		bool Get(const std::string &key, std::string &ret) {
			auto id = _entries.find(key);
			if (id != _entries.end()) {
				ret = id->second.front();
				return true;
			}
			else
				return false;
		}

		// for vectors of strings, there is no need to cast

		bool Get(const std::string &key, std::vector<std::string> &ret) {
			auto id = _entries.find(key);
			if (id != _entries.end()) {
				ret = id->second;
				return bool(ret.size());
			}
			else
				return false;
		}

		bool Get(const std::string &key) {
			return (_entries.find(key) != _entries.end());
		}


		std::vector<std::string> ListKeys() {
			std::vector<std::string> keys;
			keys.reserve(_entries.size());
			for (const auto &id : _entries)
				keys.push_back(id.first);
			return keys;
		}



	private:
		std::map<std::string, int, ic> st;
		std::map<std::string, std::vector<std::string>, ic> _entries;

		bool kVerbosity;
		std::string _cardFile;
};

#endif
