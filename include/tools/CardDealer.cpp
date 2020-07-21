#include "CardDealer.h"

CardDealer::CardDealer(std::string cardFile, bool verb) :
	_cardFile(cardFile),
	kVerbosity(verb)
{
	ReadCard(_cardFile);
}

CardDealer::CardDealer(char *filename, bool verb) :
	_cardFile(filename),
	kVerbosity(verb)
{
	ReadCard(_cardFile);
}

std::string CardDealer::CardName() {
	return _cardFile;
}

bool CardDealer::ReadCard(const std::string cardFile)
{
	key_doubles.clear();
	key_strings.clear();

	std::ifstream inFile(cardFile.c_str());
	
	if(!inFile.is_open())
	{
		std::cerr << "CardDealer::ReadCard " << cardFile << " does not exist." << std::endl; 
		return false;
	}
	else if (kVerbosity)
		std::cout << "CardDealer::ReadCard reading input from " << cardFile << std::endl;

   	unsigned int nLine = 0;
	std::string line;
	while (std::getline(inFile, line))
	{
		++nLine;
		if (kVerbosity)
			std::cout << "CardDealer::ReadCard reading line " << nLine << ": " << line << std::endl;

		std::string key, token;
		std::vector<double> numbers;
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
			std::cout << "CardDealer::ReadCard key " << key << " has empty argument" << std::endl;
			continue;
		}

		std::getline(ssl, token);
		trim(token);	//trim white spaces

		ssl.clear();
		ssl.str("");
		ssl << token;

		while (ssl)	//while not end of line
		{
			//check if start of string, should be delimited by "..." or '...'
			if (ssl.peek() == '\'' || ssl.peek() == '\"')
			{
				char end = ssl.peek();
				ssl.ignore();

				std::string word;
				std::getline(ssl, word, end);
				if (ssl.fail()) {
					std::cerr << "ERROR - CardDealer : string variable is not closed with proper delimiter (" << end << ")" << std::endl;

					break;
				}
				words.push_back(word);
			}
			//check if digit
			else if (std::isdigit(ssl.peek()) || ssl.peek() == '-'
				     || ssl.peek() == '+' || ssl.peek() == '.')
			{
				double value;
				ssl >> value;
				if (ssl.fail())
					break;
				numbers.push_back(value);
			}
			else
				ssl.ignore();
		}


		/*
		bool isString = true;
		if (ssl.peek() == '.')		//line starting with . can be (e.g. ./)  or a double (e.g. .0321 )
		{
			ssl.get();
			if (ssl.peek() == '/' || ssl.peek() == '.')	//must be a path
				isString = true;
			else				//must be a double
				isString = false;

			ssl.seekg(-1, std::ios::cur);	//go back one
		}

		if (isString && !std::isdigit(ssl.peek()) && ssl.peek() != '-' && ssl.peek() != '+')	//a string
		{
			std::getline(ssl, token, '\"');
			token.erase(token.find_last_not_of(" \t\n\r\f\v")+1);		//trim empty spaces
		}
		else						//a double (or array)
		{
			double value;
			while (ssl)
			{
				if (!std::isdigit(ssl.peek()) &&
				    ssl.peek() != '.' && ssl.peek() != '-' && ssl.peek() != '+')
					ssl.ignore();
				else if (ssl >> value)
					vArray.push_back(value);
			}
		}
		*/


		if (!key.empty())
		{
			if (kVerbosity)
				std::cout << "CardDealer::ReadCard obtained key " << key << " and ";

			if (words.size())
			{
				if (kVerbosity)
					std::cout << "strings" << std::endl;
				key_strings[key] = words;
			}
			else if (numbers.size())
			{
				if (kVerbosity)
					std::cout << "doubles" << std::endl;
				key_doubles[key] = numbers;
			}
		}
	}
}

bool CardDealer::Find(const std::string key)
{
	return (key_strings.count(key) + key_doubles.count(key)) > 0;
}


////////////// get first value

//get first int
bool CardDealer::Get(const std::string key, int &ii)
{
	double dd;
	if (Get(key, dd))
	{
		ii = int(dd);
		return true;
	}
	else
		return false;
}


//get first double
bool CardDealer::Get(const std::string key, double &dd)
{
	imd = key_doubles.find(key);

	if (imd != key_doubles.end())
	{
		std::vector<double> values = imd->second;
		if (values.size() > 1 && kVerbosity)
			std::cerr << "CardDealer::Get \"" << key << "\" not a single value!" << std::endl;

		dd = values.front();
		return true;
	}
	else
	{
		if( kVerbosity )
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}

//get first string
bool CardDealer::Get(const std::string key, std::string & ss)
{
	ims = key_strings.find(key);

	if (ims != key_strings.end())
	{
		std::vector<std::string> values = ims->second;
		if (values.size() > 1 && kVerbosity)
			std::cerr << "CardDealer::Get \"" << key << "\" not a single value!" << std::endl;
		ss = values.front();
		return true;
	}
	else
	{
		if( kVerbosity )
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		ss = "";
		return false;
	}
}

////////////// get vectors

//get vector of int
bool CardDealer::Get(const std::string key, std::vector<int>& vi)
{
	std::vector<double> vd;
	std::vector<double>::iterator id;

	if (Get(key, vd))
	{
		vi.clear();
		for (id = vd.begin(); id != vd.end(); ++id)
			vi.push_back(int(*id));
		return true;
	}
	else
		return false;
}

//get vector of double
bool CardDealer::Get(const std::string key, std::vector<double>& vd)
{
	imd = key_doubles.find(key);

	if (imd != key_doubles.end())
	{
		vd.insert(vd.end(), imd->second.begin(), imd->second.end());
		return true;
	}
	else
	{
		if( kVerbosity )
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}

//get vector of strings
bool CardDealer::Get(const std::string key, std::vector<std::string>& vs)
{
	ims = key_strings.find(key);

	if (ims != key_strings.end())
	{
		vs.insert(vs.end(), ims->second.begin(), ims->second.end());
		return true;
	}
	else
	{
		if( kVerbosity )
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}


////////////// get map of vectors

//build map of vectors of int
bool CardDealer::Get(const std::string key, std::map<std::string, std::vector<int> >& mi)
{
	mi.clear();

	if(kVerbosity)
		std::cout << "CardDealer::Get looking for substring \"" << key << "\"" << std::endl;

	for (imd = key_doubles.begin(); imd != key_doubles.end(); ++imd)
	{
		if (imd->first.find(key) != std::string::npos)	//found subkey
		{
			std::vector<int> vi;
			Get(imd->first, vi);

			std::string newKey = imd->first;
			newKey.erase(0, key.length());
			mi[newKey] = vi;

			if(kVerbosity)
			{
				std::cout << "CardDealer::Get found \"" << imd->first << "\"" << std::endl;
				std::cout << "            created new key \"" << newKey << "\"" << std::endl;
			}
		}
	}

	if (mi.size())
		return true;
	else
	{
		if(kVerbosity)
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}

//build map of vectors of double
bool CardDealer::Get(const std::string key, std::map<std::string, std::vector<double> >& mv)
{
	mv.clear();

	if(kVerbosity)
		std::cout << "CardDealer::Get looking for substring \"" << key << "\"" << std::endl;

	for (imd = key_doubles.begin(); imd != key_doubles.end(); ++imd)
	{
		if (imd->first.find(key) != std::string::npos)	//found subkey
		{
			std::string newKey = imd->first;
			newKey.erase(0, key.length());
			mv[newKey] = imd->second;

			if(kVerbosity)
			{
				std::cout << "CardDealer::Get found \"" << imd->first << "\"" << std::endl;
				std::cout << "            created new key \"" << newKey << "\"" << std::endl;
			}
		}
	}

	if (mv.size())
		return true;
	else
	{
		if(kVerbosity)
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}

//build map of vectors of strings
bool CardDealer::Get(const std::string key, std::map<std::string, std::vector<std::string> >& ms)
{
	ms.clear();

	if(kVerbosity)
		std::cout << "CardDealer::Get looking for substring \"" << key << "\"" << std::endl;

	for (ims = key_strings.begin(); ims != key_strings.end(); ++ims)
	{
		if (ims->first.find(key) != std::string::npos)	//found subkey
		{
			std::string newKey = ims->first;
			newKey.erase(0, key.length());
			ms[newKey] = ims->second;

			if(kVerbosity)
			{
				std::cout << "CardDealer::Get found \"" << ims->first << "\"" << std::endl;
				std::cout << "            created new key \"" << newKey << "\"" << std::endl;
			}
		}
	}

	if (ms.size())
		return true;
	else
	{
		if(kVerbosity)
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}


//////////// get map of first objects

//build map of first doubles
bool CardDealer::Get(const std::string key, std::map<std::string, int>& mi)
{
	mi.clear();

	if(kVerbosity)
		std::cout << "CardDealer::Get looking for substring \"" << key << "\"" << std::endl;

	for (imd = key_doubles.begin(); imd != key_doubles.end(); ++imd)
	{
		if (imd->first.find(key) != std::string::npos)	//found subkey
		{
			std::string newKey = imd->first;
			newKey.erase(0, key.length());

			int ii;
			Get(imd->first, ii);
			mi[newKey] = ii;

			if(kVerbosity)
			{
				std::cout << "CardDealer::Get found \"" << imd->first << "\"" << std::endl;
				std::cout << "            created new key \"" << newKey << "\"" << std::endl;
			}
		}
	}

	if (mi.size())
		return true;
	else
	{
		if(kVerbosity)
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}

//build map of first doubles
bool CardDealer::Get(const std::string key, std::map<std::string, double>& mv)
{
	mv.clear();

	if(kVerbosity)
		std::cout << "CardDealer::Get looking for substring \"" << key << "\"" << std::endl;

	for (imd = key_doubles.begin(); imd != key_doubles.end(); ++imd)
	{
		if (imd->first.find(key) != std::string::npos)	//found subkey
		{
			std::string newKey = imd->first;
			newKey.erase(0, key.length());
			mv[newKey] = imd->second.front();

			if(kVerbosity)
			{
				std::cout << "CardDealer::Get found \"" << imd->first << "\"" << std::endl;
				std::cout << "            created new key \"" << newKey << "\"" << std::endl;
			}
		}
	}

	if (mv.size())
		return true;
	else
	{
		if(kVerbosity)
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}

//build map of first strings
bool CardDealer::Get(const std::string key, std::map<std::string, std::string>& ms)
{
	ms.clear();

	if(kVerbosity)
		std::cout << "CardDealer::Get looking for substring \"" << key << "\"" << std::endl;

	for (ims = key_strings.begin(); ims != key_strings.end(); ++ims)
	{
		if (ims->first.find(key) != std::string::npos)	//found subkey
		{
			std::string newKey = ims->first;
			newKey.erase(0, key.length());
			ms[newKey] = ims->second.front();

			if(kVerbosity)
			{
				std::cout << "CardDealer::Get found \"" << ims->first << "\"" << std::endl;
				std::cout << "                created new key \"" << newKey << "\"" << std::endl;
			}
		}
	}

	if (ms.size())
		return true;
	else
	{
		if(kVerbosity)
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}



std::vector<std::string> CardDealer::ListStringKeys()
{
	std::vector<std::string> vKeys;

	for (ims = key_strings.begin(); ims != key_strings.end(); ++ims)
		vKeys.push_back(ims->first);

	return vKeys;
}


std::vector<std::string> CardDealer::ListDoubleKeys()
{
	std::vector<std::string> vKeys;

	for (imd = key_doubles.begin(); imd != key_doubles.end(); ++imd)
		vKeys.push_back(imd->first);

	return vKeys;
}

std::vector<std::string> CardDealer::ListKeys()
{
	std::vector<std::string> sKeys = ListStringKeys();
	std::vector<std::string> dKeys = ListDoubleKeys();

	sKeys.insert(sKeys.end(), dKeys.begin(), dKeys.end());

	return sKeys;
}
