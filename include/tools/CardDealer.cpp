#include "CardDealer.h"

CardDealer::CardDealer(std::string cardFile) :
	kVerbosity(Utils::verbosity)
{
	ReadCard(cardFile);
}

bool CardDealer::Status()
{
	return kStatus;
}

bool CardDealer::ReadCard(const std::string cardFile)
{
	mKeyDouble.clear();
	mKeyString.clear();

	std::ifstream inFile(cardFile.c_str());
	
	if(!inFile.is_open())
	{
		std::cerr << "CardDealer::ReadCard " << cardFile << " does not exist." << std::endl; 
		kStatus = false;
	}
	else
		std::cout << "CardDealer::ReadCard reading input from " << cardFile << std::endl;

   	unsigned int nLine = 0;
	std::string line;
	while (std::getline(inFile, line))
	{
		if (kVerbosity)
			std::cout << "CardDealer::ReadCard reading line " << nLine << ": " << line << std::endl;
		std::string key, token;
		std::vector<double> vArray;
		++nLine;

		//removes all comments (anything after #)
		if (line.find_first_of('#') != std::string::npos)
			line.erase(line.find_first_of('#'));

		if (line.empty())
			continue;

		std::stringstream ssl(line);
		
		if (std::isalpha(ssl.peek()))	//first character of key must be a letter
			ssl >> key;
		else
		{
			std::cout << char(ssl.peek()) << "\t" << ssl.str() << std::endl;
			std::cerr << "CardDealer::ReadCard line " << nLine << ": irregular first character. Skipping." << std::endl;
			kStatus = false;
			continue;
		}


		while (std::isspace(ssl.peek()))	//or isblank? move to first non-blank character
			ssl.ignore();

		if (ssl.peek() == '{' || ssl.peek() == '\"')	//accepted characters for strings and arrays
			ssl.ignore();

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

		if (kVerbosity)
			std::cout << "CardDealer::ReadCard obtained key \"" << key << "\"" << std::endl;

		if (!key.empty())
		{
			if (!token.empty())
				mKeyString[key] = token;
			else if (vArray.size())
				mKeyDouble[key] = vArray;
		}
	}

	kStatus = true;
}

bool CardDealer::Find(const std::string key)
{
	return (mKeyString.count(key) + mKeyDouble.count(key)) > 0;
}

bool CardDealer::Get(const std::string key, std::vector<double>& vv)
{
	imd = mKeyDouble.find(key);

	if (imd != mKeyDouble.end())
	{
		vv = imd->second;
		return true;
	}
	else
	{
		if( kVerbosity )
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}

bool CardDealer::Get(const std::string key, double &dd)
{
	imd = mKeyDouble.find(key);

	if (imd != mKeyDouble.end())
	{
		std::vector<double> vArray = imd->second;
		if (vArray.size() > 1 && kVerbosity)
			std::cerr << "CardDealer::Get \"" << key << "\" not a single value!" << std::endl;

		dd = vArray.front();
		return true;
	}
	else
	{
		if( kVerbosity )
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}

bool CardDealer::Get(const std::string key, std::string & ss)
{
	ims = mKeyString.find(key);

	if (ims != mKeyString.end())
	{
		ss = ims->second;
		return true;
	}
	else
	{
		if( kVerbosity )
			std::cerr << "CardDealer::Get \"" << key << "\" not found" << std::endl;
		return false;
	}
}

//build list of string
bool CardDealer::Get(const std::string subKey, std::map<std::string, std::vector<double> >& mv)
{
	mv.clear();

	if(kVerbosity)
		std::cout << "CardDealer::Get looking for substring \"" << subKey << "\"" << std::endl;

	for (imd = mKeyDouble.begin(); imd != mKeyDouble.end(); ++imd)
	{
		if (imd->first.find(subKey) != std::string::npos)	//found subkey
		{
			std::string newKey = imd->first;
			newKey.erase(0, subKey.length());
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
			std::cerr << "CardDealer::Get \"" << subKey << "\" not found" << std::endl;
		return false;
	}
}

bool CardDealer::Get(const std::string subKey, std::map<std::string, double>& mv)
{
	mv.clear();

	if(kVerbosity)
		std::cout << "CardDealer::Get looking for substring \"" << subKey << "\"" << std::endl;

	for (imd = mKeyDouble.begin(); imd != mKeyDouble.end(); ++imd)
	{
		if (imd->first.find(subKey) != std::string::npos)	//found subkey
		{
			std::string newKey = imd->first;
			newKey.erase(0, subKey.length());
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
			std::cerr << "CardDealer::Get \"" << subKey << "\" not found" << std::endl;
		return false;
	}
}

bool CardDealer::Get(const std::string subKey, std::map<std::string, std::string>& ms)
{
	ms.clear();

	if(kVerbosity)
		std::cout << "CardDealer::Get looking for substring \"" << subKey << "\"" << std::endl;

	for (ims = mKeyString.begin(); ims != mKeyString.end(); ++ims)
	{
		if (ims->first.find(subKey) != std::string::npos)	//found subkey
		{
			std::string newKey = ims->first;
			newKey.erase(0, subKey.length());
			ms[newKey] = ims->second;

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
			std::cerr << "CardDealer::Get \"" << subKey << "\" not found" << std::endl;
		return false;
	}
}

std::vector<std::string> CardDealer::ListKeys()
{
	std::vector<std::string> vKeys;

	for (ims = mKeyString.begin(); ims != mKeyString.end(); ++ims)
		vKeys.push_back(ims->first);

	for (imd = mKeyDouble.begin(); imd != mKeyDouble.end(); ++imd)
		vKeys.push_back(imd->first);

	return vKeys;
}


/*

void CardReader::BuildBinMap()
{

   std::string key;
   std::string binname;
   std::string::size_type loc=0;  

   double * edges;
   int nbins;

   int counter = 0;
   for( iSMap  = StringKeys.begin() ; iSMap != StringKeys.end() ; iSMap++ )
   {
       key = iSMap->first;
       loc = key.find("bin_", 0 , 4);       

       if( loc == std::string::npos ) continue;

       // expect the prefix to be at the beginning(duh)
       binname = key.substr( loc + 4 );

       //std::cout << "Key: " << key << " " << binname << std::endl;
       nbins = ParseArray( iSMap->second , "," , &edges ); 
       //std::cout <<  binname << " " << nbins<< " " << " " << edges[0] << " " << edges[1] << " " << edges[nbins] << std::endl;

       bins   [ binname       ] = &edges[0];        
       bintype[ int(edges[1]) ] = binname;        

       counter++;
   }


}


void CardReader::BuildDimensionMap()
{


   std::string key;
   std::string binname;
   std::string::size_type loc;  

   std::vector < std::string > tokens;
   std::vector < std::string > sub_tokens;
   int ntokens;
   int nsubtokens;

   double * edges;
   int nbins;

   int counter = 0;
   for( iSMap  = StringKeys.begin() ; iSMap != StringKeys.end() ; iSMap++ )
   {
       key = iSMap->first;
       loc = key.find("dimension_");       

       if( loc == std::string::npos ) continue;

       // expect the prefix to be at the beginning(duh)
       binname = key.substr( loc + 10 );

       // split into bin_type_1 (cut_type_1): bin_type_2 (cut_type_2)
       ntokens = Tokenize( iSMap->second , ":" , tokens );

       //simple tokenize
       for( int i = 0 ; i < ntokens ; i++ )
       {
          nsubtokens = Tokenize( tokens[i],"()", sub_tokens );    
          

       }

       counter++;
   }


}

void CardReader::Erase( std::string & source, const std::string& kill, std::string::size_type l0 )
{

   std::string::size_type loc0 = source.find( kill, l0 );
  
   while( loc0 != std::string::npos )
   {
      source.erase( loc0, 1);
      loc0 = source.find(kill, loc0);
   }


}



unsigned int CardReader::Tokenize(const std::string& source,
                                        const std::string& delimiters,
                                        std::vector<std::string>& tokens)
{
    std::string::size_type prev_loc0 = 0;
    std::string::size_type loc0 = 0;
    std::string::size_type loc1 = 0;
    std::string sub;
    unsigned int ntokens = 0;

    tokens.clear();

    loc0 = source.find_first_of(delimiters, loc0);
    while (loc0 != std::string::npos)
    {

        sub = source.substr(prev_loc0, loc0 - prev_loc0);

        tokens.push_back( sub ) ;
        ntokens++;


        loc0++;
        prev_loc0 = loc0;
        loc0 = source.find_first_of(delimiters, loc0);
    }


    if (prev_loc0 < source.length())
    {
        tokens.push_back(source.substr(prev_loc0));
        ntokens++;
    }

    return ntokens;
}




double * CardReader::GetBinsForType( int a )
{

   return bins[ bintype[a] ];


}


void CardReader::MakeStringTokens( TString * line, const char * D , std::vector< TString *> & data )
{
    unsigned i;
    for( i = 0 ; i < data.size() ; i++ )
       if( data[i] )
         delete data[i]; 
    data.clear();


    TObjArray * tokens = line->Tokenize(D);
    unsigned size = tokens->GetEntries();

    TObjString * s;
    for( i = 0 ; i < size ; i++ )
    {
       s = (TObjString*) tokens->At(i);
       if( s->GetString().Length() > 0 )
          data.push_back( new TString( s->GetString().ReplaceAll(" ",1) ) );
    }

}


void CardReader::Replace( std::string & source, const std::string & kill , 
                                                const std::string & rep  , std::string::size_type l0 )
{

   std::string::size_type loc0 = source.find( kill, l0 );
  
   while( loc0 != std::string::npos )
   {
      source.replace( loc0, kill.length(), rep );
      loc0 = source.find(kill, loc0);
   }

}


TokenMap * CardReader::BuildTokenMap( const char * skey, const char * delim  )
{
   
   std::map< std::string, std::string > aMap;
   std::map< std::string, std::string >::iterator _i;
   std::vector < std::string > tokens;      
   int ntokens;

   TokenMap * tokenMap = new TokenMap( skey );

   BuildListOfStrings( skey, aMap ); 
  
   for( _i = aMap.begin() ; _i != aMap.end() ; _i++ )
   {
      tokens.clear();
      ntokens = Tokenize( _i->second , delim , tokens );   
     
      tokenMap->AddTokenList( _i->first.c_str() , ntokens, tokens ); 
   } 


  return tokenMap;
}

*/
