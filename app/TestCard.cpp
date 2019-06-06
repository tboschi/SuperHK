#include <iostream>
#include <getopt.h>

#include "tools/CardDealer.h"

void Usage(char* name);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"output", 	required_argument, 	0, 'c'},
		{"verbose", 	required_argument, 	0, 'k'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string cardFile;
	std::ofstream outFile;
	bool verb = false;
	
	while((iarg = getopt_long(argc, argv, "o:kh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'o':
				outFile.open(optarg);
				break;
			case 'k':
				verb = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	cardFile.assign(argv[optind]);

	std::ostream &Out = (outFile.is_open()) ? outFile : std::cout;

	CardDealer *cd = new CardDealer(cardFile);
	std::vector<std::string> vkeys = cd->ListKeys();

	std::vector<double> vv;
	double dd;
	std::string sv;
	for (unsigned int i = 0; i < vkeys.size(); ++i)
	{
		std::cout << vkeys.at(i) << "\t";
		if (cd->Get(vkeys.at(i), sv))
			std::cout << "\"" << sv << "\"" << std::endl;
		else if (cd->Get(vkeys.at(i), vv))
		{
			for (unsigned int j = 0; j < vv.size(); ++j)
				std::cout << vv.at(j) << "\t";
			std::cout << std::endl;
		}
	}

	std::map<std::string, std::string> mSS;
	std::map<std::string, std::string>::iterator im;

	if (cd->Get("chiave_", mSS))
	{
		for (im = mSS.begin(); im != mSS.end(); ++im)
			std::cout << im->first << "\t" << im->second << std::endl;
	}

	return 0;
}

void Usage(char *name)
{
	std::cout << "Usage : " << std::endl;
	std::cout << "\t" << name << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -1,  --option1" << std::endl;
	std::cout << "\t\tDescription" << std::endl;
	std::cout <<"\n  -2,  --option2" << std::endl;
	std::cout << "\t\tDescription" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
