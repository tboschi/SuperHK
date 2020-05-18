#include <string>
#include <fstream>

void loop(std::string pt = "./")
{
	std::string cmd = "ls " + pt + " > .tmp_list";
	system(cmd.c_str());

	std::ifstream in(".tmp_list");
	std::string line;
	while (std::getline(in, line))
		extract(line);

	in.close();
}

void extract(std::string name)
{
	TFile *inf = new TFile(name.c_str(), "OPEN");
	if (inf->IsZombie())
	{
		std::cout << "file doesn't exist\n";
		return;
	}

	name = name.replace(name.find(".root"), 5, ".dat");
	std::ofstream out(name.c_str());

	TH1D* h = static_cast<TH1D*>(inf->Get("X2minCP"));
	for (int i = 1; i < h->GetNbinsX()+1; ++i)
		out << h->GetBinCenter(i) << "\t" << h->GetBinContent(i) << "\n";

	inf->Close();
	out.close();
	delete inf;
}
