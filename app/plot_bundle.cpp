#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "physics/Oscillator.h"

#include "TF1.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TMatrixD.h"

void beginDocument(std::ofstream &out);
void endDocument(std::ofstream &out);
void addFrame(std::ofstream &out, const std::string &file, const std::string &title);
int main(int argc, char** argv)
{
	// input file should look like this
	// plot/title/title_type.tex
	//
	
	if (argc < 3) {
		std::cerr << "Plot bundle: need at least two parameters outfile and infiles\n";
		return 1;
	}

	// base = plot/title
	std::string base(argv[2]);
	base.erase(base.find_last_of('/'));
	chdir(base.c_str());	// move to plot directory

	std::map<std::string, std::string> typeMap;
	typeMap["chi2_dCP"] = "$\\chi^2$ vs $\\delta_\\text{CP}$";
	typeMap["chi2_M23"] = "$\\chi^2$ vs $\\Delta m^2_{32}$";
	typeMap["chi2_S23"] = "$\\chi^2$ vs $\\sin^2 \\theta_{23}$";
	typeMap["chi2_S13"] = "$\\chi^2$ vs $\\sin^2 2\\theta_{13}$";
	typeMap["cont_dCP_M23"] = "$\\delta_\\text{CP}$ vs $\\Delta^2 m_{32}$";
	typeMap["cont_S13_dCP"] = "$\\sin^2 2\\theta_{13}$ vs $\\delta_\\text{CP}$";
	typeMap["cont_S23_dCP"] = "$\\sin^2 \\theta_{23}$ vs $\\delta_\\text{CP}$";
	typeMap["cont_S23_M23"] = "$\\sin^2 \\theta_{23}$ vs $\\Delta^2 m_{32}$";
	typeMap["cont_S13_M23"] = "$\\sin^2 2\\theta_{13}$ vs $\\Delta^2 m_{32}$";
	typeMap["cont_S13_S23"] = "$\\sin^2 2\\theta_{13}$ vs $\\sin^2 \\theta_{23}$ vs $\\Delta^2 m_{32}$";
	typeMap["sensitivity"] = "sensitivity to $\\delta_\\text{CP}$";
	typeMap["detail"] = "sensitivity to $\\delta_\\text{CP}$ (differece)";
	// to be added
	//typeMap["octant"] = "sensitivity to $\\sin^2 \\theta_{23}$";
	//typeMap["octant"] = "sensitivity to $\\sin^2 \\theta_{23}$ (differece)";


	std::string out(argv[1]);
	if (out.find(".tex") != std::string::npos)
		out += ".tex";
	std::ofstream tout(out.c_str());
	beginDocument(tout);

	// input file should look like this
	// plot/title/title_type.tex
	for (int f = 1; f < argc; ++f) {
		std::string file(argv[f]), title = file;
		std::cout << "processing file " << file << std::endl;

		// title_type.tex
		file.erase(0, file.find_last_of('/')+1);
		std::cout << "X " << file << std::endl;

		// title
		title.erase(title.find_last_of('/'));		// plot/title
		std::cout << "X " << title << std::endl;
		title.erase(0, title.find_first_of('/')+1); 	// title
		std::cout << "X " << title << std::endl;


		// type
		std::string type = file;		// title_type.tex
		std::cout << "X " << type << std::endl;
		type.erase(type.find(".tex"));		// title_type
		std::cout << "X " << type << std::endl;
		type.erase(type.find(title), title.length()+1); // type
		std::cout << "X " << type << std::endl;

		std::cout << "grab " << file << " " << title << " " << type << std::endl;
		while (title.find_first_of('_') != std::string::npos)
			title.replace(title.find_first_of('_'), 1, " ");
		title += " --- " + typeMap[type];
		addFrame(tout, file, title);
	}

	endDocument(tout);
	tout.close();

	std::string cmd = "pdflatex " + out;
	system(cmd.c_str());

	return 0;
}

void beginDocument(std::ofstream &out)
{
	out << "\\documentclass[aspectratio=169]{beamer}\n";
	out << "\\setbeamertemplate{navigation symbols}{}\n";
	out << "\\setbeamerfont{frametitle}{size=\\small}\n";
	//out << "\\setbeamertemplate{frametitle} {\n"
	//    << "\\nointerlineskip\n"
	//    << "\\begin{beamercolorbox}[wd=\\paperwidth,ht=2.0ex,dp=0.6ex]{frametitle}\n"
	//    << "\\hspace*{1ex}\\insertframetitle\n\\end{beamercolorbox} }\n";
	//out << "\\setbeamertemplate{footline}{}\n";
	out <<"\n\\begin{document}\n\n";
}

void addFrame(std::ofstream &out, const std::string &file, const std::string &title)
{
	out << "\\begin{frame}\n"
	    << "\t\\frametitle{" << title << "}\n"
	    << "\t\\small\n"
	    << "\t\\begin{center}\n";

	out << "\t\t\\resizebox{0.70\\linewidth}{!}{\\input{" << file << "}}\n";

	out << "\t\\end{center}\n"
	    << "\\end{frame}\n\n";
}

void endDocument(std::ofstream &out)
{
	out <<"\n\\end{document}\n\n";
}
