#include <iostream>

#include "CardDealer.h"
#include "FileRecord.h"
#include "EventParser.h"
#include "Oscillator.h"

#include "Event.h"

int main(int argc, char** argv)
{

	std::vector<double> vVariables;
	double weight;
	for (unsigned int i = 0; i < numberofevents; ++i)
	{
		GetEntry(i);
		Event::Type et = Parser->Parse(datamanager, vVariables, w);
		if (w < 0)
			DataStacker->Stack(et, vVariables);
		else
			MCStacker->Stack(et, vVariables);
	}

	TTree *parmTree = new TTree ("parameters", "parameters");

	parmTree->Branch("CHI", &chi, "fchi/D");
	parmTree->Branch("dCP", &dCP, "fdCP/D");
	parmTree->Branch("M21", &M21, "fM21/D");
	parmTree->Branch("M31", &M31, "fM31/D");
	parmTree->Branch("dCP", &dCP, "fdCP/D");
	parmTree->Branch("dCP", &dCP, "fdCP/D");

	for (unsigned int i = 0; i < oscparspace; ++i)
	{
		MCStacker->ChangeStacks(Oscillator.at(i));
		dCP = Oscillator.at(i)->GetCP();
		chi = ChiTest(DataStacker, MCStacker);
	}

	return 0;
}
