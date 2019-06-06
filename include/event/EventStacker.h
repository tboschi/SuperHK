/* EventStacker
 * creates stacks (histograms) and fill them
 */

#ifndef EventStacker_H
#define EventStacker_H

#include <map>
#include <vector>
#include <string>

#include "tools/CardDealer.h"
#include "event/Event.h"
#include "physics/Oscillator.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

class Oscillator;

class EventStacker
{
	public:
		EventStacker(CardDealer *& binCard);
		~EventStacker();

		void CreateStacks(CardDealer * &binCard);
		void Stack(Event::Type et, std::vector<double> &vars, double w);

		void ChangeStacks(Oscillator *&l);
		void ResetStacks();

		unsigned int Bins();
		double GetBin(unsigned int b);

		std::vector<double> Vectorise();



	private:
		void InitMap();
		void CreateStack(std::string name, Event::Type et, std::map<char, std::vector<double> > &mAx);
		unsigned int FillStack(Event::Type et, std::vector<double> &vars, double w);

		void ChangeWeights(Event::Type et, Oscillator *&l);
		void ResetWeights(Event::Type et);

		std::map<Event::Type, std::string> mEvent;
		std::map<Event::Type, std::string>::iterator ime;

		std::map<Event::Type, TH1*> mStack;
		std::map<Event::Type, TH1*>::iterator ims;

		std::map<Event::Type, std::vector<Eigen::Triplet<double> > > trWeight;
		std::map<Event::Type, std::vector<Eigen::Triplet<double> > > trEnergy;
		std::map<Event::Type, unsigned int> mColumnIdx;

		bool kVerbosity;
};

#endif
