#ifndef EventType_H
#define EventType_H

namespace Event
{
	enum Type
	{
		undefined,

		SubGeV_elike_0dcy,
		SubGeV_elike_1dcy,
		SubGeV_SingleRing_pi0like,
		SubGeV_mulike_0dcy,
		SubGeV_mulike_1dcy,
		SubGeV_mulike_2dcy,
		SubGeV_pi0like,
		MultiGeV_elike_nue,
		MultiGeV_elike_nuebar,
		MultiGeV_mulike,
		MultiRing_elike_nue,
		MultiRing_elike_nuebar,
		MultiRing_mulike,
		MultiRingOther_1,
		PCStop,
		PCThru,
		UpStop_mu,
		UpThruNonShower_mu,
		UpThruShower_mu,
		
		T2Knumu,
		T2Knumu1,
		T2Knumu2,
		T2Knumu3,
		T2Knumu4,
		T2Knumu5,
		T2Knue,
		T2Knue1,
		T2Knue2,
		T2Knue3,
		T2Knue4,
		T2Knue5,
		T2Knumubar,
		T2Knumubar1,
		T2Knumubar2,
		T2Knumubar3,
		T2Knumubar4,
		T2Knumubar5,
		T2Knuebar,
		T2Knuebar1,
		T2Knuebar2,
		T2Knuebar3,
		T2Knuebar4,
		T2Knuebar5,
		
		MINOSnumu,
		MINOSnue
	};
};

#endif
