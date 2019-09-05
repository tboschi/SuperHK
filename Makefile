.PHONY: clean

INCDIR =	include
ProbINC =	../T2HKK/Prob3++
APPDIR =	app
BINDIR =	bin
LIBDIR =	lib

ROOTLIB		= $(shell root-config --glibs)
ROOTCXX		= $(shell root-config --cflags)

LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIB) -L$(LIBDIR)
CXXFLAGS := $(CXXFLAGS) -fPIC -std=c++11 -O3 -march=native -mavx $(ROOTCXX) -I$(INCDIR) -I$(ProbINC)


#apps and exctuables
CPP =	purifysystematic	\
	addmatrix	\
	addpenalty	\
	oscillator	\
	#pmns	\
	penalty	\
	simple	\
	#exclusion
	#PredictEvents	\
	ValorEvents	\
	GetNevents	\
	GetSpectra	\
	TestOsc
	#TestCard	\

#header folders
HPP =	tools/Const		\
	tools/Integration	\
	tools/CardDealer	\
	tools/DataManager	\
	event/Flux	\
	event/Reco	\
	physics/Oscillator	
	#event/Event		\
	event/EventStacker	\
	event/EventParser	
HH =	EarthDensity		\
	mosc3		\
	mosc		\
	BargerPropagator	\

#main target
TGT :=	$(CPP:%=$(APPDIR)/%)

#dependencies of target
INCDEP := $(HPP:%=$(INCDIR)/%.cpp)
INccc := $(HH:%=$(ProbINC)/%.o)
#INCccc := $(HH:%=$(INCDIR)/%.cpp)
DEP :=	$(patsubst %.cpp,%.o,$(wildcard $(INCDEP)))
#Dcc :=	$(patsubst %,%.o,$(wildcard INC$(HH)))

all: $(TGT)
	@mkdir -p $(BINDIR)
	@mkdir -p $(LIBDIR)
	@echo "Moving stuff..."
	@cp $(TGT) $(BINDIR)
	@cp $(DEP) $(LIBDIR)
	@echo "Done!"

$(TGT): $(DEP) $(INccc)

clean:
	find $(INCDIR) -mindepth 1 -name "*.o" -delete
	find $(INCDIR) -mindepth 1 -name "*~"  -delete
	find $(APPDIR) -mindepth 1 -name "*~"  -delete
	find $(BINDIR) -mindepth 1 -name "*"   -delete
