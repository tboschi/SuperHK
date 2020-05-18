INCDIR =	include
ProbINC =	../T2HKK/Prob3++
APPDIR =	app
BINDIR =	bin
LIBDIR =	lib

ROOTLIB		= $(shell root-config --glibs)
ROOTCXX		= $(shell root-config --cflags)

OSC3LIB = /data/tboschi/HKsens/OscAna/Osc3++/lib
OSC3INC = /data/tboschi/HKsens/OscAna/Osc3++/src

LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIB) -L$(LIBDIR) -L$(OSC3LIB)
LDLIBS   := -losc3pp
CXXFLAGS := $(CXXFLAGS) -fPIC -std=c++11 -O3 -march=native $(ROOTCXX) -I$(INCDIR) -I$(ProbINC) -I$(OSC3INC)


#apps and exctuables
CPP =	get_hessian	\
	escale 	\
	escale_other 	\
	escale_deriv
	#TestCard	\
	vis		\
	escale		\
	#fitter		\
	testfitter	\
	newfitter	\
	addpenalty	\
	buildcontours	\
	dropchi2	\
	dropsens	\
	exclusion	\
	getepsilon
	
	#atmofitter	\
	convertSKsyst	\
	predictions	\
	addatmo		\
	getmin
	#baseball	\
	plotosc		\
	#momspaghetti	\
	testchi2	\
	convertpoint	\
	test
	#nevents	\
	validation	\
	addmatrix	\
	exclusion_filter	\
	purifysystematic	\
	oscillator	\
	BuildContourPlots \
	runchi2		\
	compare		\
	filter 		\
	viewchi2	\
	baseball	\

	#pmns	\
	penalty	\
	simple	\
	#PredictEvents	\
	ValorEvents	\
	GetNevents	\
	GetSpectra	\
	TestOsc

#header folders
HPP =	tools/Const		\
	tools/Integration	\
	tools/CardDealer	\
	event/Flux	\
	event/Reco	\
	event/ChiSquared	\
	physics/Oscillator	\
	physics/ParameterSpace	\
	event/NewChiSquared	\
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

$(TGT): $(DEP)
#$(INccc)

clean:
	find $(INCDIR) -mindepth 1 -name "*.o" -delete
	find $(INCDIR) -mindepth 1 -name "*~"  -delete
	find $(APPDIR) -mindepth 1 -name "*~"  -delete
	find $(BINDIR) -mindepth 1 -name "*"   -delete

.PHONY: all clean

