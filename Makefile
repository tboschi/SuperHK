INCDIR = include
APPDIR = app
BINDIR = bin
DOCDIR = doc

## root
ROOTLIB	= $(shell root-config --glibs)
ROOTCXX	= $(shell root-config --cflags)

## Eigen matrix library
#TARGETS := $(if $(APP), $(APPDIR)/$(APP), $(CPP:.cpp=))
EIGENINC = $(if $(EIGEN), -I$(EIGEN), )

#optimization
ARCH ?= -march=native

WARNING := -Wall
LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIB)
CXXFLAGS := $(DEBUG) $(WARNING) -fPIC -std=c++11 -O3 $(ARCH) $(ROOTCXX) -I$(INCDIR) $(EIGENINC)


#apps and exctuables
TARGETS := $(shell find $(APPDIR) -maxdepth 1 -name '*.cpp')
SOURCES := $(shell find $(INCDIR) -maxdepth 2 -name '*.cpp')
HEADERS := $(shell find $(INCDIR) -maxdepth 2 -name '*.h')

OBJECTS := $(SOURCES:.cpp=.o)
DEPENDS := $(SOURCES:.cpp=.d)
TARGETS := $(if $(APP), $(APPDIR)/$(APP), $(TARGETS:.cpp=))


##main target
#OBJECTS := $(SOURCES:.cpp=.o)
#DEPEND := $(SRC:.cpp=.o)
#DEPEND := $(DEPEND:.cc=.o)
#DEPEND := $(DEPEND:.c=.o)

all: $(TARGETS)
	@mkdir -p $(BINDIR)
	@cp $^ $(BINDIR)
	@echo "Done!"


doc: 
	$(MAKE) -C $(DOCDIR)

help:
	@echo Targets found are $(TARGETS)
	@echo Sources found are $(SOURCES)
	@echo Headers found are $(HEADERS)
	@echo "If you need to build just one file, do make APP=name"
	@echo "or if you need to specify an architecture, do make ARCH=arch"
	@echo "To build documentation, make doc"
	@echo "Enjoy your compilation"

$(TARGETS): $(OBJECTS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPENDS)

#$(OBJECT): $(HEADERS)


clean:
	$(RM) $(TARGETS)
	$(RM) $(OBJECTS)
	$(RM) $(DEPENDS)
	$(MAKE) -C $(DOCDIR) clean

#$(RM) -r $(BINDIR)



.PHONY: all doc help clean
