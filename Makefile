INCDIR =	include
APPDIR =	app
BINDIR =	bin
LIBDIR =	lib

## root
ROOTLIB		= $(shell root-config --glibs)
ROOTCXX		= $(shell root-config --cflags)

## Eigen matrix library
EIGENINC = $(EIGEN)

#optimization
ARCH ?= -march=native

LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIB) -L$(LIBDIR)
#LDLIBS   := -losc3pp
CXXFLAGS := $(CXXFLAGS) $(DEBUG) -fPIC -std=c++11 -O3 $(ARCH) $(ROOTCXX) -I$(INCDIR) -I$(EIGENINC)




#apps and exctuables
CPP := $(shell find $(APPDIR) -maxdepth 1 -name '*.cpp')
SRC := $(shell find $(INCDIR) -maxdepth 2 -path $(INCDIR)/ignore -prune -o -name '*.c*' -not -name '.*' -print)


#main target
TARGET := $(if $(APP), $(APPDIR)/$(APP), $(CPP:.cpp=))
#TARGET := $(CPP:.cpp=)
DEPEND := $(SRC:.cpp=.o)
DEPEND := $(DEPEND:.cc=.o)
DEPEND := $(DEPEND:.c=.o)

all: welcome $(TARGET)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(BINDIR)
	@echo "Cleaning up..."
	@cp $(DEPEND) $(LIBDIR)
	@cp $(TARGET) $(BINDIR)
	@echo "Done!"

welcome:
	@echo "If you need to build just one file, do make APP=name"
	@echo "or if you need to specify an architecture, do make ARCH=arch"
	@echo "Enjoy your compilation"


$(TARGET): $(DEPEND)

include:
	$(eval DEPEND := $(shell find $(LIBDIR) -maxdepth 1 -name '*.o'))


clean:
	-find $(INCDIR) -name "*.o"  -delete
	-find $(LIBDIR) -maxdepth 1 -type f -name "*"   -delete
	-find $(APPDIR) -maxdepth 1 -type f -name "*~"  -delete
	-find $(BINDIR) -maxdepth 1 -type f -name "*"   -delete


.PHONY: all clean
