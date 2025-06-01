CXX     = g++

CXXFLAGS= $(shell root-config --cflags)
LIBS    = $(shell root-config --libs) 


SOURCES = HGCNtupleVariables_D1.cc AnalyzeHGCMuons_D1.cc 
HEADERS = HGCNtupleVariables_D1.h AnalyzeHGCMuons_D1.h 
OBJECTS = $(SOURCES:.cc=.o)

EXECUTABLE = analyzeHGCMuons_D1

all: $(SOURCES) $(EXECUTABLE)

%.o: %.cc
	@echo Compiling $<...
	$(CXX) $(CXXFLAGS) -c -o $@ $< 


$(EXECUTABLE): $(OBJECTS)
	@echo "Linking $(PROGRAM) ..."
	@echo "@$(CXX) $(LIBS) $(OBJECTS) -o $@"
	@$(CXX) -o $@ $^ $(LIBS)
	@echo "done"


# Specifying the object files as intermediates deletes them automatically after the build process.
.INTERMEDIATE:  $(OBJECTS)

# The default target, which gives instructions, can be called regardless of whether or not files need to be updated.
.PHONY : clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

###
HGCNtupleVariables_D1.o: HGCNtupleVariables_D1.h
AnalyzeHGCMuons_D1.o:HGCNtupleVariables_D1.h AnalyzeHGCMuons_D1.h
