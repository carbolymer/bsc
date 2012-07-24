DIR_HPP=./inc/
DIR_CPP=./src/
DIR_OBJ=./tmp/
DIR_BIN=./
DIR_TMP=./tmp/

CXXFLAGS=`root-config --cflags` -I$(DIR_HPP) -O3 -g 
LFLAGS=`root-config --libs` -lboost_regex -g
OBJS=fit1dcould.o fitshanalyticaaabackshdircovcoulpars.o merger.o siniukow2therminator.o

# search paths
vpath %.xx $(DIR_HPP)
vpath %.cxx $(DIR_CPP)

all: fitsh merger fit1d plotter

fitsh: fitshanalyticaaabackshdircovcoulpars.o
	echo -e "\033[01;33m[make]\033[00;32m Generating fitsh..."
	echo -e "\033[01;33m[make]\033[01;36m $(addprefix $(DIR_OBJ), $^) \t\033[00;31m$(LFLAGS)\033[00m"
	$(CXX) $(LFLAGS) $(addprefix $(DIR_OBJ), $^) -o $(DIR_BIN)fitsh	

merger: merger.o
	echo -e "\033[01;33m[make]\033[00;32m Generating merger..."
	echo -e "\033[01;33m[make]\033[01;36m $(addprefix $(DIR_OBJ), $^) \t\033[00;31m$(LFLAGS)\033[00m"
	$(CXX) $(LFLAGS) $(addprefix $(DIR_OBJ), $^) -o $(DIR_BIN)merger

fit1d: fit1dcould.o
	echo -e "\033[01;33m[make]\033[00;32m Generating fit1d..."	
	echo -e "\033[01;33m[make]\033[01;36m $(addprefix $(DIR_OBJ), $^) \t\033[00;31m$(LFLAGS)\033[00m"
	$(CXX) $(LFLAGS) $(addprefix $(DIR_OBJ), $^) -o $(DIR_BIN)fit1d

plots:
	echo -e "\033[01;33m[make]\033[00;32m Making plots..."
	root -l -b -q macros/plotter.C
	root -l -b -q macros/plotter_lcms.C

s2t: siniukow2therminator.o
	echo -e "\033[01;33m[make]\033[00;32m Generating s2t..."	
	echo -e "\033[01;33m[make]\033[01;36m $(addprefix $(DIR_OBJ), $^) \t\033[00;31m$(LFLAGS)\033[00m"
	$(CXX) $(LFLAGS) $(addprefix $(DIR_OBJ), $^) /home/mgalazyn/workspace/Therminator2/build/obj/ParticleCoor.o -o $(DIR_BIN)s2t

$(OBJS): %.o: %.cxx
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	@[ -d $(DIR_TMP) ] || mkdir -p $(DIR_TMP)
	echo -e "\033[01;33m[make]\033[01;36m $< \t\033[00;31m$(CXXFLAGS)\033[00m"
	$(CXX) $(CXXFLAGS) $< -o $(DIR_OBJ)$@ -c

clean:
	rm -f $(DIR_OBJ)*.o
	rm -f $(DIR_OBJ)*.root
	rm -f ./fitsh
	rm -f ./fit1d
	rm -f ./merger
	rm -f ./plotter
	rm -f ./s2t
	echo -e "\033[01;33m[make]\033[00;32m All *.o and binary files removed.\033[00m"

.PHONY: all clean
.SILENT :
