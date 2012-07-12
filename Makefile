DIR_HPP=./inc/
DIR_CPP=./src/
DIR_OBJ=./tmp/
DIR_BIN=./
DIR_TMP=./tmp/

CXXFLAGS=`root-config --cflags` -I$(DIR_HPP) -O0 -g 
LFLAGS=`root-config --libs` -lboost_regex -g
OBJS=fit1dcould.o fitshanalyticaaabackshdircovcoulpars.o main.o merger.o

# search paths
vpath %.xx $(DIR_HPP)
vpath %.cxx $(DIR_CPP)

main: $(OBJS)
	echo -e "\033[01;33m[make]\033[00;32m Linking all files..."
	echo -e "\033[01;33m[make]\033[01;36m $(addprefix $(DIR_OBJ), $^) \t\033[00;31m$(LFLAGS)\033[00m"
	$(CXX) $(LFLAGS) $(addprefix $(DIR_OBJ), $^) -o $(DIR_BIN)main
	echo -e "\033[01;33m[make]\033[01;36m $(DIR_BIN)main \033[00;32m has been built successfully. \033[00m"

$(OBJS): %.o: %.cxx
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	@[ -d $(DIR_TMP) ] || mkdir -p $(DIR_TMP)
	echo -e "\033[01;33m[make]\033[01;36m $< \t\033[00;31m$(CXXFLAGS)\033[00m"
	$(CXX) $(CXXFLAGS) $< -o $(DIR_OBJ)$@ -c

clean:
	rm -f $(DIR_OBJ)*.o
	rm -f $(DIR_OBJ)*.root
	rm -f ./main
	echo -e "\033[01;33m[make]\033[00;32m All *.o and binary files removed.\033[00m"

.PHONY: all clean
.SILENT :
