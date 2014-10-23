COMPILER = g++

CFLAGS = -g -fPIC -Wall -O3 `root-config --cflags` -Iinclude
LDLIBS = `root-config --libs`
LDFLAGS = `root-config --glibs`
ROOT_INC = `root-config --incdir`

SOURCES = vikar_core.cpp detectors.cpp materials.cpp vikar.cpp
TOOLS = vikarFront angleConvert kinematics kindist dump
OBJECTS = $(addprefix $(C_OBJ_DIR)/,$(SOURCES:.cpp=.o))

TOP_LEVEL = $(shell pwd)
INCLUDE_DIR = $(shell pwd)/include
TOOL_DIR = $(shell pwd)/tools
SOURCE_DIR = $(shell pwd)/source
C_OBJ_DIR = $(shell pwd)/obj
DICT_DIR = $(shell pwd)/dict
DICT_OBJ_DIR = $(DICT_DIR)/obj

# ROOT dictionary stuff
DICT_SOURCE = RootDict
STRUCT_FILE = structures

ROOTOBJ = $(DICT_OBJ_DIR)/$(DICT_SOURCE).o 
ROOTOBJ += $(C_OBJ_DIR)/$(STRUCT_FILE).o
SFLAGS = $(addprefix -l,$(DICT_SOURCE))

PROG = vikar

all: $(PROG)
	@echo " Finished Compilation"

dictionary: $(DICT_OBJ_DIR) $(DICT_OBJ_DIR)/$(DICT_SOURCE).so
#	Create root dictionary objects

directory: $(DICT_OBJ_DIR) $(C_OBJ_DIR)

libs: core detectors main
	@echo " Done making libs"

.PHONY: clean tidy directory

.SECONDARY: $(DICT_DIR)/$(DICT_SOURCE).cpp $(ROOTOBJ)
#	Want to keep the source files created by rootcint after compilation
#	as well as keeping the object file made from those source files

#####################################################################

$(C_OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
#	Compile C++ source files
	$(COMPILER) -c $(CFLAGS) $< -o $@

$(C_OBJ_DIR)/%.o: $(TOOL_DIR)/%.cpp
#	Compile C++ source files
	$(COMPILER) -c $(CFLAGS) $< -o $@

#####################################################################

$(DICT_OBJ_DIR):
#	Make root dictionary object file directory
	mkdir $(DICT_OBJ_DIR)

$(C_OBJ_DIR):
#	Make c++ object file directory
	mkdir $(C_OBJ_DIR)

#####################################################################

$(DICT_OBJ_DIR)/%.o: $(DICT_DIR)/%.cpp
#	Compile rootcint source files
	$(COMPILER) -c $(CFLAGS) $< -o $@

$(DICT_OBJ_DIR)/%.so: $(C_OBJ_DIR)/$(STRUCT_FILE).o $(DICT_OBJ_DIR)/$(DICT_SOURCE).o
#	Generate the root shared library (.so) for the dictionary
	$(COMPILER) -g -shared -Wl,-soname,lib$(DICT_SOURCE).so -o $(DICT_OBJ_DIR)/lib$(DICT_SOURCE).so $(C_OBJ_DIR)/$(STRUCT_FILE).o $(DICT_OBJ_DIR)/$(DICT_SOURCE).o -lc

$(DICT_DIR)/%.cpp: $(INCLUDE_DIR)/$(STRUCT_FILE).h $(DICT_DIR)/LinkDef.h
#	Generate the dictionary source files using rootcint
	@cd $(DICT_DIR); rootcint -f $@ -c $(INCLUDE_DIR)/$(STRUCT_FILE).h $(DICT_DIR)/LinkDef.h

#####################################################################

vikarFront: $(SOURCE_DIR)/vikarFront.cpp
	$(COMPILER) $(CFLAGS) -o $@ $(SOURCE_DIR)/vikarFront.cpp
	@echo " Done making "$@

angleConvert: $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/angleConvert.o
	$(COMPILER) $(CFLAGS) $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/angleConvert.o -o $@
	@echo " Done making "$@

kinematics: $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/Kinematics.o
	$(COMPILER) $(CFLAGS) $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/Kinematics.o -o $@
	@echo " Done making "$@

kindist: $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/Kindist.o
	$(COMPILER) $(CFLAGS) $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/Kindist.o -o $@
	@echo " Done making "$@	

dump: $(TOOL_DIR)/dump.cpp
	$(COMPILER) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)
	@echo " Done making "$@	

test: $(C_OBJ_DIR)/vikar_core.o $(TOOL_DIR)/test.cpp
	$(COMPILER) $(CFLAGS) $(C_OBJ_DIR)/vikar_core.o -o $@ $(TOOL_DIR)/test.cpp

$(PROG): dictionary $(OBJECTS)
	$(COMPILER) $(LDFLAGS) $(OBJECTS) $(ROOTOBJ) -L$(DICT_OBJ_DIR) $(SFLAGS) -o $@ $(LDLIBS)
	@echo " Done making "$(PROG)

clean:
	rm -f obj/*.o
	rm -f $(PROG) $(TOOLS)
	rm -f dict/obj/* dict/RootDict.*
