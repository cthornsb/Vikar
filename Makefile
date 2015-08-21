#####################################################################

# Set the RootPixieScan directory
PIXIE_SCAN_DIR = /home/cthorns/RootPixieScan

#####################################################################

COMPILER = g++

CFLAGS = -g -fPIC -Wall -O3 `root-config --cflags` -Iinclude
LDLIBS = `root-config --libs`
LDFLAGS = `root-config --glibs`
ROOT_INC = `root-config --incdir`

SOURCES = vikar_core.cpp detectors.cpp materials.cpp vikar.cpp
SOURCES2 = vikar_core.cpp detectors.cpp materials.cpp
OBJECTS = $(addprefix $(OBJ_DIR)/,$(SOURCES:.cpp=.o))
OBJECTS2 = $(addprefix $(OBJ_DIR)/,$(SOURCES2:.cpp=.o))

TOP_LEVEL = $(shell pwd)
INCLUDE_DIR = $(shell pwd)/include
SOURCE_DIR = $(shell pwd)/source
OBJ_DIR = $(shell pwd)/obj

TOOL_DIR = $(shell pwd)/tools
TOOL_SRC_DIR = $(TOOL_DIR)/source
TOOL_OBJ_DIR = $(TOOL_DIR)/obj

DICT_DIR = $(shell pwd)/dict
DICT_OBJ_DIR = $(DICT_DIR)/obj

# ROOT dictionary stuff
DICT_SOURCE = RootDict
STRUCT_FILE = Structures

ROOTOBJ = $(DICT_OBJ_DIR)/$(DICT_SOURCE).o 
ROOTOBJ += $(OBJ_DIR)/$(STRUCT_FILE).o
SFLAGS = $(addprefix -l,$(DICT_SOURCE))

INSTALL_DIR = ~/bin

# Tools
VIKAR_FRONT_EXE = $(TOOL_DIR)/vikarFront
VIKAR_FRONT_SRC = $(TOOL_SRC_DIR)/vikarFront.cpp
ANGLE_CONVERT_EXE = $(TOOL_DIR)/angleConvert
ANGLE_CONVERT_SRC = $(TOOL_SRC_DIR)/angleConvert.cpp
KINEMATICS_EXE = $(TOOL_DIR)/kinematics
KINEMATICS_SRC = $(TOOL_SRC_DIR)/Kinematics.cpp
KINDIST_EXE = $(TOOL_DIR)/kindist
KINDIST_SRC = $(TOOL_SRC_DIR)/Kindist.cpp
ROOT2RAW_EXE = $(TOOL_DIR)/root2raw
ROOT2RAW_SRC = $(TOOL_SRC_DIR)/root2raw.cpp
INTEGRATOR_EXE = $(TOOL_DIR)/integrator
INTEGRATOR_SRC = $(TOOL_SRC_DIR)/Integrator.cpp
RANGE_EXE = $(TOOL_DIR)/range
RANGE_OBJ = $(TOOL_OBJ_DIR)/range.o
RANGE_SRC = $(TOOL_SRC_DIR)/range.cpp

TOOLS = $(VIKAR_FRONT_EXE) $(ANGLE_CONVERT_EXE) $(KINEMATICS_EXE) $(KINDIST_EXE) $(ROOT2RAW_EXE) $(INTEGRATOR_EXE) $(RANGE_EXE)
TOOLS_OBJ = $(RANGE_OBJ)

# Main executable
EXECUTABLE = vikar

########################################################################

all: directory dictionary $(EXECUTABLE)
	@echo " Finished Compilation"

dictionary:
#	Create root dictionary objects
	@$(PIXIE_SCAN_DIR)/tools/rcbuild.sh -t $(PIXIE_SCAN_DIR)/tools -d $(DICT_DIR) -s $(SOURCE_DIR) -i $(INCLUDE_DIR) -o $(OBJ_DIR)

directory: $(OBJ_DIR) $(TOOL_OBJ_DIR)
#	Make all directories

libs: $(OBJECTS)
#	Make only library objects
	@echo " Done making libs"

tools: $(TOOLS)
#	Make all support tools
	@echo " Done making tools"

#####################################################################

$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
#	Compile C++ source files
	$(COMPILER) -c $(CFLAGS) $< -o $@

$(TOOL_OBJ_DIR)/%.o: $(TOOL_SRC_DIR)/%.cpp
#	Compile C++ source files
	$(COMPILER) -c $(CFLAGS) $< -o $@

#####################################################################

$(OBJ_DIR):
#	Make the object file directory
	@if [ ! -d $@ ]; then \
		echo "Making directory: "$@; \
		mkdir $@; \
	fi
	
$(TOOL_OBJ_DIR):
#	Make the object file directory
	@if [ ! -d $@ ]; then \
		echo "Making directory: "$@; \
		mkdir $@; \
	fi

#####################################################################

$(VIKAR_FRONT_EXE): $(VIKAR_FRONT_SRC)
	$(COMPILER) -g -fPIC -Wall -O3 -Iinclude -o $@ $^

$(ANGLE_CONVERT_EXE): $(OBJ_DIR)/vikar_core.o $(ANGLE_CONVERT_SRC)
	$(COMPILER) -g -fPIC -Wall -O3 -Iinclude $^ -o $@

$(KINEMATICS_EXE): $(OBJ_DIR)/vikar_core.o $(KINEMATICS_SRC)
	$(COMPILER) -g -fPIC -Wall -O3 -Iinclude $^ -o $@

$(KINDIST_EXE): $(OBJ_DIR)/vikar_core.o $(KINDIST_SRC)
	$(COMPILER) -g -fPIC -Wall -O3 -Iinclude $^ -o $@

$(ROOT2RAW_EXE): $(ROOT2RAW_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@ $(LDLIBS)

$(INTEGRATOR_EXE): $(INTEGRATOR_SRC)
	$(COMPILER) -g -fPIC -Wall -O3 -Iinclude $^ -o $@

$(RANGE_EXE): $(OBJECTS2) $(ROOTOBJ) $(RANGE_OBJ)
	$(COMPILER) $(LDFLAGS) $^ -o $@ $(LDLIBS)

########################################################################

$(EXECUTABLE): $(OBJECTS)
	$(COMPILER) $(LDFLAGS) $(OBJECTS) $(ROOTOBJ) -L$(DICT_OBJ_DIR) $(SFLAGS) -o $@ $(LDLIBS)
	@echo " Done making "$(EXECUTABLE)

########################################################################

install: tools
#	Install tools into the install directory
	@echo "Installing tools to "$(INSTALL_DIR)
	@ln -s -f $(KINEMATICS_EXE) $(INSTALL_DIR)/kinematics
	@ln -s -f $(KINDIST_EXE) $(INSTALL_DIR)/kindist
	@ln -s -f $(ROOT2RAW_EXE) $(INSTALL_DIR)/root2raw

#####################################################################

tidy: clean_obj clean_tools clean_dict

clean: clean_obj

clean_obj:
	@echo "Cleaning up..."
	@rm -f $(OBJ_DIR)/*.o
	@rm -f $(EXECUTABLE)
	
clean_dict:
	@echo "Removing ROOT dictionaries..."
	@rm -f $(DICT_DIR)/$(DICT_SOURCE).cpp $(DICT_DIR)/$(DICT_SOURCE).h $(DICT_OBJ_DIR)/*.so
	@rm -f $(ROOTOBJ) $(SOURCE_DIR)/Structures.cpp $(INCLUDE_DIR)/Structures.h $(DICT_DIR)/LinkDef.h

clean_tools:
	@echo "Removing tools..."
	@rm -f $(TOOLS) $(TOOLS_OBJ)
