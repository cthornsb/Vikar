#####################################################################

# Set the RootPixieScan directory
PIXIE_SCAN_DIR = $(HOME)/RootPixieScan

#####################################################################

COMPILER = g++

CFLAGS = -g -fPIC -Wall -O3 -Iinclude
RFLAGS = -g -fPIC -Wall -O3 `root-config --cflags` -Iinclude
LDLIBS = `root-config --libs`
LDFLAGS = `root-config --glibs`
ROOT_INC = `root-config --incdir`

SOURCES = vikar_core.cpp detectors.cpp materials.cpp
OBJECTS = $(addprefix $(OBJ_DIR)/,$(SOURCES:.cpp=.o))

TOP_LEVEL = $(shell pwd)
INCLUDE_DIR = $(shell pwd)/include
SOURCE_DIR = $(shell pwd)/source
OBJ_DIR = $(shell pwd)/obj

TOOL_DIR = $(shell pwd)/tools
TOOL_SRC_DIR = $(TOOL_DIR)/source

DICT_DIR = $(shell pwd)/dict
DICT_OBJ_DIR = $(DICT_DIR)/obj

# ROOT dictionary stuff
DICT_SOURCE = RootDict
STRUCT_FILE = Structures

DEFINITION_FILE = $(DICT_DIR)/def.struct
SHARED_OBJ_FILE = $(DICT_OBJ_DIR)/libRootDict.so

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
RANGE_SRC = $(TOOL_SRC_DIR)/range.cpp
ELOSS_EXE = $(TOOL_DIR)/eloss
ELOSS_SRC = $(TOOL_SRC_DIR)/eloss.cpp
TEST_SETUP_EXE = $(TOOL_DIR)/testSetup
TEST_SETUP_SRC = $(TOOL_SRC_DIR)/testSetup.cpp
TEST_VIEWER_EXE = $(TOOL_DIR)/testViewer
TEST_VIEWER_SRC = $(TOOL_SRC_DIR)/testViewer.cpp

TOOLS = $(VIKAR_FRONT_EXE) $(ANGLE_CONVERT_EXE) $(KINEMATICS_EXE) $(KINDIST_EXE) $(ROOT2RAW_EXE) \
        $(INTEGRATOR_EXE) $(RANGE_EXE) $(ELOSS_EXE) $(TEST_SETUP_EXE) $(TEST_VIEWER_EXE)

# Main executable
MAIN_SRC = $(SOURCE_DIR)/vikar.cpp
MAIN_OBJ = $(OBJ_DIR)/vikar.o

EXECUTABLE = vikar

########################################################################

$(EXECUTABLE): $(OBJ_DIR) $(SHARED_OBJ_FILE) $(OBJECTS) $(MAIN_OBJ)
#	Link the executable
	$(COMPILER) $(LDFLAGS) $(OBJECTS) $(MAIN_OBJ) $(ROOTOBJ) -L$(DICT_OBJ_DIR) $(SFLAGS) -o $@ $(LDLIBS)
	@echo " Done making "$(EXECUTABLE)

#	Make all directories and setup the default config files
	@if [ ! -e ./input/default.in ] || [ ! -e ./detectors/default.det ] || [ ! -e ./xsections/default.xsect ]; then \
		echo " Unpacking default config files"; \
		tar -xf $(TOP_LEVEL)/config.tar; \
	fi

########################################################################

$(SHARED_OBJ_FILE): $(OBJ_DIR) $(DEFINITION_FILE)
#	Create root dictionary objects
	@$(PIXIE_SCAN_DIR)/tools/rcbuild.sh -t $(PIXIE_SCAN_DIR)/tools -d $(DICT_DIR) -s $(SOURCE_DIR) -i $(INCLUDE_DIR) -o $(OBJ_DIR)

libs: $(OBJ_DIR) $(OBJECTS)
#	Make only library objects
	@echo " Done making libs"

tools: $(TOOLS)
#	Make all support tools
	@echo " Done making tools"
	
.PHONY: libs tools install tidy clean clean_obj clean_dict clean_tools

#####################################################################

$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
#	Compile C++ source files
	$(COMPILER) -c $(CFLAGS) $< -o $@

$(MAIN_OBJ): $(MAIN_SRC)
#	Compile the main
	$(COMPILER) -c $(RFLAGS) $< -o $@

#####################################################################

$(OBJ_DIR):
#	Make the object file directory
	mkdir -p $@

#####################################################################

$(VIKAR_FRONT_EXE): $(VIKAR_FRONT_SRC)
	$(COMPILER) $(CFLAGS) -o $@ $^

$(ANGLE_CONVERT_EXE): $(OBJ_DIR)/vikar_core.o $(ANGLE_CONVERT_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@

$(KINEMATICS_EXE): $(OBJ_DIR)/vikar_core.o $(KINEMATICS_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@

$(KINDIST_EXE): $(OBJ_DIR)/vikar_core.o $(KINDIST_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@

$(ROOT2RAW_EXE): $(ROOT2RAW_SRC)
	$(COMPILER) $(RFLAGS) $^ -o $@ $(LDLIBS)

$(INTEGRATOR_EXE): $(INTEGRATOR_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@

$(RANGE_EXE): $(OBJECTS) $(RANGE_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@

$(ELOSS_EXE): $(OBJECTS) $(ELOSS_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@
	
$(TEST_SETUP_EXE): $(SHARED_OBJ_FILE) $(OBJECTS) $(TEST_SETUP_SRC)
	$(COMPILER) $(RFLAGS) $(OBJECTS) $(ROOTOBJ) -L$(DICT_OBJ_DIR) $(SFLAGS) $(TEST_SETUP_SRC) -o $@ $(LDLIBS)
	
$(TEST_VIEWER_EXE): $(TEST_VIEWER_SRC)
	$(COMPILER) $(RFLAGS) $^ -o $@ $(LDLIBS)

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
	@rm -f $(TOOLS)
