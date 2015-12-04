#####################################################################

# Set the RootPixieScan directory
SIMPLE_SCAN_DIR = $(HOME)/Research/Pixie16/SimplePixieScan

#####################################################################

COMPILER = g++

CFLAGS = -Wall -O3 -Iinclude
RFLAGS = -Wall -O3 `root-config --cflags` -Iinclude
LDLIBS = `root-config --libs`
LDFLAGS = `root-config --glibs`
ROOT_INC = `root-config --incdir`

SOURCES = vikar_core.cpp kindeux.cpp geometry.cpp detectors.cpp materials.cpp
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

RCBUILD_DIR = $(SIMPLE_SCAN_DIR)/rcbuild

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
RELATIVISTIC_EXE = $(TOOL_DIR)/relativistic
RELATIVISTIC_SRC = $(TOOL_SRC_DIR)/relativistic.cpp
LIGHT_EXE = $(TOOL_DIR)/light
LIGHT_SRC = $(TOOL_SRC_DIR)/light.cpp

TOOLS = $(VIKAR_FRONT_EXE) $(ANGLE_CONVERT_EXE) $(KINEMATICS_EXE) $(KINDIST_EXE) $(ROOT2RAW_EXE) \
        $(INTEGRATOR_EXE) $(RANGE_EXE) $(ELOSS_EXE) $(TEST_SETUP_EXE) $(TEST_VIEWER_EXE) $(RELATIVISTIC_EXE) \
        $(LIGHT_EXE)

RENDERER_DIR = $(TOOL_DIR)/qtfiles
RENDERER_EXE = renderer

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
	@$(RCBUILD_DIR)/rcbuild.sh -c $(SIMPLE_SCAN_DIR) -d $(DICT_DIR) -s $(SOURCE_DIR) -i $(INCLUDE_DIR) -o $(OBJ_DIR)

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

$(ANGLE_CONVERT_EXE): $(OBJ_DIR)/vikar_core.o $(OBJ_DIR)/kindeux.o $(ANGLE_CONVERT_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@

$(KINEMATICS_EXE): $(OBJ_DIR)/vikar_core.o $(OBJ_DIR)/kindeux.o $(KINEMATICS_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@

$(KINDIST_EXE): $(OBJ_DIR)/vikar_core.o $(OBJ_DIR)/kindeux.o $(KINDIST_SRC)
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

$(RELATIVISTIC_EXE): $(OBJECTS) $(RELATIVISTIC_SRC)
	$(COMPILER) $(CFLAGS) $^ -o $@

$(LIGHT_EXE): $(OBJECTS) $(LIGHT_SRC)
	$(COMPILER) $(RFLAGS) $^ -o $@ $(LDLIBS)
	
$(RENDERER_EXE): $(OBJECTS)
#	Use qmake to generate the renderer Makefile, if it doesn't exist.
	@if [ ! -e $(RENDERER_DIR)/Makefile ]; then \
		echo " Calling qmake to generate renderer Makefile"; \
		cd $(RENDERER_DIR) && qmake; \
	fi
	
#	Call the renderer Makefile.
	@cd $(RENDERER_DIR) && $(MAKE)

########################################################################

install: tools
#	Install tools into the install directory
	@echo "Installing tools to "$(INSTALL_DIR)
	@ln -s -f $(ROOT2RAW_EXE) $(INSTALL_DIR)/root2raw
	@ln -s -f $(RELATIVISTIC_EXE) $(INSTALL_DIR)/relativistic

#####################################################################

tidy: clean_obj clean_tools clean_dict clean_renderer
	@rm -f $(EXECUTABLE)

clean: clean_obj

clean_obj:
	@echo "Cleaning up object files..."
	@rm -f $(OBJ_DIR)/*.o
	@rm -f $(TOOL_DIR)/obj/*
	@rm -f $(RENDERER_DIR)/obj/*
	@rm -f $(ROOTOBJ)
	
clean_dict:
	@echo "Removing ROOT dictionaries..."
	@rm -f $(DICT_DIR)/$(DICT_SOURCE).cpp $(DICT_DIR)/$(DICT_SOURCE).h $(DICT_OBJ_DIR)/*.so
	@rm -f $(SOURCE_DIR)/Structures.cpp $(INCLUDE_DIR)/Structures.h $(DICT_DIR)/LinkDef.h

clean_tools:
	@echo "Removing tools..."
	@rm -f $(TOOLS)

clean_renderer:
	@rm -f $(RENDERER_DIR)/moc/*
	@rm -f $(RENDERER_DIR)/ui/*
	#@cd $(RENDERER_DIR) && $(MAKE) clean
