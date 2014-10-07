# Makefile for VIKAR 4.0 (c++ version)
# Based on version 3.11 by Steve Pain

COMPILER = g++

CFLAGS = -g -fPIC -Wall -O3 `root-config --cflags` -Iinclude
LDLIBS = `root-config --libs`
LDFLAGS = `root-config --glibs`
ROOT_INC = `root-config --incdir`

SOURCES = vikar_core.cpp planar.cpp vikar.cpp
OBJECTS = $(addprefix $(C_OBJ_DIR)/,$(SOURCES:.cpp=.o))

TOOL_DIR = ./tools
SOURCE_DIR = ./source
C_OBJ_DIR = ./obj

PROG = vikar

all: $(PROG)
	@echo " Finished Compilation"

libs: core planar main
	@echo " Done making libs"

$(C_OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
#	Compile C++ source files
	$(COMPILER) -c $(CFLAGS) $< -o $@

$(C_OBJ_DIR)/%.o: $(TOOL_DIR)/%.cpp
#	Compile C++ source files
	$(COMPILER) -c $(CFLAGS) $< -o $@

vikarFront: $(SOURCE_DIR)/vikarFront.cpp
	$(COMPILER) $(FFLAGS) -o $@ $(SOURCE_DIR)/vikarFront.cpp
	@echo " Done making "$@

angleConvert: $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/angleConvert.o
	$(COMPILER) $(FFLAGS) $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/angleConvert.o -o $@
	@echo " Done making "$@

Kinematics: $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/Kinematics.o
	$(COMPILER) $(FFLAGS) $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/Kinematics.o -o $@
	@echo " Done making "$@

Kindist: $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/Kindist.o
	$(COMPILER) $(FFLAGS) $(C_OBJ_DIR)/vikar_core.o $(C_OBJ_DIR)/Kindist.o -o $@
	@echo " Done making "$@	

$(PROG): $(OBJECTS)
	$(COMPILER) $(LDFLAGS) $(OBJECTS) -o $@ $(LDLIBS)
	@echo " Done making "$(PROG)

clean:
	rm -f obj/*.o
	rm -f vikar
