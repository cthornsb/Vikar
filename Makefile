# Makefile for VIKAR 4.0 (c++ version)
# Based on version 3.11 by Steve Pain

COMPILER = g++

CFLAGS = -g -fPIC -Wall -O3 `root-config --cflags`
LDLIBS = `root-config --libs`
LDFLAGS = `root-config --glibs`
ROOT_INC = `root-config --incdir`

SOURCES = vikar_core.cpp planar.cpp vikar.cpp
OBJECTS = $(addprefix $(C_OBJ_DIR)/,$(SOURCES:.cpp=.o))

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

front:
	$(COMPILER) $(FFLAGS) -o vikarFront ./source/vikarFront.cpp
	@echo " Done making vikarFront"

angle:
	$(COMPILER) $(FFLAGS) -c -o ./tools/angleConvert.o ./tools/angleConvert.cpp
	$(COMPILER) $(FFLAGS) vikar_core.o ./tools/angleConvert.o -o angleConvert
	@echo " Done making angleConvert"

kind:
	$(COMPILER) $(FFLAGS) -c -o ./tools/Kinematics.o ./tools/Kinematics.cpp
	$(COMPILER) $(FFLAGS) ./obj/vikar_core.o ./tools/Kinematics.o -o Kinematics
	@echo " Done making Kinematics"

$(PROG): $(OBJECTS)
	$(COMPILER) $(LDFLAGS) $(OBJECTS) -o $@ $(LDLIBS)
	@echo " Done making "$(PROG)

clean:
	rm -f obj/*.o
	rm -f vikar
