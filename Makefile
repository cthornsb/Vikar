# Makefile for VIKAR 4.0 (c++ version)
# Based on version 3.11 by Steve Pain

COMPILER = g++
FFLAGS = -Wall

SOURCES = ./obj/vikar_core.o ./obj/planar.o ./obj/vikar.o

all: libs vikar
	@echo " Finished Compilation"

libs: core planar main
	@echo " Done making libs"

core:
	$(COMPILER) $(FFLAGS) -c -o ./obj/vikar_core.o ./source/vikar_core.cpp

planar:
	$(COMPILER) $(FFLAGS) -c -o ./obj/planar.o ./source/planar.cpp

main:
	$(COMPILER) $(FFLAGS) -c -o ./obj/vikar.o ./source/vikar.cpp
	
vikar:
	$(COMPILER) $(FFLAGS) -o vikar $(SOURCES)
	@echo " Done making vikar"

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

$(PROG): libs
	$(COMPILER) $(FFLAGS) -o $(PROG) $(SOURCES)
	@echo " Done making"$(PROG)

clean:
	rm -f obj/*.o
	rm -f vikar
