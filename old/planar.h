// vikar_lib.h
// Cory Thornsberry

#ifndef PLANAR_H
#define PLANAR_H

#include <string>
#include <cmath>
#include <fstream>
#include <vector>

void DetAng_plan(double, double, unsigned short, double, double, double, bool, bool, double*, double*, double*, double*);

void DetHit_plan(double**, double**, double**, double**, unsigned short, unsigned short*, double*, double, double, unsigned short&, unsigned short&, unsigned short&);

void DetSet_plan_read(const char*, unsigned short, unsigned short*, double*, double*, double*, double*, double*, bool*, bool*, double*, double*, double*, double*, double*, double*, double*, double*);

void det_thick_plan(double, double, double, double&);

void resolution_plan(double, double, double, double, double, double, double, bool, bool, double, double, double, double, unsigned short, unsigned short, double, double, double);

void strag_dE_plan(double, double, double, double, double, double, double, double&, double&, double, double);

#endif
