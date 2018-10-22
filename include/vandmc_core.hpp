/** \file vandmc_core.cpp
 * \brief Handles low level mathematics and calculations.
 *
 * This file contains classes and functions which perform a very
 * specific function for the vandmc program (i.e. matrix and
 * vector calculations) but are not overly useful on their own.
 * Most of the functions in this file are used in other classes.
 *
 * \author C. R. Thornsberry
 * \date Feb. 26th, 2016
 */
#ifndef VIKAR_LIB_H
#define VIKAR_LIB_H

#include <cmath>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "Vector3.hpp"

/////////////////////////////////////////////////////////////////////
// Globals
/////////////////////////////////////////////////////////////////////

extern const double c, pi, deg2rad;
extern const double rad2deg, LN2;

/////////////////////////////////////////////////////////////////////
// Classes
/////////////////////////////////////////////////////////////////////

// To handle the circular dependencies.
class Ray;
class Line;
class Target;

class Ray{
  public:
	Vector3 pos; /// The originating point of the 2d ray.
	Vector3 dir; /// The direction of the ray in 2d space.

	/// Default constructor.
	Ray(){}

	/** Construct a ray by supplying its starting point (x1_, y1_) and
	  * a point through which it passes (x2_, y2_).
	  */	
	Ray(const double &x1_, const double &y1_, const double &x2_, const double &y2_);

	/** Construct a ray by supplying its starting point pos_ and
	  * and its direction (dx_, dy_).
	  */
	Ray(const Vector3 &pos_, const double &dx_, const double &dy_);

	/** Construct a ray by supplying its starting point (x_, y_) and
	  * and its direction dir_.
	  */
	Ray(const double &x_, const double &y_, const Vector3 &dir_);

	/** Construct a ray by supplying its starting point pos_ and
	  * and its direction dir_.
	  */
	Ray(const Vector3 &pos_, const Vector3 &dir_);

	/// Construct a ray from a line segment.
	Ray(const Line &line_);
	
	/// Assignment operator.
	const Ray& operator = (const Ray&);
	
	/// Return true if this ray intersects another ray in 2d space.
	bool Intersect(const Ray &other_, Vector3 &p);

	/// Return true if this ray intersects a line segment in 2d space.
	bool Intersect(const Line &line_, Vector3 &p);
};

class Line{
  public:
 	Vector3 p1; /// The originating point of the 2d line segment.
	Vector3 p2; /// The terminating point of the 2d line segment.
	Vector3 dir; /// The direction vector of the line segment.
	double length; /// The length of the line segment.

	/// Default constructor.
	Line(){ length = 0.0; }

	/** Construct a line segment by supplying its starting point (x1_, y1_) and
	  * its ending point (x2_, y2_).
	  */
	Line(const double &x1_, const double &y1_, const double &x2_, const double &y2_);

	/** Construct a line segment by supplying its starting point pos_ its
	  * direction (dx_, dy_) and its length_.
	  */
	Line(const Vector3 &pos_, const double &dx_, const double &dy_, const double &length_);

	/** Construct a line segment by supplying its starting point (x_, y_) its
	  * direction dir_ and its length_.
	  */
	Line(const double &x_, const double &y_, const Vector3 &dir_, const double &length_);

	/** Construct a line segment by supplying its starting point pos_ its
	  * direction dir_ and its length_.
	  */
	Line(const Vector3 &pos_, const Vector3 &dir_, const double &length_);
	
	/// Construct a line segment from a ray by specifying its length_.
	Line(const Ray &ray_, const double &length_);
	
	/// Assignment operator.
	const Line& operator = (const Line&);

	/// Return true if this line segment intersects a ray in 2d space.
	bool Intersect(const Line &other_, Vector3 &p);

	/// Return true if this line segment intersects another line segment in 2d space.
	bool Intersect(const Ray &ray_, Vector3 &p);
};

class AngularDist{
  private:
	double *com_theta; /// Array for storing the center of mass angle of the distribution (rad).
	double *dsigma_domega; /// Array for storing the differential cross section of the distribution (mb/Sr).
	double *integral; /// Array for storing the integral of the distribution.
	double reaction_xsection; /// The total reaction cross section (mb).
	double rate; /// The expected reaction rate.
	unsigned int num_points; /// The number of entries in the distribution arrays.
	bool init; /// Set to true if the distribution arrays have been initialized.
	
  public:
  	/// Default constructor.
	AngularDist();
	
	/// Destructor;
	~AngularDist();
	
	/// Setup the angular distribution by reading it from a file.
	bool Initialize(const char* fname, const double &beam_intensity=0, Target *targ_=NULL);
	
	/// Setup the angular distribution using arrays.
	bool Initialize(const unsigned int &num_points_, double *angle_, double *xsection_, const double &beam_intensity=0, Target *targ_=NULL);
	
	/// Setup the angular distribution using an isotropic distribution.
	bool Initialize(const double &xsection_);
	
	/// Return the reaction rate (pps).
	double GetRate(){ return rate; }
	
	/// Return the number of entries in the distribution.
	unsigned int GetNumPoints(){ return num_points; }
	
	/// Return the total reaction cross section (mb).
	double GetReactionXsection(){ return reaction_xsection; }
	
	/// Return a random angle sampled from the distribution (rad).
	double Sample();
};

/////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////

bool IsInVector(const std::string &input_, const std::vector<std::string> &str_vector_);
void RandomGauss(double fwhm_, double offset_, Vector3 &beam);
void RandomCircle(double radius_, double offset_, Vector3 &beam);
void RandomHalo(double fwhm_, double offset_, Vector3 &beam);
bool SetBool(std::string input_, std::string text_, bool &output);
bool SetBool(std::string input_, bool &output);
bool Prompt(std::string prompt_);
std::string Parse(std::string input);
void Parse(std::string input, float *arr, unsigned int num_values);
double Order(double input_);
double Dist3d(const Vector3&, const Vector3&);
double dabs(double);
double min(double, double);
double max(double, double);
double frand();
double frand(double, double);
void UnitSphereRandom(Vector3&);
void UnitSphereRandom(double&, double&);
double UnitCircleRandom();
double WrapValue(double, double, double);
unsigned int GetLines(const char*);
double radlength(unsigned int, unsigned int);
double rndgauss0(double);
void straggleA(double&, double, double, double, double, double);
double Interpolate(double, double, double, double, double);
bool Interpolate(const double &x, double &y, double *x_, double *y_, const size_t &len_);

#endif
