// renderer.cpp
// Cory Thornsberry
// Oct. 16th, 2014
// Render a VIKAR detector setup scene in three dimensions

#include <cmath>

#include "vikar_core.h"
#include "detectors.h"

#define WINDOW_WIDTH 1024
#define WINDOW_HEIGHT 768
#define DEPTH 1.0

#define CAM_POS_X 0.0
#define CAM_POS_Y 10.0
#define CAM_POS_Z -10.0

#define CAM_POINT_X 0.0
#define CAM_POINT_Y 0.0
#define CAM_POINT_Z 0.0

#define VIEWPORT_PIXEL_SIZE 0.01

// Render the 3-dimensional scene using the Planar ray tracing algorithm (WIP)
void RenderScene(Planar *bar_array, unsigned int num_bars){
	double 
	
	Vector3 view_ray;
	for(unsigned int i = 0; i < WINDOW_HEIGHT; i++){
		for(unsigned int j = 0; j < WINDOW_WIDTH; j++){
			ray.axis[0]
		}
	}
}
