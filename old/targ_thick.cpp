// targ_thick.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:41:31 2014

void targ_thick(theta_in,phi_in,thicknessZ,depth,theta_targ,thickness){ 
	// targ_thick 1.0 written by S.D.Pain on 2/2/2005
	//
	//   theta_in = polar angle of particle in lab
	//   phi_in = azimuthal angle of particle in lab
	//   thicknessZ = total target thickness (mg/cm^2)
	//   depth = depth through target for interaction, perpendicular to target plane
	//   theta_targ = angle of target (rotated clockwise about y axis, viewed from above)
	//   thickness = thickness of material the particle must pass through
	
	double theta_in, phi_in, thicknessZ, depth; 
	double x, y, z, length, thickness, theta_targ; 
	double dummy, pi; 
	double newx, newy, newz; 
	double dummyr, dummytheta, dummyphi; 
	pi = 3.14159; 
	
	// Convert the ion's direction to cartesian coordinates
	dummy = 1.0; 
	call sphere2cart(dummy,theta_in,phi_in,x,y,z); 
	length = sqrt(pow(x, 2)+pow(y, 2)+pow(z, 2)); 
	
	// Rotate the ion's vector so it is measured wrt to the target
	newx = x*cos(theta_targ) - z*sin(theta_targ); 
	newz = z*cos(theta_targ) + x*sin(theta_targ); 
	newy = y; 
	
	// Convert the ion's vector back to spherical polars, now wrt the
	// target
	cart2sphere(newx,newy,newz,dummyr,dummytheta,dummyphi); 
	
	// Calculate the thickness seen by the ion. Bloody marvelous.
	if(dummytheta <= (0.5*pi)){ thickness = (thicknessZ-depth)/std::cos(dummytheta); }
	else{ thickness = (-depth)/std::cos(dummytheta); }
	
	// This line was added to account for an occasion where dummytheta was greater
	// than 0.5*pi (only just), but for some reason, cos(dummytheta) was positive
	// which gave a negative thickness.
	if(thickness < 0){ thickness = -thickness; }
} 



