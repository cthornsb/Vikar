# Cory Thornsberry
# 7Be(d,n)8B Experimental Setup
# December 22nd, 2013

"""
Description of model coordinate system...
	+X is defined as the beam axis
	+Z is defined as beam-right
	+Y is defined as vertically upward
	
Note! X,Z axes are swapped in newVIKAR!
"""

import wx
import math
import visual

class Vector3:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

# Returns a color on the spectrum
# Values range from (0,1]
def SmoothSpectrum(value, maximum):
	code = int((float(value)/maximum)*1195)
	if code == 0: output = (153,204,255)
	elif code < 175: output = (175-code,0,255)
	elif code < 430: output = (0,(code-175),255)
	elif code < 685: output = (0,255,255-(code-430))
	elif code < 940: output = ((code-685),255,0)
	elif code < 1195: output = (255,255-(code-940),0)
	else: output = (255,0,0)
	return (output[0]/255.0,output[1]/255.0,output[2]/255.0)

# Calculate proper bar spacing for a wall of VANDLE bars
# Leave half gaps at either edge for clearance to other walls
def BarSpacing(total_width, bar_width, num_bars):
	return (total_width-num_bars*bar_width)/num_bars

# Calculate the angular spacing between adjacent bars
def BarSpacingAngle(radius, spacing, bar_width, curved=True):
	if not curved: return 2.0*math.asin((spacing+bar_width)/(2.0*radius))
	else: return (spacing+bar_width)/radius

def GetTuple(inp):
	return (float(inp[0]), float(inp[2]), float(inp[1]))

# Convert x from inches to meters
def ConvIM(x):
	return 0.0254*x

# Convert x, y, and z from inches to meters and return the 3-tuple (x', y', z')
def ConvIM3(x, y, z):
	return(0.0254*x, 0.0254*y, 0.0254*z)

# Return the cartesian coordinates of a point in spherical space as the 3-tuple (x, y, z)
def GetCartesianCoords(r, theta, phi):
	return(r*math.sin(theta)*math.cos(phi), r*math.cos(theta), r*math.sin(theta)*math.sin(phi))

# Return the spherical coordinates of a point in cartesian space as the 3-tuple (r, theta, phi)
def GetSphericalCoords(coord, deg=False):
	r = math.sqrt(coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2])
	if r == 0: return(0.0, 0.0, 0.0)
	else: 
		if(not deg): return(r, WrapValue(math.acos(coord[1]/r),0,math.pi), WrapValue(math.atan2(coord[2], coord[0]),0,2*math.pi))
		else: return(r, WrapValue(RadToDeg(math.acos(coord[1]/r)),0,180), WrapValue(RadToDeg(math.atan2(coord[2], coord[0])),0,360))

# Return the VPython color code for an 255-RGB code as the 3-tuple (R, G, B)
def GetRGBcode(r, g, b):
	return(r/255.0, g/255.0, b/255.0)

# Convert an angle in raidans to degrees (rad must be in radians)
def RadToDeg(rad):
	return((180/math.pi)*rad)
	
# Convert an angle in degrees to radians (deg must be in degrees)
def DegToRad(deg):
	return((math.pi/180)*deg)

# Translate coords1 by a vector coords2
def Translate(coords1, coords2, factor=1):
	if(len(coords1) == len(coords2)): return(coords1[0]+factor*coords2[0],coords1[1]+factor*coords2[1],coords1[2]+factor*coords2[2])
	else: 
		print " Warning: Vectors must have the same number of components, ignoring command"
		return coords1

# Translate in the XZ plane (phi must be in radians)
def TranslateXZ(coords, d, phi):
	return(coords[0]+d*math.cos(phi),coords[1],coords[2]+d*math.sin(phi))

# Translate in 3d space by a distance d (theta and phi must be in radians)
def Translate3D(coords, d, theta, phi):
	return(coords[0]+d*math.cos(phi),coords[1]+d*math.cos(theta),coords[2]+d*math.sin(phi))

# Return a string from a coordinate tuple of arbitrary length
def GetString(coords):
	return_str = "("
	for i in range(len(coords)):
		if(i == 0): return_str += str(coords[0])
		else: return_str += ", " + str(coords[i])
	
	return_str += ")"
	return return_str

# Wrap a value between min_val and max_val
def WrapValue(value, min_val, max_val):
	if(value < min_val): return max_val-(min_val-value)
	elif(value > max_val): return min_val+(value-max_val)
	else: return value

# Class for the drawing of the VANDLE bars
class VandleBar:
	def __init__(self, target=(0,0,0), R=1, theta=0, phi=0, bar_length=1, bar_width=1, bar_height=1, include_pmt=False, bar_color=visual.color.white):
		self.target = target
		self.R = R
		self.theta = theta
		self.phi = phi
		
		self.position = Translate(self.target, GetCartesianCoords(self.R, self.theta, self.phi))
		self.body = visual.box(pos=self.position, axis=(0, bar_length, 0), width=bar_width, height=bar_height, color=bar_color, material=visual.materials.emissive)
		if(include_pmt):
			self.pmt1 = visual.cylinder(pos=(bar_length/2.0,0,0), axis=(bar_length/4.0,0,0), radius=bar_width/2.0)
			self.pmt2 = visual.cylinder(pos=(-bar_length/2.0,0,0), axis=(-bar_length/4.0,0,0), radius=bar_width/2.0)
		self.rotation = [0.0, 0.0, 0.0]
		
	def SetSpherical(self, R, theta, phi):
		self.R = R
		self.theta = theta
		self.phi = phi
		self.position = Translate(self.target, GetCartesianCoords(self.R, self.theta, self.phi))
		self.body.pos = self.position
		
	def SetPosition(self, pos):
		self.body.pos[0] = pos[0]
		self.body.pos[1] = pos[1]
		self.body.pos[2] = pos[2]
		
	def SetX(self, x): self.body.pos[0] = x
	def SetY(self, y): self.body.pos[1] = y
	def SetZ(self, z): self.body.pos[2] = z
		
	def SetRotation(self, theta, phi):
		self.rotation[0] = theta
		self.rotation[1] = phi	
		self.rotation[2] = psi
		self.body.rotate(angle=theta, axis=(0, 0, 1))
		self.body.rotate(angle=phi, axis=(1, 0, 0))
		self.body.rotate(angle=psi, axis=(0, 1, 0))
	
	def SetTheta(self, theta):
		self.rotation[0] = theta
		self.body.rotate(angle=theta, axis=(0, 1, 0))
	
	def SetPhi(self, phi):
		self.rotation[1] = phi
		self.body.rotate(angle=phi, axis=(0, 0, 1))
	
	def SetPsi(self, psi):
		self.rotation[2] = psi
		self.body.rotate(angle=psi, axis=(1, 0, 0))
	
	def SetColor(self, color_):
		self.body.color = color_
	
	# [x, y, z]
	def GetSpherical(self): return(self.R, self.theta, self.phi)
	def GetPosition(self): return self.position
	def GetX(self): return self.body.pos[0]
	def GetY(self): return self.body.pos[1]
	def GetZ(self): return self.body.pos[2]
	
	# [theta, phi]	
	def GetRotation(self): return self.rotation
	def GetTheta(self): return self.rotation[0]
	def GetPhi(self): return self.rotation[1]
	def GetPsi(self): return self.rotation[2]
	
	def GetLength(self): return self.body.length
	def GetWidth(self): return self.body.width
	def GetHeight(self): return self.body.height
		
	def GetSolidAngle(self, R):
		return 4*math.asin(1/math.sqrt((2*R/self.body.width)**2 + 1))/math.sqrt((2*R/self.body.length)**2 + 1)
		
	def GetSubtendedAngles(self, R):
		return [2*math.asin(1/math.sqrt((2*R/self.body.length)**2 + 1)), 2*math.asin(1/math.sqrt((2*R/self.body.width)**2 + 1))]
		
	def GetSubtendedTheta(self, R):
		return 2*math.asin(1/math.sqrt((2*R/self.body.length)**2 + 1))
		
	def GetSubtendedPhi(self, R):
		return 2*math.asin(1/math.sqrt((2*R/self.body.width)**2 + 1))


# Class for the drawing of the lab model
class BeExpSetup:
	def __init__(self, room_=1, beam_axis_=0, target_pos_=(0,0,0)):
		self.labScene = visual.display(title="7Be(d,n)8B Experiment", width=800, height=600, background=GetRGBcode(153,204,255))
		axisx = visual.box(pos=(0,0,0), axis=(10.0,0,0), width=0.05, height=0.05, color=visual.color.red)
		axisy = visual.box(pos=(0,0,0), axis=(0,10.0,0), width=0.05, height=0.05, color=visual.color.blue)
		axisz = visual.box(pos=(0,0,0), axis=(0,0,10.0), width=0.05, height=0.05, color=visual.color.green)
		labelx = visual.label(pos=(5.0,0,0), text="Z-Axis")
		labely = visual.label(pos=(0,5.0,0), text="Y-Axis")
		labelz = visual.label(pos=(0,0,5.0), text="X-Axis")
		self.labScene.center = target_pos_
		self.labScene.autoscale = False
		
		self.room = room_
		self.beam_axis = beam_axis_
		self.target_pos = target_pos_
		
		self.Floors = []
		self.Walls = []
		self.Columns = []
		self.Others = []		
		
		if(self.room == 1): self.BuildRoom1()
		elif(self.room == 2): self.BuildRoom2()

		chamber_radius = 0.25
		self.Beamline1 = visual.cylinder(pos=Translate(self.target_pos,GetCartesianCoords(chamber_radius, math.pi/2.0, DegToRad(180+self.beam_axis))), axis=ConvIM3(71.75,0,-71.75*math.tan(DegToRad(180-self.beam_axis))), radius=ConvIM(1.75), color=visual.color.blue) # East beamline
		self.Beamline2 = visual.cylinder(pos=Translate(self.target_pos,GetCartesianCoords(chamber_radius, math.pi/2.0, DegToRad(self.beam_axis))), axis=ConvIM3(-217.5,0,217.5*math.tan(DegToRad(180-self.beam_axis))), radius=ConvIM(1.75), color=visual.color.blue) # West beamline
		self.OneMeterChamber = visual.cylinder(pos=self.target_pos, axis=(0,chamber_radius*2,0), radius=chamber_radius, color=visual.color.blue)
		self.OneMeterChamber.pos[1] = -0.5
			
	def BuildRoom1(self):
		self.Floors.append(visual.box(pos=ConvIM3(0,-69.0,24.625), axis=ConvIM3(241.5,0,0), width=ConvIM(346.5), height=0.01)) # Main floor
		self.Floors.append(visual.box(pos=ConvIM3(144.625,-69.0,0), axis=ConvIM3(47.75,0,0), width=ConvIM(72.75), height=0.01)) # Beam pipe access floor
		self.Walls.append(visual.box(pos=ConvIM3(0.0,29.425,-148.625), axis=ConvIM3(241.5,0.0,0.0), width=0.01, height=ConvIM(196.85))) # North Wall
		self.Walls.append(visual.box(pos=ConvIM3(0.0,29.425,197.875), axis=ConvIM3(241.5,0.0,0.0), width=0.01, height=ConvIM(196.85)))	# South Wall
		self.Walls.append(visual.box(pos=ConvIM3(120.75,29.425,-92.5), axis=(0.01,0,0), width=ConvIM(112.25), height=ConvIM(196.85))) # North-East Wall
		self.Walls.append(visual.box(pos=ConvIM3(144.625,-27.125,-36.375), axis=ConvIM3(47.75,0,0), width=0.01, height=ConvIM(83.75))) # North beam pipe access wall
		self.Walls.append(visual.box(pos=ConvIM3(168.5,-27.125,0), axis=(0.01,0,0), width=ConvIM(72.75), height=ConvIM(83.75))) # East beam pipe access wall
		self.Walls.append(visual.box(pos=ConvIM3(120.75,71.3,0), axis=(0.01,0,0), width=ConvIM(72.75), height=ConvIM(113.1))) # Opening header wall
		self.Walls.append(visual.box(pos=ConvIM3(144.625,-27.125,36.375), axis=ConvIM3(47.75,0,0), width=0.01, height=ConvIM(83.75))) # South beam pipe access wall
		self.Walls.append(visual.box(pos=ConvIM3(120.75,29.425,117.125), axis=(0.01,0,0), width=ConvIM(161.5), height=ConvIM(196.85))) # South-East Wall
		self.Others.append(visual.box(pos=ConvIM3(96.75,-9.04,114.385), axis=ConvIM3(48.0,0,0), width=ConvIM(156.02), height=ConvIM(119.92))) # Shield door
		self.Columns.append(visual.box(pos=ConvIM3(-93.75,29.425,109.375), axis=ConvIM3(14.0,0,0), width=ConvIM(14), height=ConvIM(196.85))) # South Column
		self.Columns.append(visual.box(pos=ConvIM3(-93.75,29.425,-148.625), axis=ConvIM3(14.0,0,0), width=ConvIM(14), height=ConvIM(196.85))) # North Column

	def BuildRoom2(self):
		self.Floors.append(visual.box(axis=(19.1, 0, 0), width=9.5, height=0.01)) # Main floor
		self.Floors.append(visual.box(pos=(-6.85, 0, -7.8), axis=(5.4, 0, 0), width=6.1, height=0.01)) # Small northern floor
		self.Walls.append(visual.box(pos=(-9.55, 0, -3.05), axis=(0.01, 0, 0), width=15.6, height=5)) # West wall 
		self.Walls.append(visual.box(pos=(9.55, 0, 0), axis=(0.01, 0, 0), width=9.5, height=5)) # East wall
		self.Walls.append(visual.box(pos=(0, 0, 4.75), axis=(19.1, 0, 0), width=0.01, height=5)) # South wall
		self.Walls.append(visual.box(pos=(2.7, 0, -4.75), axis=(13.9, 0, 0), width=0.01, height=5)) # North Wall
		self.Walls.append(visual.box(pos=(-4.15, 0, -7.8), axis=(0.01, 0, 0), width=6.1, height=5)) # East small wall
		self.Walls.append(visual.box(pos=(-6.85, 0, -10.85), axis=(5.4, 0, 0), width=0.01, height=5)) # North small wall
		self.Others.append(visual.box(pos=(9.25, 0, -0.95), axis=(0.6, 0, 0), width=4.1, height=5)) # Shield wall
		self.Columns.append(visual.box(pos=(-7.65, 0, 1.85), axis=(0.6, 0, 0), width=0.6, height=5)) # Southernmost row column a
		self.Columns.append(visual.box(pos=(-2.05, 0, 1.85), axis=(0.6, 0, 0), width=0.6, height=5))	# Southernmost row column b
		self.Columns.append(visual.box(pos=(3.55, 0, 1.85), axis=(0.6, 0, 0), width=0.6, height=5)) # Southernmost row column c
		self.Columns.append(visual.box(pos=(-7.65, 0, -4.85), axis=(0.6, 0, 0), width=0.6, height=5)) # Middle row column a
		self.Columns.append(visual.box(pos=(-2.05, 0, -4.85), axis=(0.6, 0, 0), width=0.6, height=5)) # Middle row column b
		self.Columns.append(visual.box(pos=(3.55, 0, -4.85), axis=(0.6, 0, 0), width=0.6, height=5))	# Middle row column c
		self.Columns.append(visual.box(pos=(-7.65, 0, -10.75), axis=(0.6, 0, 0), width=0.6, height=5)) # Northernmost row column a
		
	def KeyPress(self, s):
		if(s == " "): print "Camera Position =", GetString(GetSphericalCoords(Translate(self.labScene.mouse.camera,self.labScene.center,-1), True))
		elif(s == "x"): self.labScene.forward = (1,0,0)
		elif(s == "y"): self.labScene.forward = (0,1,0)
		elif(s == "z"): self.labScene.forward = (0,0,1)
		elif(s == "t"):
			self.labScene.forward = (0,0,1)
			self.labScene.forward = (0,-1,0)
		elif(s == "ctrl+x"): self.labScene.forward = (-1,0,0)
		elif(s == "ctrl+y"): self.labScene.forward = (0,-1,0)
		elif(s == "ctrl+z"): self.labScene.forward = (0,0,-1)
		elif(s == "ctrl+t"):
			self.labScene.forward = (1,0,0)
			self.labScene.forward = (0,-1,0)
		elif(s == "f"): self.labScene.forward = (1,-1,0)
		elif(s == "c"): self.labScene.center = self.target_pos
		elif(s in ["w","a","s","d"]): 
			cam_pos = GetSphericalCoords(Translate(self.labScene.mouse.camera,self.labScene.center,-1))
			if(s == "w"): self.labScene.center = TranslateXZ(self.labScene.center, 0.1, WrapValue(cam_pos[2]+math.pi,0,2*math.pi))
			if(s == "a"): self.labScene.center = TranslateXZ(self.labScene.center, 0.1, WrapValue(cam_pos[2]+math.pi/2.0,0,2*math.pi))
			if(s == "s"): self.labScene.center = TranslateXZ(self.labScene.center, 0.1, cam_pos[2])
			if(s == "d"): self.labScene.center = TranslateXZ(self.labScene.center, 0.1, WrapValue(cam_pos[2]-math.pi/2.0,0,2*math.pi))
		elif(s == "ctrl+w"): self.labScene.center[1] += 0.1
		elif(s == "ctrl+s"): self.labScene.center[1] -= 0.1
		
	def PrintHelp(self):
		print "Help---------------------"
		print " Press \'space\' to show position about the center of the camera view"
		print " Press \'x\', \'y\', or \'z\' to look along that axis"
		print " Press \'t\' to look at the target from the top down parallel to the z-axis"
		print " Press \'f\' to look at the target along the x-axis at 45 degrees"
		print " Use \'w\', \'a\', \'s\', and \'d\' to move the center of the camera in the x-z plane"
		print " Press \'c\' to reset the camera to the center of the target\n"
	
		print " Hold \'ctrl\' and press \'x\', \'y\', or \'z\' to look down that axis"
		print " Hold \'ctrl\' and press \'t\' to look from the top down parallel to the x-axis"
		print " Hold \'ctrl\' and use \'w\' and \'s\' to move the center of the camera up or down\n"
		
		print " Press \'p\' to generate NewVIKAR detector data file for VANDLE setup"
		print " Hold \'ctrl\' and Press \'p\' to generate VIKAR detector data file\n"
		
	def ColorFloors(self, color_): 
		for i in range(len(self.Floors)): self.Floors[i].color = color_
		
	def ColorWalls(self, color_): 
		for i in range(len(self.Walls)): self.Walls[i].color = color_
		
	def ColorOthers(self, color_): 
		for i in range(len(self.Others)): self.Others[i].color = color_
		
	def ColorColumns(self, color_): 
		for i in range(len(self.Columns)): self.Columns[i].color = color_
	
	def ColorRoom(self, color_):
		for i in range(len(self.Floors)): self.Floors[i].color = color_
		for i in range(len(self.Walls)): self.Walls[i].color = color_
		for i in range(len(self.Others)): self.Others[i].color = color_
		for i in range(len(self.Columns)): self.Columns[i].color = color_


# Generate VIKAR planar detector string
# Offset is the physical position of the target
# Angle is the beam-axis angle (must be in degrees)
def VIKARdet(detector, offset=(0,0,0), angle=0):
	# Length(cm) Width(cm) Strips Z(cm) X(cm) Y(cm) Orientation dE_thickness(um) E_thickness(um) dE_res(MeV) 
	# dE_pres(at 5.8 MeV)(cm) E_res(MeV) E_pres(at 5.8 MeV)(cm) dE_pres_energy(cm) E_pres_energy(cm)

	# Rotate coordinate system
	cartesian_pos = detector.GetPosition()
	if not angle == 0:
		spherical_pos = detector.GetSpherical()
		phi = WrapValue(spherical_pos[2] - DegToRad(angle), 0, 2*math.pi)
		cartesian_pos = Translate(offset, GetCartesianCoords(spherical_pos[0], spherical_pos[1], phi))

	output_str = str(detector.GetLength()*100.0) # Length (cm)
	output_str += "\t" + str(detector.GetWidth()*100.0) # Width (cm)
	output_str += "\t1\t" + str(cartesian_pos[0]*100.0 - offset[0]*100.0) # VIKAR Z-position (cm)
	output_str += "\t" + str(cartesian_pos[2]*100.0 - offset[2]*100.0) # VIKAR X-position (cm)
	output_str += "\t" + str(cartesian_pos[1]*100.0 - offset[1]*100.0) # VIKAR Y-position (cm)
	output_str += "\ty\t50000.0\t1000.0\t0.050\t7.0\t0.030\t7.0\t5.8\t5.8\n" # Default data
	
	return output_str

# This function automatically swaps the X and Z axes
# to meet the requirements of VIKAR code
def NewVIKAR(detector, offset=(0,0,0), angle=0):
	# X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Length(m) Width(m) Depth(m)

	# Rotate coordinate system
	cartesian_pos = detector.GetPosition()
	if not angle == 0:
		spherical_pos = detector.GetSpherical()
		phi = WrapValue(spherical_pos[2] - DegToRad(angle), 0, 2*math.pi)
		cartesian_pos = Translate(offset, GetCartesianCoords(spherical_pos[0], spherical_pos[1], phi))
		
	output_str = str(cartesian_pos[2] - offset[2]) # VIKAR X-position (m)
	output_str += "\t" + str(cartesian_pos[1] - offset[1]) # VIKAR Y-position (m)
	output_str += "\t" + str(cartesian_pos[0] - offset[0]) # VIKAR Z-position (m)
	output_str += "\t" + str(-(detector.GetTheta()+DegToRad(angle))) # Rotation theta about y-axis (rad)	
	output_str += "\t" + str(detector.GetPhi()) # Rotation phi about z-axis (rad)
	output_str += "\t" + str(detector.GetPsi()) # Rotation psi about x-axis (rad)
	if(detector.GetLength() == 0.6 and detector.GetWidth() == 0.03 and detector.GetHeight() == 0.03): output_str += "\tsmall\n"
	elif(detector.GetLength() == 1.2 and detector.GetWidth() == 0.05 and detector.GetHeight() == 0.03): output_str += "\tmedium\n"
	elif(detector.GetLength() == 2.0 and detector.GetWidth() == 0.05 and detector.GetHeight() == 0.05): output_str += "\tlarge\n"
	else:
		output_str += "\t" + str(detector.GetLength()) # Length (m)
		output_str += "\t" + str(detector.GetWidth()) # Width (m)
		output_str += "\t" + str(detector.GetHeight()) + "\n" # Depth (m)

	return output_str

# +x is the beam direction
# +y is vertically upward
# +z is beam right
def main():
	# axis = dx
	# height = dy
	# width = dz
	
	# Target and beam data
	target_origin = ConvIM3(96.75,0.0,6.4308)
	#target_origin = TranslateXZ(target_origin, 1.0, DegToRad(165)) # Shift chamber 1 meter down beam
	beam_axis = 165
	#beam_axis = 0
	
	# Large bar data
	num_large_bars = 0 #45
	large_bar_R = 3.0
	large_start_angle = 20
	large_packing_coeff = 0.2
	large_packing_angle = BarSpacingAngle(large_bar_R, 0.01350, 0.05)
	
	# Small bar data
	num_small_bars = 28
	small_bar_R = 0.5
	small_start_angle = 190 #70
	small_packing_coeff = 0.2
	small_packing_angle = BarSpacingAngle(small_bar_R, 0.01488, 0.03)
	
	# Wall data
	minimum_z = ConvIM(-148.625)
	maximum_x = ConvIM(120.75)

	SmallBars = []
	LargeBars = []
	dphiSmall = 0.0
	dphiLarge = 0.0
	
	myLab = BeExpSetup(room_=1, beam_axis_=beam_axis, target_pos_=target_origin)
	
	for i in range(num_large_bars):
		if large_start_angle+RadToDeg(i*dphiLarge) >= 180: print "Warning! Large bar", i, "intersects beam pipe..."
		LargeBars.append(VandleBar(target=target_origin, R=large_bar_R, theta=math.pi/2.0, phi=DegToRad(large_start_angle+beam_axis)+i*dphiLarge, bar_length=2.0, bar_width=0.05, bar_height=0.05, bar_color=visual.color.red))
		LargeBars[i].SetTheta(-DegToRad(large_start_angle+beam_axis)-i*dphiLarge)
		if i == 0: 
			#dphiLarge = LargeBars[0].GetSubtendedPhi(large_bar_R)
			#dphiLarge += dphiLarge*large_packing_coeff
			dphiLarge = large_packing_angle
			
		#Position Check
		elif LargeBars[i].GetX() >= maximum_x: 
			print " Warning! Large bar", i+1, "intersects East wall...", "Position = (", large_bar_R, ",", RadToDeg(math.pi/2.0), ",", large_start_angle+RadToDeg(i*dphiLarge), ")"
			LargeBars[i].SetColor(visual.color.black)
		elif LargeBars[i].GetZ() <= minimum_z: 
			print " Warning! Large bar", i+1, "intersects North wall...", "Position = (", large_bar_R, ",", RadToDeg(math.pi/2.0), ",", large_start_angle+RadToDeg(i*dphiLarge), ")"
			LargeBars[i].SetColor(visual.color.black)

	for i in range(num_small_bars):	
		if small_start_angle+RadToDeg(i*dphiSmall) >= 180: print "Warning! Small bar", i, "intersects beam pipe..."
		SmallBars.append(VandleBar(target=target_origin, R=small_bar_R, theta=math.pi/2.0, phi=DegToRad(small_start_angle+beam_axis)+i*dphiSmall, bar_length=0.6, bar_width=0.03, bar_height=0.03, bar_color=visual.color.green))
		SmallBars[i].SetTheta(-DegToRad(small_start_angle+beam_axis)-i*dphiSmall)
		if i == 0: 
			#dphiSmall = SmallBars[0].GetSubtendedPhi(small_bar_R)
			#dphiSmall += dphiSmall*small_packing_coeff
			dphiSmall = small_packing_angle

		#Position Check
		if SmallBars[i].GetX() >= maximum_x: 
			print " Warning! Small bar", i+1, "intersects East wall...", "Position = (", small_bar_R, ",", RadToDeg(math.pi/2.0), ",", small_start_angle+RadToDeg(i*dphiSmall), ")"
			SmallBars[i].SetColor(visual.color.black)
		if SmallBars[i].GetZ() <= minimum_z: 
			print " Warning! Small bar", i+1, "intersects North wall...", "Position = (", small_bar_R, ",", RadToDeg(math.pi/2.0), ",", small_start_angle+RadToDeg(i*dphiSmall), ")"
			SmallBars[i].SetColor(visual.color.black)

	# Model information
	print " Target X:", target_origin[0], "m"
	print " Target Y:", target_origin[1], "m"
	print " Target Z:", target_origin[2], "m"
	print " Beam Axis:", beam_axis, "degrees\n"

	print " Number Large Bars:", num_large_bars
	print " Large Bar Radial Distance:", large_bar_R, "m"
	print " Large Bar Angular Packing: 1 bar per", RadToDeg(dphiLarge), "degrees"
	print " Large Bar Start Angle:", large_start_angle, "degrees"
	print " Large Bar End Angle:", large_start_angle+RadToDeg((num_large_bars-1)*dphiLarge), "degrees\n"
	
	print " Number Small Bars:", num_small_bars
	print " Small Bar Radial Distance:", small_bar_R, "m"
	print " Small Bar Angular Packing: 1 bar per", RadToDeg(dphiSmall), "degrees"
	print " Small Bar Start Angle:", small_start_angle, "degrees"
	print " Small Bar End Angle:", small_start_angle+RadToDeg((num_small_bars-1)*dphiSmall), "degrees\n"

	# Main drawing loop
	s_ = ""
	print "Press \'escape\' to quit"
	print "Press \'h\' for help\n"
	while True:
		visual.rate(60)
		if myLab.labScene.kb.keys:
			s_ = myLab.labScene.kb.getkey()
			
			if(s_ == "esc"): wx.Exit()
			elif(s_ == "h"): myLab.PrintHelp()
			elif(s_ == "p"): 
				f = open("/home/cory/Research/vikar311/detectors/VANDLE.det", "w")
				for i in range(len(LargeBars)): f.write(NewVIKAR(LargeBars[i], offset=target_origin, angle=beam_axis))
				for i in range(len(SmallBars)): f.write(NewVIKAR(SmallBars[i], offset=target_origin, angle=beam_axis))
				print "Wrote detector file /home/cory/Research/vikar311/detectors/VANDLE.det"
				f.close()
			elif(s_ == "ctrl+p"): 
				f = open("/home/cory/Research/vikar311/detectors/LargeBars.det", "w")
				f2 = open("/home/cory/Research/vikar311/detectors/SmallBars.det", "w")
				for i in range(len(LargeBars)): f.write(VIKARdet(LargeBars[i], offset=target_origin, angle=beam_axis))
				for i in range(len(SmallBars)): f2.write(VIKARdet(SmallBars[i], offset=target_origin, angle=beam_axis))
				print " Wrote Large bar detector file /home/cory/Research/vikar311/detectors/LargeBars.det"
				print " Wrote Small bardetector file /home/cory/Research/vikar311/detectors/SmallBars.det"
				f.close()
				f2.close()
			else: myLab.KeyPress(s_)

# Read xyz coordinate tuple data
def RawRead(fname, pt_size):
	f = open(fname,"r")
	tuples = []
	for line in f:
		line = line.strip()
		arr = line.split("\t")
		tuples.append((float(arr[2]),float(arr[1]),float(arr[0])))
	f.close()
	return visual.points(pos=tuples, size=pt_size, color=visual.color.black)

# Read a standard VIKAR detector file
def DetRead(fname, opacity=0.25):
	f = open(fname,"r")
	lines = f.readlines()
	num_det = 0
	output = []
	for i in range(len(lines)):
		if lines[i][0] == "#": continue
		arr = lines[i].strip().split("\t")
		if len(arr) >= 8:
			if arr[7] == "small": output.append(VandleBar(bar_length=0.6, bar_width=0.03, bar_height=0.03))
			elif arr[7] == "medium": output.append(VandleBar(bar_length=1.2, bar_width=0.05, bar_height=0.03))
			elif arr[7] == "large": output.append(VandleBar(bar_length=2.0, bar_width=0.05, bar_height=0.05))
			else: output.append(VandleBar(bar_length=float(arr[8]), bar_width=float(arr[9]), bar_height=float(arr[10])))
			output[num_det].SetX(float(arr[2]))
			output[num_det].SetY(float(arr[1]))
			output[num_det].SetZ(float(arr[0]))
			output[num_det].SetTheta(-1*float(arr[3]))
			output[num_det].SetPhi(float(arr[4]))
			output[num_det].SetPsi(float(arr[5]))
			output[num_det].body.opacity = opacity
			num_det += 1
	f.close()
	return output

def Hist3D(det_fname, data_fname):
	f = open(det_fname,"r")
	lines = f.readlines()
	output = []
	starts = []
	bins = []
	for i in range(len(lines)):
		arr = lines[i].strip().split("\t")
		if len(arr) >= 7:
			temp_array = []
			if arr[6] == "small": 
				starts.append(float(arr[1])-0.3)
				bins.append(6)
				for j in range(6):
					temp_array.append(VandleBar(bar_length=0.1, bar_width=0.03, bar_height=0.03))
					temp_array[j].SetY(starts[i] + 0.05 + j*0.1)
					temp_array[j].SetX(float(arr[2]))
					temp_array[j].SetZ(float(arr[0]))
					temp_array[j].SetTheta(-1*float(arr[3]))
					temp_array[j].SetPhi(float(arr[4]))
					temp_array[j].SetPsi(float(arr[5]))
			elif arr[6] == "medium":
				starts.append(float(arr[1])-0.6)
				bins.append(12)
				for j in range(12):
					temp_array.append(VandleBar(bar_length=0.1, bar_width=0.05, bar_height=0.03))
					temp_array[j].SetY(starts[i] + 0.05 + j*0.1)
					temp_array[j].SetX(float(arr[2]))
					temp_array[j].SetZ(float(arr[0]))
					temp_array[j].SetTheta(-1*float(arr[3]))
					temp_array[j].SetPhi(float(arr[4]))
					temp_array[j].SetPsi(float(arr[5]))
			elif arr[6] == "large": 
				starts.append(float(arr[1])-1.0)
				bins.append(20)
				for j in range(20):
					temp_array.append(VandleBar(bar_length=0.1, bar_width=0.05, bar_height=0.05))
					temp_array[j].SetY(starts[i] + 0.05 + j*0.1)
					temp_array[j].SetX(float(arr[2]))
					temp_array[j].SetZ(float(arr[0]))
					temp_array[j].SetTheta(-1*float(arr[3]))
					temp_array[j].SetPhi(float(arr[4]))
					temp_array[j].SetPsi(float(arr[5]))
			
			output.append(temp_array)
	
	print " Found",len(output),"detectors in file"	
	f.close()
	
	f = open(data_fname,"r")
	lines = f.readlines()
	data = []
	
	for i in range(len(output)):
		data.append([])
		for j in range(bins[i]):
			data[i].append(0)
	
	for i in range(len(lines)):
		arr = lines[i].strip().split("\t")
		if len(arr) >= 12:
			# bar=7 y=10
			index = (float(arr[10])-starts[int(arr[7])])*10
			data[int(arr[7])][int(index)] += 1
	
	maximum = 0
	for i in range(len(data)):
		for j in range(len(data[i])):
			if data[i][j] > maximum: maximum = data[i][j]
			
	print " Max Bin =",maximum
	f.close()
	
	for i in range(len(data)):
		for j in range(len(data[i])):
			output[i][j].SetColor(SmoothSpectrum(data[i][j], maximum))
	
	return output

def main2():  
	labScene = visual.display(title="7Be(d,n)8B Experiment", width=800, height=600, background=GetRGBcode(153,204,255))
	#axisx = visual.box(pos=(0,0,0), axis=(10.0,0,0), width=0.005, height=0.005, color=visual.color.red, opacity=0.5)
	#axisy = visual.box(pos=(0,0,0), axis=(0,10.0,0), width=0.005, height=0.005, color=visual.color.blue, opacity=0.5)
	#axisz = visual.box(pos=(0,0,0), axis=(0,0,10.0), width=0.005, height=0.005, color=visual.color.green, opacity=0.5)
	#labelx = visual.label(pos=(5.0,0,0), text="Z-Axis")
	#labely = visual.label(pos=(0,5.0,0), text="Y-Axis")
	#labelz = visual.label(pos=(0,0,5.0), text="X-Axis")
	#histogram = Hist3D("/home/cory/Research/vikar311/detectors/7BeExp.det","/home/cory/Research/vikar311/source/rewrite/VIKARoutput.dat")
	#RawRead("/home/cory/Research/VANDLE/Vikar/xyz.dat", 2)
	RawRead("/home/cory/Research/VANDLE/Vikar/dump.out", 2)
	#RawRead("/home/cory/Research/VANDLE/Vikar/test.out", 2)
	VBars = DetRead("/home/cory/Research/VANDLE/Vikar/detectors/Phoswich.det")

	while True:
		visual.rate(60)
		if labScene.kb.keys:
			s_ = labScene.kb.getkey()
		
			if(s_ == "esc"): wx.Exit()

main2() 
#main()
