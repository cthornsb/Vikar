import pygame
import math
from pygame import display, draw, event, mouse, Surface
from pygame.locals import *

# Return the cartesian coordinates of a point in spherical space as the 3-tuple (x, y, z)
def GetCartesianCoords(r, theta, phi):
	return(r*math.sin(theta)*math.cos(phi), r*math.sin(theta)*math.sin(phi), r*math.cos(theta))

class vector2:
	def __init__(self, x, y):
		self.x = x
		self.y = y
		
	def __mul__(self, scalar):
		return vector2(self.x*scalar, self.y*scalar)
		
	def __neg__(self):
		return vector2(-self.x, -self.y)
		
	def __add__(self, other):
		return vector2(self.x+other.x, self.y+other.y)
		
	def __sub__(self, other):
		return vector2(self.x-other.x, self.y-other.y)
		
	def __print__(self):
		print self.x, self.y
		
	def get_dist(self, other):
		return math.sqrt((self.x-other.x)**2+(self.y-other.y)**2)
		
	def get_tuple(self):
		return (self.x, self.y)
		
	def dump(self):
		print self.x, self.y

class line:
	def __init__(self, p1, p2):
		self.p1 = p1
		self.p2 = p2
		self.length = self.p1.get_dist(self.p2)
		self.center = self.p1+((self.p2-self.p1)*(self.p1.get_dist(self.p2)/2.0))
		
	def rotate(self, theta_):
		unit_x = vector2(math.cos(theta_), math.sin(theta_))
		unit_y = vector2(-math.sin(theta_), math.cos(theta_))
		
	def normalize(self, win_width, win_height, Bars_):
		maxx = Bars_.maximumx
		maxy = Bars_.maximumy
		minx = Bars_.minimumx
		miny = Bars_.minimumy
		self.p1.x *= (maxx-minx)/win_width
		self.p2.x *= (maxx-minx)/win_width
		self.p1.y *= (maxy-miny)/win_height
		self.p2.y *= (maxy-miny)/win_height
		
	def draw(self, screen_):
		pygame.draw.polygon(screen_, (0,0,0), [self.p1.get_tuple(),self.p2.get_tuple()], 1)

class box:
	def __init__(self, center, width, height):
		self.center = center
		self.width = width
		self.height = height
		self.unit_x = vector2(1,0)
		self.unit_y = vector2(0,1)
		self._calculate()
		
	def _calculate(self):
		self.p1 = self.center + (self.unit_x*(self.width/2.0) + self.unit_y*(self.height/2.0))
		self.p2 = self.center + (self.unit_x*(self.width/2.0) - self.unit_y*(self.height/2.0))
		self.p3 = self.center + (-self.unit_x*(self.width/2.0) - self.unit_y*(self.height/2.0))
		self.p4 = self.center + (-self.unit_x*(self.width/2.0) + self.unit_y*(self.height/2.0))
		
	def rotate(self, theta_):
		self.unit_x = vector2(math.cos(theta_), math.sin(theta_))
		self.unit_y = vector2(-math.sin(theta_), math.cos(theta_))
		self._calculate()
		
	def normalize(self, win_width, win_height, minx, maxx, miny, maxy):
		self.width *= (win_width/(maxx-minx))
		self.height *= (win_height/(maxy-miny))
		self.center.x = (self.center.x - minx)*(win_width/(maxx-minx))
		self.center.y = (self.center.y - miny)*(win_height/(maxy-miny))
		self._calculate()
		
	def draw(self, screen_):
		pygame.draw.polygon(screen_, (0,0,0), [self.p1.get_tuple(),self.p2.get_tuple(),self.p3.get_tuple(),self.p4.get_tuple()], 1)

class boxArray:
	def __init__(self):
		self.array = []
		self.maximumx = -1E22
		self.maximumy = -1E22
		self.minimumx = 1E22
		self.minimumy = 1E22
		
	def append(self, inBox):
		self.array.append(inBox)
		
	def normalize(self, margins=0.1):
		for i in range(len(self.array)):
			if float(self.array[i].p1.x) > self.maximumx: self.maximumx = float(self.array[i].p1.x)
			elif float(self.array[i].p1.x) < self.minimumx: self.minimumx = float(self.array[i].p1.x)
			if float(self.array[i].p1.y) > self.maximumy: self.maximumy = float(self.array[i].p1.y)
			elif float(self.array[i].p1.y) < self.minimumy: self.minimumy = float(self.array[i].p1.y)
		
			if float(self.array[i].p2.x) > self.maximumx: self.maximumx = float(self.array[i].p2.x)
			elif float(self.array[i].p2.x) < self.minimumx: self.minimumx = float(self.array[i].p2.x)
			if float(self.array[i].p2.y) > self.maximumy: self.maximumy = float(self.array[i].p2.y)
			elif float(self.array[i].p2.y) < self.minimumy: self.minimumy = float(self.array[i].p2.y)
		
			if float(self.array[i].p3.x) > self.maximumx: self.maximumx = float(self.array[i].p3.x)
			elif float(self.array[i].p3.x) < self.minimumx: self.minimumx = float(self.array[i].p3.x)
			if float(self.array[i].p3.y) > self.maximumy: self.maximumy = float(self.array[i].p3.y)
			elif float(self.array[i].p3.y) < self.minimumy: self.minimumy = float(self.array[i].p3.y)
		
			if float(self.array[i].p4.x) > self.maximumx: self.maximumx = float(self.array[i].p4.x)
			elif float(self.array[i].p4.x) < self.minimumx: self.minimumx = float(self.array[i].p4.x)
			if float(self.array[i].p4.y) > self.maximumy: self.maximumy = float(self.array[i].p4.y)
			elif float(self.array[i].p4.y) < self.minimumy: self.minimumy = float(self.array[i].p4.y)
		
		for i in range(len(self.array)):
			self.array[i].normalize(720, 720, self.minimumx*(1+margins), self.maximumx*(1+margins), self.minimumy*(1+margins), self.maximumy*(1+margins))

# Read a standard VIKAR detector file
def DetRead(fname):
	f = open(fname,"r")
	lines = f.readlines()
	output = boxArray()
	for i in range(len(lines)):
		arr = lines[i].strip().split("\t")
		if len(arr) >= 7:
			if arr[6] == "small": output.append(box(vector2(float(arr[2]),float(arr[0])), 0.03, 0.03))
			elif arr[6] == "medium": output.append(box(vector2(float(arr[2]),float(arr[0])), 0.03, 0.05))
			elif arr[6] == "large": output.append(box(vector2(float(arr[2]),float(arr[0])), 0.05, 0.05))
			else: output.append(box(vector2(float(arr[2]),float(arr[0])), float(arr[7]), float(arr[8])))
			output.array[i].rotate(float(arr[3]))	
	
	output.normalize()
		
	f.close()
	return output

def main():		
	pygame.init()
	screen = pygame.display.set_mode((720,720))
	pygame.display.set_caption('Schematic')
	clock = pygame.time.Clock()
	
	background = pygame.Surface(screen.get_size())
	background = background.convert()
	background.fill((255, 255, 255))

	font = pygame.font.Font(None, 36)
	screen.blit(background, (0, 0))
	pygame.display.flip()
	myline = line(vector2(0,380), vector2(360,380))
	#VBars = DetRead("/home/cory/Research/VANDLE/vikar2/detectors/7BeExp.det")
	VBars = DetRead("/home/cory/Research/VANDLE/vikar2/detectors/MichArray.det")
	#myline.normalize(720, 720, VBars)

	while 1:
		clock.tick(10)
		myline.draw(screen)
		for i in range(len(VBars.array)): VBars.array[i].draw(screen)
		for event in pygame.event.get():
			if event.type == QUIT: return
			
		display.update()

main()
"""temp = [268, 284, 301, 316, 331, 353]
for i in range(len(temp)):
	angle = temp[i]*math.pi/180.0 #rad
	cart = GetCartesianCoords(0.5, angle, 0.0)
	string = str(cart[0]) + "\t" + str(cart[1]) + "\t" + str(cart[2]) + "\t" + str(angle-math.pi/2.0) + "\t"
	string += "0.0\t0.0\t0.0\t0.1016\t0.1524"
	print string"""
