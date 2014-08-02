import random
import math

f = open("pizza.out","w")

for i in range(10000):
	zeta = random.random()
	#time = (1.0/math.log(2.0))*math.log(1.0/(1.0-zeta))
	time = math.sqrt(zeta)
	f.write(str(zeta)+"\t"+str(time)+"\n")
	
f.close()
	
