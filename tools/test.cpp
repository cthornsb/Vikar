#include <iostream>

#include "vikar_core.h"

int main(){
	/*Matrix3 pizza;
	pizza.SetRow1(1, 0, 0);
	pizza.SetRow2(0, 0, 1);
	pizza.SetRow3(0, 1, 0);
	pizza.Dump();*/
	
	Vector3 taco(0, 1, 0);
	Vector3 burrito(1, 0, 0);
	
	taco.TransformSpherical(0, 0, 0);
	burrito.TransformSpherical(3.14159/2.0, 3.14159, 0);
	
	std::cout << taco.Dump() << std::endl;
	std::cout << burrito.Dump() << std::endl;

	return 0;
}
