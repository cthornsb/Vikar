#include <QApplication>

#include "camera.h"
#include "vikar_core.h"
#include "detectors.h"

void help(char * prog_name_){
	std::cout << "  SYNTAX: " << prog_name_ << " [filename] <options>\n";
	std::cout << "   Available options:\n";
}

int main(int argc, char* argv[]){
	if(argc < 2){
		std::cout << " Error: Invalid number of arguments to " << argv[0] << ". Expected 1, received " << argc-1 << ".\n";
		help(argv[0]);
		return 1;
	}
	
    QApplication app(argc, argv);

    Camera cam(0.0, 0.0, -1.0);
    cam.SetApp(&app);
    cam.show();

    ReadDetFile(argv[1], *cam.GetPrimitives());
    std::cout << "Found " << cam.GetPrimitives()->size() << " objects.\n";

    return app.exec();
}
