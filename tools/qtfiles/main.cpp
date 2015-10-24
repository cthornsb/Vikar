#include <QApplication>

#include "camera.h"
#include "ui_camera.h"
#include "vikar_core.h"
#include "detectors.h"

int main(int argc, char* argv[]){
    QApplication app(argc, argv);

    Camera cam(0.0, 0.0, -1.0);
    cam.SetApp(&app);
    cam.show();

	if(argc >= 2){    
    	cam.GetUI()->lineEdit_2->setText(QString(argv[1]));
    	cam.LoadDetectors();
	}    

    return app.exec();
}
