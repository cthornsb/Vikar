#include <QApplication>

#include "camera.h"
#include "vikar_core.h"
#include "detectors.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    Camera cam(0.0, 0.0, -1.0);
    cam.SetApp(&app);
    cam.show();

    ReadDetFile("/home/cory/Research/Vikar/detectors/temp.det", *cam.GetPrimitives());
    std::cout << "Found " << cam.GetPrimitives()->size() << " objects.\n";

    return app.exec();
}
