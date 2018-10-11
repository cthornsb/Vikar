#include <QPainter>
#include <QPixmap>

#include "camera.h"
#include "ui_camera.h"

#include "detectors.hpp"

RGBcolor::RGBcolor(int color_/*=0xFFFFFF*/){
	SetColor(color_);
}

void RGBcolor::SetColor(int color_){
	unsigned char temparr[3];
	memcpy(temparr, &color_, 3);
	b = temparr[0];
	g = temparr[1];
	r = temparr[2];
}

void RGBcolor::GetColor(unsigned char &red, unsigned char &green, unsigned char &blue, float luminosity_/*=1.0*/){
	red = (unsigned char)(r * luminosity_);
	green = (unsigned char)(g * luminosity_);
	blue = (unsigned char)(b * luminosity_);
}

void Camera::set_pixel_size(){
    pixelX = 2.0*std::tan(fov/(2.0*scaling));
    pixelY = 2.0*std::tan(fov/(2.0*scaling));
    center = pos + dir;
    origin = center + screenX*(-pixelX/2.0) + screenY*(-pixelY/2.0);
    pixelX = pixelX / sizeX;
    pixelY = pixelY / sizeY;
}

void Camera::set_rotation(){
    dir = Vector3(0.0, 0.0, 1.0);
    screenX = Vector3(1, 0, 0);
    screenY = Vector3(0, -1, 0);

    rot.SetRotationMatrixSphere(theta, phi, psi);
    rot.Transform(dir);
    rot.Transform(screenX);
    rot.Transform(screenY);
    rotated = false;
}

void Camera::set_colors(){
	int rgbcodes[12] = {0xFF0000, 0xFF7F00, 0xFFFF00, 0x7FFF00,
	                    0x00FF00, 0x00FF7F, 0x00FFFF, 0x007FFF,
	                    0x0000FF, 0x7F00FF, 0xFF00FF, 0xFF007F};
	
	for(int i = 0; i < 12; i++){
		colors[i].SetColor(rgbcodes[i]);    
	}
}

void Camera::update_detSelect(){
    if(!currDetector) return;
    ui->doubleSpinBox_xPos->setValue(currDetector->GetX());
    ui->doubleSpinBox_yPos->setValue(currDetector->GetY());
    ui->doubleSpinBox_zPos->setValue(currDetector->GetZ());
    ui->doubleSpinBox_theta->setValue(currDetector->GetTheta()*rad2deg);
    ui->doubleSpinBox_phi->setValue(currDetector->GetPhi()*rad2deg);
    ui->doubleSpinBox_psi->setValue(currDetector->GetPsi()*rad2deg);
    ui->doubleSpinBox_length->setValue(currDetector->GetLength());
    ui->doubleSpinBox_width->setValue(currDetector->GetWidth());
    ui->doubleSpinBox_depth->setValue(currDetector->GetDepth());
    ui->lineEdit_type->setText(QString::fromStdString(currDetector->GetType()));
    ui->lineEdit_subtype->setText(QString::fromStdString(currDetector->GetSubtype()));
    ui->lineEdit_material->setText(QString::fromStdString(currDetector->GetMaterialName()));
}

void Camera::disable_all(){
    ui->pushButton->setDisabled(true);
    ui->pushButton_2->setDisabled(true);
    ui->doubleSpinBox_xPos->setDisabled(true);
    ui->doubleSpinBox_yPos->setDisabled(true);
    ui->doubleSpinBox_zPos->setDisabled(true);
    ui->doubleSpinBox_theta->setDisabled(true);
    ui->doubleSpinBox_phi->setDisabled(true);
    ui->doubleSpinBox_psi->setDisabled(true);
    ui->doubleSpinBox_length->setDisabled(true);
    ui->doubleSpinBox_width->setDisabled(true);
    ui->doubleSpinBox_depth->setDisabled(true);
    ui->spinBox_detSelect->setDisabled(true);
    ui->lineEdit_type->setDisabled(true);
    ui->lineEdit_subtype->setDisabled(true);
    ui->lineEdit_material->setDisabled(true);
}

void Camera::enable_all(){
    ui->pushButton->setEnabled(true);
    ui->pushButton_2->setEnabled(true);
    ui->doubleSpinBox_xPos->setEnabled(true);
    ui->doubleSpinBox_yPos->setEnabled(true);
    ui->doubleSpinBox_zPos->setEnabled(true);
    ui->doubleSpinBox_theta->setEnabled(true);
    ui->doubleSpinBox_phi->setEnabled(true);
    ui->doubleSpinBox_psi->setEnabled(true);
    ui->doubleSpinBox_length->setEnabled(true);
    ui->doubleSpinBox_width->setEnabled(true);
    ui->doubleSpinBox_depth->setEnabled(true);
    ui->spinBox_detSelect->setEnabled(true);
    ui->lineEdit_type->setEnabled(true);
    ui->lineEdit_subtype->setEnabled(true);
    ui->lineEdit_material->setEnabled(true);
}

Camera::Camera(QWidget* parent): QMainWindow(parent), ui(new Ui::Camera) {
    fov = pi/2.0;
    sizeX = 240;
    sizeY = 240;
    scaling = 1;

    SetRotation(0.0, 0.0, 0.0);

    ui->setupUi(this);

    scene = new QGraphicsScene(this);
    ui->graphicsView->setScene(scene);

    pixmap = new QPixmap(sizeX, sizeY);
    pixmap->fill(Qt::black);

    scene->addPixmap(*pixmap);
    the_app = NULL;
    
    currDetector = NULL;

    set_colors();
}

Camera::Camera(double x_, double y_, double z_, double theta_, double phi_, QWidget* parent): QMainWindow(parent), ui(new Ui::Camera) {
    fov = pi/2.0;
    sizeX = 240;
    sizeY = 240;
    scaling = 1;

    pos = Vector3(x_, y_, z_);
    SetRotation(theta_, phi_, 0.0);

    ui->setupUi(this);

    scene = new QGraphicsScene(this);
    ui->graphicsView->setScene(scene);

    pixmap = new QPixmap(sizeX, sizeY);
    pixmap->fill(Qt::black);

    scene->addPixmap(*pixmap);
    the_app = NULL;
    
    disable_all();

    set_colors();
}

Camera::~Camera(){
    for(std::vector<Primitive*>::iterator iter = primitives.begin(); iter != primitives.end(); iter++){
        delete (*iter);
    }
    primitives.clear();
    delete ui;
    delete scene;
    delete pixmap;
}

void Camera::SetPosition(double x_, double y_, double z_){
    pos = Vector3(x_, y_, z_);
}

void Camera::SetRotation(double theta_, double phi_, double psi_){
    theta = theta_;
    phi = phi_;
    psi = psi_;
    rotated = true;
}

void Camera::PointCamera(const Vector3 &g_){
    dir = g_ - pos;
    dir.Normalize();
}

void Camera::SetFOV(double fov_){
    fov = fov_*deg2rad;
}

void Camera::SetScaling(int scaling_){
    scaling = scaling_;
}

void Camera::SetSizeX(int sizeX_){
    sizeX = sizeX_;
}

void Camera::SetSizeY(int sizeY_){
    sizeY = sizeY_;
}

void Camera::AddPrimitive(Primitive *prim_){
    primitives.push_back(prim_);
}

void Camera::Render() {
    if(rotated){ set_rotation(); }
    set_pixel_size();

    Vector3 p1;

    pixmap->fill(Qt::black);
    QPainter pen(pixmap);
    QColor pen_color;

    Vector3 ray;
    Vector3 normal;

    Vector3 face_hit1;
    Vector3 face_hit2;

    unsigned char red, green, blue;

    float luminosity; // Intensity of the color of the drawing pen.
    int index; // Current detector index in the Primitive vector.
    int detector = 0; // Index of visible detector in the Primitive vector.
    double currentX; // Current pixel along the x-axis.
    double currentY; // Current pixel along the y-axis.
    double depth; // Distance of the "pixel" from the viewer.
    double t1; // "Distance" parameter such that P1 = position + ray*t1.
    double t2; // Not used.
    for(int i = 0; i < sizeY; i++){
    	ui->progressBar->setValue((float(i)*100/sizeY));
        currentY = (pixelY/2.0) + i*pixelY;
        for(int j = 0; j < sizeX; j++){ // Loop over all rows of pixels.
            currentX = (pixelX/2.0) + j*pixelX;
            ray = origin + screenX*currentX + screenY*currentY - pos;
            ray.Normalize(); // Normalize the direction vector.
            
            // Set the initial conditions.
            depth = 9999;
            luminosity = 0;
            
            // Loop over all pixels in this row.
           index = 0;
            for(std::vector<Primitive*>::iterator iter = primitives.begin(); iter != primitives.end(); iter++){
                if((*iter)->IntersectPrimitive(pos, ray, p1, normal, t1, t2) && t1 < depth){
			        luminosity = fabs(ray.CosAngle(normal));
			        if(luminosity < 0.0 || luminosity > 1.0){ continue; }
			        detector = index;
                	depth = t1;
                }
                index++;
            }
            
            // Only draw the pixel if we found a valid intersection.
            if(depth != 9999){ 
        		while(detector >= 12){ // Wrap the detector id around so we don't go outside the bounds of the color array.
        			detector -= 12;
        		}
            	colors[detector].GetColor(red, green, blue, luminosity);
                pen_color.setRgb(red, green, blue);
                pen.setPen(pen_color);
                pen.drawPoint(j, i);
            }
        }
    }

    ui->progressBar->setValue(100);

    scene->addPixmap(*pixmap);
}

void Camera::paintEvent(QPaintEvent*) {
}

void Camera::on_pushButton_clicked(){
    this->Render();
}

void Camera::on_spinBox_detSelect_valueChanged(int i){
    currDetector = ((unsigned int)i < primitives.size() ? primitives.at(i) : NULL);
    update_detSelect();
}

void Camera::on_doubleSpinBox_valueChanged(double arg1){
    pos.axis[0] = arg1;
}

void Camera::on_doubleSpinBox_2_valueChanged(double arg1){
    pos.axis[1] = arg1;
}

void Camera::on_doubleSpinBox_3_valueChanged(double arg1){
    pos.axis[2] = arg1;
}

void Camera::on_doubleSpinBox_5_valueChanged(double arg1){
    theta = -arg1*deg2rad;
    rotated = true;
}

void Camera::on_doubleSpinBox_6_valueChanged(double arg1){
    phi = -arg1*deg2rad;
    rotated = true;
}

void Camera::on_doubleSpinBox_7_valueChanged(double arg1){
    psi = arg1*deg2rad;
    rotated = true;
}

void Camera::on_doubleSpinBox_xPos_valueChanged(double arg1){
    currDetector->SetPosition(arg1, currDetector->GetY(), currDetector->GetZ());
}

void Camera::on_doubleSpinBox_yPos_valueChanged(double arg1){
    currDetector->SetPosition(currDetector->GetX(), arg1, currDetector->GetZ());
}

void Camera::on_doubleSpinBox_zPos_valueChanged(double arg1){
    currDetector->SetPosition(currDetector->GetX(), currDetector->GetY(), arg1);
}

void Camera::on_doubleSpinBox_theta_valueChanged(double arg1){
    currDetector->SetRotation(arg1*deg2rad, currDetector->GetPhi(), currDetector->GetPsi());
}

void Camera::on_doubleSpinBox_phi_valueChanged(double arg1){
    currDetector->SetRotation(currDetector->GetTheta(), arg1*deg2rad, currDetector->GetPsi());
}

void Camera::on_doubleSpinBox_psi_valueChanged(double arg1){
    currDetector->SetRotation(currDetector->GetTheta(), currDetector->GetPhi(), arg1*deg2rad);
}

void Camera::on_doubleSpinBox_length_valueChanged(double arg1){
    currDetector->SetSize(arg1, currDetector->GetWidth(), currDetector->GetDepth());
}

void Camera::on_doubleSpinBox_width_valueChanged(double arg1){
    currDetector->SetSize(currDetector->GetLength(), arg1, currDetector->GetDepth());
}

void Camera::on_doubleSpinBox_depth_valueChanged(double arg1){
    currDetector->SetSize(currDetector->GetLength(), currDetector->GetWidth(), arg1);
}

void Camera::on_spinBox_valueChanged(int arg1){
    fov = arg1*deg2rad;
}

void Camera::on_spinBox_2_valueChanged(int arg1){
    scaling = arg1;
}

void Camera::on_actionScreenshot_triggered(){
    pixmap->save(ui->lineEdit->text());
    std::cout << "Screenshot saved to " << ui->lineEdit->text().toStdString() << std::endl;
}

void Camera::on_pushButton_2_clicked(){
    pixmap->save(ui->lineEdit->text());
    std::cout << "Screenshot saved to " << ui->lineEdit->text().toStdString() << std::endl;
}

void Camera::on_pushButton_3_clicked(){
    for(std::vector<Primitive*>::iterator iter = primitives.begin(); iter != primitives.end(); iter++){
        delete (*iter);
    }
    primitives.clear();
    ReadDetFile(ui->lineEdit_2->text().toStdString().c_str(), primitives);
    std::cout << "Loaded detector file " << ui->lineEdit_2->text().toStdString() << std::endl;
    std::cout << " Found " << primitives.size() << " objects.\n";

    if(!primitives.empty()){ // Set position diagnostic tools to the position of the first detector read.
        int count = 1;
        for(std::vector<Primitive*>::iterator iter = primitives.begin(); iter != primitives.end(); iter++){
            std::stringstream stream;
            stream << count++;
        }
	ui->spinBox_detSelect->setMaximum(count-1);
	ui->lcdNumber->display(count-1);
        currDetector = primitives.front();
	enable_all();
        update_detSelect();
    }
    else{
        ui->lcdNumber->display(0);
	disable_all();
        currDetector = NULL;
    }
}

void Camera::on_actionExit_triggered(){
    if(the_app){ the_app->quit(); }
}
