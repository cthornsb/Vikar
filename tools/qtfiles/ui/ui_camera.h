/********************************************************************************
** Form generated from reading UI file 'camera.ui'
**
** Created by: Qt User Interface Compiler version 5.0.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CAMERA_H
#define UI_CAMERA_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Camera
{
public:
    QAction *actionScreenshot;
    QAction *actionExit;
    QWidget *centralWidget;
    QGraphicsView *graphicsView;
    QPushButton *pushButton;
    QDoubleSpinBox *doubleSpinBox;
    QDoubleSpinBox *doubleSpinBox_2;
    QDoubleSpinBox *doubleSpinBox_3;
    QSpinBox *spinBox;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_5;
    QDoubleSpinBox *doubleSpinBox_5;
    QDoubleSpinBox *doubleSpinBox_6;
    QLabel *label_6;
    QLabel *label_7;
    QLabel *label_8;
    QDoubleSpinBox *doubleSpinBox_7;
    QSpinBox *spinBox_2;
    QLineEdit *lineEdit;
    QPushButton *pushButton_2;
    QLabel *label_9;
    QProgressBar *progressBar;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *Camera)
    {
        if (Camera->objectName().isEmpty())
            Camera->setObjectName(QStringLiteral("Camera"));
        Camera->resize(537, 363);
        actionScreenshot = new QAction(Camera);
        actionScreenshot->setObjectName(QStringLiteral("actionScreenshot"));
        actionExit = new QAction(Camera);
        actionExit->setObjectName(QStringLiteral("actionExit"));
        centralWidget = new QWidget(Camera);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        graphicsView = new QGraphicsView(centralWidget);
        graphicsView->setObjectName(QStringLiteral("graphicsView"));
        graphicsView->setGeometry(QRect(10, 10, 246, 246));
        pushButton = new QPushButton(centralWidget);
        pushButton->setObjectName(QStringLiteral("pushButton"));
        pushButton->setGeometry(QRect(80, 260, 87, 27));
        doubleSpinBox = new QDoubleSpinBox(centralWidget);
        doubleSpinBox->setObjectName(QStringLiteral("doubleSpinBox"));
        doubleSpinBox->setGeometry(QRect(270, 30, 62, 25));
        doubleSpinBox->setMinimum(-50);
        doubleSpinBox->setMaximum(50);
        doubleSpinBox->setSingleStep(0.1);
        doubleSpinBox_2 = new QDoubleSpinBox(centralWidget);
        doubleSpinBox_2->setObjectName(QStringLiteral("doubleSpinBox_2"));
        doubleSpinBox_2->setGeometry(QRect(350, 30, 62, 25));
        doubleSpinBox_2->setMinimum(-50);
        doubleSpinBox_2->setMaximum(50);
        doubleSpinBox_2->setSingleStep(0.1);
        doubleSpinBox_3 = new QDoubleSpinBox(centralWidget);
        doubleSpinBox_3->setObjectName(QStringLiteral("doubleSpinBox_3"));
        doubleSpinBox_3->setGeometry(QRect(430, 30, 62, 25));
        doubleSpinBox_3->setMinimum(-50);
        doubleSpinBox_3->setMaximum(50);
        doubleSpinBox_3->setSingleStep(0.1);
        doubleSpinBox_3->setValue(-1);
        spinBox = new QSpinBox(centralWidget);
        spinBox->setObjectName(QStringLiteral("spinBox"));
        spinBox->setGeometry(QRect(310, 140, 52, 25));
        spinBox->setMaximum(180);
        spinBox->setSingleStep(5);
        spinBox->setValue(90);
        label = new QLabel(centralWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(300, 120, 81, 16));
        label_2 = new QLabel(centralWidget);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(410, 120, 41, 16));
        label_3 = new QLabel(centralWidget);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(270, 10, 71, 16));
        label_4 = new QLabel(centralWidget);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(350, 10, 71, 16));
        label_5 = new QLabel(centralWidget);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(430, 10, 71, 16));
        doubleSpinBox_5 = new QDoubleSpinBox(centralWidget);
        doubleSpinBox_5->setObjectName(QStringLiteral("doubleSpinBox_5"));
        doubleSpinBox_5->setGeometry(QRect(430, 90, 62, 25));
        doubleSpinBox_5->setMinimum(-180);
        doubleSpinBox_5->setMaximum(180);
        doubleSpinBox_5->setSingleStep(5);
        doubleSpinBox_6 = new QDoubleSpinBox(centralWidget);
        doubleSpinBox_6->setObjectName(QStringLiteral("doubleSpinBox_6"));
        doubleSpinBox_6->setGeometry(QRect(350, 90, 62, 25));
        doubleSpinBox_6->setMinimum(-180);
        doubleSpinBox_6->setMaximum(180);
        doubleSpinBox_6->setSingleStep(5);
        label_6 = new QLabel(centralWidget);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setGeometry(QRect(450, 70, 31, 16));
        label_7 = new QLabel(centralWidget);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setGeometry(QRect(370, 70, 41, 16));
        label_8 = new QLabel(centralWidget);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(280, 70, 41, 16));
        doubleSpinBox_7 = new QDoubleSpinBox(centralWidget);
        doubleSpinBox_7->setObjectName(QStringLiteral("doubleSpinBox_7"));
        doubleSpinBox_7->setGeometry(QRect(270, 90, 62, 25));
        doubleSpinBox_7->setMinimum(-180);
        doubleSpinBox_7->setMaximum(180);
        doubleSpinBox_7->setSingleStep(5);
        spinBox_2 = new QSpinBox(centralWidget);
        spinBox_2->setObjectName(QStringLiteral("spinBox_2"));
        spinBox_2->setGeometry(QRect(400, 140, 52, 25));
        spinBox_2->setMinimum(1);
        spinBox_2->setMaximum(20);
        spinBox_2->setSingleStep(1);
        spinBox_2->setValue(1);
        lineEdit = new QLineEdit(centralWidget);
        lineEdit->setObjectName(QStringLiteral("lineEdit"));
        lineEdit->setGeometry(QRect(270, 200, 221, 25));
        lineEdit->setDragEnabled(true);
        pushButton_2 = new QPushButton(centralWidget);
        pushButton_2->setObjectName(QStringLiteral("pushButton_2"));
        pushButton_2->setGeometry(QRect(340, 230, 87, 27));
        label_9 = new QLabel(centralWidget);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setGeometry(QRect(320, 180, 131, 16));
        progressBar = new QProgressBar(centralWidget);
        progressBar->setObjectName(QStringLiteral("progressBar"));
        progressBar->setGeometry(QRect(317, 270, 131, 23));
        progressBar->setValue(0);
        Camera->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(Camera);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 537, 25));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        Camera->setMenuBar(menuBar);
        mainToolBar = new QToolBar(Camera);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        Camera->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(Camera);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        Camera->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuFile->addAction(actionScreenshot);
        menuFile->addSeparator();
        menuFile->addAction(actionExit);

        retranslateUi(Camera);

        QMetaObject::connectSlotsByName(Camera);
    } // setupUi

    void retranslateUi(QMainWindow *Camera)
    {
        Camera->setWindowTitle(QApplication::translate("Camera", "Camera", 0));
        actionScreenshot->setText(QApplication::translate("Camera", "Screenshot", 0));
        actionExit->setText(QApplication::translate("Camera", "Exit", 0));
        pushButton->setText(QApplication::translate("Camera", "Render", 0));
        label->setText(QApplication::translate("Camera", "Field of View", 0));
        label_2->setText(QApplication::translate("Camera", "Zoom", 0));
        label_3->setText(QApplication::translate("Camera", "X-Position", 0));
        label_4->setText(QApplication::translate("Camera", "Y-Position", 0));
        label_5->setText(QApplication::translate("Camera", "Z-Position", 0));
        label_6->setText(QApplication::translate("Camera", "Yaw", 0));
        label_7->setText(QApplication::translate("Camera", "Roll", 0));
        label_8->setText(QApplication::translate("Camera", "Pitch", 0));
        lineEdit->setText(QApplication::translate("Camera", "./renderer.png", 0));
        pushButton_2->setText(QApplication::translate("Camera", "Snapshot", 0));
        label_9->setText(QApplication::translate("Camera", "Snapshot FIlename", 0));
        menuFile->setTitle(QApplication::translate("Camera", "File", 0));
    } // retranslateUi

};

namespace Ui {
    class Camera: public Ui_Camera {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CAMERA_H
