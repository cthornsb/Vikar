#-------------------------------------------------
#
# Project created by QtCreator 2015-10-18T17:58:37
#
#-------------------------------------------------

QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = renderer
TEMPLATE = app

OBJECTS_DIR = ./obj
INCLUDEPATH += ../../include
DESTDIR = ../
UI_DIR = ./ui
MOC_DIR = ./moc

HEADERS += camera.h

SOURCES += main.cpp camera.cpp

OBJECTS += ../../obj/geometry.o \
           ../../obj/detectors.o \
           ../../obj/vikar_core.o

FORMS += camera.ui
