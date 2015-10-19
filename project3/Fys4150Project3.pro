QT += core
QT -= gui

TARGET = Fys4150Project2
CONFIG += console
CONFIG -= app_bundle

INCLUDEPATH += C:\Qt\include
# LIBS += -L C:\Qt
LIBS += -lliblapack -llibblas

TEMPLATE = app

SOURCES += main.cpp
HEADERS +=
