include(../defaults.pri)
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
TEMPLATE = lib
TARGET = myapp
SOURCES += tridiag.cpp \
    Closed_form.cpp \
    Explicit.cpp \
    Implicit.cpp \
    Mon-Car.cpp \
    Random.cpp
QMAKE_CXXFLAGS += -std=c++0x
QMAKE_CXXFLAGS+= -openmp
QMAKE_LFLAGS +=  -openmp
HEADERS += tridiag.h \
    Closed_form.h \
    Explicit.h \
    Implicit.h \
    Mon-Car.h \
    Random.h

