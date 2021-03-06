#-------------------------------------------------
#
# Project created by QtCreator 2016-11-04T16:14:30
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = ChosDistribution
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    chosdistribution.cpp \
    qcustomplot/qcustomplot.cpp \
    embededtest.cpp \
    FileProcessing/filereader.cpp \
    DataModel/distributiondata.cpp \
    Algorithm/algorithms.cpp \
    descentprogressview.cpp

HEADERS  += mainwindow.h \
    chosdistribution.h \
    qcustomplot/qcustomplot.h \
    embededtest.h \
    FileProcessing/filereader.h \
    DataModel/distributiondata.h \
    Algorithm/algorithms.h \
    descentprogressview.h

FORMS    += mainwindow.ui \
    embededtest.ui \
    descentprogressview.ui

CONFIG += c++14

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
