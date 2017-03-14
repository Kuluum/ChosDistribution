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
    alglib/alglibinternal.cpp \
    alglib/alglibmisc.cpp \
    alglib/ap.cpp \
    alglib/dataanalysis.cpp \
    alglib/diffequations.cpp \
    alglib/fasttransforms.cpp \
    alglib/integration.cpp \
    alglib/interpolation.cpp \
    alglib/linalg.cpp \
    alglib/optimization.cpp \
    alglib/solvers.cpp \
    alglib/specialfunctions.cpp \
    alglib/statistics.cpp \
    Algorithm/algorithms.cpp

HEADERS  += mainwindow.h \
    chosdistribution.h \
    qcustomplot/qcustomplot.h \
    embededtest.h \
    FileProcessing/filereader.h \
    DataModel/distributiondata.h \
    alglib/alglibinternal.h \
    alglib/alglibmisc.h \
    alglib/ap.h \
    alglib/dataanalysis.h \
    alglib/diffequations.h \
    alglib/fasttransforms.h \
    alglib/integration.h \
    alglib/interpolation.h \
    alglib/linalg.h \
    alglib/optimization.h \
    alglib/solvers.h \
    alglib/specialfunctions.h \
    alglib/statistics.h \
    alglib/stdafx.h \
    Algorithm/algorithms.h

FORMS    += mainwindow.ui \
    embededtest.ui

CONFIG += c++14

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
