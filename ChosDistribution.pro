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
    alglib/interpolation.cpp \
    alglib/integration.cpp \
    alglib/fasttransforms.cpp \
    alglib/ap.cpp \
    alglib/statistics.cpp \
    alglib/specialfunctions.cpp \
    alglib/solvers.cpp \
    alglib/optimization.cpp \
    alglib/linalg.cpp \
    alglib/diffequations.cpp \
    alglib/dataanalysis.cpp \
    alglib/alglibmisc.cpp \
    alglib/alglibinternal.cpp

HEADERS  += mainwindow.h \
    chosdistribution.h \
    qcustomplot/qcustomplot.h \
    alglib/stdafx.h \
    alglib/interpolation.h \
    alglib/integration.h \
    alglib/fasttransforms.h \
    alglib/ap.h \
    alglib/statistics.h \
    alglib/specialfunctions.h \
    alglib/solvers.h \
    alglib/optimization.h \
    alglib/linalg.h \
    alglib/diffequations.h \
    alglib/dataanalysis.h \
    alglib/alglibmisc.h \
    alglib/alglibinternal.h

FORMS    += mainwindow.ui

CONFIG += c++14
