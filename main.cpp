#include "mainwindow.h"
#include "embededtest.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    //MainWindow w;
    //w.show();
    embededTest e;
    e.show();
    return a.exec();
}
