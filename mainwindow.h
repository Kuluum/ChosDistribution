#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot/qcustomplot.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void setupQuadraticDemo(QCustomPlot *customPlot);

private slots:
    void on_drawButton_clicked();

    void spinboxValueChanged();

private:
    Ui::MainWindow *ui;
    void setupSlotConnection();
    QCPBars* createBar(QCustomPlot *plotParent, double x, double height, double width);
    void fit();
};

#endif // MAINWINDOW_H
