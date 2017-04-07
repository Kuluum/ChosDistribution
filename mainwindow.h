#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot/qcustomplot.h"
#include "DataModel/distributiondata.h"
#include <vector>

class DistributionData;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void drawChos(QCustomPlot *customPlot);

private slots:
    void on_drawButton_clicked();

    void spinboxValueChanged();

    void on_actionOpen_triggered();

    void on_fitButton_clicked();

    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;
    DistributionData *data;
    std::vector<double> descentProgress;

    void setupSlotConnection();
    QCPBars* createBar(QCustomPlot *plotParent, double x, double height, double width);

    double getChosValue(double x);
    DisVector sqrDiffVector();
    void drawDiffChos(QCustomPlot *customPlot);

    double getChosValue(double x, double x1, double x2,
                                    double m1, double a1, double beta1, double ny1,
                                    double m2, double a2, double beta2, double ny2,
                                    double m3, double a3, double beta3, double ny3);

    double getChosValueWithDistrParams(double x, double mean, double sig, double as, double ex);

};

#endif // MAINWINDOW_H
