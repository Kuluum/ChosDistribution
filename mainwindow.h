#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot/qcustomplot.h"
#include "DataModel/distributiondata.h"
#include "chosdistribution.h"
#include <vector>
#include <stdlib.h>

class DistributionData;

namespace Ui {
class MainWindow;
}

struct fitResults
{
    double rss = 0;
    int quantil1 = 0;
    int quantil2 = 0;
    vector<vector<double>> params;
};


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void drawChos(QCustomPlot *customPlot);

private slots:
    void spinboxValueChanged();

    void on_actionOpen_triggered();

    void on_fitButton_clicked();

    void on_pushButton_clicked();

    void on_sigAffectionButton_clicked();

    void on_checkRssButton_clicked();

private:
    Ui::MainWindow *ui;
    DistributionData *data;
    std::vector<double> descentProgress;
    std::vector<double> fitParams;
    std::vector<ChosDistribution> chosVector;

    void setupSlotConnection();
    QCPBars* createBar(QCustomPlot *plotParent, double x, double height, double width);

    double getChosValue(double x);
    double getChosValueWithDistrParams(double x, double mean, double sig, double as, double ex);


    fitResults calcFit2(PointsVector points, int quantilElem, int shakeCount);
    fitResults calcFit(PointsVector points, int quantilElem, int shakeCount);

};

#endif // MAINWINDOW_H
