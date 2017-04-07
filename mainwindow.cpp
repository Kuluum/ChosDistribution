#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "chosdistribution.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "FileProcessing/filereader.h"

//UI
#include "descentprogressview.h"

//TODO: remove
#include <numeric>
#include <omp.h>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    drawChos(ui->plotWidget);
    ui->plotWidget->setInteraction(QCP::iRangeDrag, true);
    ui->plotWidget->setInteraction(QCP::iRangeZoom, true);

    ui->plotWidget->xAxis->setLabel("x");
    ui->plotWidget->yAxis->setLabel("y");

    ui->plotWidget->xAxis->setRange(-5, 5);
    ui->plotWidget->yAxis->setRange(0, 0.5);
    setupSlotConnection();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::drawChos(QCustomPlot *customPlot)
{
  QVector<double> x, y, dx, dm, da, dbeta, dny;

  double a = ui->aDoubleSpinBox->value();
  double b = ui->bDoubleSpinBox->value();

  double mean1 = ui->mean1DoubleSpinBox->value();
  double sig1 = ui->sig1DoubleSpinBox->value();
  double as1 = ui->as1DoubleSpinBox->value();
  double ex1 = ui->ex1DoubleSpinBox->value();

  double a1 = mean1;
  double m1 = 2 / (ex1 - as1*as1);
  double ny1 = m1*sig1*as1/2;
  double beta1 = sqrt(m1*sig1*sig1 - ny1*ny1);

  double x1 = ui->x1DoubleSpinBox->value();
/*
  double mean2 = ui->mean2DoubleSpinBox->value();
  double sig2 = ui->sig2DoubleSpinBox->value();
  double as2 = ui->as2DoubleSpinBox->value();
  double ex2 = ui->ex2DoubleSpinBox->value();

  double a2 = mean2;
  double m2 = 2 / (ex2 - as2*as2);
  double ny2 = m2*sig2*as2/2;
  double beta2 = sqrt(m2*sig2*sig2 - ny2*ny2);

  double x2 = ui->x2DoubleSpinBox->value();

  double mean3 = ui->mean3DoubleSpinBox->value();
  double sig3 = ui->sig3DoubleSpinBox->value();
  double as3 = ui->as3DoubleSpinBox->value();
  double ex3 = ui->ex3DoubleSpinBox->value();

  double a3 = mean3;
  double m3 = 2 / (ex3 - as3*as3);
  double ny3 = m3*sig3*as3/2;
  double beta3 = sqrt(m3*sig3*sig3 - ny3*ny3);

  qDebug() << "a1 = " << a1 << " m1 = " << m1 << " ny1 = " << ny1 << " beta1 = " << beta1;
  qDebug() << "a2 = " << a2 << " m2 = " << m2 << " ny2 = " << ny2 << " beta2 = " << beta2;
  qDebug() << "a3 = " << a3 << " m3 = " << m3 << " ny3 = " << ny3 << " beta3 = " << beta3;
*/
  for (double i=-10; i<10; i+=0.001)
  {

    x.append(i);

//    if (i < x1) {
        y.append(ChosDistribution::value(i, m1, a1, beta1, ny1));
//    }
//    else if (i < x2) {
//        y.append(ChosDistribution::value(i, m2, a2, beta2, ny2));
//    }
//    else {
//        y.append(ChosDistribution::value(i, m3, a3, beta3, ny3));
//    }
  }

  customPlot->addGraph();
  customPlot->graph(0)->setData(x, y);
}

//void MainWindow::drawDiffChos(QCustomPlot *customPlot)
//{
//    QVector<double> x, y;

//    auto sqrDiff = sqrDiffVector();
//    for (auto &p : sqrDiff)
//    {
//        x.append(p.first);
//        y.append(p.second);
//    }

//    customPlot->addGraph();
//    customPlot->graph(1)->setData(x, y);
//}

DisVector MainWindow::sqrDiffVector() {
    DisVector sqrDiffVector;

    for (auto &p : data->getStepRelativePoints()) {
        double x = p.first;
        double chos = getChosValue(x);
        double y = (p.second - chos) * (p.second - chos);
        sqrDiffVector.push_back(make_pair(x, y));
    }
    return sqrDiffVector;
}

double MainWindow::getChosValue(double x)
{

  double mean1 = ui->mean1DoubleSpinBox->value();
  double sig1 = ui->sig1DoubleSpinBox->value();
  double as1 = ui->as1DoubleSpinBox->value();
  double ex1 = ui->ex1DoubleSpinBox->value();

  double a1 = mean1;
  double m1 = 2 / (ex1 - as1*as1);
  double ny1 = m1*sig1*as1/2;
  double beta1 = sqrt(m1*sig1*sig1 - ny1*ny1);

  double x1 = ui->x1DoubleSpinBox->value();

  double mean2 = ui->mean2DoubleSpinBox->value();
  double sig2 = ui->sig2DoubleSpinBox->value();
  double as2 = ui->as2DoubleSpinBox->value();
  double ex2 = ui->ex2DoubleSpinBox->value();

  double a2 = mean2;
  double m2 = 2 / (ex2 - as2*as2);
  double ny2 = m2*sig2*as2/2;
  double beta2 = sqrt(m2*sig2*sig2 - ny2*ny2);

  double x2 = ui->x2DoubleSpinBox->value();

  double mean3 = ui->mean3DoubleSpinBox->value();
  double sig3 = ui->sig3DoubleSpinBox->value();
  double as3 = ui->as3DoubleSpinBox->value();
  double ex3 = ui->ex3DoubleSpinBox->value();

  double a3 = mean3;
  double m3 = 2 / (ex3 - as3*as3);
  double ny3 = m3*sig3*as3/2;
  double beta3 = sqrt(m3*sig3*sig3 - ny3*ny3);


//    if (x < x1) {
        return getChosValueWithDistrParams(x, mean1, sig1, as1, ex1);
//    }
//    else if (x < x2) {
//        return ChosDistribution::value(x, m2, a2, beta2, ny2);
//    }
//    else {
//        return ChosDistribution::value(x, m3, a3, beta3, ny3);
//    }

}

double MainWindow::getChosValue(double x, double x1, double x2,
                                double m1, double a1, double beta1, double ny1,
                                double m2, double a2, double beta2, double ny2,
                                double m3, double a3, double beta3, double ny3)
{
    if (x < x1) {
        return ChosDistribution::value(x, m1, a1, beta1, ny1);
    }
    else if (x < x2) {
        return ChosDistribution::value(x, m2, a2, beta2, ny2);
    }
    else {
        return ChosDistribution::value(x, m3, a3, beta3, ny3);
    }
}

double MainWindow::getChosValueWithDistrParams(double x, double mean, double sig, double as, double ex) {
    double a = mean;
    double m = 2 / (ex - as*as);
    double ny = m*sig*as/2;
    double beta = sqrt(m*sig*sig - ny*ny);

    return ChosDistribution::value(x, m, a, beta, ny);
}

void MainWindow::on_drawButton_clicked()
{
    ui->plotWidget->clearGraphs();
    drawChos(ui->plotWidget);
    //drawDiffChos(ui->plotWidget);
    ui->plotWidget->replot();
    ChosDistribution d;
    d.setDistribution(data);

    double mean1 = ui->mean1DoubleSpinBox->value();
    double sig1 = ui->sig1DoubleSpinBox->value();
    double as1 = ui->as1DoubleSpinBox->value();
    double ex1 = ui->ex1DoubleSpinBox->value();
}

void MainWindow::setupSlotConnection() {
    connect(ui->mean1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->sig1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->as1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->ex1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));

    connect(ui->x1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));

    connect(ui->mean2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->sig2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->as2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->ex2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));

    connect(ui->x2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));

    connect(ui->mean3DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->sig3DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->as3DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->ex3DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
}

void MainWindow::spinboxValueChanged()
{
    ui->plotWidget->clearGraphs();
    drawChos(ui->plotWidget);
    ui->plotWidget->replot();
}

QCPBars* MainWindow::createBar(QCustomPlot *plotParent, double x, double height, double width) {
    QCPBars *bar = new QCPBars(plotParent->xAxis, plotParent->yAxis);
    bar->setPen(QPen(QColor(0, 168, 140).lighter(130)));
    bar->addData(x, height);
    bar->setWidth(width);
    return bar;
}

void MainWindow::on_actionOpen_triggered()
{
    FileReader f;
    QStringList stringList = f.readFileWithDialog(this);
    QVector<QPair<double, double>> points;
    data = new DistributionData();
    for (QString str : stringList) {
        QStringList pointCoords = str.split(" ");
        if (pointCoords.size() == 2) {
            QString xStr = pointCoords.at(0);
            double x = xStr.toDouble();
            double y = pointCoords.at(1).toDouble();
            points.append(qMakePair<double, double>(x, y));
            data->addPoint(x, y);
        }
    }
    ui->plotWidget->clearPlottables();
    for (auto &p : data->getStepRelativePoints()) {
        createBar(ui->plotWidget, p.first, p.second, 0.5);
    }
    ui->plotWidget->replot();
}



void MainWindow::on_fitButton_clicked()
{
       if (data != nullptr)
       {
           ChosDistribution distr;
           distr.setDistribution(data);
           auto params = data->getDistributionParameters(0, 200);
           distr.setInitialParams({params[0], params[1], 0, 0.3});
           auto matchParams = distr.gradDescent();
           descentProgress = distr.descentProgress;
           ui->mean1DoubleSpinBox->setValue(matchParams[0]);
           ui->sig1DoubleSpinBox->setValue(matchParams[1]);
           ui->as1DoubleSpinBox->setValue(matchParams[2]);
           ui->ex1DoubleSpinBox->setValue(matchParams[3]);
           drawChos(ui->plotWidget);
       }
}

void MainWindow::on_pushButton_clicked()
{
    DescentProgressView  *descentProgressView = new DescentProgressView();
    descentProgressView->setDescentProgress(this->descentProgress);
    descentProgressView -> show();
}
