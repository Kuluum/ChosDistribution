#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "chosdistribution.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "FileProcessing/filereader.h"


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
  for (double i=a; i<x1; i+=0.001)
  {
    //auto grad = ChosDistribution::functionGradient(i, mean1, sig1, as1, ex1);
    x.append(i);
//    y.append();
//    if (i < x1) {
        y.append(ChosDistribution::value(i, m1, a1, beta1, ny1));
       // std::vector<double> grad = ChosDistribution::grad(i, m1, a1, beta1, ny1);
       // dx.append(grad[0]);
       // dm.append(grad[1]);
       // da.append(grad[2]);
       // dbeta.append(grad[3]);
       // dny.append(grad[4]);
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
 // customPlot->addGraph();
 // customPlot->graph(1)->setData(x, dny);
}

void MainWindow::drawDiffChos(QCustomPlot *customPlot)
{
    QVector<double> x, y;

    auto sqrDiff = sqrDiffVector();
    for (auto &p : sqrDiff)
    {
        x.append(p.first);
        y.append(p.second);
    }

    customPlot->addGraph();
    customPlot->graph(1)->setData(x, y);
}

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
    drawDiffChos(ui->plotWidget);
    ui->plotWidget->replot();
    ChosDistribution d;
    d.setDistribution(data);

    double mean1 = ui->mean1DoubleSpinBox->value();
    double sig1 = ui->sig1DoubleSpinBox->value();
    double as1 = ui->as1DoubleSpinBox->value();
    double ex1 = ui->ex1DoubleSpinBox->value();

    //qDebug() << d.RSS(mean1, sig1, as1, ex1);
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
   // auto relPoints = data->getStepRelativePoints();
  //  qDebug() << relPoints;
   // qDebug() << "point sum =" << std::accumulate(relPoints.begin(), relPoints.end(), 0.0, [](double a, QPair<double, double>p) -> double {
   //     return a + p.second;
   // });
    //qDebug() << data->getDistributionParameters(0, 6);
    //qDebug() << data->getDistributionParameters(7, 13);
}



//#include "alglib/interpolation.h"

//using namespace alglib;
//void function_cx_1_func(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
//{
//    // this callback calculates f(c,x)=exp(-c0*sqr(x0))
//    // where x is a position on X-axis and c is adjustable parameter

//    if (c[1] < 0.4 || c[1]>10) {
//        func = 1e+13;
//    }
//    if (fabs(c[2]) > 2.5) {
//        func = 1e+13;
//    }
//    if (c[3] <= 1.5*c[2]*c[2]) {
//        func = 1e+13;
//    }
//    double a = c[0];
//    double m = 2 / (c[3] - c[2]*c[2]);
//    double ny = m*c[1]*c[2]/2;
//    double beta = sqrt(m*c[1]*c[1] - ny*ny);

//    func = ChosDistribution::value(x[0], m, a, beta, ny);

//    //func = exp(-c[0]*pow(x[0],2));
//}

void MainWindow::on_fitButton_clicked()
{
  //  double a = data->getStepRelativePoints().first().first - 0.25;
   // double b = data->getStepRelativePoints().last().first + 0.25;
    //double x1 = b;

  //  double minDiff = 10000000000.0;
  //  double minMean = 0;
  //  double minSig = 0;
  //  double minAs = 0;
  //  double minEx = 0;

//    std::vector<QPair<double, double>> stdVect = data->getStepRelativePoints().toStdVector();

    //for (auto &p : data->getStepRelativePoints()) {
    //    std::pair pair = make_pair(p.first, p.second);
    //    stdVect
    //}
/*
    int numThreads = ui->numThreadsSpinBox->value();

    {
        std::cout << "num threads = " << numThreads << " ";
        PROFILE_BLOCK("time");

  //  for (int mean = -5; mean < 5; mean += 1 ) {
       int mean = 0;
#pragma omp parallel for num_threads(numThreads)
        for (int sig10 = 0; sig10 < 95; sig10++) {
            std::cout << sig10/95. * 100.0 << "%"<<std::endl;
            for (int as10 = 0; as10 < 50; as10++) {
                for (int ex10 = 0; ex10 < 100; ex10++) {

                    double sig = sig10/10.0 + 0.4;
                    double as = as10/10.0 - 2.5;
                    double ex = as*as + ex10/10.0 + 0.1;

                    //for (int x10 = 0; x10 <100; x10++) {
                       // std::cout << "sig = " << sig10 << " as = " << as10 << " ex = " << ex10 << " x = " << x10;
                        double x = 5;
                        double value = getChosValueWithDistrParams(x, 0, sig, as, ex);
                        //chosValuesArr[sig10][as10][ex10][x10] = value;
                    //}
                }
            }
        }
    //}
    }

    */

    //
    // In this example we demonstrate exponential fitting
    // by f(x) = exp(-c*x^2)
    // using function value and gradient (with respect to c).
    //
    //
       // In this example we demonstrate exponential fitting
       // by f(x) = exp(-c*x^2)
       // using function value only.
       //
       // Gradient is estimated using combination of numerical differences
       // and secant updates. diffstep variable stores differentiation step
       // (we have to tell algorithm what step to use).
       //
//       real_2d_array x = "[[-1.5],[-1.0],[-0.5],[0.0],[0.5],[1.0],[1.5]]";
//       real_1d_array y = "[0.129, 0.241, 0.352, 0.398, 0.352, 0.241, 0.129]";
//       real_1d_array c = "[-0.5, 0.15, 0.1, 0.2]";
//       real_1d_array bndl = "[-1.0, 0.05, -2.5, 0.1]";
//       real_1d_array bndu = "[1.0, 10.0, 2.5, 10.0]";
//       double epsf = 1;
//       double epsx = 0.00001;
//       ae_int_t maxits = 0;
//       ae_int_t info;
//       lsfitstate state;
//       lsfitreport rep;
//       double diffstep = 0.0001;

//       //
//       // Fitting without weights
//       //
//       lsfitcreatef(x, y, c, diffstep, state);
//       lsfitsetbc(state, bndl, bndu);
//       lsfitsetcond(state, epsf, epsx, maxits);
//       alglib::lsfitfit(state, function_cx_1_func);
//       lsfitresults(state, info, c, rep);
//       printf("%d\n", int(info)); // EXPECTED: 2
//       std::string str = c.tostring(4);
//       qDebug() << ""; // EXPECTED: [1.5]

       if (data != nullptr)
       {
           ChosDistribution distr;
           distr.setDistribution(data);
           auto params = data->getDistributionParameters(0, 100);
           distr.setInitialParams({params[0], params[1], 0, 0.3});
//           distr.iterate();
           distr.gradDescent();
           //distr.gradDescentQuadr();
       }


       //
       // Fitting with weights
       // (you can change weights and see how it changes result)
       //
     //  real_1d_array w = "[1,1,1,1,1,1,1,1,1,1,1]";
     //  lsfitcreatewf(x, y, w, c, diffstep, state);
     //  lsfitsetcond(state, epsf, epsx, maxits);
     //  alglib::lsfitfit(state, function_cx_1_func);
     //  lsfitresults(state, info, c, rep);
     //  printf("%d\n", int(info)); // EXPECTED: 2
     //  printf("%s\n", c.tostring(1).c_str()); // EXPECTED: [1.5]
}
