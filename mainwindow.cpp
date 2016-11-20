#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "chosdistribution.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "alglib/interpolation.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QCustomPlot *customPlot = ui->plotWidget;
    setupQuadraticDemo(ui->plotWidget);
    ui->plotWidget->setInteraction(QCP::iRangeDrag, true);
    ui->plotWidget->setInteraction(QCP::iRangeZoom, true);

    ui->plotWidget->xAxis->setLabel("x");
    ui->plotWidget->yAxis->setLabel("y");

    ui->plotWidget->xAxis->setRange(-5, 5);
    ui->plotWidget->yAxis->setRange(0, 0.5);
    setupSlotConnection();

    createBar(customPlot, 0, 0.2, 0.5);
    createBar(customPlot, 0.5, 0.3, 0.2);


}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setupQuadraticDemo(QCustomPlot *customPlot)
{
 // demoName = "Quadratic Demo";
  // generate some data:
  QVector<double> x, y;
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

  qDebug() << "a1 = " << a1 << " m1 = " << m1 << " ny1 = " << ny1 << " beta1 = " << beta1;
  qDebug() << "a2 = " << a2 << " m2 = " << m2 << " ny2 = " << ny2 << " beta2 = " << beta2;
  qDebug() << "a3 = " << a3 << " m3 = " << m3 << " ny3 = " << ny3 << " beta3 = " << beta3;

  for (double i=-10; i<=10; i+=0.01)
  {
    x.append(i);
    if (i < x1) {
        y.append(ChosDistribution::value(i, m1, a1, beta1, ny1));
    }
    else if (i < x2) {
        y.append(ChosDistribution::value(i, m2, a2, beta2, ny2));
    }
    else {
        y.append(ChosDistribution::value(i, m3, a3, beta3, ny3));
    }
  }
  // create graph and assign data to it:
  customPlot->addGraph();
  customPlot->graph(0)->setData(x, y);
}

void MainWindow::on_drawButton_clicked()
{
    ui->plotWidget->clearGraphs();
    setupQuadraticDemo(ui->plotWidget);
    ui->plotWidget->replot();
    fit();
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
    setupQuadraticDemo(ui->plotWidget);
    ui->plotWidget->replot();
}

QCPBars* MainWindow::createBar(QCustomPlot *plotParent, double x, double height, double width) {
    QCPBars *bar = new QCPBars(plotParent->xAxis, plotParent->yAxis);
    bar->setPen(QPen(QColor(0, 168, 140).lighter(130)));
    bar->addData(x, height);
    bar->setWidth(width);
    return bar;
}

using namespace alglib;
void function_cx_1_func(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    // this callback calculates f(c,x)=exp(-c0*sqr(x0))
    // where x is a position on X-axis and c is adjustable parameter
    //func = c[0]*x[0] + 3*c[1];
    func = ChosDistribution::value(x[0], c[0], c[1], c[2], c[3]);
}

void MainWindow::fit() {
    //
       // In this example we demonstrate exponential fitting by
       //     f(x) = exp(-c*x^2)
       // subject to bound constraints
       //     0.0 <= c <= 1.0
       // using function value only.
       //
       // Gradient is estimated using combination of numerical differences
       // and secant updates. diffstep variable stores differentiation step
       // (we have to tell algorithm what step to use).
       //
       // Unconstrained solution is c=1.5, but because of constraints we should
       // get c=1.0 (at the boundary).
       //
       real_2d_array x = "[[1.58],[2.75],[3.92], [5.08], [6.25], [7.4]]";
       real_1d_array y = "[0.12, 0.16, 0.4, 0.12, 0.16, 0.04]";
       real_1d_array c = "[0, 0, 0, 0]";
       real_1d_array bndl = "[-20, 0, 0, 0]";
       real_1d_array bndu = "[20, 10, 10, 10]";
       double epsf = 0;
       double epsx = 0.000001;
       ae_int_t maxits = 0;
       ae_int_t info;
       lsfitstate state;
       lsfitreport rep;
       double diffstep = 0.0001;

       lsfitcreatef(x, y, c, diffstep, state);
       lsfitsetbc(state, bndl, bndu);
       lsfitsetcond(state, epsf, epsx, maxits);
       alglib::lsfitfit(state, function_cx_1_func);
       lsfitresults(state, info, c, rep);
       qDebug() << ("%s\n", c.tostring(1).c_str()); // EXPECTED: [1.0]
}
