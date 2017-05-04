#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "chosdistribution.h"


#include <stdio.h>
#include <math.h>

#include "FileProcessing/filereader.h"

//UI
#include "descentprogressview.h"
#include "graphwindow.h"

// std
#include <future>
#include <thread>

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


  double x1 = ui->x1DoubleSpinBox->value();

  double mean2 = ui->mean2DoubleSpinBox->value();
  double sig2 = ui->sig2DoubleSpinBox->value();
  double as2 = ui->as2DoubleSpinBox->value();
  double ex2 = ui->ex2DoubleSpinBox->value();



  double x2 = ui->x2DoubleSpinBox->value();

  double mean3 = ui->mean3DoubleSpinBox->value();
  double sig3 = ui->sig3DoubleSpinBox->value();
  double as3 = ui->as3DoubleSpinBox->value();
  double ex3 = ui->ex3DoubleSpinBox->value();
/*
  double a3 = mean3;
  double m3 = 2 / (ex3 - as3*as3);
  double ny3 = m3*sig3*as3/2;
  double beta3 = sqrt(m3*sig3*sig3 - ny3*ny3);

  qDebug() << "a1 = " << a1 << " m1 = " << m1 << " ny1 = " << ny1 << " beta1 = " << beta1;
  qDebug() << "a2 = " << a2 << " m2 = " << m2 << " ny2 = " << ny2 << " beta2 = " << beta2;
  qDebug() << "a3 = " << a3 << " m3 = " << m3 << " ny3 = " << ny3 << " beta3 = " << beta3;
*/
  for (double i=a; i<=b; i+=0.001)
  {

    x.append(i);

    if (i < x1) {
        y.append(ChosDistribution::valueWithDistrParams(i, mean1, sig1, as1, ex1));
    }
    else if (i < x2) {
        y.append(ChosDistribution::valueWithDistrParams(i, mean2, sig2, as2, ex2));
    }
    else {
        y.append(ChosDistribution::valueWithDistrParams(i, mean3, sig3, as3, ex3));
    }
  }

  customPlot->addGraph();
  customPlot->graph(0)->setData(x, y);
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

double MainWindow::getChosValueWithDistrParams(double x, double mean, double sig, double as, double ex) {
    double a = mean;
    double m = 2 / (ex - as*as);
    double ny = m*sig*as/2;
    double beta = sqrt(m*sig*sig - ny*ny);

    return ChosDistribution::value(x, m, a, beta, ny);
}


void MainWindow::setupSlotConnection() {
    connect(ui->aDoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));

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
    if (stringList.size() == 0) {
        return;
    }
//    QVector<QPair<double, double>> points;
    PointsVector points;
    data = new DistributionData();
    for (QString str : stringList) {
        QStringList pointCoords = str.split(" ");
        if (pointCoords.size() == 2) {
            QString xStr = pointCoords.at(0);
            double x = xStr.toDouble();
            double y = pointCoords.at(1).toDouble();
            points.push_back(make_pair(x, y));
        }
    }
    data->setPoints(points);
    double step = points[1].first - points[0].first;
    data->setStep(step);

    ui->quantilSpinBox->setValue(points.size()/3);
    ui->quant2SpinBox->setValue(points.size() * 2 / 3);

    ui->plotWidget->clearPlottables();
    for (auto &p : data->getStepRelativePoints()) {
        createBar(ui->plotWidget, p.first, p.second, step);
    }
    ui->plotWidget->replot();
}


fitResults MainWindow::calcFit(PointsVector points, int quantilElem, int shakeCount)
{
    fitResults results;

    ChosDistribution distr1;
    PointsVector distr1Points(points.begin(), points.begin()+quantilElem+1);
    auto params1 = data->getDistributionParameters(0, quantilElem+1);
    distr1.setPoints(distr1Points);
    distr1.setInitialParams(params1);
//    qDebug() << "first";
    auto matchParams1 = distr1.gradDescent(shakeCount);

    fitResults best2 = calcFit2(points, quantilElem, shakeCount);

    double rss = distr1.currentRss() + best2.rss;
    results.rss = rss;
    results.quantil1 = quantilElem;
    results.quantil2 = best2.quantil2;
    vector<vector<double>> v({matchParams1, best2.params[0], best2.params[1]});
    results.params = v;



    return results;
}


fitResults MainWindow::calcFit2(PointsVector points, int quantilElem, int shakeCount)
{
    fitResults result;

    int size = points.size();
    double bestRss = 100000000;
    vector<vector<double>> bestParams;
    int bestQ2 = -1;

//    for (int quantil2Elem = quantilElem + 4; quantil2Elem < size-4; quantil2Elem+=4)
    {
        int quantil2Elem = ui->quant2SpinBox->value();

        ChosDistribution distr2;
        ChosDistribution distr3;

//        qDebug() << "q1:" << quantilElem << " q2:" << quantil2Elem;

        PointsVector distr2Points(points.begin()+ quantilElem, points.begin()+quantil2Elem+1);
        PointsVector distr3Points(points.begin()+ quantil2Elem, points.end());

        auto params2 = data->getDistributionParameters(quantilElem, quantil2Elem+1);
        auto params3 = data->getDistributionParameters(quantil2Elem, size);

        distr2.setInitialParams({params2[0], params2[1], 0, 0.1});
        distr3.setInitialParams({params3[0], params3[1], 0, 0.1});


        distr2.setPoints(distr2Points);
        distr3.setPoints(distr3Points);

//        qDebug() << "second";
        auto matchParams2 = distr2.gradDescent(shakeCount);
//        qDebug() << "third";
        auto matchParams3 = distr3.gradDescent(shakeCount);

        double rss = distr2.currentRss() + distr3.currentRss();
        if (rss < bestRss) {
            bestParams = {matchParams2, matchParams3};
            bestQ2 = quantil2Elem;
            bestRss = rss;
        }
    }

    result.rss = bestRss;
    result.params = bestParams;
    result.quantil1 = quantilElem;
    result.quantil2 = bestQ2;

    return result;
}



void MainWindow::on_fitButton_clicked()
{
    if (data != nullptr)
    {

        chosVector.clear();


        auto points = data->getStepRelativePoints();
        int size = points.size();
        vector<vector<double>> bestParams;
        double bestRss = 10000000;

        vector<shared_future<fitResults>> fitFutureVector;

        fitResults bestResults;
        int shakeCount = ui->shakeSpinBox->value();
//        for (int q1 = 4; q1 < size-8; q1+=4)
        int q1 = ui->quantilSpinBox->value();
        {
            MainWindow *that = this;
            shared_future<fitResults> fitFuture = std::async(std::launch::async, [&points, q1, shakeCount, that]() -> fitResults
            {
                qDebug() << "start fit for q1 = " << q1;
                fitResults res = that->calcFit(points, q1, shakeCount);
                qDebug() << "end fit for q1 = " << q1;
                return res;
            });
//            shared_future<fitResults> fitFuture = std::async(calcFit, this, points, q1, shakeCount);
            fitFutureVector.push_back(std::move(fitFuture));
            //            fitResults res =

        }

//        for(future<fitResults> fitFuture : fitFutureVector)
//        {
//            fitFuture.wait();
//        }

        for_each(fitFutureVector.begin(), fitFutureVector.end(), [&](shared_future<fitResults> fitFuture)
        {
            fitResults res = fitFuture.get();
            if (res.rss < bestRss)
            {
                bestResults = res;
                bestRss = res.rss;
            }
        });

        auto matchParams1 = bestResults.params[0];
        auto matchParams2 = bestResults.params[1];
        auto matchParams3 = bestResults.params[2];


        ui->mean1DoubleSpinBox->setValue(matchParams1[0]);
        ui->sig1DoubleSpinBox->setValue(matchParams1[1]);
        ui->as1DoubleSpinBox->setValue(matchParams1[2]);
        ui->ex1DoubleSpinBox->setValue(matchParams1[3]);

        ui->mean2DoubleSpinBox->setValue(matchParams2[0]);
        ui->sig2DoubleSpinBox->setValue(matchParams2[1]);
        ui->as2DoubleSpinBox->setValue(matchParams2[2]);
        ui->ex2DoubleSpinBox->setValue(matchParams2[3]);

        ui->mean3DoubleSpinBox->setValue(matchParams3[0]);
        ui->sig3DoubleSpinBox->setValue(matchParams3[1]);
        ui->as3DoubleSpinBox->setValue(matchParams3[2]);
        ui->ex3DoubleSpinBox->setValue(matchParams3[3]);

        ui->aDoubleSpinBox->setValue(points[0].first);
        ui->x1DoubleSpinBox->setValue(points[bestResults.quantil1].first);
        ui->x2DoubleSpinBox->setValue(points[bestResults.quantil2].first);
        ui->bDoubleSpinBox->setValue(points.back().first);

//        ui->rss1Label->setText(QString::number(distr1.currentRss()));
//        ui->rss2Label->setText(QString::number(distr2.currentRss()));
//        ui->rss3Label->setText(QString::number(distr3.currentRss()));

        drawChos(ui->plotWidget);
   }
}



void MainWindow::on_pushButton_clicked()
{
    DescentProgressView  *descentProgressView = new DescentProgressView();
    descentProgressView->setDescentProgress(this->descentProgress);
    descentProgressView->show();
}

void MainWindow::on_sigAffectionButton_clicked()
{
    GraphWindow *graphWindow = new GraphWindow();

    vector<double> params = fitParams;

    ChosDistribution distr;
    distr.setPoints(data->getStepRelativePoints());
    distr.setInitialParams(fitParams);
    vector<pair<double, double>> points;
    for(double d = 0.1; d < 5; d+=0.01)
    {
        params[1] = d;
        points.push_back(make_pair(d, distr.RSS(data->getStepRelativePoints(), params[0], params[1], params[2], params[3])));
    }

    graphWindow->setPointsVector(points);
    graphWindow->show();
}

void MainWindow::on_checkRssButton_clicked()
{
    auto points = data->getStepRelativePoints();


    int quantilElem = ui->quantilSpinBox->value();
    int quantil2Elem = ui->quant2SpinBox->value();
    int size = points.size();

    //348
    PointsVector distr1Points(points.begin(), points.begin()+quantilElem+1);
    PointsVector distr2Points(points.begin()+ quantilElem, points.begin()+quantil2Elem+1);
    PointsVector distr3Points(points.begin()+ quantil2Elem, points.end());

    double m1 = ui->mean1DoubleSpinBox->value();
    double s1 = ui->sig1DoubleSpinBox->value();
    double a1 = ui->as1DoubleSpinBox->value();
    double e1 = ui->ex1DoubleSpinBox->value();

    ChosDistribution::RSS(distr1Points, m1, s1, a1, e1);

}
