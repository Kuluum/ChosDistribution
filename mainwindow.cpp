#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "chosdistribution.h"
#include "tabledelegate.h"

#include <stdio.h>
#include <math.h>

#include "FileProcessing/filereader.h"
#include "Algorithm/algorithms.h"

//UI
#include "descentprogressview.h"
#include "graphwindow.h"

// std
#include <future>
#include <thread>

#include <QLocale>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
//    drawChos(ui->plotWidget);
    ui->plotWidget->setInteraction(QCP::iRangeDrag, true);
    ui->plotWidget->setInteraction(QCP::iRangeZoom, true);
    ui->plotWidget->setInteraction(QCP::iSelectAxes, true);

    ui->plotWidget->xAxis->setLabel("x");
    ui->plotWidget->yAxis->setLabel("y");

    ui->plotWidget->xAxis->setRange(-5, 5);
    ui->plotWidget->yAxis->setRange(0, 0.5);

    model = new QStandardItemModel(0, 2, this);
    connect(model, SIGNAL(dataChanged(QModelIndex,QModelIndex)), this, SLOT(tableDataChanged(QModelIndex,QModelIndex)));
    connect(ui->plotWidget, SIGNAL(mouseWheel(QWheelEvent*)), this, SLOT(mouseWheel()));
    connect(ui->plotWidget, SIGNAL(selectionChangedByUser()), this, SLOT(selectionChanged()));
    model->setHorizontalHeaderItem(0, new QStandardItem(QString("x")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("y")));
    TableDelegate *delegate = new TableDelegate;
    ui->pointsTable->setItemDelegate(delegate);
    ui->pointsTable->setModel(model);
    ui->pointsTable->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    ui->pointsTable->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    data = nullptr;
    setupSlotConnection();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::selectionChanged()
{
  /*
   normally, axis base line, axis tick labels and axis labels are selectable separately, but we want
   the user only to be able to select the axis as a whole, so we tie the selected states of the tick labels
   and the axis base line together. However, the axis label shall be selectable individually.

   The selection state of the left and right axes shall be synchronized as well as the state of the
   bottom and top axes.

   Further, we want to synchronize the selection of the graphs with the selection state of the respective
   legend item belonging to that graph. So the user can select a graph by either clicking on the graph itself
   or on its legend item.
  */

  // make top and bottom axes be selected synchronously, and handle axis and tick labels as one selectable object:
  if (ui->plotWidget->xAxis->selectedParts().testFlag(QCPAxis::spAxis) || ui->plotWidget->xAxis->selectedParts().testFlag(QCPAxis::spTickLabels) ||
      ui->plotWidget->xAxis2->selectedParts().testFlag(QCPAxis::spAxis) || ui->plotWidget->xAxis2->selectedParts().testFlag(QCPAxis::spTickLabels))
  {
    ui->plotWidget->xAxis2->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
    ui->plotWidget->xAxis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
  }
  // make left and right axes be selected synchronously, and handle axis and tick labels as one selectable object:
  if (ui->plotWidget->yAxis->selectedParts().testFlag(QCPAxis::spAxis) || ui->plotWidget->yAxis->selectedParts().testFlag(QCPAxis::spTickLabels) ||
      ui->plotWidget->yAxis2->selectedParts().testFlag(QCPAxis::spAxis) || ui->plotWidget->yAxis2->selectedParts().testFlag(QCPAxis::spTickLabels))
  {
    ui->plotWidget->yAxis2->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
    ui->plotWidget->yAxis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
  }

  // synchronize selection of graphs with selection of corresponding legend items:
  for (int i=0; i<ui->plotWidget->graphCount(); ++i)
  {
    QCPGraph *graph = ui->plotWidget->graph(i);
    QCPPlottableLegendItem *item = ui->plotWidget->legend->itemWithPlottable(graph);
    if (item->selected() || graph->selected())
    {
      item->setSelected(true);
      graph->setSelection(QCPDataSelection(graph->data()->dataRange()));
    }
  }
}

void MainWindow::mouseWheel()
{
  // if an axis is selected, only allow the direction of that axis to be zoomed
  // if no axis is selected, both directions may be zoomed

  if (ui->plotWidget->xAxis->selectedParts().testFlag(QCPAxis::spAxis)){
    ui->plotWidget->axisRect()->setRangeZoomAxes(ui->plotWidget->xAxis,ui->plotWidget->yAxis);
    ui->plotWidget->axisRect()->setRangeZoom(ui->plotWidget->xAxis->orientation());
  }
  else if (ui->plotWidget->yAxis->selectedParts().testFlag(QCPAxis::spAxis)){
    ui->plotWidget->axisRect()->setRangeZoomAxes(ui->plotWidget->xAxis,ui->plotWidget->yAxis);
    ui->plotWidget->axisRect()->setRangeZoom(ui->plotWidget->yAxis->orientation());
  }
  else if (ui->plotWidget->xAxis2->selectedParts().testFlag(QCPAxis::spAxis)){
    ui->plotWidget->axisRect()->setRangeZoomAxes(ui->plotWidget->xAxis2,ui->plotWidget->yAxis2);
    ui->plotWidget->axisRect()->setRangeZoom(ui->plotWidget->xAxis2->orientation());
  }
  else if (ui->plotWidget->yAxis2->selectedParts().testFlag(QCPAxis::spAxis)){
    ui->plotWidget->axisRect()->setRangeZoomAxes(ui->plotWidget->xAxis2,ui->plotWidget->yAxis2);
    ui->plotWidget->axisRect()->setRangeZoom(ui->plotWidget->yAxis2->orientation());
  }
  else
    ui->plotWidget->axisRect()->setRangeZoom(Qt::Horizontal|Qt::Vertical);
}

void MainWindow::drawChos(QCustomPlot *plotWidget)
{
  QVector<double> x, y;

  double a = ui->aDoubleSpinBox->value();
  double b = ui->bDoubleSpinBox->value();

  double a1 = ui->a1DoubleSpinBox->value();
  double m1 = ui->m1DoubleSpinBox->value();
  double ny1 = ui->ny1DoubleSpinBox->value();
  double beta1 = ui->beta1DoubleSpinBox->value();


  double x1 = ui->x1DoubleSpinBox->value();

  double a2 = ui->a2DoubleSpinBox->value();
  double m2 = ui->m2DoubleSpinBox->value();
  double ny2 = ui->ny2DoubleSpinBox->value();
  double beta2 = ui->beta2DoubleSpinBox->value();

  double x2 = ui->x2DoubleSpinBox->value();

  double a3 = ui->a3DoubleSpinBox->value();
  double m3 = ui->m3DoubleSpinBox->value();
  double ny3 = ui->ny3DoubleSpinBox->value();
  double beta3 = ui->beta3DoubleSpinBox->value();

  double step = 1;
  if (data)
  {
      step = data->getStep();
  }
  for (double i=a-step/2; i<=b+step/2; i+=step/100)
  {

    x.append(i);

    if (i < x1) {
        // TODO a<->m beta<->ny
        y.append(ChosDistribution::value(i, m1, a1, beta1, ny1));
    }
    else if (i < x2) {
        y.append(ChosDistribution::value(i, m2, a2, beta2, ny2));
    }
    else {
        y.append(ChosDistribution::value(i, m3, a3, beta3, ny3));
    }
  }

  if (data) {
      double firstRss = 0.0;
      double secRss = 0.0;
      double thrRss = 0.0;

      double rss = 0.0;
      double khi = 0.0;
      double size = data ->getDistributionSize();
      double integral = 0.0;
      for(auto &p : data->getStepRelativePoints()) {
          double value = 0.0;
          double i = p.first;
          double integralValue;
          if (i < x1) {
              value = ChosDistribution::value(i, m1, a1, beta1, ny1);

              integralValue = Algorithms::Integral(i-step/2.0, i+step/2.0, [&](double d)->double{
                  double v = ChosDistribution::value(d, m1, a1, beta1, ny1);
                  return v;
              });
          }
          else if (i < x2) {
              value = ChosDistribution::value(i, m2, a2, beta2, ny2);

              integralValue = Algorithms::Integral(i-step/2.0, i+step/2.0, [&](double d)->double{
                  double v = ChosDistribution::value(d, m2, a2, beta2, ny2);
//                  qDebug() << "integer value 2" << v;
                  return v;
              });
          }
          else {
              value = ChosDistribution::value(i, m3, a3, beta3, ny3);

              integralValue = Algorithms::Integral(i-step/2.0, i+step/2.0, [&](double d)->double{
                  double v = ChosDistribution::value(d, m3, a3, beta3, ny3);
//                  qDebug() << "integer value 3 ("<<d<<", "<<m3<<", "<<a3<<", "<<beta3<<", "<<ny3<<") = " << v;
                  return v;
              });
          }
          if (isnan(value)) {
              break;
          }
          integral += integralValue;
          double r = p.second - value;
          rss += r*r;

          if (i < x1) {
            firstRss += r*r;
          }
          else if (i < x2) {
              secRss += r*r;
          }
          else {
              thrRss += r*r;
          }

          double sq = size * p.second * step;
//          qDebug() << "sq="<<sq;
          integralValue *= size;
//          qDebug() << "integral=" << integralValue;
          double chi = 0.0;
          chi += integralValue - sq;
          chi *= chi;
          chi /= integralValue;
          khi += chi;

      }
      qDebug()<<"firsRss="<<firstRss<<" secRss=" << secRss << "thrRss=" << thrRss;
      int freed = data->getPoints().size() - 14;
      if (freed <= 0) {
          freed = 1;
      }
      double chiCrit95 = Algorithms::ChiCritical(freed, 0.95);
      double chiCrit99 = Algorithms::ChiCritical(freed, 0.99);
      QString result = QString("chi = %1\nchiCrit(0.05) = %2\nchiCrit(0.01) = %3\n").arg(
                   QString::number(khi),
                  QString::number(chiCrit95),QString::number(chiCrit99)
                  );
//                  QString::number(integral));



      ui->resultTextEdit->setText(result);
  }
  plotWidget->addGraph();
  plotWidget->graph(0)->setData(x, y);
  plotWidget->replot();
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

    connect(ui->a1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->m1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->ny1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->beta1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));

    connect(ui->x1DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));

    connect(ui->a2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->m2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->ny2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->beta2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));

    connect(ui->x2DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));

    connect(ui->a3DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->m3DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->ny3DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
    connect(ui->beta3DoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(spinboxValueChanged()));
}

void MainWindow::spinboxValueChanged()
{
    QSpinBox *sender = (QSpinBox*)QObject::sender();
    if (sender->hasFocus()) {
        ui->plotWidget->clearGraphs();
        drawChos(ui->plotWidget);
        ui->plotWidget->replot();
    }
}

QCPBars* MainWindow::createBar(QCustomPlot *plotParent, QVector<double> x, QVector<double> height, double width) {
    QCPBars *bar = new QCPBars(plotParent->xAxis, plotParent->yAxis);
    bar->setPen(QPen(QColor(0, 168, 140).lighter(130)));
    bar->setData(x, height);
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
    PointsVector points;
    data = new DistributionData();
    for (QString str : stringList) {
        QStringList pointCoords = str.split(" ");
        if (pointCoords.size() >= 2) {
            QString xStr = pointCoords.at(0);
            double x = xStr.toDouble();
            double y = pointCoords.at(1).toDouble();
            points.push_back(make_pair(x, y));
        }
    }

    model->clear();
    model->setHorizontalHeaderItem(0, new QStandardItem(QString("x")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("y")));

    data->setPoints(points);
   //points[1].first - points[0].first;
    //data->setStep(step);

    ui->plotWidget->clearPlottables();

    auto stepRelPoints = data->getStepRelativePoints();
    double step = data->getStep();
    ui->plotWidget->xAxis->setRange(stepRelPoints.front().first - 2*step, stepRelPoints.back().first + 2*step);
    ui->plotWidget->yAxis->setRange(0, 2);

    QLocale locale = QLocale::system();
    for (auto &p : data->getPoints())
    {
        QStandardItem *xItem = new QStandardItem(locale.toString(p.first, 'f', 3));
        xItem->setFlags(xItem->flags() &  ~Qt::ItemIsEditable);
        QStandardItem *yItem = new QStandardItem(locale.toString(p.second, 'f', 3));
        QList<QStandardItem*> row({xItem, yItem});
        model->appendRow(row);
    }

    QVector<double>x,y;
    for (auto &p : stepRelPoints) {
        x.append(p.first);
        y.append(p.second);
    }
    ui->pointsTable->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    createBar(ui->plotWidget, x, y, step);
    ui->plotWidget->replot();
}

void MainWindow::on_fitButton_clicked()
{
    if (data != nullptr)
    {
        chosVector.clear();

        auto points = data->getStepRelativePoints();
        int size = points.size();
        double bestRss = 10000000;

        vector<shared_future<fitResults>> fitFutureVector;

        fitResults bestResults;
        int shakeCount = ui->shakeSpinBox->value();
        int from = 5;
        int to = size - 10;
        int step = 5;

        if (size <= 20) {
            from = 0;
            to = size-1;
            step = 1;
        }

        if (!ui->autoQuantilCheckBox->isChecked())
        {
            from = ui->q1SpinBox->value();
            to = from+1;
        }

        for (int q1 = from; q1 < to; q1 += step)
        {
            MainWindow *that = this;
//            shared_future<fitResults> fitFuture = std::async(std::launch::async, [&points, q1, shakeCount, that]() -> fitResults
//            {
//                qDebug() << "start fit for q1 = " << q1;
                fitResults res = that->calcFit(points, q1, shakeCount);
//                qDebug() << "end fit for q1 = " << q1;
              //  return res;
//            });
//            fitFutureVector.push_back(std::move(fitFuture));
//        }

//        for_each(fitFutureVector.begin(), fitFutureVector.end(), [&](shared_future<fitResults> fitFuture)
//        {
//            fitResults res = fitFuture.get();
            if (res.rss < bestRss)
            {
                bestResults = res;
                bestRss = res.rss;
            }
//        });
        }

        auto matchParams1 = bestResults.params[0];
        auto matchParams2 = bestResults.params[1];
        auto matchParams3 = bestResults.params[2];

        double mean1 = matchParams1[0];
        double sig1 = matchParams1[1];
        double as1 = matchParams1[2];
        double ex1 = matchParams1[3];

        double a1 = mean1;
        double m1 = 2 / (ex1 - as1*as1);
        double ny1 = m1*sig1*as1/2;
        double beta1 = sqrt(m1*sig1*sig1 - ny1*ny1);

        double mean2 = matchParams2[0];
        double sig2 = matchParams2[1];
        double as2 = matchParams2[2];
        double ex2 = matchParams2[3];

        double a2 = mean2;
        double m2 = 2 / (ex2 - as2*as2);
        double ny2 = m2*sig2*as2/2;
        double beta2 = sqrt(m2*sig2*sig2 - ny2*ny2);

        double mean3 = matchParams3[0];
        double sig3 = matchParams3[1];
        double as3 = matchParams3[2];
        double ex3 = matchParams3[3];

        double a3 = mean3;
        double m3 = 2 / (ex3 - as3*as3);
        double ny3 = m3*sig3*as3/2;
        double beta3 = sqrt(m3*sig3*sig3 - ny3*ny3);


        ui->a1DoubleSpinBox->setValue(a1);
        ui->m1DoubleSpinBox->setValue(m1);
        ui->ny1DoubleSpinBox->setValue(ny1);
        ui->beta1DoubleSpinBox->setValue(beta1);

        ui->a2DoubleSpinBox->setValue(a2);
        ui->m2DoubleSpinBox->setValue(m2);
        ui->ny2DoubleSpinBox->setValue(ny2);
        ui->beta2DoubleSpinBox->setValue(beta2);

        ui->a3DoubleSpinBox->setValue(a3);
        ui->m3DoubleSpinBox->setValue(m3);
        ui->ny3DoubleSpinBox->setValue(ny3);
        ui->beta3DoubleSpinBox->setValue(beta3);

        ui->aDoubleSpinBox->setValue(points[0].first);
        ui->x1DoubleSpinBox->setValue(points[bestResults.quantil1].first);
        ui->x2DoubleSpinBox->setValue(points[bestResults.quantil2].first);
        ui->bDoubleSpinBox->setValue(points.back().first);


        drawChos(ui->plotWidget);


//        double x1 = points[bestResults.quantil1].first;
//        double x2 = points[bestResults.quantil2].first;
//        double firstRss = 0.0;
//        double secRss = 0.0;
//        double thrRss = 0.0;

//        double rss = 0.0;
//        double khi = 0.0;
//        //double size = data ->getDistributionSize();
//        double integral = 0.0;
//        for(auto &p : data->getStepRelativePoints()) {
//            double value = 0.0;
//            double i = p.first;
//            double integralValue;
//            if (i < x1) {
//                value = ChosDistribution::value(i, m1, a1, beta1, ny1);

//                integralValue = Algorithms::Integral(i-step/2.0, i+step/2.0, [&](double d)->double{
//                    double v = ChosDistribution::value(d, m1, a1, beta1, ny1);
//                    return v;
//                });
//            }
//            else if (i < x2) {
//                value = ChosDistribution::value(i, m2, a2, beta2, ny2);

//                integralValue = Algorithms::Integral(i-step/2.0, i+step/2.0, [&](double d)->double{
//                    double v = ChosDistribution::value(d, m2, a2, beta2, ny2);
//  //                  qDebug() << "integer value 2" << v;
//                    return v;
//                });
//            }
//            else {
//                value = ChosDistribution::value(i, m3, a3, beta3, ny3);

//                integralValue = Algorithms::Integral(i-step/2.0, i+step/2.0, [&](double d)->double{
//                    double v = ChosDistribution::value(d, m3, a3, beta3, ny3);
//  //                  qDebug() << "integer value 3 ("<<d<<", "<<m3<<", "<<a3<<", "<<beta3<<", "<<ny3<<") = " << v;
//                    return v;
//                });
//            }
//            if (isnan(value)) {
//                break;
//            }
//            integral += integralValue;
//            double r = p.second - value;
//            rss += r*r;

//            if (i < x1) {
//              firstRss += r*r;
//            }
//            else if (i < x2) {
//                secRss += r*r;
//            }
//            else {
//                thrRss += r*r;
//            }

//            double sq = size * p.second * step;
//  //          qDebug() << "sq="<<sq;
//            integralValue *= size;
//  //          qDebug() << "integral=" << integralValue;
//            double chi = 0.0;
//            chi += integralValue - sq;
//            chi *= chi;
//            chi /= integralValue;
//            khi += chi;
//        }
       // fitDensity(firstRss, secRss, thrRss, integral, bestResults);

    }
}

fitResults MainWindow::fitDensity(double rss1, double rss2, double rss3, double realDensity, fitResults prevRes)
{
    auto points = data->getStepRelativePoints();
    double rssSum = rss1 + rss2 + rss3;

    double per1 = rss1/rssSum;
    double per2 = rss2/rssSum;
    double per3 = rss3/rssSum;

    double densityShift = 1 - realDensity;

    double shift1 = densityShift * per1;
    double shift2 = densityShift * per2;
    double shift3 = densityShift * per3;

    int q1 = prevRes.quantil1;
    int q2 = prevRes.quantil2;

    ChosDistribution distr1;
    ChosDistribution distr2;
    ChosDistribution distr3;
    PointsVector distr1Points(points.begin(), points.begin()+q1+1);
    PointsVector distr2Points(points.begin() + q1, points.begin() + q2+1);
    PointsVector distr3Points(points.begin() + q2, points.end());
    distr1.setPoints(distr1Points);
    distr2.setPoints(distr2Points);
    distr3.setPoints(distr3Points);


  //  std::vector<double> newRes1 = distr1.gradLinDens(prevRes.params[0], shift1, points[0].first, points[prevRes.quantil1].first );
  //  std::vector<double> newRes2 = distr2.gradLinDens(prevRes.params[0], shift2, points[prevRes.quantil1].first,  points[prevRes.quantil2].first);
  //  std::vector<double> newRes3 = distr3.gradLinDens(prevRes.params[0], shift3, points[prevRes.quantil2].first,  points.back().first);

//    qDebug()<<newRes1;
  //  qDebug()<<newRes2;
  //  qDebug()<<newRes3;
}

fitResults MainWindow::calcFit(PointsVector points, int quantilElem, int shakeCount)
{
    fitResults results;

    ChosDistribution distr1;
    PointsVector distr1Points(points.begin(), points.begin()+quantilElem+1);
    auto params1 = data->getDistributionParameters(0, quantilElem+1);
    distr1.setPoints(distr1Points);
    distr1.setInitialParams(params1);
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

    int from = quantilElem + 5;
    int to = size - 5;
    int step = 5;
    if (size <= 20)
    {
        from = quantilElem+1;
        to = size;
        step = 1;
    }
    if (!ui->autoQuantilCheckBox->isChecked())
    {
        from = ui->q2SpinBox->value();
        to = from+1;
    }

    for (int quantil2Elem = from; quantil2Elem < to; quantil2Elem += step)
    {
        ChosDistribution distr2;
        ChosDistribution distr3;

        PointsVector distr2Points(points.begin() + quantilElem, points.begin() + quantil2Elem+1);
        PointsVector distr3Points(points.begin() + quantil2Elem, points.end());

        auto params2 = data->getDistributionParameters(quantilElem, quantil2Elem+1);
        auto params3 = data->getDistributionParameters(quantil2Elem, size);

        distr2.setInitialParams({params2[0], params2[1], 0, 0.1});
        distr3.setInitialParams({params3[0], params3[1], 0, 0.1});

        distr2.setPoints(distr2Points);
        distr3.setPoints(distr3Points);

        auto matchParams2 = distr2.gradDescent(shakeCount);
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

void MainWindow::tableDataChanged(const QModelIndex &from, const QModelIndex &to) {
    tableDataChanged();
}

void MainWindow::tableDataChanged()
{
    PointsVector points;

    ui->plotWidget->clearPlottables();
    QLocale locale(QLocale::system());
    for (int i = 0; i < model->rowCount(); i++)
    {
        auto xItem = model->item(i, 0);
        auto yItem = model->item(i, 1);

        double x = locale.toDouble(xItem->text());
        double y = locale.toDouble(yItem->text());
        points.push_back(make_pair(x, y));
    }

    if (!data)
    {
        data = new DistributionData();
    }

    data->setPoints(points);
    QVector<double>x,y;
    auto d = data->getStepRelativePoints();
    for (Point p : d)
    {
        x.append(p.first);
        y.append(p.second);
    }
    createBar(ui->plotWidget, x, y, data->getStep());

    ui->plotWidget->replot();
    ui->plotWidget->replot();

}

void MainWindow::on_pushButton_clicked()
{
    auto selectionModel = ui->pointsTable->selectionModel();

    QStandardItem *xItem = new QStandardItem("");
    QStandardItem *yItem = new QStandardItem("");

    QList<QStandardItem*> itemList({xItem, yItem});

    if (selectionModel->hasSelection()) {
        int index = 0;
        QModelIndexList indexes = selectionModel->selectedIndexes();
        if (indexes.size() > 0)
        {
            index = indexes.last().row();
        }
        model->insertRow(index, itemList);
    }
    else {
        model->appendRow(itemList);
    }
}

void MainWindow::on_pushButton_2_clicked()
{
    auto selectionModel = ui->pointsTable->selectionModel();
    if (selectionModel->hasSelection()) {
        int index = 0;
        QModelIndexList indexes = selectionModel->selectedIndexes();
        if (indexes.size() > 0)
        {
            index = indexes.last().row();
        }
        model->removeRow(index);
    }
    tableDataChanged();
}

