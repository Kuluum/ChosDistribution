#include "graphwindow.h"
#include "ui_graphwindow.h"

GraphWindow::GraphWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::GraphWindow)
{
    ui->setupUi(this);
    ui->plotWidget->setInteraction(QCP::iRangeDrag, true);
    ui->plotWidget->setInteraction(QCP::iRangeZoom, true);

    ui->plotWidget->xAxis->setLabel("x");
    ui->plotWidget->yAxis->setLabel("y");
}

GraphWindow::~GraphWindow()
{
    delete ui;
}

void GraphWindow::setPointsVector(std::vector<std::pair<double, double> > points)
{
    ui->plotWidget->xAxis->setRange(points[0].first, points.back().first);
    ui->plotWidget->yAxis->setRange(0, 5);

    QVector<double> x, y;
    for (auto p: points)
    {
        x.append(p.first);
        y.append(p.second);
    }

    ui->plotWidget->addGraph();
    ui->plotWidget->graph(0)->setData(x, y);
}
