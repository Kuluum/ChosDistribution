#include "descentprogressview.h"
#include "ui_descentprogressview.h"

DescentProgressView::DescentProgressView(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::DescentProgressView)
{
    ui->setupUi(this);
    ui->plotWidget->setInteraction(QCP::iRangeDrag, true);
    ui->plotWidget->setInteraction(QCP::iRangeZoom, true);

    ui->plotWidget->xAxis->setLabel("x");
    ui->plotWidget->yAxis->setLabel("y");

//    ui->plotWidget->xAxis->setRange(0, 5);

}

DescentProgressView::~DescentProgressView()
{
    delete ui;
}

void DescentProgressView::setDescentProgress(std::vector<double> descentProgress)
{
    this->descenntProgress = descentProgress;
    drawProgress();
}

void DescentProgressView::drawProgress()
{
    ui->plotWidget->xAxis->setRange(0, descenntProgress.size());
    ui->plotWidget->yAxis->setRange(0, descenntProgress[0]);

    QVector<double> x, y;
    for (int i = 0; i < descenntProgress.size(); i++)
    {
        x.append(i);
        y.append(descenntProgress[i]);
    }

    ui->plotWidget->addGraph();
    ui->plotWidget->graph(0)->setData(x, y);
}
