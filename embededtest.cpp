#include "embededtest.h"
#include "ui_embededtest.h"
#include "mainwindow.h"
embededTest::embededTest(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::embededTest)
{
    ui->setupUi(this);

    MainWindow *m = new MainWindow(this);
    ui->gridLayout->addWidget(m);
}

embededTest::~embededTest()
{
    delete ui;
}

void embededTest::on_pushButton_clicked()
{
    QList<double> list1({1, 2, 3, 4, 5});
    QList<double> list2({2, 3, 4, 5, 6});

    double summ1 = 0.0;
    foreach (double e1, list1) {
        summ1 += e1;
    }

    double mean1 = summ1 / list1.size();

    double summ2 = 0.0;
    foreach (double e2, list2) {
        summ2 += e2;
    }

    double mean2 = summ2 / list2.size();


    double covSumm = 0.0;
    for (int i = 0; i < list1.size(); ++i) {
        covSumm += (list1[i] - mean1) * (list2[i] - mean2);
    }
    double cov = covSumm / list1.size();
    qDebug() << "cov = " << cov;
}
