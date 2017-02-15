#ifndef EMBEDEDTEST_H
#define EMBEDEDTEST_H

#include <QWidget>

namespace Ui {
class embededTest;
}

class embededTest : public QWidget
{
    Q_OBJECT

public:
    explicit embededTest(QWidget *parent = 0);
    ~embededTest();

private slots:
    void on_pushButton_clicked();

private:
    Ui::embededTest *ui;
};

#endif // EMBEDEDTEST_H
