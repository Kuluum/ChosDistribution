#ifndef DESCENTPROGRESSVIEW_H
#define DESCENTPROGRESSVIEW_H

#include <QMainWindow>
#include <vector>

namespace Ui {
class DescentProgressView;
}

class DescentProgressView : public QMainWindow
{
    Q_OBJECT

public:
    explicit DescentProgressView(QWidget *parent = 0);
    void setDescentProgress(std::vector<double> descentProgress);
    ~DescentProgressView();

private:
    Ui::DescentProgressView *ui;
    std::vector<double> descenntProgress;

    void drawProgress();
};

#endif // DESCENTPROGRESSVIEW_H
