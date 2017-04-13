#ifndef GRAPHWINDOW_H
#define GRAPHWINDOW_H

#include <QWidget>

#include <vector>

namespace Ui {
class GraphWindow;
}

class GraphWindow : public QWidget
{
    Q_OBJECT

public:
    explicit GraphWindow(QWidget *parent = 0);
    ~GraphWindow();

    void setPointsVector(std::vector<std::pair<double, double>> points);
private:
    Ui::GraphWindow *ui;


};

#endif // GRAPHWINDOW_H
