#ifndef BARCHARTWINDOW_HPP
#define BARCHARTWINDOW_HPP

#include <QMainWindow>
#include <QSpinBox>

#include "../eigen3/Eigen/Dense"

namespace Ui {
class BarChartWindow;
}

class BarChartWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit BarChartWindow(QWidget *parent = nullptr, Eigen::MatrixXd* _pi_matrix_h = nullptr);
    ~BarChartWindow();

private slots:
    void updateChart(int iter);
    void backButtonTriggered();
    void nextButtonTriggered();

private:
    Ui::BarChartWindow *ui;
    Eigen::MatrixXd* pi_matrix_h;
    QSpinBox* spinbox;
};

#endif // BARCHARTWINDOW_HPP
