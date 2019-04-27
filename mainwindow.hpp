#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>

#include "eigen3/Eigen/Dense"
#include <cmath>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void loadFile();
    void on_spinBox_valueChanged(int arg1);
    void on_setButton_clicked();
    void on_solveButton_clicked();
    void on_spinBox_2_valueChanged(int arg1);
/*
    void on_pi_table_cellChanged(int row, int column);
    void on_probs_table_cellChanged(int row, int column);
*/

    void on_comboBox_currentIndexChanged(const QString &arg1);

private:
    double poisson(double n, double lambda);
    Ui::MainWindow *ui;
    Eigen::MatrixXd probs_matrix;
    Eigen::MatrixXd q_matrix;
    Eigen::MatrixXd pi_f;
    Eigen::MatrixXd pi_vector;
    Eigen::MatrixXd pi_matrix_h;
};

#endif // MAINWINDOW_HPP
