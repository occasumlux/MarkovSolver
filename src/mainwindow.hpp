#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>

#include "../eigen3/Eigen/Dense"
#include <cmath>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void loadFile();
    void on_spinBox_valueChanged(int arg1);
    void on_setButton_clicked();
    void on_solveButton_clicked();
    //void on_spinBox_2_valueChanged(int arg1);
/*
    void on_pi_table_cellChanged(int row, int column);
    void on_probs_table_cellChanged(int row, int column);
*/

    void on_comboBox_currentIndexChanged(const QString &arg1);

    void on_actionM_M_m_triggered();

    void on_action_M_Er_1_triggered();

    void on_actionM_G_1_triggered();

private:
    double poisson(double n, double lambda);
    int n_for_error(double error, double lambda_t);
    int factorial(int n);
    Ui::MainWindow *ui;
    Eigen::MatrixXd probs_matrix;
    Eigen::MatrixXd q_matrix;
    Eigen::MatrixXd pi_vector;
};

#endif // MAINWINDOW_HPP
