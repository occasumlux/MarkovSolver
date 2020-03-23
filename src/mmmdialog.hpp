#ifndef MMMDIALOG_HPP
#define MMMDIALOG_HPP

#include <QDialog>
#include <QPushButton>

namespace Ui {
class MMmDialog;
}

class MMmDialog : public QDialog
{
    Q_OBJECT

public:
    explicit MMmDialog(QWidget *parent = nullptr);
    ~MMmDialog();
    int getM();
    int getK();
    double getLambda();
    double getMu();

private slots:
    void on_mSpinBox_valueChanged(int arg1);

    void on_kSpinBox_valueChanged(int arg1);

    void on_lambdaSpinBox_valueChanged(double arg1);

    void on_muSpinBox_valueChanged(double arg1);

private:
    Ui::MMmDialog *ui;
    int m;
    int k;
    double lambda;
    double mu;
};

#endif // MMMDIALOG_HPP
