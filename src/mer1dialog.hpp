#ifndef MER1DIALOG_HPP
#define MER1DIALOG_HPP

#include <QDialog>
#include <QPushButton>

namespace Ui {
class MEr1Dialog;
}

class MEr1Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit MEr1Dialog(QWidget *parent = nullptr);
    ~MEr1Dialog();
    int getR();
    int getK();
    double getLambda();
    double getMu();

private slots:
    void on_rSpinBox_valueChanged(int arg1);

    void on_kSpinBox_valueChanged(int arg1);

    void on_lambdaSpinBox_valueChanged(double arg1);

    void on_muSpinBox_valueChanged(double arg1);

private:
    Ui::MEr1Dialog *ui;
    int r;
    int k;
    double lambda;
    double mu;
};

#endif // MER1DIALOG_HPP
