#include "mmmdialog.hpp"
#include "ui_mmmdialog.h"

MMmDialog::MMmDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MMmDialog),
    m(1), k(1), lambda(1), mu(1)
{
    ui->setupUi(this);
    this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);


}

MMmDialog::~MMmDialog()
{
    delete ui;
}

int MMmDialog::getM()
{
    return m;
}

int MMmDialog::getK()
{
    return k;
}

double MMmDialog::getLambda()
{
    return lambda;
}

double MMmDialog::getMu()
{
    return mu;
}

void MMmDialog::on_mSpinBox_valueChanged(int arg1)
{
    if ( arg1 != m) {
        m = arg1;
    }
}

void MMmDialog::on_kSpinBox_valueChanged(int arg1)
{
    if ( arg1 != k) {
        k = arg1;
    }
}

void MMmDialog::on_lambdaSpinBox_valueChanged(double arg1)
{
    lambda = arg1;
    if (lambda >= mu * m) {
        this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);
    }
    else {
        this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);
    }
}

void MMmDialog::on_muSpinBox_valueChanged(double arg1)
{
    mu = arg1;
    if (lambda >= mu * m) {
        this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);
    }
    else {
        this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);
    }
}
