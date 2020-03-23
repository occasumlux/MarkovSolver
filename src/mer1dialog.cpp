#include "mer1dialog.hpp"
#include "ui_mer1dialog.h"

MEr1Dialog::MEr1Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MEr1Dialog),
    r(1), k(1), lambda(1), mu(1)
{
    ui->setupUi(this);
    this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);

}

MEr1Dialog::~MEr1Dialog()
{
    delete ui;
}

int MEr1Dialog::getR()
{
    return r;
}

int MEr1Dialog::getK()
{
    return k;
}

double MEr1Dialog::getLambda()
{
    return lambda;
}

double MEr1Dialog::getMu()
{
    return mu;
}

void MEr1Dialog::on_rSpinBox_valueChanged(int arg1)
{
    if ( arg1 != r) {
        r = arg1;
    }
}

void MEr1Dialog::on_kSpinBox_valueChanged(int arg1)
{
    if ( arg1 != k) {
        k = arg1;
    }
}

void MEr1Dialog::on_lambdaSpinBox_valueChanged(double arg1)
{
    lambda = arg1;
    if (lambda >= mu) {
        this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);
    }
    else {
        this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);
    }
}

void MEr1Dialog::on_muSpinBox_valueChanged(double arg1)
{
    mu = arg1;
    if (lambda >= mu) {
        this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);
    }
    else {
        this->ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);
    }
}
