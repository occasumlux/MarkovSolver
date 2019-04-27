#include "barchartwindow.hpp"
#include "ui_barchartwindow.h"

#include <QVector>

BarChartWindow::BarChartWindow(QWidget *parent, Eigen::MatrixXd *_pi_matrix_h) :
    QMainWindow(parent),
    ui(new Ui::BarChartWindow),
    pi_matrix_h(_pi_matrix_h),
    spinbox(new QSpinBox(this))
{
    ui->setupUi(this);

    if (pi_matrix_h != nullptr) {

        ui->toolBar->addAction("Back", this, &BarChartWindow::backButtonTriggered);

        //QSpinBox* spinbox = new QSpinBox(this);
        spinbox->setMinimum(0);
        spinbox->setMaximum(pi_matrix_h->rows() - 1);
        spinbox->setValue(0);
        ui->toolBar->addWidget(spinbox);
        QObject::connect(spinbox, QOverload<int>::of(&QSpinBox::valueChanged),
                         this, &BarChartWindow::updateChart);


        ui->toolBar->addAction("Next", this, &BarChartWindow::nextButtonTriggered);


        ui->barChart->clearPlottables();
        ui->barChart->xAxis->setRange(-0.5, pi_matrix_h->cols() - 0.5);
        ui->barChart->yAxis->setRange(0,1);

        //ui->barChart->setInteraction(QCP::iRangeDrag, true);
        //ui->barChart->setInteraction(QCP::iRangeZoom, true);

        //ui->barChart->axisRect()->setra

        QCPBars* prob_bars = new QCPBars(ui->barChart->xAxis, ui->barChart->yAxis);
        prob_bars->setName("Probabilities for each state");
        QVector<double> keyData;
        QVector<double> valueData;

        for (int i = 0; i < pi_matrix_h->cols(); ++i)
            keyData.append(i);
        for (int i = 0; i < pi_matrix_h->cols(); ++i)
            valueData.append(pi_matrix_h->coeff(0,i));

        prob_bars->setData(keyData, valueData, true);

        //ui->barChart->rescaleAxes();
        ui->barChart->replot();
    }
}

BarChartWindow::~BarChartWindow()
{
    delete ui;
}

void BarChartWindow::updateChart(int iter)
{
    ui->barChart->clearPlottables();
    ui->barChart->xAxis->setRange(-0.5, pi_matrix_h->cols() - 0.5);
    ui->barChart->yAxis->setRange(0,1);

    QCPBars* prob_bars = new QCPBars(ui->barChart->xAxis, ui->barChart->yAxis);
    prob_bars->setName("Probabilities for each state");
    QVector<double> keyData;
    QVector<double> valueData;

    for (int i = 0; i < pi_matrix_h->cols(); ++i)
        keyData.append(i);
    for (int i = 0; i < pi_matrix_h->cols(); ++i)
        valueData.append(pi_matrix_h->coeff(iter,i));

    prob_bars->setData(keyData, valueData, true);

    //ui->barChart->rescaleAxes();
    ui->barChart->replot();

}

void BarChartWindow::backButtonTriggered()
{
    int current_iter = spinbox->value();
    int next_iter = current_iter - 1;
    if (next_iter > -1) {
        spinbox->setValue(next_iter);
        updateChart(next_iter);
    }
}

void BarChartWindow::nextButtonTriggered()
{
    int current_iter = spinbox->value();
    int next_iter = current_iter + 1;
    if (next_iter < pi_matrix_h->rows()) {
        spinbox->setValue(next_iter);
        updateChart(next_iter);
    }
}
