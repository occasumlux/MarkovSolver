#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "tableview.hpp"
#include "barchartwindow.hpp"

#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QTableWidgetItem>
#include <QMessageBox>

#include <iostream>
#include <iomanip>
#include <limits>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->comboBox->addItem("CMPD");
    ui->comboBox->addItem("CMPC");

    ui->pi_table->setRowCount(1);
    ui->pi_table->setColumnCount(1);
    QTableWidgetItem* pi_item = new QTableWidgetItem("0");
    ui->pi_table->setItem(0, 0, pi_item);

    ui->probs_table->setRowCount(1);
    ui->probs_table->setColumnCount(1);
    QTableWidgetItem* probs_item = new QTableWidgetItem("0");
    ui->probs_table->setItem(0, 0, probs_item);

    probs_matrix = Eigen::MatrixXd::Zero(1,1);
    pi_vector    = Eigen::MatrixXd::Zero(1,1);
    pi_matrix_h  = Eigen::MatrixXd::Zero(1,1);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::loadFile() {
    QString filename = QFileDialog::getOpenFileName(this, tr("Open CSV file"), "./", tr("CSV (*.csv);; Text file (*.txt)"));
    if (filename == "")
        return;

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) {
        std::cerr << file.errorString().toStdString() << std::endl;
        return;
    }
    QTextStream stream(&file);
    int line_count = -2;
    int real_count = 0;
    int dims = 0;
    int iters = 0;
    while (!stream.atEnd() && line_count <= dims) {
        ++real_count;
        QString line = stream.readLine();
        if (line.startsWith("#") || line == "")
            continue;
        QStringList splitted = line.split(',');
        if (line_count == -2) {
            if (line.compare("CMPD", Qt::CaseInsensitive) == 0) {
                ui->comboBox->setCurrentText("CMPD");
            }
            else if (line.compare("CMPC", Qt::CaseInsensitive) == 0) {
                ui->comboBox->setCurrentText("CMPC");
            }
            else {
                std::cerr << "Invalid mode!!!" << std::endl;
                QMessageBox::information(this, tr("Invalid Mode"),
                                         "Check Mode in line #" +
                                         QString::number(real_count));
                return;
            }
        }
        else if (line_count == -1) {
            //parse dim and iters
            if (splitted.length() >= 2) {
                bool ok = false;
                dims = splitted.at(0).toInt(&ok);
                if (ok) {
                    ui->spinBox->setValue(dims);
                }
                else {
                    QMessageBox::information(this, tr("Non number found"),
                                             "Check number of states in line #" +
                                             QString::number(real_count));
                    return;
                }

                ok = false;
                iters = splitted.at(1).toInt(&ok);
                if (ok) {
                    ui->spinBox_2->setValue(iters);
                }
                else {
                    QMessageBox::information(this, tr("Non number found"),
                                             "Check iterations in line #" +
                                             QString::number(real_count));
                    return;
                }
                if (ui->comboBox->currentText() == "CMPC" && splitted.length() >= 3) {
                    //Parse t
                    ok = false;
                    double t = splitted.at(2).toDouble(&ok);
                    if (ok) {
                        ui->t_SpinBox->setValue(t);
                    }
                    else {
                        QMessageBox::information(this, tr("Non number found"),
                                                 "Check t in line #" +
                                                 QString::number(real_count));
                        return;
                    }
                }
            }
        }
        else if (line_count == 0) {
            if (splitted.length() < dims) {
                std::cerr << "Error parsing file" << std::endl;
                QMessageBox::information(this, tr("Missing column(s)"),
                                         "Check data in line #" +
                                         QString::number(real_count));
                return;
            }
            for (int i = 0; i < dims; ++i) {
                bool ok = false;
                double value = splitted.at(i).toDouble(&ok);

                //std::cout << value << std::endl;

                if (ok) {
                    pi_vector(0,i) = value;
                    QTableWidgetItem* item = new QTableWidgetItem(splitted.at(i));
                    ui->pi_table->setItem(0, i, item);
                }
                else {
                    std::cerr << "Non double found!" << std::endl;
                    QMessageBox::information(this, tr("Non double found"),
                                             "Check numbers in line #" +
                                             QString::number(real_count));
                    return;
                }
            }
            pi_matrix_h.row(0) = pi_vector;
        }

        else {
            if (splitted.length() < dims) {
                std::cerr << "Error parsing file" << std::endl;
                QMessageBox::information(this, tr("Missing column(s)"),
                                         "Check data in line #" +
                                         QString::number(real_count));
                return;
            }
            for (int i = 0; i < dims; ++i) {
                bool ok = false;
                double value = splitted.at(i).toDouble(&ok);

                if (ok) {
                    if (ui->comboBox->currentText() == "CMPD") {
                        probs_matrix(line_count - 1, i) = value;
                    }
                    else {
                        // Fill Q Matrix
                    }

                    QTableWidgetItem* item = new QTableWidgetItem(splitted.at(i));
                    //ui->probs_table->blockSignals(true);
                    ui->probs_table->setItem(line_count - 1, i, item);
                    //ui->probs_table->blockSignals(false);
                }
                else {
                    std::cerr << "Non double found!" << std::endl;
                    QMessageBox::information(this, tr("Non double found"),
                                             "Check numbers in line #" +
                                             QString::number(real_count));
                    return;
                }
            }
        }
        //wordList.append(line.split(',').first());
        ++line_count;
    }
    ui->solveButton->setEnabled(true);

}

void MainWindow::on_spinBox_valueChanged(int arg1)
{
    int iters = ui->spinBox_2->value();

    ui->pi_table->setColumnCount(arg1);

    ui->probs_table->setRowCount(arg1);
    ui->probs_table->setColumnCount(arg1);

    pi_vector    = Eigen::MatrixXd::Zero(1, arg1);
    pi_matrix_h  = Eigen::MatrixXd::Zero(iters, arg1);
    probs_matrix = Eigen::MatrixXd::Zero(arg1, arg1);
}

void MainWindow::on_setButton_clicked()
{
    int dims = ui->spinBox->value();
    QString current_mode = ui->comboBox->currentText();

    for (int col = 0; col < dims; ++col) {
        for (int row = 0; row < dims; ++row) {
            //Iterate Qt Table
            bool ok = false;
            QTableWidgetItem* item = ui->probs_table->item(row, col);
            if (item != nullptr) {
                QString text = item->text();
                double value = text.toDouble(&ok);
                if (ok){
                    if (current_mode == "CMPD") {
                        probs_matrix(row, col) = value;
                    }
                    else {
                        // Fill Q Matrix
                    }
                }
                else {
                    ui->solveButton->setEnabled(false);
                    QMessageBox::information(this, tr("Non number found"),
                                             "Check (Row, Column): (" +
                                             QString::number(row + 1) + "," +
                                             QString::number(col + 1) + ")");
                    return;
                }
            }
            else {
                ui->solveButton->setEnabled(false);
                QMessageBox::information(this, tr("Missing data"),
                                         "Check (Row, Column): (" +
                                         QString::number(row + 1) + "," +
                                         QString::number(col + 1) + ")");
                return;
            }
        }
        bool ok = false;
        QTableWidgetItem* pi_item = ui->pi_table->item(0, col);
        if (pi_item != nullptr) {
            QString text = pi_item->text();
            double value = text.toDouble(&ok);
            if (ok){
                pi_vector(0, col) = value;
            }
            else {
                ui->solveButton->setEnabled(false);
                QMessageBox::information(this, tr("Non number found"),
                                         "Check (Row, Column): (" +
                                         QString::number(0 + 1) + "," +
                                         QString::number(col + 1) + ")");
                return;
            }
        }
        else {
            ui->solveButton->setEnabled(false);
            QMessageBox::information(this, tr("Missing data"),
                                     "Check (Row, Column): (" +
                                     QString::number(0 + 1) + "," +
                                     QString::number(col + 1) + ")");
            return;
        }
    }

    pi_matrix_h.row(0) = pi_vector;

    ui->solveButton->setEnabled(true);
}

void MainWindow::on_solveButton_clicked()
{

    int iters = ui->spinBox_2->value();
    QString current_mode = ui->comboBox->currentText();

    if (current_mode == "CMPC") {
        // Find max in diag
        // Then populate probs_matrix
    }

    for (int i = 1; i < iters; ++i) {
        pi_matrix_h.row(i) = pi_matrix_h.row(i - 1) * probs_matrix;
    }

    TableView* tab = new TableView(this, pi_matrix_h);
    tab->setWindowTitle("CMPD Transition Values");

    BarChartWindow* chart = new BarChartWindow(this, &pi_matrix_h);
    chart->setWindowTitle("CMPD Transition Charts");

    if (current_mode == "CMPC") {
        // Iterate over pi_matrix_h and compute Sum[Poisson(n, lambda) * Pi(n)]
        // Then show the final result
    }
    tab->show();
    chart->show();
}

void MainWindow::on_spinBox_2_valueChanged(int arg1)
{
    int dims = ui->spinBox->value();
    pi_matrix_h = Eigen::MatrixXd::Zero(arg1, dims);
    pi_matrix_h.row(0) = pi_vector;
}

void MainWindow::on_comboBox_currentIndexChanged(const QString &arg1)
{
    ui->solveButton->setDisabled(true);
    if (arg1 == "CMPC") {
        ui->t_SpinBox->setEnabled(true);
    }
    else {
        ui->t_SpinBox->setDisabled(true);
    }
}
