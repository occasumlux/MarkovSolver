#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "tableview.hpp"
#include "barchartwindow.hpp"
#include "mmmdialog.hpp"
#include "mer1dialog.hpp"

#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QTableWidgetItem>
#include <QMessageBox>

#include <iostream>
//#include <iomanip> No longer used
//#include <limits> No longer used

#include<vector> // Used in M/G/1


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

    probs_matrix    = Eigen::MatrixXd::Zero(1,1);
    q_matrix        = Eigen::MatrixXd::Zero(1,1);
    pi_vector       = Eigen::MatrixXd::Zero(1,1);
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

                //std::cout<< std::setprecision(std::numeric_limits<double>::digits10 + 1) << value << std::endl;

                if (ok) {
                    if (ui->comboBox->currentText() == "CMPD") {
                        probs_matrix(line_count - 1, i) = value;
                    }
                    else {
                        q_matrix(line_count - 1, i) = value;
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
    ui->pi_table->setColumnCount(arg1);

    ui->probs_table->setRowCount(arg1);
    ui->probs_table->setColumnCount(arg1);

    pi_vector    = Eigen::MatrixXd::Zero(1, arg1);
    probs_matrix = Eigen::MatrixXd::Zero(arg1, arg1);
    q_matrix     = Eigen::MatrixXd::Zero(arg1, arg1);

    //std::cout << rows << "," << cols << std::endl;
}

void MainWindow::on_setButton_clicked()
{
    int dims = ui->spinBox->value();
    QString current_mode = ui->comboBox->currentText();
    //Eigen::MatrixXd probs_matrix(dims, dims);
    //Eigen::MatrixXd pi_vector(1, dims);

    //std::cout << probs_matrix.rows() << "x" << probs_matrix.cols() << std::endl;
    //std::cout << pi_vector.rows() << "x" << pi_vector.cols() << std::endl;

    for (int col = 0; col < dims; ++col) {
        for (int row = 0; row < dims; ++row) {
            //std::cout << row << "," << col << std::endl;
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
                        q_matrix(row, col) = value;
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

    ui->solveButton->setEnabled(true);
    /*
    std::cout << "Pi Vector:" << std::endl;
    std::cout << pi_vector << std::endl;

    std::cout << "Probs Matrix:" << std::endl;
    std::cout << probs_matrix << std::endl;

    std::cout << "Next Pi:" << std::endl;
    std::cout << pi_vector * probs_matrix << std::endl;
*/
}

void MainWindow::on_solveButton_clicked()
{
    // TODO: Include partition size in GUI
    // Or a way to calculate more than one pi(t)
    int partitions = 10;
    int iters = ui->spinBox_2->value();
    int dims = ui->spinBox->value();
    QString current_mode = ui->comboBox->currentText();
    double t = ui->t_SpinBox->value();
    double lambda_max = 1;

    Eigen::MatrixXd* pi_matrix_h = new Eigen::MatrixXd(iters, dims);
    *pi_matrix_h = Eigen::MatrixXd::Zero(iters, dims);
    pi_matrix_h->row(0) = pi_vector;

    if (current_mode == "CMPC") {
        //Find max in diag
        lambda_max = q_matrix.diagonal().array().abs().maxCoeff();
        probs_matrix = Eigen::MatrixXd::Identity(dims, dims) + q_matrix / lambda_max;
    }

    double lambda_t = lambda_max * t;
    if (current_mode == "CMPC") {
        // Compute minimum n for desired error level
        // TODO: Put error level in GUI
        int n = n_for_error(0.000001, lambda_t);
        if (n > iters) {
            iters = n;
            ui->spinBox_2->setValue(iters);
        }
        //iters = (n > iters) ? n : iters;
    }

    // C.K eq
    for (int i = 1; i < iters; ++i) {
        pi_matrix_h->row(i) = pi_matrix_h->row(i - 1) * probs_matrix;
    }

    //std::cout << pi_matrix_h << std::endl;

    //std::cout << ui->comboBox->currentText().toStdString() << std::endl;

    TableView* tab = new TableView(this, pi_matrix_h);
    tab->setWindowTitle("CMPD Transition Values");

    BarChartWindow* chart = new BarChartWindow(this, pi_matrix_h);
    chart->setWindowTitle("CMPD Transition Charts");

    if (current_mode == "CMPC") {
        Eigen::MatrixXd* pi_f = new Eigen::MatrixXd(partitions, dims);
        *pi_f = Eigen::MatrixXd::Zero(partitions, dims);
        for (int i = 0; i < iters; ++i) {
            for (int k = 0; k < partitions; ++k) {
                pi_f->row(k) += pi_matrix_h->row(i) * poisson(i, lambda_t * (k + 1)/partitions);
            }
        }
        //std::cout << *pi_f << std::endl;

        TableView* tab_2 = new TableView(this, pi_f);
        tab_2->setWindowTitle("CMPC Final Value");
        tab_2->show();
        BarChartWindow* chart_2 = new BarChartWindow(this, pi_f);
        chart_2->setWindowTitle("CMPC Result Chart");
        chart_2->show();
    }
    tab->show();
    chart->show();
    //std::cout << probs_matrix << std::endl;
}

// DEPRECATED - Now pi_matrix_h is a pointer
//void MainWindow::on_spinBox_2_valueChanged(int arg1)
//{
//    int dims = ui->spinBox->value();
//    pi_matrix_h = Eigen::MatrixXd::Zero(arg1, dims);
//    pi_matrix_h.row(0) = pi_vector;
//}

/*
void MainWindow::on_pi_table_cellChanged(int row, int column)
{
    ui->solveButton->setEnabled(false);
}

void MainWindow::on_probs_table_cellChanged(int row, int column)
{
    ui->solveButton->setEnabled(false);
}
*/

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

double MainWindow::poisson(double n, double lambda)
{
    return std::exp(n * std::log(lambda) - std::lgamma(n + 1.0) - lambda);
}

int MainWindow::n_for_error(double error, double lambda_t)
{
    int n = 0;
    double cumulative_prob = 0;
    while(cumulative_prob < 1 - error) {
        cumulative_prob += poisson(n, lambda_t);
        ++n;
    }
    return n;
}

int MainWindow::factorial(int n)
{
    return (n == 0 || n == 1) ? 1 : factorial(n - 1) * n;
}

void MainWindow::on_actionM_M_m_triggered()
{
    // Open new window
    MMmDialog dialog(this);
    if ( dialog.exec() == QDialog::Accepted) {
        int m = dialog.getM();
        int k = dialog.getK();
        double lambda = dialog.getLambda();
        double mu = dialog.getMu();

        double rho = lambda / (m * mu);
        double l_mu = lambda / mu;

        // Compute probabilities
        Eigen::MatrixXd* pi_vector_queue = new Eigen::MatrixXd(k + 1, 1);
        *pi_vector_queue = Eigen::MatrixXd::Zero(k + 1, 1);
        pi_vector_queue->coeffRef(0,0) = 0;
        // pi_0 calculation
        pi_vector_queue->coeffRef(0,0) += 1;
        for (int i = 1; i <= m - 1; ++i) {
            pi_vector_queue->coeffRef(0,0) += std::pow(m*rho, i) / factorial(i);
         }
        pi_vector_queue->coeffRef(0,0) += std::pow(m * rho, m) / ( factorial(m) * (1 - rho) );
        pi_vector_queue->coeffRef(0,0) = 1 / pi_vector_queue->coeffRef(0,0);

        // pi_k calculation
        for (int i = 1; i <= k; ++i) {
            if (i < m) {
                pi_vector_queue->coeffRef(i,0) = pi_vector_queue->coeffRef(i-1,0) * l_mu / (i);
            }
            else {
                pi_vector_queue->coeffRef(i,0) = pi_vector_queue->coeffRef(i-1,0) * rho;
            }
        }

        // Show graph (?)
        TableView* tab = new TableView(this, pi_vector_queue);
        tab->setWindowTitle("M/M/m values");
        tab->show();
        BarChartWindow* chart = new BarChartWindow(this, pi_vector_queue);
        chart->setWindowTitle("M/M/m Result Chart");
        chart->show();

   }
}

void MainWindow::on_action_M_Er_1_triggered()
{
    // Open new window
    MEr1Dialog dialog(this);
    if ( dialog.exec() == QDialog::Accepted) {
        int r = dialog.getR();
        int k = dialog.getK();
        double lambda = dialog.getLambda();
        double mu = dialog.getMu();

        double rho = lambda / (r * mu);
        double rho_0 = lambda / mu;

        // Stages to calculate
        int stages = k * r + 1;

        // Compute probabilities
        // Result vector, for number of users
        Eigen::MatrixXd* pi_vector_queue = new Eigen::MatrixXd(k + 1, 1);
        *pi_vector_queue = Eigen::MatrixXd::Zero(k + 1, 1);
        // Temporary vector, for number of remaining stages
        Eigen::MatrixXd pi_vector_queue_t = Eigen::MatrixXd::Zero(stages + 1, 1);
        // pi_0 calculation (both)
        pi_vector_queue->coeffRef(0,0) = 1 - rho_0;
        pi_vector_queue_t(0,0) = 1 - rho_0;
        // pi_i calculation (stages)
        for (int i = 1; i <= stages; ++i) {
            int j = (i - r < 0) ? 0 : i - r;
            pi_vector_queue_t(i,0) = rho * pi_vector_queue_t.block(j, 0, i-j, 1).sum(); // i-1-j (Equation) +1 (Block asks size, not position)
        }

        //pi_k calculation (users)
        for(int i = 1; i <= k; ++i) {
            pi_vector_queue->coeffRef(i,0) = pi_vector_queue_t(i*r + 1, 0) / rho;
        }

        // Show graph
        TableView* tab = new TableView(this, pi_vector_queue);
        tab->setWindowTitle("M/Er/1 values");
        tab->show();
        BarChartWindow* chart = new BarChartWindow(this, pi_vector_queue);
        chart->setWindowTitle("M/Er/1 Result Chart");
        chart->show();

   }
}

void MainWindow::on_actionM_G_1_triggered()
{
    std::vector<double> akVector;

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
    double lambda = 0;
    double e_x = 0;
    int iters = 0;

    while ( !stream.atEnd() ) {
        ++real_count;
        QString line = stream.readLine();
        if (line.startsWith("#") || line == "")
            continue;

        if (line_count == -2) {
            // Parse Poisson lambda
            bool ok = false;
            lambda = line.simplified().toDouble(&ok);

            if (!ok || lambda <= 0) {
                std::cerr << "Invalid Poisson lambda!!!" << std::endl;
                QMessageBox::information(this, tr("Invalid Poisson lambda"),
                                         "Check Lambda in line #" +
                                         QString::number(real_count));
                return;
            }
        }
        else if (line_count == -1) {
            // Parse Distribution E[x]
            bool ok = false;
            e_x = line.simplified().toDouble(&ok);

            if (!ok || e_x <= 0) {
                std::cerr << "Invalid E[x]!!!" << std::endl;
                QMessageBox::information(this, tr("Invalid E[x]"),
                                         "Check E[x] in line #" +
                                         QString::number(real_count));
                return;
            }
        }
        else if (line_count == 0) {
            // Parse number of states to compute
            bool ok = false;
            iters = line.simplified().toDouble(&ok);

            if (!ok || iters <= 0) {
                std::cerr << "Invalid number of states!!!" << std::endl;
                QMessageBox::information(this, tr("Invalid number of states"),
                                         "Check line #" +
                                         QString::number(real_count));
                return;
            }
        }

        else {
            // Parse ak
            bool ok = false;
            double ak = line.simplified().toDouble(&ok);

            if (!ok || ak <= 0) {
                std::cerr << "Invalid ak!!!" << std::endl;
                QMessageBox::information(this, tr("Invalid ak"),
                                         "Check ak in line #" +
                                         QString::number(real_count));
                return;
            }

            // Add ak to list
            akVector.push_back(ak);
        }
        ++line_count;
    }

    // Check if it's possible to compute the required state
    if (akVector.size() + 1 < iters) {
        std::cerr << "Too few ak values!!!" << std::endl;
        QMessageBox::information(this, tr("Too few ak values"),
                                 "Please add more data to the file");
        return;
    }

    // Ak calculation
    std::vector<double> AkVector(akVector.size());
    AkVector[0] = akVector[0];
    for(size_t k = 1; k < akVector.size(); ++k ) {
        AkVector[k] = AkVector[k-1] + akVector[k];
    }

    double ak_0 = akVector[0];

    //nAk calculation - Reutilization of akVector variable
    for(size_t k = 0; k < AkVector.size(); ++k ) {
        akVector[k] = 1 - AkVector[k];
    }
    AkVector.clear();

    // pi vector definition
    Eigen::MatrixXd* pi_vector_queue = new Eigen::MatrixXd(iters + 1, 1);
    *pi_vector_queue = Eigen::MatrixXd::Zero(iters + 1, 1);
    pi_vector_queue->coeffRef(0,0) = 1 - lambda * e_x;

    //pi vector calculation
    for(int n = 1; n <= iters; ++n){
        pi_vector_queue->coeffRef(n, 0) = pi_vector_queue->coeffRef(0,0)*(akVector[n-1] / ak_0);

        if (n - 1 > 0) {
            // Second part
            for(int i = 1; i <= n - 1; ++i) {
                pi_vector_queue->coeffRef(n,0) += pi_vector_queue->coeffRef(i,0) * (akVector[n-i] / ak_0);
            }
        }
    }

    // Typical graph show
    TableView* tab = new TableView(this, pi_vector_queue);
    tab->setWindowTitle("M/G/1 values");
    tab->show();
    BarChartWindow* chart = new BarChartWindow(this, pi_vector_queue);
    chart->setWindowTitle("M/G/1 Result Chart");
    chart->show();

}
