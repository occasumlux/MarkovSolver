#include "tableview.hpp"
#include "ui_tableview.h"

#include <QPushButton>
#include <QKeySequence>
#include <QFileDialog>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>

#include <iostream>
//#include <iomanip>
//#include <limits>
#include <QDebug>

TableView::TableView(QWidget *parent, const Eigen::MatrixXd* pi_matrix_h) :
    QMainWindow(parent),
    ui(new Ui::TableView)
{
    ui->setupUi(this);

    QPushButton* button = new QPushButton("Save", this);
    button->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_S));
    ui->menubar->setCornerWidget(button, Qt::TopRightCorner);
    QObject::connect(button, &QPushButton::clicked,
                     this, &TableView::saveCSV);

    int cols = pi_matrix_h->cols();
    int rows = pi_matrix_h->rows();

    //std::cout << rows << "," << cols << std::endl;

    ui->tableWidget->setColumnCount(cols);
    ui->tableWidget->setRowCount(rows);

    for (int col = 0; col < cols; ++col) {
        for (int row = 0; row < rows; ++row) {
            //QString value = QString("%1").arg(pi_matrix_h(row, col));
            QString value2 = QString::number(pi_matrix_h->coeff(row, col), 'f', 10); // std::numeric_limits<double>::digits10 + 1);
            //qDebug() << value << "\n";
            //qDebug() << value2 << "\n";
            QStringList splitted = value2.split('.');
            bool ok = false;
            int post_dot = splitted.at(1).toInt(&ok);
            if (ok && post_dot == 0) {
                value2 = splitted.at(0);
                //qDebug() << splitted << "\n";
                //qDebug() << value2 << "\n";
            }
            QTableWidgetItem* item = new QTableWidgetItem(value2);
            ui->tableWidget->setItem(row, col, item);
        }
    }
    //double va = 0.2499999999;
    //QString vv = QString::number(va, 'f', 10);
    //qDebug() << vv.split('.') << "\n";
}

TableView::~TableView()
{
    delete ui;
}

void TableView::saveCSV()
{
    std::cout << "Clicked" << std::endl;
    QString filename = QFileDialog::getSaveFileName(this, tr("Save as CSV file"), "", tr("CSV (*.csv);; Text file (*.txt)"));

    if (filename.isEmpty())
        return;
    else {
        QFile file(filename);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QMessageBox::information(this, tr("Unable to open file"),
                                     file.errorString());
            return;
        }
        QTextStream out(&file);
        out << "# Pi vector values, with each line containing one" << endl;
        out << "# Enjoy!" << endl;
        for (int row = 0; row < ui->tableWidget->rowCount(); ++row) {
            for (int col = 0; col < ui->tableWidget->columnCount(); ++col) {
                out << ui->tableWidget->item(row, col)->text() << ",";
            }
            out << endl;
        }
    }
}
