#include "tableview.hpp"
#include "ui_tableview.h"

#include <QPushButton>
#include <QKeySequence>
#include <QFileDialog>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>

TableView::TableView(QWidget *parent, Eigen::MatrixXd pi_matrix_h) :
    QMainWindow(parent),
    ui(new Ui::TableView)
{
    ui->setupUi(this);

    QPushButton* button = new QPushButton("Save", this);
    button->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_S));
    ui->menubar->setCornerWidget(button, Qt::TopRightCorner);
    QObject::connect(button, &QPushButton::clicked,
                     this, &TableView::saveCSV);

    int cols = pi_matrix_h.cols();
    int rows = pi_matrix_h.rows();

    ui->tableWidget->setColumnCount(cols);
    ui->tableWidget->setRowCount(rows);

    for (int col = 0; col < cols; ++col) {
        for (int row = 0; row < rows; ++row) {
            QString value2 = QString::number(pi_matrix_h(row, col), 'f', 10);
            QStringList splitted = value2.split('.');
            bool ok = false;
            int post_dot = splitted.at(1).toInt(&ok);
            if (ok && post_dot == 0) {
                value2 = splitted.at(0);
            }
            QTableWidgetItem* item = new QTableWidgetItem(value2);
            ui->tableWidget->setItem(row, col, item);
        }
    }
}

TableView::~TableView()
{
    delete ui;
}

void TableView::saveCSV()
{
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
