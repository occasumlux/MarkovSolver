#ifndef TABLEVIEW_HPP
#define TABLEVIEW_HPP

#include <QMainWindow>

#include "../eigen3/Eigen/Dense"

namespace Ui {
class TableView;
}

class TableView : public QMainWindow
{
    Q_OBJECT

public:
    explicit TableView(QWidget *parent = nullptr, const Eigen::MatrixXd* pi_matrix_h = nullptr);
    ~TableView();

private slots:
    void saveCSV();

private:
    Ui::TableView *ui;
};

#endif // TABLEVIEW_HPP
