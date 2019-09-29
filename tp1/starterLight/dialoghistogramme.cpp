#include "dialoghistogramme.h"


#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QLegend>
#include <QtCharts/QBarCategoryAxis>
QT_CHARTS_USE_NAMESPACE

#include <QDebug>
#include<cstring>

using namespace std;

DialogHistogramme::DialogHistogramme(QWidget *parent, vector<int> donnees) :
    QDialog(parent)
{
    setupUi(this);

    this->_donnees = donnees;

    vector<QBarSet*> sets (donnees.size());
    QBarSeries *series = new QBarSeries();

    for (unsigned i=0; i<_donnees.size(); i++)
    {
        //sets[i] = new QBarSet(to_string(donnees[i]), nullptr);
        sets[i] = new QBarSet("tmp");
        //*sets[i] << 1;
        //series->append(sets[i]);
    }

    QChart *chart = new QChart();
    chart->addSeries(series);
    chart->setTitle("my machin cool");

    /*
    QStringList categories;
    for (unsigned i=0; i<_donnees.size(); i++)
    {
        char truc[20];
        sprintf(truc, "%d ", _donnees[i]);
        categories << truc;
    }
    */

    QBarCategoryAxis *axis = new QBarCategoryAxis();
    //axis->append(categories);
    chart->createDefaultAxes();
    chart->setAxisX(axis, series);

    chart->legend()->setVisible(true);
    chart->legend()->setAlignment(Qt::AlignBottom);

    _myQChartView->setChart(chart);
    _myQChartView->setRenderHint(QPainter::Antialiasing);

}

void DialogHistogramme::display_histo()
{
    ;
}
