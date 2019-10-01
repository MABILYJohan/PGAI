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

DialogHistogramme::DialogHistogramme(QWidget *parent, vector<int> donnees, vector<char*> labels) :
    QDialog(parent)
{
    setupUi(this);

    this->_donnees = donnees;

    vector<QBarSet*> sets (donnees.size());
    QBarSeries *series = new QBarSeries();

    for (unsigned i=0; i<_donnees.size(); i++)
    {
        //qDebug() << "labels[" << i << "] = " << labels[i] ;
        sets[i] = new QBarSet(labels[i]);
        *sets[i] << donnees[i];
        series->append(sets[i]);
    }

    QChart *chart = new QChart();
    chart->addSeries(series);
    chart->setTitle("histogramme");


    QBarCategoryAxis *axisX = new QBarCategoryAxis();
    /*
    QBarCategoryAxis *axisY = new QBarCategoryAxis();
    axisY->setMin("0");
    axisY->setMax("100");
    QStringList categories;
    for (unsigned i=0; i<_donnees.size(); i++)
    {
        char truc[20];
        sprintf(truc, "%d ", i*10);
        categories << truc;
    }
    axisY->append(categories);
    */

    chart->createDefaultAxes();
    chart->setAxisX(axisX, series);
    //chart->setAxisY(axisY);

    chart->legend()->setVisible(true);
    chart->legend()->setAlignment(Qt::AlignBottom);

    _myQChartView->setChart(chart);
    _myQChartView->setRenderHint(QPainter::Antialiasing);

}

