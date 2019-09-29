#ifndef DIALOGHISTOGRAMME_H
#define DIALOGHISTOGRAMME_H

#include <vector>

#include <QtWidgets>
#include <QGraphicsWidget>



#include "ui_dialoghistogramme.h"

using namespace std;

class DialogHistogramme : public QDialog, private Ui::DialogHistogramme
{
    Q_OBJECT

private:
    vector<int> _donnees;

public:
    explicit DialogHistogramme(QWidget *parent, vector<int> _donnees);

    void display_histo();
};

#endif // DIALOGHISTOGRAMME_H
