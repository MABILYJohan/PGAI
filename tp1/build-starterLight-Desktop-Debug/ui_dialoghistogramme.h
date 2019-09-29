/********************************************************************************
** Form generated from reading UI file 'dialoghistogramme.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOGHISTOGRAMME_H
#define UI_DIALOGHISTOGRAMME_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHeaderView>

QT_BEGIN_NAMESPACE

class Ui_DialogHistogramme
{
public:
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogHistogramme)
    {
        if (DialogHistogramme->objectName().isEmpty())
            DialogHistogramme->setObjectName(QStringLiteral("DialogHistogramme"));
        DialogHistogramme->resize(664, 521);
        buttonBox = new QDialogButtonBox(DialogHistogramme);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setGeometry(QRect(310, 480, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        retranslateUi(DialogHistogramme);
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogHistogramme, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogHistogramme, SLOT(reject()));

        QMetaObject::connectSlotsByName(DialogHistogramme);
    } // setupUi

    void retranslateUi(QDialog *DialogHistogramme)
    {
        DialogHistogramme->setWindowTitle(QApplication::translate("DialogHistogramme", "Dialog", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class DialogHistogramme: public Ui_DialogHistogramme {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOGHISTOGRAMME_H
