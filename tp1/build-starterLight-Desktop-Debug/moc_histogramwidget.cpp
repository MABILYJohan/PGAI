/****************************************************************************
** Meta object code from reading C++ file 'histogramwidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../starterLight/histogramwidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#include <QtCore/QVector>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'histogramwidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_FrameProcessor_t {
    QByteArrayData data[9];
    char stringdata0[94];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_FrameProcessor_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_FrameProcessor_t qt_meta_stringdata_FrameProcessor = {
    {
QT_MOC_LITERAL(0, 0, 14), // "FrameProcessor"
QT_MOC_LITERAL(1, 15, 14), // "histogramReady"
QT_MOC_LITERAL(2, 30, 0), // ""
QT_MOC_LITERAL(3, 31, 14), // "QVector<qreal>"
QT_MOC_LITERAL(4, 46, 9), // "histogram"
QT_MOC_LITERAL(5, 56, 12), // "processFrame"
QT_MOC_LITERAL(6, 69, 11), // "QVideoFrame"
QT_MOC_LITERAL(7, 81, 5), // "frame"
QT_MOC_LITERAL(8, 87, 6) // "levels"

    },
    "FrameProcessor\0histogramReady\0\0"
    "QVector<qreal>\0histogram\0processFrame\0"
    "QVideoFrame\0frame\0levels"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_FrameProcessor[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   24,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       5,    2,   27,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    4,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 6, QMetaType::Int,    7,    8,

       0        // eod
};

void FrameProcessor::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FrameProcessor *_t = static_cast<FrameProcessor *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->histogramReady((*reinterpret_cast< QVector<qreal>(*)>(_a[1]))); break;
        case 1: _t->processFrame((*reinterpret_cast< QVideoFrame(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 0:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QVector<qreal> >(); break;
            }
            break;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            typedef void (FrameProcessor::*_t)(QVector<qreal> );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&FrameProcessor::histogramReady)) {
                *result = 0;
                return;
            }
        }
    }
}

const QMetaObject FrameProcessor::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_FrameProcessor.data,
      qt_meta_data_FrameProcessor,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *FrameProcessor::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *FrameProcessor::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_FrameProcessor.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int FrameProcessor::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void FrameProcessor::histogramReady(QVector<qreal> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
struct qt_meta_stringdata_HistogramWidget_t {
    QByteArrayData data[11];
    char stringdata0[120];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_HistogramWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_HistogramWidget_t qt_meta_stringdata_HistogramWidget = {
    {
QT_MOC_LITERAL(0, 0, 15), // "HistogramWidget"
QT_MOC_LITERAL(1, 16, 12), // "processFrame"
QT_MOC_LITERAL(2, 29, 0), // ""
QT_MOC_LITERAL(3, 30, 11), // "QVideoFrame"
QT_MOC_LITERAL(4, 42, 5), // "frame"
QT_MOC_LITERAL(5, 48, 13), // "processBuffer"
QT_MOC_LITERAL(6, 62, 12), // "QAudioBuffer"
QT_MOC_LITERAL(7, 75, 6), // "buffer"
QT_MOC_LITERAL(8, 82, 12), // "setHistogram"
QT_MOC_LITERAL(9, 95, 14), // "QVector<qreal>"
QT_MOC_LITERAL(10, 110, 9) // "histogram"

    },
    "HistogramWidget\0processFrame\0\0QVideoFrame\0"
    "frame\0processBuffer\0QAudioBuffer\0"
    "buffer\0setHistogram\0QVector<qreal>\0"
    "histogram"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_HistogramWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   29,    2, 0x0a /* Public */,
       5,    1,   32,    2, 0x0a /* Public */,
       8,    1,   35,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void, 0x80000000 | 6,    7,
    QMetaType::Void, 0x80000000 | 9,   10,

       0        // eod
};

void HistogramWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        HistogramWidget *_t = static_cast<HistogramWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->processFrame((*reinterpret_cast< QVideoFrame(*)>(_a[1]))); break;
        case 1: _t->processBuffer((*reinterpret_cast< QAudioBuffer(*)>(_a[1]))); break;
        case 2: _t->setHistogram((*reinterpret_cast< QVector<qreal>(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 2:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QVector<qreal> >(); break;
            }
            break;
        }
    }
}

const QMetaObject HistogramWidget::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_HistogramWidget.data,
      qt_meta_data_HistogramWidget,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *HistogramWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *HistogramWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_HistogramWidget.stringdata0))
        return static_cast<void*>(this);
    return QWidget::qt_metacast(_clname);
}

int HistogramWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
