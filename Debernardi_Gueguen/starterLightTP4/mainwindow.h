#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QProgressDialog>
namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // les fonctions à compléter
    void showEdgeSelection(MyMesh* _mesh);
    void collapseEdge(MyMesh* _mesh, int edgeID);
    void decimation(MyMesh* _mesh, int percent, QString method);
    void updateEdgeSelectionIHM();
    int smallestEdge(MyMesh* _mesh);
    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);

    float faceArea(MyMesh* _mesh, int faceID);
    int getSmallestAngle(MyMesh *_mesh);
    int getSmallestPlan(MyMesh *_mesh);
    int getSmallestRatio(MyMesh *_mesh);
    float getAngleFF(MyMesh *_mesh, int faceID0, int faceID1);
    int getSmallestEdge(MyMesh *_mesh);
    int getSmallestEdgeFace(MyMesh *_mesh);
private slots:

    void on_pushButton_chargement_clicked();
    void on_pushButton_edgeMoins_clicked();
    void on_pushButton_edgePlus_clicked();
    void on_pushButton_delSelEdge_clicked();

    void on_pushButton_decimate_clicked();

private:

    bool modevoisinage;

    MyMesh mesh;
    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
