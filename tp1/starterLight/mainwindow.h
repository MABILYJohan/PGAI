#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>


#include "utils.h"
#include "dialoghistogramme.h"

#include <vector>
using namespace std;

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
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
    float faceArea(MyMesh* _mesh, int faceID);
    float aire_barycentrique(MyMesh* _mesh, int vertID);
    float angleFF(MyMesh *_mesh, int faceID0, int faceID1, int vertID0, int vertID1);
    float angleEE(MyMesh* _mesh, int vertexID, int faceID);
    void H_Curv(MyMesh* _mesh);
    void K_Curv(MyMesh* _mesh);
    void Bounding_box(MyMesh* _mesh);
    void delete_bound(MyMesh* _mesh);
    // fonctions perso
    float calculateCurveOnVertex(MyMesh* _mesh, int vertexID);
    // New
    MyMesh::Point normale_sommet(MyMesh *_mesh, int vertexID);
    void frequence_aire_triangles(MyMesh *_mesh);
    float aire_maillage(MyMesh *_mesh);
    void deviation_normales(MyMesh *_mesh);
    void display_my_histogramme(MyMesh *_mesh, vector<int> v, char *title, char *labelAxe, char *valType);
    void angles_diedres(MyMesh *_mesh);
    std::vector<int> liste_valence_mesh(MyMesh* _mesh);
    int valence_circulator(MyMesh* _mesh,VertexHandle vh, int valence, int n);
    int valence_n_voisins(MyMesh* _mesh, int n);

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();
    void on_pushButton_angleArea_clicked();
    void on_pushButton_H_clicked();
    void on_pushButton_K_clicked();

    void on_pushButton_clicked();
private:

    bool modevoisinage;
    MyMesh::VertexHandle sommets[8];
    MyMesh::VertexHandle barycentre;
    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
