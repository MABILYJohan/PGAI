
// TP3 MGM

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <vector>
#include <cmath>
#include <QVector3D>
#include <iostream>


using namespace std;

float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle face_h = FaceHandle(faceID);

    // on enregistre les points de la face dans un QVector
    QVector<MyMesh::Point> points;
    for(MyMesh::FaceVertexIter curVer = _mesh->fv_iter(face_h); curVer.is_valid(); curVer++) {
        VertexHandle vertex_h = *curVer;
        points.push_back(_mesh->point(vertex_h));
    }

    return norm((points[1] - points[0]) % (points[2] - points[0])) / 2;
}

float MainWindow::angleFF(MyMesh* _mesh, int faceID0, int faceID1, int vertID0, int vertID1)
{
    /* **** à compléter ! **** */

    //qDebug() << "<" << __FUNCTION__ << ">";

    FaceHandle fh0 = _mesh->face_handle(faceID0);
    FaceHandle fh1 = _mesh->face_handle(faceID1);
    MyMesh::Point vNFace0 = _mesh->calc_face_normal(fh0);
    MyMesh::Point vNFace1 = _mesh->calc_face_normal(fh1);

    vNFace0.normalize();
    vNFace1.normalize();

    VertexHandle vh1 = _mesh->vertex_handle(vertID0);
    VertexHandle vh2 = _mesh->vertex_handle(vertID1);

    float angle=0.f;
    if ( ((vNFace0 % vNFace1) | (_mesh->point(vh2) - _mesh->point(vh1))) >=0)
    {
        //qDebug() << "\t+";
        angle = acos( (vNFace0|vNFace1) );
    }
    else {
        //qDebug() << "\t-";
        angle = -acos( (vNFace0|vNFace1) );
    }

    //qDebug() << "</" << __FUNCTION__ << ">";
    return angle;
}

float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    /* **** à compléter ! **** */

    MyMesh::Point v1;
    MyMesh::Point v2;
    FaceHandle fhId = _mesh->face_handle(faceID);
    VertexHandle vhId = _mesh->vertex_handle(vertexID);
    std::vector<VertexHandle> vh;
    for (MyMesh::FaceVertexCWIter fv_it = _mesh->fv_cwiter(fhId); fv_it.is_valid(); fv_it++)
    {
        vh.push_back(*fv_it);
    }

    MyMesh::Point A;
    MyMesh::Point B;
    MyMesh::Point C;
    for (unsigned i=0; i<vh.size(); i++)
    {
        if (vh[i] == vhId) {
            A = _mesh->point (vh[i]);
            int k=i+1;
            if (k>=vh.size()) k=0;
            B = _mesh->point (vh[k]);
            k++;
            if (k>=vh.size()) k=0;
            C = _mesh->point (vh[k]);
            break;
        }
    }
    v1 = B-A;
    v2 = C-A;
    v1.normalize();
    v2.normalize();

    float angle = acos((v1 | v2));

    return angle;
}

float MainWindow::aire_barycentrique(MyMesh* _mesh, int vertID)
{
    float sommeA = 0.f;
    VertexHandle vh = _mesh->vertex_handle(vertID);
    for (MyMesh::VertexFaceCWIter vf_it = _mesh->vf_cwiter(vh); vf_it.is_valid(); vf_it++)
    {
        FaceHandle fh = *vf_it;
        sommeA += faceArea(_mesh, fh.idx());
    }

    float aireB = (1.f/3.f) * sommeA;
    return aireB;
}

float MainWindow::aire_maillage(MyMesh *_mesh)
{
    float aireTotale=0.f;
    for (MyMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it)
    {
        FaceHandle fh = *f_it;
        aireTotale += faceArea(_mesh, fh.idx());
    }
    return aireTotale;
}

void MainWindow::test_histogramme(MyMesh *_mesh, vector<int> v, vector<char*> labels)
{
    DialogHistogramme dlh(nullptr, v, labels);
    if (dlh.exec()) {
        ;
    }
    else {
        ;
    }
}

void MainWindow::frequence_aire_triangles(MyMesh *_mesh)
{
    bool flagColor=false;
    float minAire=DBL_MAX;
    float maxAire = 0.f;
    FaceHandle faceMin, faceMax;
    vector<int> nbTriangles (10, 0); // 10 cases par tranches de 10%

    // MIN ET MAX
    for (MyMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it)
    {
        FaceHandle fh = *f_it;
        float aire = faceArea(_mesh, fh.idx());
        if (aire<=minAire) {
            minAire = aire;
            faceMin = fh;
        }
        if (aire>=maxAire) {
            maxAire = aire;
            faceMax = fh;
        }
    }
    qDebug() << "max aire = " << maxAire;
    qDebug() << "min aire = " << minAire << endl;

    // FREQUENCE TRIANGLE
    for (MyMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it)
    {
        FaceHandle fh = *f_it;
        float aire = faceArea(_mesh, fh.idx());

        float intervInf, intervSup;
        for (float i=0.f; i<nbTriangles.size(); i+=1.f)
        {
            if (i==0.f)   intervInf=minAire;
            else        intervInf = minAire + (maxAire-minAire)*(i*10.f)/100.f;
            intervSup = minAire + (maxAire-minAire)*(i+1.f)*10.f/100.f;

            if (aire>intervInf && aire<=intervSup) {
                nbTriangles[i]+=1;
                if (flagColor) {
                    _mesh->set_color(fh, MyMesh::Color(0, (10-i)*(255/10), i*(255/10)));
                }
            }
        }
    }

    // AFFICHAGE
    for (int i=0; i<nbTriangles.size(); i++)
    {
        qDebug() << " triangles entre " << i*10 << " et " << (i+1)*10 << "% de l'aire max = "
                 << nbTriangles[i] ;
    }
    qDebug() << _mesh->n_faces() << " triangles au total";
    if (flagColor) {
        _mesh->set_color(faceMin, MyMesh::Color(255, 200, 255));
        _mesh->set_color(faceMax, MyMesh::Color(0, 0, 50));
        displayMesh(_mesh);
    }

    vector<char[20]> labels(10);
    vector<char*> l(labels.size());
    for (int i=0; i<(int)nbTriangles.size(); i++)
    {
        sprintf(labels[i], "%d-%d%c", i*10, (i+1)*10, '%');
        l[i] = labels[i];
    }
    test_histogramme(_mesh, nbTriangles, l);
}

/*-------------------------------------------------------------------------
 * Renvoi la normale du sommet d'indice @vertexID
 * (pas normalisée)
 * ----------------------------------------------------------------------*/
MyMesh::Point MainWindow::normale_sommet(MyMesh *_mesh, int vertexID)
{
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    _mesh->request_face_normals();
    _mesh->request_vertex_normals();
    MyMesh::Point p = _mesh->calc_vertex_normal(vh);
    mesh.release_vertex_normals();
    return p;
}

void MainWindow::deviation_normales(MyMesh *_mesh)
{
    //qDebug() << "<" << __FUNCTION__ << ">";
    float maxAngle = 0.f;
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
    {
        VertexHandle vh = *v_it;
        MyMesh::Point normSommet = normale_sommet(_mesh, vh.idx());
        normSommet.normalize();

        vector<float> angles;
        for (MyMesh::VertexFaceIter vf_it = _mesh->vf_iter(vh); vf_it.is_valid(); vf_it++)
        {
            FaceHandle fh = *vf_it;
            MyMesh::Point nf = _mesh->calc_face_normal(fh);
            nf.normalize();
            maxAngle = acos( (normSommet | nf) );
            angles.push_back(maxAngle);
        }
        maxAngle = 0.f;
        for (unsigned i=0; i<angles.size(); i++)
        {
            if (maxAngle <= angles[i]) {
                maxAngle = angles[i];
            }
        }
        maxAngle = Utils::RadToDeg(maxAngle);
        //qDebug() << "déviation max au sommet " << vh.idx() << " = " << maxAngle << "degrés";
        _mesh->data(vh).thickness = 8;
        _mesh->set_color(vh, MyMesh::Color(0, 180-maxAngle, 0));
    }
    displayMesh(_mesh);
    //qDebug() << "</" << __FUNCTION__ << ">";
}

void MainWindow::H_Curv(MyMesh* _mesh)
{
    /* **** à compléter ! **** */

    qDebug() << "<" << __FUNCTION__ << ">";
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle vh0 = *curVert;
        float a = 1 / ( 4 * aire_barycentrique(_mesh, vh0.idx()) );

        float somme=0.f;
        for (MyMesh::VertexEdgeCWIter ve_it = _mesh->ve_cwiter(vh0); ve_it.is_valid(); ve_it++)
        {
            EdgeHandle eh = *ve_it;
            float edgeLength = _mesh->calc_edge_length(eh);

            MyMesh::HalfedgeHandle hehc1 = _mesh->halfedge_handle(eh, 0);
            MyMesh::HalfedgeHandle hehc2 = _mesh->halfedge_handle(eh, 1);
            VertexHandle vhhe1 = _mesh->to_vertex_handle(hehc1);
            VertexHandle vhhe2 = _mesh->to_vertex_handle(hehc2);
            VertexHandle vh1;
            if (vhhe1 != vh0)
                vh1 = vhhe1;
            else
                vh1 = vhhe2;

            FaceHandle fh0;
            FaceHandle fh1;
            if (vh1 == vhhe1) {
                fh0 = _mesh->face_handle(hehc1);
                fh1 = _mesh->face_handle(hehc2);
            }
            else {
                fh0 = _mesh->face_handle(hehc2);
                fh1 = _mesh->face_handle(hehc1);
            }

            float gamma = angleFF(_mesh, fh0.idx(), fh1.idx(), vh0.idx(), vh1.idx());

            somme += gamma*edgeLength;
            //qDebug() << "vh0:" << vh0.idx() << "vh1:" << vh1.idx();
        }
        //qDebug() << "";
        float H = a*somme;
        _mesh->data(vh0).value = H;
    }
    qDebug() << "</" << __FUNCTION__ << ">";

}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    /* **** à compléter ! **** */

    qDebug() << "<" << __FUNCTION__ << ">";
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle vh = *curVert;
        float a = 1 / aire_barycentrique(_mesh, vh.idx());
        //qDebug() << "a = " << a;
        float theta = 0.f;
        for (MyMesh::VertexFaceCWIter vf_it = _mesh->vf_cwiter(vh); vf_it.is_valid(); vf_it++)
        {
            FaceHandle fh = *vf_it;
            theta += angleEE(_mesh, vh.idx(), fh.idx());
        }
        float b = 2*M_PI - theta;
        //qDebug() << "b = " << b;
        float K = a*b;
        //qDebug() << "K = " << K;

        _mesh->data(vh).value = K;
    }
    qDebug() << "</" << __FUNCTION__ << ">";
}
/* **** fin de la partie à compléter **** */



/* **** début de la partie boutons et IHM **** */
void MainWindow::on_pushButton_H_clicked()
{
    H_Curv(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

void MainWindow::on_pushButton_K_clicked()
{
    K_Curv(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
    //displayMesh(&mesh, false); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}



/*-------------------------------------------------------------------------------
 * Cette fonction est à utiliser UNIQUEMENT avec le fichier testAngleArea.obj
 * Elle est appelée par le bouton "Test angles/aires"
 *
 * Elle permet de vérifier les fonctions faceArea, angleFF et angleEE.
 * Elle doit afficher :
 *
 * Aire de la face 0 : 2
 * Aire de la face 1 : 2
 * Angle entre les faces 0 et 1 : 1.5708
 * Angle entre les faces 1 et 0 : -1.5708
 * Angle au sommet 1 sur la face 0 : 0.785398
 *----------------------------------------------------------------------------*/
void MainWindow::on_pushButton_angleArea_clicked()
{
    qDebug() << "Aire de la face 0 :" << faceArea(&mesh, 0);
    qDebug() << "Aire de la face 1 :" << faceArea(&mesh, 1);

    qDebug() << "Angle entre les faces 0 et 1 :" << angleFF(&mesh, 0, 1, 1, 2);
    qDebug() << "Angle entre les faces 1 et 0 :" << angleFF(&mesh, 1, 0, 1, 2);

    qDebug() << "Angle au sommet 1 sur la face 0 :" << angleEE(&mesh, 1, 0);
    qDebug() << "Angle au sommet 3 sur la face 1 :" << angleEE(&mesh, 3, 1);

    /*
    // TEST NORMALE SOMMET
    int sommet=1;
    MyMesh::Point p = normale_sommet(&mesh, sommet);
    qDebug() << "\nnormale du sommet " << sommet
            << " x" << p[0] << "  y" << p[1] << " z" << p[2] << endl;

    // TEST AIRE TOTALE
    float aireTotale = aire_maillage(&mesh);
    qDebug() << "aire totale" << aireTotale;

    // TEST DEVIATIONS NORMALES
    deviation_normales(&mesh);

    */
    // TEST FREQUENCE AIRE TRIANGLES
    frequence_aire_triangles(&mesh);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}
/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->color(*fIt)[0]==-1)
            {
                triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
                triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
                triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
                triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++;
            }
        }
    }

    if(i>0) ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

bool checked = true;
void MainWindow::on_pushButton_clicked()
{
    if (checked)
        Bounding_box(&mesh);
    else{
        delete_bound(&mesh);
    }
    checked = !checked;
}
// exemple pour construire un mesh face par face
void MainWindow::Bounding_box(MyMesh* _mesh)
{
    float Xmin = 5000;
    float Ymin = 5000;
    float Zmin = 5000;
    float Xmax = -5000;
    float Ymax = -5000;
    float Zmax = -5000;
    MyMesh::Point tmpBary;
    tmpBary[0] = 0;
    tmpBary[1] = 0;
    tmpBary[2] = 0;
    int valence = 0;
    // calcul des points aux extremes du mesh
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        MyMesh::VertexHandle vh = (*vit);
        MyMesh::Point point = _mesh->point (vh);
        tmpBary[0] += point[0];
        tmpBary[1] += point[1];
        tmpBary[2] += point[2];
        valence ++;
        if(point[0] < Xmin)
            Xmin = point[0];
        if(point[1] < Ymin)
            Ymin = point[1];
        if(point[2] < Zmin)
            Zmin = point[2];
        if(point[0] > Xmax)
            Xmax = point[0];
        if(point[1] > Ymax)
            Ymax = point[1];
        if(point[2] > Zmax)
            Zmax = point[2];
    }
    tmpBary[0] /= valence;
    tmpBary[1] /= valence;
    tmpBary[2] /= valence;

    qDebug() << "Barycentre : x :" << tmpBary[0] << " y : " << tmpBary[1] << " z : " << tmpBary[2];
    qDebug() << "Xmin : " << Xmin << " Ymin : " << Ymin << " Zmin : " << Zmin;
    qDebug() << "Xmax : " << Xmax << " Ymax : " << Ymax << " Zmax : " << Zmax;

    //on recupere le mesh global auquel nous ajouterons les sommets de la bound box
    MyMesh *mesh = _mesh;

    sommets[0] = _mesh->add_vertex(MyMesh::Point(Xmin, Ymin, Zmin));
    sommets[1] = _mesh->add_vertex(MyMesh::Point(Xmin, Ymin, Zmax));
    sommets[2] = _mesh->add_vertex(MyMesh::Point(Xmin, Ymax, Zmin));
    sommets[3] = _mesh->add_vertex(MyMesh::Point(Xmin, Ymax, Zmax));
    sommets[4] = _mesh->add_vertex(MyMesh::Point(Xmax, Ymin, Zmin));
    sommets[5] = _mesh->add_vertex(MyMesh::Point(Xmax, Ymin, Zmax));
    sommets[6] = _mesh->add_vertex(MyMesh::Point(Xmax, Ymax, Zmin));
    sommets[7] = _mesh->add_vertex(MyMesh::Point(Xmax, Ymax, Zmax));
    barycentre = _mesh->add_vertex(tmpBary);

    qDebug () << "Rayon de la sphere englobante : " << norm(_mesh->point(sommets[0]) - _mesh->point(barycentre)) << endl;

    _mesh->update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(_mesh);

    //on agrandis les points de la boundbox et on les met en rouge pour qu'ils soient plus visiblent
    for(int i = 0 ; i <  8; i++){
        _mesh->data(sommets[i]).thickness = 3;
        _mesh->set_color(sommets[i], MyMesh::Color(255, 0, 0));
    }
    _mesh->data(barycentre).thickness = 4;
    _mesh->set_color(barycentre, MyMesh::Color(0, 255, 0));

    // on affiche le maillage
    displayMesh(_mesh);

}

void MainWindow::delete_bound(MyMesh *_mesh)
{
    MyMesh *mesh = _mesh;

    _mesh->request_vertex_status();
    _mesh->delete_vertex(sommets[0], false);
    _mesh->delete_vertex(sommets[1], false);
    _mesh->delete_vertex(sommets[2], false);
    _mesh->delete_vertex(sommets[3], false);
    _mesh->delete_vertex(sommets[4], false);
    _mesh->delete_vertex(sommets[5], false);
    _mesh->delete_vertex(sommets[6], false);
    _mesh->delete_vertex(sommets[7], false);
    _mesh->delete_vertex(barycentre, false);
    _mesh->garbage_collection();
    _mesh->update_normals();

    displayMesh(_mesh);

    qDebug() << "Boundbox deleted";
}
