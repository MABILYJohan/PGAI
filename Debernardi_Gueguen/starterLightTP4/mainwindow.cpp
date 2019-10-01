#include "mainwindow.h"
#include "ui_mainwindow.h"
//TP4
/* **** début de la partie à compléter **** */
void MainWindow::showEdgeSelection(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    /* **** à compléter ! (Partie 1) ****
     * cette fonction utilise la variables de sélection edgeSelection qui est l'ID de
     * l'arête sélectionnée et qui est égale à -1 si la sélection est vide
     */

    if(edgeSelection >= 0){
        _mesh->set_color(_mesh->edge_handle(edgeSelection), MyMesh::Color(0, 255, 0));
        _mesh->data(_mesh->edge_handle(edgeSelection)).thickness = 12;
    }


    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

int MainWindow::getSmallestEdge(MyMesh* _mesh){

    float smallest = INT_MAX;
    int id_smallest;
    // parcours de toutes les arêtes
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        EdgeHandle eh = EdgeHandle((*curEdge).idx());
        HalfedgeHandle heh = _mesh->halfedge_handle(eh,0);
        //si la taille de l'edge courrante est plus petite que la taille min, on actualise la taille du minimum ainsi que son id
        if(_mesh->calc_edge_length(eh) < smallest && _mesh->is_collapse_ok(heh)){
            smallest = _mesh->calc_edge_length(eh);
            id_smallest = eh.idx();
        }
    }
    return id_smallest;
}

//cette fonction est utilisé pour le deuxieme critere personnalisé (taille edge * air des deux faces adjacentes)
float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    //calcul de l'aire de la face "faceID"
    FaceHandle face_h = FaceHandle(faceID);
    QVector<MyMesh::Point> points;

    //parcours de tous les sommets de la face courante
    for(MyMesh::FaceVertexIter curVer = _mesh->fv_iter(face_h); curVer.is_valid(); curVer++) {
        VertexHandle vertex_h = *curVer;
        points.push_back(_mesh->point(vertex_h));
    }

    return norm((points[1] - points[0]) % (points[2] - points[0])) / 2;
}

int MainWindow::getSmallestEdgeFace(MyMesh* _mesh){

    //cette fonction retour le plus petit rapport entre la taille de l'edge et l'aire des faces adjacentes
    float smallest = INT_MAX;
    int id_smallest;
    // parcours des arêtes
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        EdgeHandle eh = EdgeHandle((*curEdge).idx());
        HalfedgeHandle heh = _mesh->halfedge_handle(eh,0);
        HalfedgeHandle oppoheh = _mesh->opposite_halfedge_handle(heh);
        FaceHandle fh1 = _mesh->face_handle(heh);
        FaceHandle fh2 = _mesh->face_handle(oppoheh);

        //calcul de la moyenne de l'aire des deux faces adjacentes
        float f1area = faceArea(_mesh, fh1.idx());
        float f2area = faceArea(_mesh, fh2.idx());
        float moyArea = f1area + f2area;
        moyArea = moyArea/2;
        if(_mesh->calc_edge_length(eh)*moyArea < smallest && _mesh->is_collapse_ok(heh)){
            smallest = _mesh->calc_edge_length(eh)*moyArea;
            id_smallest = eh.idx();
        }
    }
    return id_smallest;
}

void MainWindow::collapseEdge(MyMesh* _mesh, int edgeID)
{
    EdgeHandle eh = EdgeHandle(edgeID);
    HalfedgeHandle heh = _mesh->halfedge_handle(eh, 0);
    //avant de supprimer une edge on vérifie qu'elle peut effectivement etre supprimée
    if (_mesh->is_collapse_ok(heh)){
        //calcul de la position du point au milieu de l'arrete que l'on souhaite supprimer
        MyMesh::Point midpoint =  _mesh->calc_edge_midpoint(heh);
        VertexHandle Destination = _mesh->to_vertex_handle(heh);

        _mesh->set_point(Destination, midpoint);
        _mesh->collapse(heh);
    }


    _mesh->garbage_collection();
}

// fonction pratique pour faire des tirages aléatoires
int randInt(int low, int high){return qrand() % ((high + 1) - low) + low;}

int MainWindow::getSmallestAngle(MyMesh* _mesh){
    //cette fonction permet de recupérer le plus petit angle entre les normales des deux faces adjacentes a une edge
    float smallest = INT_MAX;
    int id_smallest;
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
        if (! _mesh->is_boundary(*curEdge) ){
            HalfedgeHandle heh = _mesh->halfedge_handle(*curEdge, 0);
            HalfedgeHandle oppoheh = _mesh->opposite_halfedge_handle(heh);
            FaceHandle fh0 = _mesh->face_handle(heh);
            FaceHandle fh1 = _mesh->face_handle(oppoheh);

            float ffangle = getAngleFF(_mesh, fh0.idx(), fh1.idx());
            //calcul de l'angle entre les deux normales
            if(ffangle < smallest && _mesh->is_collapse_ok(heh)){
                smallest = ffangle;

                id_smallest = (*curEdge).idx();
            }
        }
    }
    return id_smallest;

}

int MainWindow::getSmallestPlan(MyMesh *_mesh)
{
    //cette fonction calcul les normale du 1-ring de chaque vertex de chaque edges
    //il fait ensuite une moyenne de ces normal et retourne la plus petite des moyennes optenue
    float smallest = INT_MAX;
    int id_smallest;
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
        if (! _mesh->is_boundary(*curEdge) ){
            HalfedgeHandle heh = _mesh->halfedge_handle(*curEdge, 0);

            //il n'y as pas besoin de faire les calcules si l'edge courante ne peut pas etre supprimée
            if(_mesh->is_collapse_ok(heh)){
                VertexHandle vertex0 = _mesh->to_vertex_handle(heh);
                VertexHandle vertex1 = _mesh->from_vertex_handle(heh);
                float moyF0 = 0;
                float moyF1 = 0;
                QVector<int> faceV0;
                QVector<int> faceV1;

                //parcours de toutes les faces du vertex courrant
                for(MyMesh::VertexFaceIter vfit = _mesh->vf_iter(vertex0); vfit.is_valid(); vfit ++){
                    faceV0.push_back((*vfit).idx());
                }
                //parcours de toutes les faces du vertex opposé au vertex courrant
                for(MyMesh::VertexFaceIter vfit = _mesh->vf_iter(vertex1); vfit.is_valid(); vfit ++){
                    faceV1.push_back((*vfit).idx());
                }

                for (int i=0; i<faceV0.size(); i++){
                    moyF0 += getAngleFF(_mesh, faceV0[i], faceV0[(i+1) % faceV0.size()]);
                }

                for (int i=0; i<faceV1.size(); i++){
                    moyF1 += getAngleFF(_mesh, faceV1[i], faceV1[(i+1) % faceV1.size()]);
                }

                moyF0 = moyF0 / faceV0.size();
                moyF1 = moyF1 / faceV1.size();

                float moy = (moyF0 + moyF1)/2;
                if(moy < smallest){
                    smallest = moy;
                    id_smallest = (*curEdge).idx();
                }
            }


        }

    }
    return id_smallest;
}

int MainWindow::getSmallestRatio(MyMesh *_mesh)
{
    float smallest = INT_MAX;
    int id_smallest;
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
        if (! _mesh->is_boundary(*curEdge) ){
            HalfedgeHandle heh = _mesh->halfedge_handle(*curEdge, 0);

            if(_mesh->is_collapse_ok(heh)){
                VertexHandle vertex0 = _mesh->to_vertex_handle(heh);
                VertexHandle vertex1 = _mesh->from_vertex_handle(heh);
                float moyF0 = 0;
                float moyF1 = 0;
                QVector<int> faceV0;
                QVector<int> faceV1;

                for(MyMesh::VertexFaceIter vfit = _mesh->vf_iter(vertex0); vfit.is_valid(); vfit ++){
                    faceV0.push_back((*vfit).idx());
                }
                for(MyMesh::VertexFaceIter vfit = _mesh->vf_iter(vertex1); vfit.is_valid(); vfit ++){
                    faceV1.push_back((*vfit).idx());
                }

                for (int i=0; i<faceV0.size(); i++){
                    moyF0 += getAngleFF(_mesh, faceV0[i], faceV0[(i+1) % faceV0.size()]);
                }

                for (int i=0; i<faceV1.size(); i++){
                    moyF1 += getAngleFF(_mesh, faceV1[i], faceV1[(i+1) % faceV1.size()]);
                }

                moyF0 = moyF0 / faceV0.size();
                moyF1 = moyF1 / faceV1.size();
                float moy = (moyF0 + moyF1)/2;
                float size = _mesh->calc_edge_length(*curEdge);
                if(moy * size< smallest){
                    smallest = moy * size;
                    id_smallest = (*curEdge).idx();
                }
            }


        }

    }
    return id_smallest;
}


float MainWindow::getAngleFF(MyMesh* _mesh, int faceID0, int faceID1)
{
    FaceHandle fh0 = FaceHandle(faceID0);
    FaceHandle fh1 = FaceHandle(faceID1);
    MyMesh::Normal Nfaceh0 = _mesh->calc_face_normal(fh0);
    MyMesh::Normal Nfaceh1 = _mesh->calc_face_normal(fh1);

    float angle = acos(Nfaceh0 | Nfaceh1);


    return angle;
}


void MainWindow::decimation(MyMesh* _mesh, int percent, QString method)
{
    /* **** à compléter ! (Partie 2 et 3) ****
     * Cette fonction supprime des arêtes jusqu'à atteindre un pourcentage d'arêtes restantes, selon un critère donné
     * percent : pourcentage de l'objet à garder
     * method  : la méthode à utiliser parmis : "Aléatoire", "Par taille", "Par angle", "Par planéité"
     */
    int edgeCount = _mesh->n_edges();
    int finalEgdeCount = (edgeCount * percent)/100;

    if(method == "Aléatoire"){
        int edge = randInt(0, _mesh->n_edges());
        while(_mesh->n_edges() > finalEgdeCount){
            collapseEdge(_mesh, edge);
            edge = randInt(0, _mesh->n_edges());
        }

    }
    else if(method == "Par taille"){

        while(_mesh->n_edges() > finalEgdeCount){
            collapseEdge(_mesh, getSmallestEdge(_mesh));
        }
    }
    else if(method == "Par angle"){

        while(_mesh->n_edges() > finalEgdeCount){
            collapseEdge(_mesh, getSmallestAngle(_mesh));
        }

    }
    else if(method == "Par planéité"){
        while(_mesh->n_edges() > finalEgdeCount){
            collapseEdge(_mesh, getSmallestPlan(_mesh));
        }
    } else if (method == "Personnalisé"){
        while(_mesh->n_edges() > finalEgdeCount){
            collapseEdge(_mesh, getSmallestRatio(_mesh));
        }
    }else if (method == "Personnalisé2"){
        while(_mesh->n_edges() > finalEgdeCount){
            collapseEdge(_mesh, getSmallestEdgeFace(_mesh));
        }
    }
    else {
        qDebug() << "Méthode inconnue !!!";
    }
}

// LE PLUS OPTI EST DE REGARDER PLANEITE ET LA TAILLE

/* **** début de la partie boutons et IHM **** */
void MainWindow::updateEdgeSelectionIHM()
{
    /* **** à compléter ! (Partie 3) ****
     * Cette fonction met à jour l'interface, les critères pourrons être affichés dans la zone de texte pour les vérifier
     */

    QString infos = "";
    infos = infos + "Surface : " + QString::number(0) + "\n";
    infos = infos + "C1 : " + QString::number(0) + "\n";
    infos = infos + "C2 : " + QString::number(0) + "\n";
    infos = infos + "C3 : " + QString::number(0) + "\n";
    ui->infoEdgeSelection->setPlainText(infos);

    ui->labelEdge->setText(QString::number(edgeSelection));

    // on montre la nouvelle sélection
    showEdgeSelection(&mesh);
}
/* **** fin de la partie à compléter **** */

void MainWindow::on_pushButton_edgeMoins_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection - 1;
    updateEdgeSelectionIHM();
}

void MainWindow::on_pushButton_edgePlus_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection + 1;
    updateEdgeSelectionIHM();
}

void MainWindow::on_pushButton_delSelEdge_clicked()
{
    // on supprime l'arête d'indice edgeSelection
    collapseEdge(&mesh, edgeSelection);

    // on actualise la sélection
    showEdgeSelection(&mesh);
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

void MainWindow::on_pushButton_decimate_clicked()
{
    decimation(&mesh, ui->horizontalSlider->value(), ui->comboBox->currentText());
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


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

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
