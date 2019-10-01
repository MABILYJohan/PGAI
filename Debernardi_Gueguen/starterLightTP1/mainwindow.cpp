#include "mainwindow.h"
#include "ui_mainwindow.h"

/* **** début de la partie à compléter **** */

void MainWindow::colorVertex(MyMesh* _mesh, int vertex, MyMesh::Color color, int thickness)
{
    if (vertex >= 0) {
        _mesh->set_color(_mesh->vertex_handle(vertex), color);
        _mesh->data(_mesh->vertex_handle(vertex)).thickness = thickness;
    }
}

void MainWindow::showSelections(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    /* **** à compléter ! ****
     * cette fonction utilise les variables de sélection vertexSelection, edgeSelection et faceSelection
     * qui sont les ID des élements sélectionnés et qui sont égales à -1 si la sélection est vide
     */

    // VERTEX
    if (vertexSelection >= 0) {
        colorVertex(_mesh, vertexSelection, MyMesh::Color(255, 0, 0));
    }

    // EDGE
    if (edgeSelection >= 0) {
        _mesh->set_color(_mesh->edge_handle(edgeSelection), MyMesh::Color(0, 255, 0));
        _mesh->data(_mesh->edge_handle(edgeSelection)).thickness = 3;

        // on récupère une des demies arêtes pour colorier les sommets bornant l'arête
        EdgeHandle edge_h = EdgeHandle(edgeSelection);
        HalfedgeHandle half_edge = _mesh->halfedge_handle(edge_h, 0);

        // colorier les deux sommets
        colorVertex(_mesh, _mesh->to_vertex_handle(half_edge).idx(), MyMesh::Color(0, 255, 0));
        colorVertex(_mesh, _mesh->from_vertex_handle(half_edge).idx(), MyMesh::Color(0, 255, 0));
    }

    // FACE
    if (faceSelection >= 0) {
        MyMesh::FaceHandle face = _mesh->face_handle(faceSelection);
        _mesh->set_color(face, MyMesh::Color(0, 0, 255));

        // colorier les trois sommets de la face triangulaire
        for (MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(face); curVertex.is_valid(); curVertex++) {
            colorVertex(_mesh, (*curVertex).idx(), MyMesh::Color(0, 0, 255));
        }
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}


void MainWindow::showSelectionsNeighborhood(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    // VERTEX
    if (vertexSelection >= 0) {
        colorVertex(_mesh, vertexSelection, MyMesh::Color(255, 0, 0));

        // colorier les arêtes adjacentes
        for (MyMesh::VertexEdgeIter ve_it = _mesh->ve_iter(_mesh->vertex_handle(vertexSelection)); ve_it.is_valid(); ++ve_it) {
            EdgeHandle eh = *ve_it;
            _mesh->set_color(eh, MyMesh::Color(255, 0, 0));
        }
    }

    // EDGE
    if (edgeSelection >= 0) {
        // colorier l'arête
        EdgeHandle edge = _mesh->edge_handle(edgeSelection);
        _mesh->set_color(edge, MyMesh::Color(0, 255, 0));
        _mesh->data(edge).thickness = 3;

        // colorier la première face
        HalfedgeHandle half0 = _mesh->halfedge_handle(edge, 0);
        _mesh->set_color(_mesh->face_handle(half0), MyMesh::Color(153, 255, 153));

        // si on n'est pas sur une bordure, alors on colorie la seconde face
        if ( ! _mesh->is_boundary(edge)) {
            HalfedgeHandle half1 = _mesh->halfedge_handle(edge, 1);
            _mesh->set_color(_mesh->face_handle(half1), MyMesh::Color(153, 255, 153));
        }
    }

    // FACE
    if (faceSelection >= 0) {
        // colorier la face
        FaceHandle face = _mesh->face_handle(faceSelection);
        _mesh->set_color(face, MyMesh::Color(0, 0, 255));

        // colorier les faces adjacentes
        for (MyMesh::FaceFaceIter curFace = _mesh->ff_iter(face); curFace.is_valid(); ++curFace) {
            FaceHandle f = *curFace;
            _mesh->set_color(f, MyMesh::Color(153, 204, 255));
        }
    }

    //affiche de nouveau le maillage
    displayMesh(_mesh);
}



void MainWindow::showBorder(MyMesh* _mesh)
{
    // on réinitialise l'affichage
    resetAllColorsAndThickness(_mesh);

    // pour toutes les arêtes on va vérifier si elles sont en bordure
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++) {
        EdgeHandle edge = *curEdge;

        // pour afficher la bordure on utilise la fonction fournie par openMesh
        if (_mesh->is_boundary(edge)) {
            _mesh->set_color(edge, MyMesh::Color(255, 0, 255));
            _mesh->data(edge).thickness = 3;
        }
    }
    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

int MainWindow::findMinDistance(MyMesh* _mesh, float* distances, bool* marked)
{
    int min_dist = INFINITY;
    int closest_vertex;

    // pour tous les sommets
    // on vérifie que le sommet n'a pas encore été visité
    // et qui a une distance de la source minimale
    for(int i = 0; i < _mesh->n_vertices(); i++) {
        if( ! marked[i] && (distances[i] < min_dist)) {
            min_dist = distances[i];
            closest_vertex = i;
        }
    }

    return closest_vertex;
}

int MainWindow::findEdge(MyMesh* _mesh, int source, int destination)
{
    // trouver une arête à partir de ses deux sommets
    for(MyMesh::VertexEdgeIter edge = _mesh->ve_iter(VertexHandle(source)); edge.is_valid(); edge++) {
        EdgeHandle eh = *edge;
        HalfedgeHandle heh = _mesh->halfedge_handle(eh, 0);
        VertexHandle v1 = _mesh->to_vertex_handle(heh);
        VertexHandle v2 = _mesh->from_vertex_handle(heh);

        if(v1.idx() == destination || v2.idx() == destination) {
            return eh.idx();
        }
    }

    //qDebug() << "error " << source << " : " << destination << endl;
    return -1;
}

void MainWindow::showPath(MyMesh* _mesh, int v1, int v2)
{
    // on réinitialise l'affichage
    resetAllColorsAndThickness(_mesh);

    // point de départ et point d'arrivée en vert et en gros
    _mesh->set_color(_mesh->vertex_handle(v1), MyMesh::Color(255, 0, 0));
    _mesh->set_color(_mesh->vertex_handle(v2), MyMesh::Color(255, 255, 255));
    _mesh->data(_mesh->vertex_handle(v1)).thickness = 12;
    _mesh->data(_mesh->vertex_handle(v2)).thickness = 12;

    // on utilise l'algorithme de Dijkstra
    // on conserve donc les distances des sommets visités à la source (v1)
    // les prédecesseurs ainsi qu'un tableau pour marquer les sommets déjà visités
    int vertexCount = _mesh->n_vertices();
    float distances[vertexCount];
    int predecessors[vertexCount];
    bool marked[vertexCount];

    // initialisation des tableaux
    for (MyMesh::VertexIter curVer = _mesh->vertices_begin(); curVer != _mesh->vertices_end(); curVer++) {
        VertexHandle vertex = *curVer;
        distances[vertex.idx()] = INFINITY;
        predecessors[vertex.idx()] = -1;
    }

    // la source est à distance 0
    distances[v1] = 0.0;
    int current_vertex = v1;

    // tant qu'il reste des sommets à visiter
    int count = 1;
    while(count != _mesh->n_vertices()) {
        // on repart du sommet 'le plus prometteur' c'est à dire le moins éloigné de la source
        current_vertex = findMinDistance(_mesh, distances, marked);

        // si on est arrivé à destination, il est inutile de continuer
        // l'algorithme de Dijkstra assure que l'on possède déjà le chemin le plus court
        if(current_vertex == v2) {
            break;
        }

        // on itère sur toutes les arêtes adjacentes au sommet courant
        // ceci est utile pour le critère de la taille de l'arête
        // et pas seulement le nombre d'arêtes
        //
        // il est peut être plus rapide d'itérer uniquement sur les vertex...
        for(MyMesh::VertexEdgeIter edge = _mesh->ve_iter(VertexHandle(current_vertex)); edge.is_valid(); edge++) {
            _mesh->set_color(edge, MyMesh::Color(0, 255, 0));
            _mesh->data(edge).thickness = 3;

            // pour récupérer les deux vertex de l'arête
            EdgeHandle eh = *edge;
            HalfedgeHandle heh = _mesh->halfedge_handle(eh, 0);
            VertexHandle vh1 = _mesh->to_vertex_handle(heh);
            VertexHandle vh2 = _mesh->from_vertex_handle(heh);

            // on calcule la taille de l'arête
            float edgeSize = _mesh->calc_edge_length(eh);
            //qDebug() << "edgeID : " << eh.idx() << " - " << edgeSize;

            // le critère choisi ici est le nombre d'arête
            float temp_dist = distances[current_vertex] + 1; //edgeSize;

            // quel sommet est la destination ?
            int v_neighbor = vh1.idx() == current_vertex ? vh2.idx() : vh1.idx();

            // si on a trouvé un chemin plus, court, mettons le à jour
            if(temp_dist < distances[v_neighbor]) {
                distances[v_neighbor] = temp_dist;
                predecessors[v_neighbor] = current_vertex;
            }
        }

        // le sommet est visité
        marked[current_vertex] = true;
        count++;
    }

    // retrouver le chemin entre la destination et la source
    int cible = v2;
    while(predecessors[cible] != -1) {
        // colorier l'arête entre le sommet courant et son prédécesseur
        int edge = findEdge(_mesh, cible, predecessors[cible]);
        if(edge >= 0) {
            EdgeHandle eh = _mesh->edge_handle(edge);
            _mesh->set_color(eh, MyMesh::Color(255, 255, 255));
            _mesh->data(eh).thickness = 3;
        }

        cible = predecessors[cible];
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

/* **** fin de la partie à compléter **** */


/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_bordure_clicked()
{
    showBorder(&mesh);
}

void MainWindow::on_pushButton_voisinage_clicked()
{
    // changement de mode entre avec et sans voisinage
    if(modevoisinage)
    {
        ui->pushButton_voisinage->setText("Repasser en mode normal");
        modevoisinage = false;
    }
    else
    {
        ui->pushButton_voisinage->setText("Passer en mode voisinage");
        modevoisinage = true;
    }

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}


void MainWindow::on_pushButton_vertexMoins_clicked()
{
    // mise à jour de l'interface
    vertexSelection = vertexSelection - 1;
    ui->labelVertex->setText(QString::number(vertexSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_vertexPlus_clicked()
{
    // mise à jour de l'interface
    vertexSelection = vertexSelection + 1;
    ui->labelVertex->setText(QString::number(vertexSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_edgeMoins_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection - 1;
    ui->labelEdge->setText(QString::number(edgeSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_edgePlus_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection + 1;
    ui->labelEdge->setText(QString::number(edgeSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_faceMoins_clicked()
{
    // mise à jour de l'interface
    faceSelection = faceSelection - 1;
    ui->labelFace->setText(QString::number(faceSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_facePlus_clicked()
{
    // mise à jour de l'interface
    faceSelection = faceSelection + 1;
    ui->labelFace->setText(QString::number(faceSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_afficherChemin_clicked()
{
    // on récupère les sommets de départ et d'arrivée
    int indexV1 = ui->spinBox_v1_chemin->value();
    int indexV2 = ui->spinBox_v2_chemin->value();

    showPath(&mesh, indexV1, indexV2);
}


void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

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
void MainWindow::displayMesh(MyMesh* _mesh)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
    MyMesh::ConstFaceVertexIter fvIt;
    int i = 0;
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

