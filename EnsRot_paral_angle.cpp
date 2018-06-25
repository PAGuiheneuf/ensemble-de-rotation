/*
Se compile avec:
  g++ EnsRot.cpp -std=c++11 -fopenmp -Ofast -o rotationEnsemble
Avant de lancer le code, spécifier le nombre de threads OpenMP à utiliser:
  export OMP_NUM_THREADS=4
Lancer le code avec:
  ./rotationEnsemble
*/

#include "EnsRot_paral_angle.h"

//Paramètres globaux
int  N, d, Tps, M; // En ligne de commande
bool v      = true;   // Verbose
int nParams = 100;
int shape   = 1; //Choix de la forme: 0=simple, 1=sinus, 2=affine

//Méthodes de Node
Node::Node(int _x, int _y){
  this->x = _x;
  this->y = _y;
  //std::cout << "initial size" << this->edges2.size() << std::endl;
}
void Node::add_edge1(Edge1* _edge){
  this->edges1.push_back(_edge);
}
void Node::add_edge2(Edge2* _edge){
  this->edges2.push_back(_edge);
}
std::string Node::get_name(){
  std::ostringstream ss;
  ss << "n" << this->x << "_" << this->y;
  return ss.str();
}

//Constructeurs de Edge1 et Edge2
Edge1::Edge1(Node& _start, Node& _end, Edge2& _edge2, double _et0, double _et1){
  this->et0   = _et0;
  this->et1   = _et1;
  this->start = &_start;
  this->end   = &_end;

  this->start->add_edge1(this);
  this->edge2 = &_edge2;
}
Edge2::Edge2(Node& _start, Node& _end, double _et0, double _et1, double _et2){
  this->et0   = std::vector<double>(M+1, _et0);
  this->et1   = std::vector<double>(M+1, _et1);
  this->et2   = std::vector<double>(M+1, _et2);
  this->start = &_start;
  this->end   = &_end;

  this->start->add_edge2(this);
}

//Création du graphe
Graph::Graph(double param){
  this->param = param;
  for (int x = 0 ; x < N ; x++){
    for (int y = 0 ; y < N ; y++){
      if (shape == 0)
        this->trans(x,y,N,d, param);
      else
        this->new_trans(x,y,N,d, param);
    }
  }
}

//Ancienne méthode
void Graph::trans(int x, int y, int N, int d, double param){
  for(int i = -d-1 ; i < d+1 ; i++){
    int j = -d-1;
    int* temp = get_temp(d, param, x, y, i, j);
    this->do_the_harlem_shake(x, y, N, temp);
    delete temp;
  }
  for(int i = -d-1 ; i < d+1 ; i++){
    int j = +d+1;
    int* temp = get_temp(d, param, x, y, i, j);
    this->do_the_harlem_shake(x, y, N, temp);
    delete temp;
  }
  for(int j = -d-1 ; j < d+1 ; j++){
    int i = -d-1;
    int* temp = get_temp(d, param, x, y, i, j);
    this->do_the_harlem_shake(x, y, N, temp);
    delete temp;
  }
  for(int j = -d-1 ; j < d+1 ; j++){
    int i = +d+1;
    int* temp = get_temp(d, param, x, y, i, j);
    this->do_the_harlem_shake(x, y, N, temp);
    delete temp;
  }
}
int* get_temp(int d, double param, int x, int y, int i, int j){
  double coords[2] = {
    double(x+double(i)/(2*d))/N,
    double(y+double(j)/(2*d))/N
  };
  double* tmp  = map(coords, param);
  int*    temp = new int[2];
  temp[0] = int(round(tmp[0]*N));
  temp[1] = int(round(tmp[1]*N));
  delete tmp;
  return temp;
}
double* map(double* coord, double param){
  double* cp = new double[2];
  cp[0] = double( coord[0] + param*sin( 2*PI*coord[1] ) );
  cp[1] = double( coord[1] + param*sin( 2*PI*cp[0] ) );
  return cp;
}
void Graph::do_the_harlem_shake(int x, int y, int N, int* temp){
  int temp2[2] = {
    modulo(temp[0], N),
    modulo(temp[1], N)
  };

  temp[0] = int((temp[0]-temp2[0])/N);
  temp[1] = int((temp[1]-temp2[1])/N);

  Node* start = this->get_node(x,y);
  Node* end   = this->get_node(temp2[0],temp2[1]);

  Edge2* edg2 = NULL;
  for(int k = 0 ; k < start->edges2.size() ; k++){
    if(start->edges2[k]->end == end){
      edg2 = start->edges2[k];
      break;
    }
  }
  if(edg2==NULL){
    edg2 = this->add_edge2(x,y,temp2[0],temp2[1], 1., 0., 0.);
  }


  bool create_edge1 = true;
  for(int k = 0 ; k < start->edges1.size() ; k++){
    if( start->edges1[k]->end == end ){
      if( (int(start->edges1[k]->et0) == int(temp[0])) && (int(start->edges1[k]->et1) == int(temp[1])) ){
        create_edge1 = false;
        break;
      }
    }
  }
  if(create_edge1){
    this->add_edge1(x,y,temp2[0], temp2[1], *edg2, temp[0], temp[1]);
  }


  //std::cout << "did one full harlem shake" << std::endl;
}

//Nouvelle méthode
void Graph::new_trans(int x, int y, int N, int d, double param){
  double Lip = (shape == 1) ? Lipschitz_sin(param) : Lipschitz_aff(param);
  for(int i = -d-1 ; i < d+1 ; i++){
    int j = -d-1;
    this->wrap_the_harlem_shake(x, y, N, d, i, j, param, Lip);
  }
  for(int i = -d-1 ; i < d+1 ; i++){
    int j = -d-1;
    this->wrap_the_harlem_shake(x, y, N, d, i, j, param, Lip);
  }
  for(int j = -d-1 ; j < d+1 ; j++){
    int i = -d-1;
    this->wrap_the_harlem_shake(x, y, N, d, i, j, param, Lip);
  }
  for(int j = -d-1 ; j < d+1 ; j++){
    int i = +d+1;
    this->wrap_the_harlem_shake(x, y, N, d, i, j, param, Lip);
  }
}
double Lipschitz_sin(double param){
  double p = 2 * PI * param;
  double t = 2 + pow(p, 2);
  return (t + pow((pow(t,2) - 4), 0.5 )) / 2.;
}
double Lipschitz_aff(double param){
  double p = 4 * param;
  double t = 2 + pow(p, 2);
  return (t + pow((pow(t,2) - 4), 0.5 )) / 2.;
}
void Graph::wrap_the_harlem_shake(int x, int y, int N, int d, int i, int j, double param, double Lip){
  int* range = get_range(x, y, i, j, N, d, param, Lip);
  //std::cout << param << " " << i << " " << j << " " << range[0] << " " << range[1] << " " << range[2] << " " << range[3] << std::endl;
  int* temp = new int[2];
  for( int xt = range[0]; xt < range[1]; xt++){
    for( int yt = range[2]; yt < range[3]; yt++){
      temp[0] = xt;
      temp[1] = yt;
      this->do_the_harlem_shake(x, y, N, temp);
    }
  }
  delete temp;
  delete range;
}
int* get_range(int x, int y, int i, int j, int N, int d, double param, double Lip){
  int* range = new int[4];
  double coords[2] = {
    double(x+double(i)/(2*d))/N,
    double(y+double(j)/(2*d))/N
  };
  double* img = (shape == 1) ? map_cos(coords, param) : map_aff(coords, param);
  range[0] = int(round(( img[0]-Lip/(2*d*N))*N ));
  range[1] = int(round(( img[0]+Lip/(2*d*N))*N ))+1;
  range[2] = int(round(( img[1]-Lip/(2*d*N))*N ));
  range[3] = int(round(( img[1]+Lip/(2*d*N))*N ))+1;
  delete img;
  return range;
}
double* map_cos(double* coord, double param){
  double* cp = new double[2];
  cp[0] = double( coord[0] + param*cos( 2*PI*coord[1] ) );
  cp[1] = double( coord[1] + param*cos( 2*PI*cp[0] ) );
  return cp;
}
double* map_aff(double* coord, double param){
  double* cp = new double[2];
  double tmp = std::fmod(coord[1], 1.);
  cp[0] = ( tmp < 0.5 ) ? coord[0] + param * ( 1 - 4*tmp ) : coord[0] + param * ( -3 + 4*tmp);
  tmp = std::fmod(cp[0], 1.);
  cp[1] = ( tmp < 0.5 ) ? coord[1] + param * ( 1 - 4*tmp ) : coord[1] + param * ( -3 + 4*tmp);
  return cp;
}

//Calcul de l'ensemble de rotation
std::vector<double> Graph::compute_the_rotation_ensemble(){

  double rot[M+1];

  double begin = omp_get_wtime();
  for (int l = 0 ; l < this->edges2.size() ; l++){
    //Por chaque edge, on redimensionne les etiquettes pour la parallelisation
    this->edges2[l]->et0 = std::vector<double>(M+1, 0.);
    this->edges2[l]->et1 = std::vector<double>(M+1, 0.);
    this->edges2[l]->et2 = std::vector<double>(M+1, 1.);
    /*
    this->edges2[l]->et1.resize(M+1);
    this->edges2[l]->et2.resize(M+1);
    for(int m = 0 ; m < M+1 ; m++){
      this->edges2[l]->et0[m] = 0;
      this->edges2[l]->et1[m] = 0;
      this->edges2[l]->et2[m] = 1;
    }
    */
  }
  std::cout << omp_get_wtime() - begin << " s" << std::endl;

#pragma omp parallel for schedule(static)
  for(int k = 0 ; k < M+1 ; k++){

    //int currentThread = omp_get_thread_num();

    double theta = double(PI*k/(4*M));
    double t = cos(theta);
    double u = sin(theta);

    std::vector<double> rotTemp;
    for(int n = 1 ; n < Tps+1 ; n++){
      double maxE = this->path_bary(n, t, u, k);
      rotTemp.push_back(maxE);
    }

    //Linear regression
    std::vector<double> X;
    std::vector<double> Y;
    int size = static_cast<int16_t>(rotTemp.size());
    for(int i = size - 30 ; i < size + 1 ; i++){
      X.push_back((double)(i));
    }
    for(int i = size - 31 ; i < size ; i++){
      Y.push_back(rotTemp[i]);
    }

    double fit = slope( X, Y );

    //#pragma omp critical
    rot[k] = fit;

    if(v){
      std::cout << "k=" << k << std::endl;
    }

  }

  std::vector<double> v(rot, rot + sizeof rot / sizeof rot[0]);

  return v;

}
double Graph::path_bary(int n, double t, double u, int currentThread){
  int m = currentThread;
  double maxE = 0;

  for(int i = 0 ; i < this->edges2.size() ; i++){

    Edge2* e = this->edges2[i];
    if (e->et2[m] >= n){

      double eti  = 0;
      if (e->et2[m] == n){
        eti = e->et1[m];
      }
      else if (e->et2[m] == n+1){
        eti = e->et0[m];
      }

      for(int j = 0 ; j < e->end->edges1.size() ; j++){
        Edge1* f  = e->end->edges1[j];
        Edge2* e2 = f->edge2;
        if(e2->et2[m] == n+1){
          e2->et1[m] = std::max( e2->et1[m] , eti + t*f->et0 + u*f->et1);
          maxE = std::max(maxE, e2->et1[m]);
        }
        else if(e2->et2[m] == n){
          e2->et0[m] = e2->et1[m];
          e2->et1[m] = eti + t*f->et0 + u*f->et1;
          e2->et2[m] = n+1;
          maxE = std::max(maxE, e2->et1[m]);
        }
      }
    }


  }
  return maxE;
}

//Get/Set méthodes
Node*  Graph::get_node(int _x, int _y){
  //std::cout << this->nodes.size() << " " << _x << std::endl;
  if(this->nodes.size() <= _x){
    while (this->nodes.size() <= _x){
      this->nodes.push_back(std::vector<Node*>{});
    }
  }
  if(this->nodes[_x].size() <= _y){
    while (this->nodes[_x].size() <= _y){
      this->nodes[_x].push_back(NULL);
    }
  }

  if(this->nodes[_x][_y] == NULL){
    Node* node = new Node(_x, _y);
    this->nodes[_x][_y] = node;
    return node;
  }
  else{
    return this->nodes[_x][_y];
  }

}
Edge1* Graph::add_edge1(int _startx, int _starty, int _endx, int _endy, Edge2& _edgelink, double _et0, double _et1){
  Node* start  = this->get_node(_startx, _starty);
  Node* end    = this->get_node(_endx, _endy);
  Edge1* edge1 = new Edge1(*start, *end, _edgelink, _et0, _et1);
  this->edges1.push_back(edge1);
  return edge1;
}
Edge2* Graph::add_edge2(int _startx, int _starty, int _endx, int _endy, double _et0, double _et1, double _et2){
  Node* start  = this->get_node(_startx, _starty);
  Node* end    = this->get_node(_endx, _endy);
  Edge2* edge2 = new Edge2(*start, *end, _et0, _et1, _et2);
  this->edges2.push_back(edge2);
  return edge2;
}

//Fonctions utilitaires
int modulo(int a, int b){
  return ((a % b) + b) % b;
}
double slope(const std::vector<double>& x, const std::vector<double>& y) {
    const auto n    = x.size();
    const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    return a;
}


//main program
int main(int argc, char *argv[]){

  //Parse the command line:
  if ( (argc!=5) && (argc!=6) ){
    std::cout << "Need to have 4 or 5 arguments: N, d, Tps, M (nParams)" << std::endl;
    std::exit(1);
  }
  if(argc==5){
    N   = atoi(argv[1]);
    d   = atoi(argv[2]);
    Tps = atoi(argv[3]);
    M   = atoi(argv[4]);
  }
  else if(argc==6){
    N       = atoi(argv[1]);
    d       = atoi(argv[2]);
    Tps     = atoi(argv[3]);
    M       = atoi(argv[4]);
    nParams = atoi(argv[5]);
  }


  double veryBeginning = omp_get_wtime();

  //Tableau des ensembles de rotation
  double ROTS[nParams][M+1];

//Zone parallèle
  for(int p = 0 ; p<nParams ; p++){
      //Paramètre utilisé
      double param = double(nParams-p)/nParams;
      std::cout << "param = " << param << std::endl;

      //Création du graphe
      double begin = omp_get_wtime();
      Graph* g = new Graph(param);
      double end = omp_get_wtime();
      //Et infos
      if(v){
        std::cout << "Graphe a " << g->nodes.size() << " noeuds, " << g->edges1.size() << " edges1 et " << g->edges2.size() << " edges2." << std::endl;
        std::cout << "Construction du graphe en :           " << double(end - begin) << " s" << std::endl;
      }

      //Calcul de l'ensemble de rotation
      begin = omp_get_wtime();
      std::vector<double> rot = g->compute_the_rotation_ensemble();
      end = omp_get_wtime();
      //Et infos
      if(v){
        std::cout << "Calcul de l'ensemble de rotation en : " << double(end - begin) << " s" << std::endl;
        std::cout << "Ensemble de rotation = [ ";
        for(int i = 0 ; i < rot.size() ; i++){
          std::cout << rot[i] << " ";
        }
        std::cout << " ]" << std::endl;

        //On passe une ligne pour séparer les résultats
        std::cout << std::endl;
      }

      //On écrit dans ROTS
      for(int i = 0 ; i < rot.size() ; i++){
        ROTS[p][i] = rot[i];
      }
  }

  //On printe les ensembles de rotation
  for(int i = 0 ; i < nParams ; i++){
    for(int j = 0 ; j < M+1 ; j++){
      std::cout << ROTS[i][j] << " ";
    }
    std::cout << std::endl;
  }

  //Et on les écrit dans un fichier
  //La première ligne contient les valeurs de N, d, Tps et M
  //La seconde ligne contient les valeurs des paramètres param (de 0 à 1 donc)
  //Chaque ligne suivante contient l'ensemble de rotation pour le paramètre correspondant
  //Le fichier est à parser avec le code python d'affichage et de plot
  std::ofstream myfile("rotationEnsemble.txt");
  if (myfile.is_open()){

    //Ecriture des variables globales
    myfile << N << " " << d << " " << Tps << " " << M << std::endl;

    //Ecriture des valeurs de paramètres
    for(int p = 0 ; p<nParams ; p++){
      double param = double(nParams-p)/nParams;
      myfile << param << " " ;
    }
    myfile << std::endl;

    //Ecriture des ensembles de rotation
    for(int i = 0 ; i < nParams ; i++){
      for(int j = 0 ; j < M+1 ; j++){
        myfile << ROTS[i][j] << " ";
      }
      myfile << std::endl;
    }

    myfile.close();
  }
  else
    std::cout << "Unable to open the result file!";

  double veryEnd = omp_get_wtime();
  std::cout << "Temps total de calcul = " << double(veryEnd - veryBeginning) << " s" << std::endl;

  return 0;
}
