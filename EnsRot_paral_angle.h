#ifndef DEF_ENSROT
#define DEF_ENSROT

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <numeric>
#include <exception>
#include "omp.h"
#include <stdlib.h>
#define PI 3.14159265358979224

#include <ctime>

class Node;
class Edge1;
class Edge2;
class Graph;

int modulo(int a, int b);
double slope(const std::vector<double>& x, const std::vector<double>& y);

double* map(double* coord, double param);
double* map_sin(double* coord, double param);
double* map_cos(double* coord, double param);
double* map_aff(double* coord, double param);
double Lipschitz_sin(double param);
double Lipschitz_aff(double param);

int* get_temp(int d, double param, int x, int y, int i, int j);
int* get_range(int x, int y, int i, int j, int N, int d, double param, double Lip);

class Node{
public:
  int x,y;
  std::vector<Edge1*> edges1;
  std::vector<Edge2*> edges2;

  Node(){};
  Node(int _x, int _y);
  void add_edge1(Edge1* _edge);
  void add_edge2(Edge2* _edge);
  std::string get_name();
protected:
private:
};

class Edge1{
public:
  double et0, et1;
  Node *start, *end;
  Edge2 *edge2;
  //Edge1(){};
  Edge1(Node& _start, Node& _end, Edge2& _edge2, double _et0, double _et1);
};

class Edge2{
public:
  std::vector<double> et0, et1, et2;
  Node *start, *end;
  //Edge2(){};
  Edge2(Node& _start, Node& _end, double _et0, double _et1, double _et2);
  //Edge2( const Edge2 &otherEdge);
};

class Graph{
public:
  std::vector< std::vector<Node*> > nodes;
  std::vector<Edge1*> edges1{};
  std::vector<Edge2*> edges2{};
  double param;
  Graph(){};
  Graph(double param);
  Node*  get_node(int x, int y);
  Edge1* add_edge1(int startx, int starty, int endx, int endy, Edge2& edgelink, double _et0, double _et1);
  Edge2* add_edge2(int startx, int starty, int endx, int endy, double _et0, double _et1, double _et2);

  //Loop functions
  void  do_the_harlem_shake(int x, int y, int N, int* temp);
  void  trans(int x, int y, int N, int d, double param);
  double path_bary(int n, double t, double u, int currentThread);

  //New functions
  void wrap_the_harlem_shake(int x, int y, int N, int d, int i, int j, double param, double Lip);
  void new_trans(int x, int y, int N, int d, double param);

  //Wrappers
  std::vector<double> compute_the_rotation_ensemble();

};

#endif
