#ifndef NODE_H
#define NODE_H

#include "../../gmlib/modules/trianglesystem/gmtrianglesystem.h"
using namespace GMlib;

class Node {

public:
  TSVertex<float> *_v;

public:
  Node();

  Node(TSVertex<float> *ve);

  Array<TSTriangle<float> *> getTriangles();

  TSEdge<float> *getNeighbor(Node &n);

  bool isThis(TSVertex<float> *v);

  void setZ(float z);
};

#endif // NODE_H
