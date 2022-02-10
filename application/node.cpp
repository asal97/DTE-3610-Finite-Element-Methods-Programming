#include "node.h"

Node::Node() {}

Node::Node(TSVertex<float> *ve) { this->_v = ve; }

Array<TSTriangle<float> *> Node::getTriangles() { return _v->getTriangles(); }

TSEdge<float> *Node::getNeighbor(Node &n) {

  Array<TSEdge<float> *> edg = _v->getEdges();
  for (int i = 0; i < edg.size(); i++) {

    if (n.isThis(edg[i]->getOtherVertex(*_v)))
      return edg[i];
  }
  return NULL;
}

bool Node::isThis(TSVertex<float> *v) { return (v == _v); }

void Node::setZ(float z) { _v->setZ(z); }
