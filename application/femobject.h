#ifndef FEMOBJECT_H
#define FEMOBJECT_H

#include "../../gmlib/modules/trianglesystem/gmtrianglesystem.h"

#include "node.h"

using namespace GMlib;

class FEMObject : public TriangleFacets<float> {

private:
  ArrayLX<Node> _nodes; // array of nodes
  DMatrix<float> _A;    // stiffness matrix
  DVector<float> _b;    // load vector
  float _f;             // external force
  DVector<float> _h;    // height for all the nodes

public:
  FEMObject();

  void regularTriangulation(int n, int m, float r);

  Vector<Vector<float, 2>, 3> findVectors(TSEdge<float> *e);

  Vector<Vector<float, 2>, 3> findVectors(TSTriangle<float> *tr, Node *n);

  void fillInternalNodes();

  void computation();
  void UpdateHeight(float f);
  void setForce(float f);
};

#endif // FEMOBJECT_H
