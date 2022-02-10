#include "femobject.h"
#include <cmath>
#include <cstdlib>

using namespace GMlib;

FEMObject::FEMObject() {}

void FEMObject::regularTriangulation(int n, int m, float r) {
  this->insertAlways(TSVertex<float>(Point<float, 2>(0.0, 0.0)));

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n * (j + 1); i++) {
      auto const alpha = (i * M_2PI) / (n * (j + 1));

      auto const u = Vector<float, 2>(std::cos(alpha), -std::sin(alpha));
      auto const v = Vector<float, 2>(std::sin(alpha), std::cos(alpha));

      Matrix<float, 2, 2> R = Matrix<float, 2, 2>(u, v);

      Point<float, 2> p = R * Vector<float, 2>((r * (j + 1)) / m, 0);
      this->insertAlways(TSVertex<float>(p));
    }
  }
  this->triangulateDelaunay();
}

Vector<Vector<float, 2>, 3> FEMObject::findVectors(TSEdge<float> *e) {
  Array<TSTriangle<float> *> tr = e->getTriangle();

  Array<TSVertex<float> *> v1 = tr[0]->getVertices();
  Array<TSVertex<float> *> v2 = tr[1]->getVertices();

  Point<float, 2> p0, p1, p2, p3;

  p0 = e->getFirstVertex()->getParameter();
  p1 = e->getLastVertex()->getParameter();

  for (int i = 0; i < 3; i++) {
    if (v1[i] != e->getFirstVertex() && v1[i] != e->getLastVertex())
      p2 = v1[i]->getParameter();
    if (v2[i] != e->getFirstVertex() && v2[i] != e->getLastVertex())
      p3 = v2[i]->getParameter();
  }

  Vector<Vector<float, 2>, 3> d{};
  d[0] = p1 - p0;
  d[1] = p2 - p0;
  d[2] = p3 - p0;

  return d;
}

Vector<Vector<float, 2>, 3> FEMObject::findVectors(TSTriangle<float> *tr,
                                                   Node *n) {
  Point<float, 2> p0, p1, p2;
  Array<TSVertex<float> *> v = tr->getVertices();

  if (n->isThis(v[1])) {
    std::swap(v[0], v[1]);
    std::swap(v[1], v[2]);
  }
  if (n->isThis(v[2])) {
    std::swap(v[0], v[2]);
    std::swap(v[1], v[2]);
  }

  p0 = v[0]->getParameter();
  p1 = v[1]->getParameter();
  p2 = v[2]->getParameter();

  Vector<Vector<float, 2>, 3> d{};
  d[0] = p1 - p0;
  d[1] = p2 - p0;
  d[2] = p2 - p1;

  return d;
}

void FEMObject::computation() {

  for (int i = 0; i < this->size(); i++) {
    if (!getVertex(i)->boundary())
      _nodes += Node(getVertex(i));
  }

  // initiliazing the stiffnes matrix
  _A = DMatrix<float>(_nodes.size(), _nodes.size());

  for (int i = 0; i < _A.getDim1(); i++) {
    for (int j = 0; j < _A.getDim2(); j++) {
      _A[i][j] = 0;
    }
  }

  // initilizing b vector
  _b = DVector<float>(_nodes.size());

  //  Assembly of the stiffness matrix
  for (int i = 0; i < _nodes.size(); i++) {
    for (int j = 0; j < _nodes.size(); j++) {

      //      std::cout << "A" << i << "," << j << " " << _A[i][j] <<
      //      std::endl;

      TSEdge<float> *edge = _nodes[i].getNeighbor(_nodes[j]);

      if (edge != NULL) {

        // compute non-diagonal element of the stiffness matrix

        auto const vec = findVectors(edge);

        auto const dd = 1 / (vec[0] * vec[0]);
        auto const dh1 = dd * vec[1] * vec[0];
        auto const dh2 = dd * vec[2] * vec[0];

        auto const area1 = std::abs(vec[0] ^ vec[1]);
        auto const area2 = std::abs(vec[0] ^ vec[2]);

        auto const h1 = dd * area1 * area1;
        auto const h2 = dd * area2 * area2;

        _A[i][j] = _A[j][i] = ((dh1 * (1 - dh1) / h1 - dd) * area1 / 2) +
                              ((dh2 * (1 - dh2) / h2 - dd) * area2 / 2);
      }
    }

    Array<TSTriangle<float> *> triangles = _nodes[i].getTriangles();

    float diagonal_element = 0;

    for (int k = 0; k < triangles.size(); k++) {

      auto const vectors = findVectors(triangles[k], &_nodes[i]);
      diagonal_element +=
          (vectors[2] * vectors[2]) / (2 * std::abs(vectors[0] ^ vectors[1]));
    }
    _A[i][i] = diagonal_element;

    float stk = 0;
    for (int k = 0; k < triangles.size(); k++) {
      stk += triangles[k]->getArea2D();
    }
    _b[i] = _f * 1. / 3 * stk;
  }
}
void FEMObject::UpdateHeight(float f) {
  for (int i = 0; i < _nodes.size(); i++) {
    _nodes[i].setZ(_h[i] * f);
  }
}

void FEMObject::setForce(float f) { _f = f; }
