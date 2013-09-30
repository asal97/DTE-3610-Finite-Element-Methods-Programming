#include "glscenerenderer.h"

// local
#include "window.h"

// Qt
#include <QSGSimpleTextureNode>

GLSceneRenderer::GLSceneRenderer(QQuickItem *parent)
  : QQuickItem(parent),
    _tex_size(1,1)
{
  setFlag(ItemHasContents);
  setSmooth(false);
}

void GLSceneRenderer::doUpdate() {
  update();
}

QSGNode* GLSceneRenderer::updatePaintNode(QSGNode* old_node, QQuickItem::UpdatePaintNodeData*) {

  if (width() <= 0 || height() <= 0) {

    delete old_node;
    return 0x0;
  }

  QSGSimpleTextureNode *node = static_cast<QSGSimpleTextureNode *>(old_node);
  if( !node ) {
    node = new QSGSimpleTextureNode;
    _tex = GMlib::GL::Texture("display_render_target", GL_TEXTURE_2D );
  }

  const QRectF r = boundingRect();
  node->setRect( QRectF(r.left(), r.top()+r.height(), r.width(), -r.height() ) );
  node->setTexture( window()->createTextureFromId( _tex.getId(), _tex_size ) );

  return node;
}

void GLSceneRenderer::itemChange(QQuickItem::ItemChange change, const QQuickItem::ItemChangeData& value) {

  if( change == QQuickItem::ItemSceneChange ) {

    if( value.window ) {

      Window *w = qobject_cast<Window*>(window());
      connect( w,     &Window::signFrameReady,
               this,  &GLSceneRenderer::doUpdate );
      connect( this,  &GLSceneRenderer::signRenderGeometryChanged,
               w,     &Window::signSceneRenderGeometryChanged );
    }
  }
}

void GLSceneRenderer::geometryChanged(const QRectF& new_geometry, const QRectF& /*oldGeometry*/) {

  _tex_size = new_geometry.toRect().size();
  emit signRenderGeometryChanged(new_geometry);
}