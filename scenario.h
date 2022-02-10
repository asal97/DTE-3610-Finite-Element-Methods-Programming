#ifndef SCENARIO_H
#define SCENARIO_H

#include "application/gmlibwrapper.h"

// qt
#include <QObject>
class FEMObject;

class Scenario : public GMlibWrapper {
  Q_OBJECT
public:
  using GMlibWrapper::GMlibWrapper;

  void initializeScenario() override;
  void cleanupScenario() override;

public slots:
  void callDefferedGL();

private:
  FEMObject *fem{nullptr};
};

#endif // SCENARIO_H
