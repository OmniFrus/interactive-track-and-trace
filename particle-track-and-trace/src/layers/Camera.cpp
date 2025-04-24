#include "Camera.h"

#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>

using namespace std;

double Camera::getScreenShakeOffset() {
  double maxOffset = screenShakeProgress * maxScreenShake / maxScreenShakeDuration;
  return maxOffset * dis(gen);
}

double Camera::getZoomOffset() {
  if (zoomIn)
    return 1 + zoomProgress * zoomProgress * (maxZoom - 1) / (maxZoomDuration * maxZoomDuration);
  else
    return 1 + zoomProgress * (maxZoom - 1) / maxZoomDuration;
}

void Camera::clampCamera(double pos[3]) {
  pos[0] += getScreenShakeOffset();
  pos[1] += getScreenShakeOffset();

  auto cam = renderer->GetActiveCamera();
  double ogpos[3];
  cam->GetPosition(ogpos);

  cam->SetParallelScale(zoomLevel * getZoomOffset());
  double scale = cam->GetParallelScale();

  // Only check the x, y coords of the camera; we don't care about z
  for (int i = 0; i < 2; i++) {
    if (abs(pos[i]) + scale >= 1.0) {
      pos[i] = pos[i] >= 0 ? 1 - scale : scale - 1;
    }
  }

  cam->SetPosition(pos[0], pos[1], ogpos[2]);
  cam->SetFocalPoint(pos[0], pos[1], 0);

  // Handle background and rendering issues
  if (renderer->GetRenderWindow()) {
    int* size = renderer->GetRenderWindow()->GetSize();
    double aspect = static_cast<double>(size[0]) / size[1];
    cam->SetParallelScale(scale * std::max(1.0, 1.0 / aspect));

    // Clear background properly
    renderer->SetBackground(0.0, 0.0, 0.0);  // Ensure black background
    renderer->EraseOn();                     // Ensure background is erased
    renderer->GetRenderWindow()->SetErase(1); // Enable erase
    renderer->GetRenderWindow()->SetMultiSamples(0); // Disable multisampling

    // Force redraw if zooming
    if (zoomProgress > 0 || zoomIn) {
      renderer->GetRenderWindow()->Render();
    }
  }
}

Camera::Camera() : gen(rd()), dis(0.0, 1.0) {}

void Camera::updateData(int t) {
  if (screenShakeProgress > 0) screenShakeProgress--;
  if (zoomIn) zoomProgress += 3;

  if (zoomProgress > maxZoomDuration) {
    zoomProgress = maxZoomDuration;
    zoomIn = false;
  }

  if (zoomProgress > 0 && !zoomIn) {
    zoomProgress -= 1;
  }

  if (zoomProgress < 0) zoomProgress = 0;

  // Always ensure a clean frame on update
  if (renderer && renderer->GetRenderWindow()) {
    renderer->SetBackground(0.0, 0.0, 0.0);
    renderer->EraseOn();
    renderer->GetRenderWindow()->SetErase(1);
    renderer->GetRenderWindow()->SetMultiSamples(0);
    renderer->GetRenderWindow()->Render();
  }
}

void Camera::shakeScreen() {
  screenShakeProgress = maxScreenShakeDuration;
}

void Camera::zoomScreen() {
  zoomIn = true;
}

void Camera::handleGameOver() {
  screenShakeProgress = 0;
  zoomProgress = 0;
  zoomIn = false;

  if (renderer && renderer->GetRenderWindow()) {
    renderer->SetBackground(0.0, 0.0, 0.0);
    renderer->EraseOn();
    renderer->GetRenderWindow()->SetErase(1);
    renderer->GetRenderWindow()->Render();
  }
}