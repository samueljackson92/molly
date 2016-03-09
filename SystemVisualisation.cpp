//
// Created by Samuel Jackson on 09/03/2016.
//

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkVertex.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCubeAxesActor2D.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleTrackball.h>
#include <vtkAxisActor2D.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkCommand.h>

#include <Eigen/Dense>
#include <vtkDoubleArray.h>
#include "SystemVisualisation.h"
#include "Simulation.h"

using namespace Molly;

class vtkTimerCallback : public vtkCommand
{
public:
    static vtkTimerCallback *New()
    {
        vtkTimerCallback *cb = new vtkTimerCallback;
        cb->TimerCount = 0;
        return cb;
    }

    void set_sim(Simulation* sim, SystemVisualisation* viz) {
        this->sim = sim;
        this->viz = viz;
    }

    virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId,
                         void *vtkNotUsed(callData))
    {
        if (vtkCommand::TimerEvent == eventId)
        {
            sim->single_step(this->TimerCount);
            viz->display(sim);
            ++this->TimerCount;
        }
    }

private:
    Simulation* sim;
    SystemVisualisation* viz;
    int TimerCount;

};

SystemVisualisation::SystemVisualisation()
        : renderer(vtkSmartPointer<vtkRenderer>::New()),
          renderWindowInteractor(vtkSmartPointer<vtkRenderWindowInteractor>::New()),
          points(vtkSmartPointer<vtkPoints>::New())
{
}

void SystemVisualisation::set_up(Simulation* sim) {
    // create collection of points
    for (auto iter = sim->get_mols().begin(); iter != sim->get_mols().end(); ++iter) {
        points->InsertNextPoint(iter->r.data());
    }

    // create PolyData from list of points
    polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);

    // create sphere sources for each points
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetRadius(0.3);
    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();

    glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
    glyph3D->SetInputData(polyData);
    glyph3D->Update();

    // Visualize
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyph3D->GetOutputPort());
    mapper->ImmediateModeRenderingOn();

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(600, 600);

    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderer->AddActor(actor);
    renderer->SetBackground(.3, .3, .3);
    renderWindow->Render();

    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    renderWindowInteractor->SetInteractorStyle(style);

//    draw_axes();
    renderWindowInteractor->Initialize();

    // Sign up to receive TimerEvent
    vtkSmartPointer<vtkTimerCallback> cb =
            vtkSmartPointer<vtkTimerCallback>::New();
    cb->set_sim(sim, this);

    renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, cb);
    renderWindowInteractor->CreateRepeatingTimer(1);
    renderWindowInteractor->Start();
}

void SystemVisualisation::display(Simulation* sim) {

    int i = 0;
    for (auto iter = sim->get_mols().cbegin(); iter != sim->get_mols().cend(); ++iter)
    {
        points->SetPoint(i, iter->r.data());
        ++i;
    }

    points->Modified();
    mapper->Update();
    renderWindowInteractor->GetRenderWindow()->Render();
}

void SystemVisualisation::draw_axes() {
    // add & render CubeAxes
    vtkSmartPointer<vtkCubeAxesActor2D> axes = vtkSmartPointer<vtkCubeAxesActor2D>::New();
    axes->SetInputData(polyData);
    axes->SetFontFactor(3.0);
    axes->SetFlyModeToNone();
    axes->SetCamera(renderer->GetActiveCamera());

    vtkSmartPointer<vtkAxisActor2D> xAxis = axes->GetXAxisActor2D();
    xAxis->SetAdjustLabels(1);

    renderer->AddViewProp(axes);
}