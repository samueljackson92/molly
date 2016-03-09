//
// Created by Samuel Jackson on 09/03/2016.
//

#ifndef MOLLY_SYSTEMVISUALISATION_H
#define MOLLY_SYSTEMVISUALISATION_H

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include "Molecule.h"
#include "Simulation.h"

namespace Molly {
    class SystemVisualisation {
    public:
        SystemVisualisation();
        void set_up(Simulation* sim);
        void display(Simulation* sim);
        void draw_axes();

    private:
        vtkSmartPointer<vtkPoints> points;
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkRenderer> renderer;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
    };
}

#endif //MOLLY_SYSTEMVISUALISATION_H
