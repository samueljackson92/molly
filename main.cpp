#include <iostream>

#include "SystemParameters.h"
#include "Simulation.h"
#include "SystemVisualisation.h"

using namespace std;
using namespace Molly;

int main(int argc, char** argv) {
    SystemParameters params;
    params.density = .8;
    params.temperature = 1.0;
    params.delta_time = .005;
    params.step_average = 100;
    params.r_cutoff = pow(2., 1./6.);

    Simulation sim(params);
    sim.set_up();
//    sim.run();

    SystemVisualisation viz;
    viz.set_up(&sim);
    return 0;
}