//
// Created by Samuel Jackson on 08/03/2016.
//

#ifndef MOLLY_SIMULATION_H
#define MOLLY_SIMULATION_H

#include <vector>
#include <Eigen/Dense>

#include "SystemParameters.h"
#include "Molecule.h"

namespace Molly {

    class Simulation {

    public:
        Simulation(SystemParameters params);

        void set_up();
        void run();
        void single_step(int step);
        std::vector<Molecule>& get_mols() { return mols; };

    private:
        const SystemParameters params;

        unsigned int num_mols;
        std::vector<Molecule> mols;
        double initial_magnitude;
        double time_now;
        const int step_limit;
        int step_count;

        Eigen::Vector3d unit_cell;
        Eigen::Vector3d region;

        void init_props();
        void init_coordinates();

        void integrate(int part);

        void apply_boundary_check(Eigen::Vector3d* vector);
        void update_bound(double& value, const double& bound);
        void apply_boundary_checks();
        void compute_forces();

        double potential(double rr);
    };

}


#endif //MOLLY_SIMULATION_H
