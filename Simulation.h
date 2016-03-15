//
// Created by Samuel Jackson on 08/03/2016.
//

#ifndef MOLLY_SIMULATION_H
#define MOLLY_SIMULATION_H

#include <vector>
#include <Eigen/Dense>

#include "SystemParameters.h"
#include "Molecule.h"
#include "Property.h"
#include "Cell.h"
#include "CellMatrix.h"

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

        Property total_energy;
        Property kinetic_energy;
        Property pressure;

        unsigned int num_mols;
        std::vector<Molecule> mols;

        double initial_magnitude;
        double time_now;
        const int step_limit;
        int step_count;
        double vir_sum;
        double u_sum;

        Eigen::Vector3d v_sum;
        Eigen::Vector3d unit_cell;
        Eigen::Vector3d region;
        Eigen::Vector3d cell_size;

        CellMatrix cell_matrix;

        void init_props();
        void init_cells();
        void init_coordinates();
        void print_header();
        void integrate(int part);
        void apply_boundary_check(Eigen::Vector3d& vector);
        void update_bound(double& value, const double& bound);
        void apply_boundary_checks();
        void compute_forces();
        double potential(double rr);
        void evaluate_properties();
        void average_properties();
        void print_summary();
        void clear_properties();
    };

}


#endif //MOLLY_SIMULATION_H
