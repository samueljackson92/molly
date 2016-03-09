//
// Created by Samuel Jackson on 08/03/2016.
//

#include <iostream>
#include <Eigen/Dense>
#include "SystemParameters.h"
#include "Simulation.h"

using namespace Eigen;
using namespace Molly;

Simulation::Simulation(SystemParameters params)
        : params(params), step_limit(10), step_count(0) {

}

void Simulation::set_up() {
    init_props();
    init_coordinates();
}

void Simulation::init_props() {
    // create unit cell of volume
    unit_cell << 10, 10, 10;
    // scale the region size down by density
    region = unit_cell * (1. / sqrt(params.density));

    num_mols = 1000; // 10^3 molecules
    mols.resize(num_mols);

    initial_magnitude = sqrt(3 * (1. - 1./num_mols) * params.temperature);
}

void Simulation::init_coordinates() {
    Vector3d gap = region.cwiseQuotient(unit_cell);

    int n = 0;
    for (int nx = 0; nx < unit_cell(0); ++nx) {
        for (int ny = 0; ny < unit_cell(1); ++ny) {
            for(int nz = 0; nz < unit_cell(2); ++nz) {
                // set up position
                mols[n].r(0) = nx+0.5;
                mols[n].r(1) = ny+0.5;
                mols[n].r(2) = nz+0.5;

                mols[n].r = mols[n].r.cwiseProduct(gap);
                mols[n].r += -0.5 * region;

                // set up initial velocity
                mols[n].rv = Vector3d::Random();
                mols[n].rv *= initial_magnitude;
                ++n;
            }
        }
    }
}

void Simulation::run() {
    int n = 0;
    while(n < step_limit) {
        single_step(n);
        ++n;
    }

}

void Simulation::single_step(int step) {
    step_count = step;
    time_now = step_count * params.delta_time;

    integrate(1);
    apply_boundary_checks();
    integrate(2);
    compute_forces();
//    evaluate_measurements();
//    accumulate_measurements();

    std::cout << "Iteration " << step_count << std::endl;

    if (step_count % params.step_average == 0) {
//        average_measurements();
//        print_measurements();
//        clear_measurements();
    }
}

void Simulation::compute_forces() {
    Vector3d dr;
    double fc_value, rr;
    double uSum = 0.;
    double virSum = 0.;
    double norm;
    double cut_off = params.r_cutoff * params.r_cutoff;

    // reset acceleration vectors
    for(auto iter = mols.begin(); iter != mols.end(); ++iter) {
        iter->ra.fill(0);
    }

    // compute interactions between molecules
    for (int i=0; i < mols.size()-1; ++i) {
        for(int j=i+1; j < mols.size(); ++j) {
            // find difference between mols i and j
            dr = mols[i].r - mols[j].r;
//            // apply boundary checks to position difference
//            apply_boundary_check(&dr);
//
//            norm = dr.norm();
//            rr = norm*norm;
//
//            if(rr < cut_off) {
//                fc_value = potential(rr);
//                dr *= fc_value;
//                mols[i].ra += dr;
//                mols[j].ra -= dr;
//            }
        }
    }
}

/**
 * Lennard Jones potential function
 */
double Simulation::potential(double rr) {
    double rri = 1./rr;
    double rri3 = rri * rri * rri;
    double fc_val = 48. * rri3 * (rri3 - 0.5) * rri;
//    uSum += 4. * rri3 * (rri3 - 1.) + 1.;
//    virSum += fc_val * rr;
    return fc_val;
}

/**
 * Uses leapfrog integration method.
 *
 */
void Simulation::integrate(int part) {

    for (auto iter = mols.begin(); iter != mols.end(); ++iter) {
        Molecule mol = *iter;
        switch(part) {
            case 1:
                // update velocity to by half a time step.
                // update position one whole step.
                mol.rv += 0.5 * params.delta_time * mol.ra;
                mol.r += params.delta_time * mol.rv;
                *iter = mol;
                break;
            case 2:
                // just update velocity by half a time step.
                mol.rv += 0.5 * params.delta_time * mol.ra;
                *iter = mol;
                break;
            default:
                throw std::runtime_error("Invalid part operation for integrate");
        }
    }
}

void Simulation::apply_boundary_checks() {
    for(auto iter = mols.begin(); iter != mols.end(); ++iter) {
        apply_boundary_check(&(iter->r));
    }
}

void Simulation::apply_boundary_check(Vector3d* vector) {
    update_bound((*vector)(0), region(0));
    update_bound((*vector)(1), region(1));
    update_bound((*vector)(2), region(2));
}

void Simulation::update_bound(double& value, const double& bound) {
    value = (value >= 0.5 * bound) ? value - bound : value;
    value = (value < -0.5 * bound) ? value + bound : value;
}
