//
// Created by Samuel Jackson on 08/03/2016.
//

#include <memory>
#include <iostream>
#include <Eigen/Dense>

#include "Molecule.h"
#include "SystemParameters.h"
#include "Simulation.h"

using namespace Eigen;
using namespace Molly;

Simulation::Simulation(SystemParameters params)
        : params(params), step_limit(10), step_count(0), vir_sum(0), u_sum(0) {

}

void Simulation::set_up() {

    printf("MOLLY\n--------------------------------------------------------------------------------------------\n");
    printf("%5s %8s %7s %7s %7s %7s %7s %7s %7s\n",
           "step", "time", "sum", "tot_sum", "tot_sq", "kin_sum", "kin_sq", "press_sum", "press_sq");
    init_props();
    init_cells();
    init_coordinates();
}

void Simulation::init_props() {
    // create unit cell of volume x * y * z
    unit_cell << 10, 10, 1;
    // scale the unit cell size by density to create the bound box
    region = unit_cell * (1. / sqrt(params.density));
    // one cell for each unit of volume
    num_mols = unit_cell.prod();
    mols.resize(num_mols);

    initial_magnitude = std::sqrt(3 * (1. - 1./num_mols) * params.temperature);
}

void Simulation::init_cells() {
    // create cells proportional to cell cut off distance.
    cell_size = region / params.r_cutoff ;
    cell_size(0) = std::ceil(cell_size(0));
    cell_size(1) = std::ceil(cell_size(1));
    cell_size(2) = std::ceil(cell_size(2));
    cell_matrix.resize(cell_size(0), cell_size(1), cell_size(2));
}

void Simulation::init_coordinates() {
    Vector3d gap = region.cwiseQuotient(unit_cell);
    std::vector<Molecule_ptr> mols;
    mols.resize(num_mols);

    int n = 0;
    for (int nx = 0; nx < unit_cell(0); ++nx) {
        for (int ny = 0; ny < unit_cell(1); ++ny) {
            for(int nz = 0; nz < unit_cell(2); ++nz) {
                // set up position
                mols[n] = std::make_shared<Molecule>();
                mols[n]->r(0) = nx+0.5;
                mols[n]->r(1) = ny+0.5;
                mols[n]->r(2) = nz+0.5;

                mols[n]->r = mols[n]->r.cwiseProduct(gap);
                mols[n]->r += -0.5 * region;

                // set up initial velocity
                mols[n]->rv = Vector3d::Random();
                mols[n]->rv.normalize();
                mols[n]->rv *= initial_magnitude;
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
    compute_forces();
    integrate(2);
    evaluate_properties();

    if (step_count % params.step_average == 0) {
        average_properties();
        print_summary();
        clear_properties();
    }
}

void Simulation::compute_forces() {
    Vector3d dr = Vector3d::Zero();
    double rr = 0.;
    double cut_off = params.r_cutoff * params.r_cutoff;

    // reset sums for this iteration
    u_sum = 0.;
    vir_sum = 0.;

    // reset acceleration vectors
    for(auto iter = mols.begin(); iter != mols.end(); ++iter) {
        iter->ra.fill(0);
    }

    // compute interactions between molecules
    for (int i=0; i < mols.size()-1; ++i) {
        for(int j=i+1; j < mols.size(); ++j) {
            // find difference between mols i and j
            dr = mols[i].r - mols[j].r;
            // apply boundary checks to position difference
            apply_boundary_check(dr);

            rr = dr.squaredNorm();

            if(rr < cut_off) {
                dr*= potential(rr);
                mols[i].ra += dr;
                mols[j].ra -= dr;
            }
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
    u_sum += 4. * rri3 * (rri3 - 1.) + 1.;
    vir_sum += fc_val * rr;
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
        apply_boundary_check(iter->r);
    }
}

void Simulation::apply_boundary_check(Vector3d& vector) {
    update_bound(vector(0), region(0));
    update_bound(vector(1), region(1));
    update_bound(vector(2), region(2));
}

void Simulation::update_bound(double& value, const double& bound) {
    value = (value >= 0.5 * bound) ? value - bound : value;
    value = (value < -0.5 * bound) ? value + bound : value;
}

void Simulation::evaluate_properties() {
    v_sum = Vector3d::Zero();
    double v_sqrd_sum = 0.;

    for(auto iter = mols.cbegin(); iter != mols.cend(); ++iter) {
        v_sum += iter->rv;
        v_sqrd_sum += iter->rv.squaredNorm();
    }

    double kin_energy = 0.5 * v_sqrd_sum / num_mols;
    kinetic_energy.set_value(kin_energy);
    total_energy.set_value(kin_energy + (u_sum / num_mols));
    pressure.set_value(params.density * (v_sqrd_sum + vir_sum) / (num_mols / 3));
}

void Simulation::average_properties() {
    kinetic_energy.average(params.step_average);
    total_energy.average(params.step_average);
    pressure.average(params.step_average);
}

void Simulation::print_summary() {
    printf("%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
           step_count,
           time_now,
           v_sum.sum() / num_mols,
           total_energy.get_sum(),
           total_energy.get_sum_sqrd(),
           kinetic_energy.get_sum(),
           kinetic_energy.get_sum_sqrd(),
           pressure.get_sum(),
           pressure.get_sum_sqrd());
}

void Simulation::clear_properties() {
    kinetic_energy.clear();
    total_energy.clear();
    pressure.clear();
}
