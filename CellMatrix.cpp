//
// Created by Samuel Jackson on 11/03/2016.
//

#include <Eigen/Dense>
#include <set>
#include "CellMatrix.h"

using namespace Eigen;
using namespace Molly;

// offset indicies for neighbour cells for each cell
// we only need to look at half the cells at one time.
// includes 13 neighbours (half the total)
const size_t NUM_OFFSETS = 13;
const int OFFSETS[][3] {
    {1, 0, 0},
    {1, 1, 0},
    {0, 1, 0},
    {-1, 1, 0},
    {0, 0, 1},
    {1, 0, 1},
    {1, 1, 1},
    {0, 1, 1},
    {-1, 1, 1},
    {-1, 0, 1},
    {-1, -1, 1},
    {0, -1, 1},
    {1, -1, 1}
};


CellMatrix::CellMatrix(const double cell_width, size_t d1, size_t d2, size_t d3)
        : cell_width(cell_width) {
    resize(d1, d2, d3);
}

void CellMatrix::resize(size_t d1, size_t d2, size_t d3) {
    data.clear();
    data.resize(d1*d2*d3);

    for (auto iter = data.begin(); iter != data.end(); ++iter) {
        *iter = std::make_shared<Cell>();
    }

    // update the dimensions
    this->d1 = d1;
    this->d2 = d2;
    this->d3 = d3;

    // reinitialise neighbours
    create_neighbours();
}

Cell_ptr CellMatrix::operator()(size_t i, size_t j, size_t k) {
    const size_t index = i*d2*d3 + j*d3 + k;

    if (index >= data.size()) {
        std::stringstream ss;
        ss << "Index is out of points for matrix of size " << d1 << ", " << d2 << ", " << d3;
        throw std::runtime_error(ss.str());
    }

    return data[index];
}

Cell_ptr const CellMatrix::operator()(size_t i, size_t j, size_t k) const {
    const size_t index = i*d2*d3 + j*d3 + k;

    if (index >= data.size()) {
        std::stringstream ss;
        ss << "Index is out of points for matrix of size " << d1 << ", " << d2 << ", " << d3;
        throw std::runtime_error(ss.str());
    }

    return data[index];
}

void CellMatrix::create_neighbours() {
    // loop over each dimension of the cell matrix
    for(size_t i = 0; i < d1; ++i) {
        for(size_t j = 0; j < d2; ++j) {
            for(size_t k = 0; k < d3; ++k) {
                // add pointers to each neighbour of a cell in the cell matrix
                Cell_ptr cell = (*this)(i, j, k);
                for (size_t offset_num = 0; offset_num < NUM_OFFSETS; ++offset_num) {
                    // get offset cell indices
                    const int x_offset = OFFSETS[offset_num][0];
                    const int y_offset = OFFSETS[offset_num][1];
                    const int z_offset = OFFSETS[offset_num][2];

                    // correct the offsets to wrap around
                    size_t x = (i+x_offset) % d1;
                    size_t y = (j+y_offset) % d2;
                    size_t z = (k+z_offset) % d3;

                    //set pointer to neighbour cell
                    Cell_ptr neighbour = (*this)(x, y, z);
                    cell->add_neighbour(neighbour);
                }
            }
        }
    }
}

void CellMatrix::add_molecule(Molecule_ptr mol) {
    Eigen::Vector3i index;
    convert_vector_to_index(mol->r, index);
    (*this)(index(0), index(1), index(2))->add_molecule(mol);
}

void CellMatrix::convert_vector_to_index(const Vector3d& vec, Vector3i& index) const {
    index = (vec / cell_width).cast<int>();
    index(0) = index(0) % d1;
    index(1) = index(1) % d2;
    index(2) = index(2) % d3;
}

void CellMatrix::wrap_cells() {
    Vector3i index;
    int i = 0;

    // loop over every cell in the matrix
    for (auto iter = data.begin(); iter != data.end(); ++iter) {
        auto& molecules = (*iter)->get_molecules();

        // check the molecules in each cell to see if they need to move
        // to a different cell
        for(auto mol = molecules.begin(); mol != molecules.end();) {
            convert_vector_to_index((*mol)->r, index);
            // check if this molecule is in the same cell
            if(i != (index(0)*d2*d3 + index(1)*d3 + index(2))) {
                // if it isn't then move it to the correct cell mark for removal from the current.
                (*this)(index(0), index(1), index(2))->add_molecule(*mol);
                // note: must update iterator after removal
                mol = molecules.erase(mol);
            } else {
                // otherwise just skip to next molecule
                ++mol;
            }
        }
        ++i;
    }
}
