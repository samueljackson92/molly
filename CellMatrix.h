//
// Created by Samuel Jackson on 11/03/2016.
//

#ifndef MOLLY_CELLMATRIX_H
#define MOLLY_CELLMATRIX_H

#include <vector>
#include "Cell.h"

namespace Molly {
    class CellMatrix {
    public:
        CellMatrix(size_t d1=0, size_t d2=0, size_t d3=0);
        void resize(size_t d1, size_t d2, size_t d3);
        void add_molecule(Molecule_ptr mol);
        void wrap_cells();
        Cell_ptr operator()(size_t i, size_t j, size_t k);
        Cell_ptr const operator()(size_t i, size_t j, size_t k) const;

    private:
        void create_neighbours();
        void convert_vector_to_index(const Eigen::Vector3d& vec, Eigen::Vector3i& index) const;

        size_t d1, d2, d3;
        std::vector<Cell_ptr> data;
    };
}


#endif //MOLLY_CELLMATRIX_H
