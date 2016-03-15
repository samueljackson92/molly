//
// Created by Samuel Jackson on 11/03/2016.
//

#ifndef MOLLY_CELL_H
#define MOLLY_CELL_H

#include <list>
#include <Eigen/Dense>
#include <vector>
#include <set>

#include "Molecule.h"
#include "SystemParameters.h"

namespace Molly {
    class Cell; // forward declare class for typedef

    typedef std::shared_ptr<Cell> Cell_ptr;

    class Cell {

    public:
        void add_neighbour(Cell_ptr neighbour);
        void add_molecule(Molecule_ptr mol);
        std::list<Molecule_ptr>& get_molecules();
        std::vector<Cell_ptr> get_neighbours();

    private:
        std::vector<Cell_ptr> neighbours;
        std::list<Molecule_ptr> cell_mols;
    };
}



#endif //MOLLY_CELL_H
