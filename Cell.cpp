//
// Created by Samuel Jackson on 11/03/2016.
//

#include "Cell.h"

using namespace Molly;

void Cell::add_molecule(Molecule_ptr mol) {
    cell_mols.push_back(mol);
}

std::list<Molecule_ptr>& Cell::get_molecules() {
    return cell_mols;
}

void Cell::add_neighbour(Cell_ptr neighbour) {
    neighbours.push_back(neighbour);
}

std::vector<Cell_ptr> Cell::get_neighbours() {
    return neighbours;
}