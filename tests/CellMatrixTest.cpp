//
// Created by Samuel Jackson on 14/03/2016.
//

#include "../Cell.h"
#include "../CellMatrix.h"
#include "catch.hpp"

using namespace Molly;

TEST_CASE( "Create Cell Matrix", "[cell matrix]" ) {
    CellMatrix matrix(10, 10, 10);
    Cell_ptr cell = matrix(0, 0, 0);

    REQUIRE( cell->get_molecules().size() == 0 );
    REQUIRE( cell->get_neighbours().size() == 13);
}

TEST_CASE( "Check matrix indexing", "[cell matrix]") {
    CellMatrix matrix(10, 10, 10);

    // valid indexing operations
    REQUIRE_NOTHROW( matrix(9, 9, 9) );
    REQUIRE_NOTHROW( matrix(0, 0, 0) );

    // check invalid ranges
    REQUIRE_THROWS( matrix(10, 20, 10) );
    REQUIRE_THROWS( matrix(10, 10, 10) );
}

TEST_CASE( "Resize matrix", "[cell matrix]") {
    CellMatrix matrix(10, 10, 10);

    // valid indexing operations
    REQUIRE_NOTHROW( matrix(9, 9, 9) );
    REQUIRE_NOTHROW( matrix(0, 0, 0) );
    // check invalid range
    REQUIRE_THROWS( matrix(10, 10, 10) );

    matrix.resize(20, 20, 20);

    // valid indexing operations
    REQUIRE_NOTHROW( matrix(19, 19, 19) );
    REQUIRE_NOTHROW( matrix(0, 0, 0) );
    // check invalid range
    REQUIRE_THROWS( matrix(20, 20, 20) );
}

TEST_CASE( "Add Molecule", "[cell matrix]") {
    CellMatrix matrix(10, 10, 10);

    Molecule_ptr mol = std::make_shared<Molecule>();
    mol->r(0) = 0;
    mol->r(1) = 0;
    mol->r(2) = 1.13;

    REQUIRE_NOTHROW( matrix.add_molecule(mol) );
    REQUIRE( matrix(0, 0, 1)->get_molecules().size() == 1 );
    REQUIRE( matrix(0, 0, 0)->get_molecules().size() == 0 );
}

TEST_CASE( "Move Molecule", "[cell matrix]") {
    CellMatrix matrix(10, 10, 10);

    Molecule_ptr mol = std::make_shared<Molecule>();
    mol->r(0) = 0;
    mol->r(1) = 0;
    mol->r(2) = 1.13;

    REQUIRE_NOTHROW( matrix.add_molecule(mol) );

    // check it got added to the correct place
    REQUIRE( matrix(0, 0, 1)->get_molecules().size() == 1 );

    // update the postion of the molecule
    // this should cause it to move to another cell
    mol->r(0) = 0;
    mol->r(1) = 0;
    mol->r(2) = 1.12;

    // wrap cells shifts molecules to correct place
    matrix.wrap_cells();

    // check it was moved between cells
    REQUIRE( matrix(0, 0, 1)->get_molecules().size() == 0 );
    REQUIRE( matrix(0, 0, 0)->get_molecules().size() == 1 );
}



