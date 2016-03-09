//
// Created by Samuel Jackson on 08/03/2016.
//
#include "CatchRequire.h"
#include "../Molecule.h"

using namespace Molly;

TEST_CASE( "Create Molecule", "[molecule]" ) {
    Molecule m;
    m.r << 10, 10, 10;
    m.rv << 0.1, 0.2, 0.3;
    m.ra << 0.01, 0.02, 0.03;
}