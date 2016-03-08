//
// Created by Samuel Jackson on 08/03/2016.
//

#ifndef MOLLY_MOLECULE_H
#define MOLLY_MOLECULE_H

#include <Eigen/Dense>

class Molecule {

private:
    Eigen::Vector3d r;
    Eigen::Vector3d rv;
    Eigen::Vector3d ra;
};


#endif //MOLLY_MOLECULE_H