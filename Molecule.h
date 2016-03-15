//
// Created by Samuel Jackson on 08/03/2016.
//

#ifndef MOLLY_MOLECULE_H
#define MOLLY_MOLECULE_H

#include <Eigen/Dense>

namespace Molly {

    class Molecule {

    public:
        Eigen::Vector3d r;
        Eigen::Vector3d rv;
        Eigen::Vector3d ra;
    };

    typedef std::shared_ptr<Molecule> Molecule_ptr;
}

#endif //MOLLY_MOLECULE_H
