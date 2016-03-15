//
// Created by Samuel Jackson on 09/03/2016.
//

#ifndef MOLLY_SYSTEMPARAMETERS_H
#define MOLLY_SYSTEMPARAMETERS_H

#include <memory>

namespace Molly {
    class SystemParameters {
    public:
        double temperature;
        double density;
        double delta_time;
        double r_cutoff;

        int step_average;
    };
}

#endif //MOLLY_SYSTEMPARAMETERS_H
