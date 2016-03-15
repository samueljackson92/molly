//
// Created by Samuel Jackson on 10/03/2016.
//

#ifndef MOLLY_PROPERTY_H
#define MOLLY_PROPERTY_H

#include <cmath>
#include <algorithm>

class Property {

public:
    Property() : sum(0), sum_sqrd(0) {}

    void set_value(double value) {
        this->sum += value;
        this->sum_sqrd += value * value;
    }

    double get_sum() const { return  sum; }
    double get_sum_sqrd() const { return sum_sqrd; }

    void clear() {
        sum = 0;
        sum_sqrd = 0;
    }

    void average(int n) {
        sum /= n;
        sum_sqrd = std::sqrt(std::max<double>(sum_sqrd / n - (sum * sum), 0));
    }

private:
    double sum;
    double sum_sqrd;
};


#endif //MOLLY_PROPERTY_H
