#include <iostream>
#include "bsplines.hpp"

int main() {

    int p = 4;
    int k = 0;

    int num_knots = 20;
    
    double knot_vector[20] = {0, 0, 0, 0, 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5};

    double N[15];

    BSPLINES::BsplineInterface* Bspline = BSPLINES::Factory<>::create(p, k);

    Bspline->evaluate(0.5, num_knots, knot_vector, N);

    for (int i = 0; i < num_knots - p - 1; ++i) {
        std::cout << N[i]  << std::endl;
    }

    return 0;
}
