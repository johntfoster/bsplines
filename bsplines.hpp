#include <math.h>
#include <iostream>
#include <memory>

const unsigned END_VAL = 10;

namespace BSPLINES {

class BsplineInterface {
    public:
    virtual void evaluate(const double, const unsigned, const double*, double*) = 0;
};

template<unsigned p, unsigned k> 
//p = interpolation order, k = derivative of basis function, k = 0 returns basis
//function itself
class Bspline : public BsplineInterface  {
    public:
    virtual void evaluate(const double x, const unsigned num_knots, const double* knot_vector, double* N) override {
        evaluate_impl(x, num_knots, knot_vector, N);
    }

    static void evaluate_impl(const double x, const unsigned num_knots, const double* knot_vector, double* N) {

        auto Npm1 = std::make_unique<double[]>(num_knots - p + 1);

        Bspline<p - 1, k - 1>::evaluate_impl(x, num_knots, knot_vector, Npm1.get());

        for (int i = 0; i < num_knots - p - 1; ++i) {

            if (fabs(knot_vector[i + p] - knot_vector[i]) > 1e-15)
                N[i] = p / (knot_vector[i + p] - knot_vector[i]) * Npm1[i];
            else
                N[i] = 0.0;

            if (fabs(knot_vector[i + p + 1] - knot_vector[i + 1]) > 1e-15)
                N[i] -= p / (knot_vector[i + p + 1] - knot_vector[i + 1]) * Npm1[i + 1];
        }

        return;
    }
};

template<unsigned p> 
class Bspline<p, 0> : public BsplineInterface {
    public:
    virtual void evaluate(const double x, const unsigned num_knots, const double* knot_vector, double* N) override {
        evaluate_impl(x, num_knots, knot_vector, N);
    }

    static void evaluate_impl(const double x, const unsigned int num_knots, const double* knot_vector, double* N) {

        auto Npm1 = std::make_unique<double[]>(num_knots - p + 1);

        Bspline<p - 1, 0>::evaluate_impl(x, num_knots, knot_vector, Npm1.get());

        for (int i = 0; i < num_knots - p - 1; ++i) {

            if (fabs(knot_vector[i + p] - knot_vector[i]) > 1e-15)
                N[i] = (x - knot_vector[i]) / (knot_vector[i + p] - knot_vector[i]) * Npm1[i];
            else
                N[i] = 0.0;

            if (fabs(knot_vector[i + p + 1] - knot_vector[i + 1]) > 1e-15)
                N[i] += (knot_vector[i + p + 1] - x) / 
                        (knot_vector[i + p + 1] - knot_vector[i + 1]) * Npm1[i + 1];
            
        }
    }
};

template<> 
class Bspline<0, 0> : public BsplineInterface {
    public:
    virtual void evaluate(const double x, const unsigned num_knots, const double* knot_vector, double* N) override {
        evaluate_impl(x, num_knots, knot_vector, N);
    }

    static void evaluate_impl(const double x, const int num_knots, const double* knot_vector, double* N) {

        for (int i = 0; i < num_knots; ++i) {

            if (x > knot_vector[i] && x < knot_vector[i + 1])
                N[i] = 1.;
            else
                N[i] = 0.;
        }
    }
};

struct ThrowError
{
    static inline BsplineInterface* create (int, int)
    {
        throw std::runtime_error ("Could not create Bspline class");
    }
};

template<int DEPTH = 0, int N1 = 1, int N2 = 0>
struct Factory : ThrowError {};

template<int N2>
struct Factory<0, END_VAL, N2> : ThrowError {};

template<int N1>
struct Factory<1, N1, N1> : ThrowError {};

template<int N1, int N2>
struct Factory<0, N1, N2>
{
    static inline BsplineInterface* create(int p, int k)
    {
        if (p == N1)
        {
            return Factory<1, N1, 0>::create(p, k);
        }
        else
            return Factory<0, N1 + 1, N2>::create(p, k);
    }
};

template<int N1, int N2>
struct Factory<1, N1, N2>
{
    static inline BsplineInterface* create(int p, int k)
    {
        if (k == N2)
        {
            return new Bspline<N1, N2> ();
        }
        else
            return Factory<1, N1, N2 + 1>::create(p, k);
    }
};


} //end BSPLINES namespace
