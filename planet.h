#include <iostream>
#include <cmath>
#include "constants.h"

using namespace std;

// class to make planet objects containing data of each planet: mass, x0, y0, vx0, vy0
class Planet {
public:
    double X0;
    double Y0;
    double vx0;
    double vy0;

    double M;

    // constructor
    Planet(double, double, double, double, double);
};


