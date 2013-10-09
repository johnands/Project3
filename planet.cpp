#include "planet.h"

// constructor
Planet::Planet(double dX0, double dY0, double dvx0, double dvy0, double dM) {
    X0 = dX0;
    Y0 = dY0;
    vx0 = dvx0;
    vy0 = dvy0;

    M = dM/MSun;
}

