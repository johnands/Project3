#include <iostream>
#include <armadillo>
#include <cmath>
#include "constants.h"
#include "planet.h"
#include <fstream>

using namespace std;
using namespace arma;

// function declarations
vec f(vec A, Planet **B);
double kin_energy(vec A, Planet **B);
double pot_energy(vec A, Planet **B);
double ang_mom(vec A, Planet **B);

// number of objects
int n = 10;

int main()
{
    // initialization
    double t_min = 0.0;
    double t_max = 500.0;
    int N = 100000;
    double h = (t_max - t_min)/N;
    double h2 = h/2.0;

    // make pointers to planet objects
    // data from: http://nssdc.gsfc.nasa.gov/planetary/factsheet/
    Planet *Sun = new Planet(0.0, 0.0, 0.0, 0.0, MSun);
    Planet *Mercury = new Planet(0.387, 0.0, 0.0, 3.215*pi, 3.3e23);
    Planet *Venus = new Planet(0.723, 0.0, 0.0, 2.345*pi, 4.87e24);
    Planet *Earth = new Planet(1.0, 0.0, 0.0, 2.0*pi, 5.97e24);
    Planet *Mars = new Planet(1.523, 0.0, 0.0, 1.617*pi, 6.42e23);
    Planet *Jupiter = new Planet(5.205, 0.0, 0.0, 0.879*pi, 1.898e27);
    Planet *Saturn = new Planet(9.582, 0.0, 0.0, 0.651*pi, 5.68e26);
    Planet *Uranus = new Planet(19.20, 0.0, 0.0, 0.456*pi, 8.68e25);
    Planet *Neptune = new Planet(30.05, 0.0, 0.0, 0.362*pi, 1.02e26);
    Planet *Pluto = new Planet(39.24, 0.0, 0.0, 0.315*pi, 1.31e22);

    // declaring vectors
    // A is the system state vector
    vec A(4*n);

    //  B is an array of pointers to planet objects
    Planet **B = new Planet*[n];

    switch(n) {
    case 2:
        B[0] = Sun;
        B[1] = Earth;
        break;
    case 3:
        B[0] = Sun;
        B[1] = Earth;
        B[2] = Jupiter;
        break;
    case 10:
        B[0] = Sun;
        B[1] = Mercury;
        B[2] = Venus;
        B[3] = Earth;
        B[4] = Mars;
        B[5] = Jupiter;
        B[6] = Saturn;
        B[7] = Uranus;
        B[8] = Neptune;
        B[9] = Pluto;
        break;
    }

    // filling A, has to be the same dimension as dAdt (dim = 4*n)
    for (int i=0; i<n; i++) {
        A(i*4) = B[i]->X0;
        A(i*4+1) = B[i]->Y0;
        A(i*4+2) = B[i]->vx0;
        A(i*4+3) = B[i]->vy0;
    }

    // calculate initial energy and angular momentum
    double K, U, Lz;
    K = kin_energy(A, B);
    U = pot_energy(A, B);
    Lz = ang_mom(A, B);

    // write inital positions of all objects and Earth's energymom to file
    fstream outFile1, outFile2;
    outFile1.open("positions.dat", ios::out);
    outFile2.open("energymom.dat", ios::out);

    outFile1 << n << endl;
    for (int i=0; i<n; i++) outFile1 << A(4*i) << " " << A(4*i+1) << " ";
    outFile1 << endl;

    outFile2 << K << " " << U << " " << Lz << endl;

    // RK4 vectorized
    vec k1(4*n), k2(4*n), k3(4*n), k4(4*n);
    for (int i=0; i<N; i++) {
        k1 = f(A, B);
        k2 = f(A + k1*h2, B);
        k3 = f(A + k2*h2, B);
        k4 = f(A + k3*h, B);
        A += (h/6.0)*(k1 + 2*k2 + 2*k3 + k4);

        // calculate energy and angular momentum
        K = kin_energy(A, B);
        U = pot_energy(A, B);
        Lz = ang_mom(A, B);

        // write position values to file
        for (int j=0; j<n; j++) outFile1 << A(4*j) << " " << A(4*j+1) << " ";
        outFile1 << endl;

        // write energymom of Earth to file
        outFile2 << K << " " << U << " " << Lz << endl;

    }
    outFile1.close();
    outFile2.close();

    return 0;

}

// function that computes derivate of vector A
vec f(vec A, Planet **B)
{
    vec dAdt(4*n);  // derivative vector
    // setting dxdt and dydt for each planet
    // which in each case is just the velocity of that planet
    for (int i=0; i<n; i++) {
        dAdt(4*i) = A(4*i+2);   // dx/dt = vx
        dAdt(4*i+1) = A(4*i+3); // dy/dt = vy

        // finding the accelerations
        double ax, ay;
        ax = ay = 0.0;
        // i is current planet, j is the other planets
        // must find the force on planet i from the other planets j
        for (int j=0; j<n; j++) {
            if (i != j) {
                double r = sqrt(pow(A(4*j)-A(4*i),2) + pow(A(4*j+1)-A(4*i+1),2));
                ax += -(G*B[j]->M*(A(4*i)-A(4*j)))/(r*r*r);
                ay += -(G*B[j]->M*(A(4*i+1)-A(4*j+1)))/(r*r*r);

            }
        }
        dAdt(4*i+2) = ax;
        dAdt(4*i+3) = ay;

    }

return dAdt;
}

// function to find Earth's kinetic energy
double kin_energy(vec A, Planet **B)
{
    // finding the energies and angular momentum of only Earth first, write the general thing later maybe
    double K, vx, vy;
    if (n == 2 || n == 3) {
        vx = A(6); vy = A(7);
    }
    else {
        vx = A(14); vy = A(15);
    }
    K = 0.5*B[1]->M*(pow(vx,2) + pow(vy,2));
    return K;
}

// function to find Earth's potential energy
double pot_energy(vec A, Planet **B)
{
    double r, U;
    if (n == 2 || n == 3) {
        r = sqrt(pow(A(4)-A(0),2) + pow(A(5)-A(1),2));
        U = -(G*B[0]->M*B[1]->M)/r;
    }
    else {
        r = sqrt(pow(A(12)-A(0),2) + pow(A(13)-A(1),2));
        U = -(G*B[0]->M*B[3]->M)/r;
    }
    return U;
}

// function to find Earth's angular momentum in z-direction
double ang_mom(vec A, Planet **B)
{
    double Lz;
    if (n == 2 || n == 3) Lz = B[1]->M*(A(4)*A(7) - A(5)*A(6));
    else Lz = B[3]->M*(A(12)*A(15) - A(13)*A(14));
    return Lz;

}








