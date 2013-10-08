#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

// Program for solving the n-body kepler problem!
// We use the following units: length~ AU, time~ years, mass ~ solar masses.
// In short this program solves the differential equations (in vector form):
// dA/dt = RHS(A), where RHS is a vector of the functions of the right hand side of the differential equation for every variable in A.

const double pi = 4*atan(1.0);
const double G = 4 * pi * pi; //G given in our units


//class for stellar objects
class Object
{
public:
    Object();
    char* name;
    double mass;
    vec position;
    vec velocity;
    Object(vec pos, vec vel, double m, char* nm);
};


Object::Object() {}

Object::Object(vec pos, vec vel, double m, char* nm):
    name (nm),
    mass (m),
    position (pos),
    velocity (vel) {
}

//class for a system of stellar objects, takes in initial conditions and evolvs the system in time
class System
{
private:
    double dt;
    vec    A; //the vector that stores all the positions and velocities

public:
    int    NumberOfObjects; //the number of objects (changes as objects are added)
    int    n;               //the number of objects in a given simulation, set when you call the initializer.
    Object ObjectList[100];

    System();
    void   AddObject(vec pos, vec vel, double m, char* name);
//    double CalculateDistance(Object ob1,Object ob2);
    void   initializeSolver();
    vec    CalculateForces(vec);
    void   Solve(double);
    void   RK4();
};

System::System():
    dt (0.001),
    NumberOfObjects (0)
{
}

void System::AddObject(vec pos, vec vel, double m, char* name) // adds an object to our system
{
    Object x = Object(pos,vel,m, name);
    ObjectList[NumberOfObjects] = x;
    this->NumberOfObjects+=1;
}



void System::initializeSolver() {   //sets up the system to be ready for computation
    this->n = this->NumberOfObjects;
    this->A = zeros(4*n);

    int c = 4;
    for (int i = 0; i < n; i++) {
        this->A(i*c + 0) = this->ObjectList[i].position(0);
        this->A(i*c + 1) = this->ObjectList[i].position(1);
        this->A(i*c + 2) = this->ObjectList[i].velocity(0);
        this->A(i*c + 3) = this->ObjectList[i].velocity(1);
    }
    cout << A << endl;
}

vec System::CalculateForces(vec B) { //Computes the right hand side of the differential equations

    vec Forces = zeros(2 * n);
    vec ki      = zeros(4 * n);



    // Compute r and F between bodies.
    for (int j = 0; j < n; j++) {
        for (int k = j+1; k < n; k++) {

            double x = B[4*j]   - B[4*k];
            double y = B[4*j+1] - B[4*k+1];
            double r = sqrt(x*x + y*y);

            double f = -(G * ObjectList[j].mass * ObjectList[k].mass) / (r * r * r);
            Forces[2*j+0] += f * x;      // x-component of the force.
            Forces[2*j+1] += f * y;      // y-component of the force.
            Forces[2*k+0] -= Forces[2*j+0];    // Newtons third law, bitches.
            Forces[2*k+1] -= Forces[2*j+1];    // Newtons third law, bitches.
        }
    }


    for (int j = 0; j < n; j++) {
        double m = ObjectList[j].mass;
        ki[4*j+2] = Forces[2*j+0] / m;   // Compute dv.
        ki[4*j+3] = Forces[2*j+1] / m;

        ki[4*j+0] = B[4*j+2];   // Compute dx.
        ki[4*j+1] = B[4*j+3];
    }

    return ki;
}


void System::Solve(double numberOfYears) { //takes in the number of years to do simulation

    fstream outFile;
    outFile.open("planetData.dat", ios::out);

    double t = 0; //time (in years)
    int index = 0;
    while (t < numberOfYears) {

        RK4();

        t+=dt;
        index++;


        // Writes the position of all objects to file.
        for (int i=0;i<n; i++) {
            if (index % 10 == 0) {
                outFile << this->A(4*i+0) << " " << this->A(4*i+1) << " ";
            }
        }
        outFile << endl;
    }

    outFile.close();
}



void System::RK4() { //implements the Runge Kutta-4 method to advance one timestep

    vec k1(4*n), k2(4*n), k3(4*n), k4(4*n);

    k1 = CalculateForces(A) * dt;
    k2 = CalculateForces(A + 0.5 * k1) * dt;
    k3 = CalculateForces(A + 0.5 * k2) * dt;
    k4 = CalculateForces(A + k3) * dt;
    A += (1.0/6) * (k1 + 2 * (k2 + k3) + k4);
}


int main()
{
    System test;
    vec pos1 = zeros(2);
    pos1(0) = 0;
    pos1(1) = 0;

    vec pos2 = zeros(2);
    pos2(0) = 1;
    pos2(1) = 0;

    vec vel1 = zeros(2);
    vel1(0) = 0;
    vel1(1) = 0;

    vec vel2 = zeros(2);
    vel2(0) = 0;
    vel2(1) = 2*pi;

    test.AddObject(pos1, -vel2/2.0, 1, "Sun");
    test.AddObject(pos2, vel2/2.0, 1, "Sun2");
    test.AddObject(pos2*4, vel2/1.5, 3e-10, "Planet");
    test.initializeSolver();
    test.Solve(10);

    return 0;
}

//double System::CalculateDistance(Object ob1, Object ob2)
//{
//    double distance = 0;
//    double distx = ob1.position[0]-ob2.position[0];
//    double disty = ob1.position[1]-ob2.position[1];
//    distance += distx*distx;
//    distance += disty*disty;
//    return sqrt(distance);
//}


//    Object earth(pos1,pos1,1, "test");
//    Object moon(pos2,vel2,0.1, "test");

//    cout << earth.position[0]<< endl;
//    cout << test.CalculateDistance(earth,moon) << endl;
//    cout << test.NumberOfObjects << endl;

//    cout << test.NumberOfObjects << endl;
//    cout << test.ObjectList[0].mass << endl;
