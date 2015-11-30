#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

void eulerforward(const int N, const double dt, double* euf);
void eulerbackward(const int N, const double dt, double* eub);

int main()
{
    const double dt = M_PI/10;
    const int N = 200; //--> von t=0 bis t=20pi
    double* euf = new double[2*N];
    double* eub = new double[2*N];
    //Startbedingungen: x(0)=1, x'(0)=0
    euf[0] = 1;
    euf[N] = 0;
    eub[0] = 1;
    eub[N] = 0;
    
    eulerforward(N, dt, euf);
    eulerbackward(N, dt, eub);
    //analytische loesung: x(t)=cos(t)
    
    ofstream out("euler.txt");
    for (int i = 0; i < N; i++)
    {
        out << i*dt << "\t" << euf[i] << "\t" << eub[i]<< "\t" << cos(i*dt) << endl;
    }
    out.close();
    
    delete[] euf;
    delete[] eub;
    return 0;
    
}
    
void eulerforward(const int N, const double dt, double* euf) 
{
    for (int i = 0; i < N; i++) 
    {
        euf[i+1] = euf[i] + dt*euf[N+i]; //x(1)=x(0)+dt*x'(0)
        euf[N+i+1] = euf[N+i] - dt*euf[i]; //x'(1)=x'(0)+dt*x''(0), wobei x''=-x
    }
}

void eulerbackward(const int N, const double dt, double* eub)
{
    for (int i = 0; i < N; i++) 
    {
        eub[i+1] = (eub[i] + dt*eub[N+i]) / (1 + dt*dt);
        eub[N+i+1] = (eub[N+i] - dt*eub[i]) / (1 + dt*dt);
    }
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    