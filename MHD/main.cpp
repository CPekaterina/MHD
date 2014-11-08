#include <iostream>
#include <MHD.h>
#include <fstream>
#include <iomanip>
using namespace std;

void write(double *z, double *y, int n, char *file);

int main()
{
    //plasma parcel measurements and step lengths

    int nt=5;
    int nx=5;
    double dx= 0.1;
    double dt=0.001;

    //boundary conditions

    double NU=0.002;
    double *RHO, *VX,*H,*BX;
    RHO= new double [nx];
    VX= new double [nx];
    H= new double [nx];
    BX= new double [nx];

    for (int i=0; i<nx;i++)
    {
        RHO[i]=1.;
        BX[i]=2.;
        H[i]=16.;
        if(i < nx/2.)
        {
            VX[i]=10.;
        }
        else
        {
            VX[i]=-10.;
        }

    }
    VX[0]=-10.; //periodic conditions
    VX[nx-1]=10.;

    //definition and setup with boundary conditions

    MHD pparcel(nt, nx, dx, dt, NU);
    pparcel.setup(RHO,H,VX,BX);
    pparcel.exp();


    return 0;

}

void write(double *z, double *y, int n, char *file)
{
ofstream resout;
resout.open(file);
for (int i=0; i<n; i++)
{
resout << setprecision(15) << setw(19) << z[i] << " " << setprecision(15) << setw(19) << y[i] << endl;
}
resout.close();
}
