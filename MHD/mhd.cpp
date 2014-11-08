#include <MHD.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
MHD::MHD()
{

}

MHD::MHD(int NT, int NX, double DX, double DT, double NU)
{
     nu = NU;
     nt = NT;
     nx = NX;
     dx = DX;
     dt = DT;

}
void MHD::setup(double *RHO,double *H,double *VX,double *BX)
{

    rho = new double*[nx];
    h = new double* [nx];
    vx = new double* [nx];
    Bx = new double*[nx];
    rhovx=new double*[nx];


    for (int i = 0; i < nx; i++)
    {
        rho[i] = new double[nt];
        h[i] = new double[nt];
        vx[i] = new double[nt];
        Bx[i] = new double[nt];
        rhovx[i] = new double[nt];
    }
    for (int i=0; i<nx; i++)
    {
        rho[i][0]=RHO[i];
        h[i][0]=H[i];
        vx[i][0]=VX[i];
        Bx[i][0]=BX[i];
        rhovx[i][0]=rho[i][0]*vx[i][0];
    }
}

void MHD::exp()
{
    double a=dt/(2*dx);
    double b=dt/(dx*dx);
    p = new double[nx];
    double P;


   for(int j=0;j<nt;j++)
    {
        for(int i=0;i<nx;i++)
        {
            p[i]=2.*pow(h[i][j],gamma); //values of pressure for current time step
        }

        for(int i=1;i<nx-1;i++) //explicit scheme
        {
            rho[i][j+1]=rho[i][j]*(1-b*(vx[i+1][j]-vx[i-1][j])-2*nu*a)+rho[i+1][j]*(-b*vx[i][j]+nu*a)+rho[i-1][j]*(b*vx[i][j]+nu*a);

            P=(p[i+1]-p[i-1])/(4*dx);
            rhovx[i][j+1]=rhovx[i][j]*(1-2*nu*a)+rhovx[i+1][j]*(-b*vx[i+1][j]+nu*a)+rhovx[i-1][j]*(+b*vx[i-1][j]+nu*a)-P;

            h[i][j+1]=h[i][j]*(1-b*(vx[i+1][j]-vx[i-1][j])-2*nu*a)+h[i+1][j]*(-b*vx[i][j]+nu*a)+h[i-1][j]*(b*vx[i][j]+nu*a);

            vx[i][j+1]=rhovx[i][j+1]/rho[i][j+1];

        }
    }
}
void MHD::write(double* x, int n, char*file)
{
ofstream resout;
resout.open(file);
for (int i=0; i<n; i++)
{
resout << setprecision(15) << setw(19) << x[i] << endl;
}
resout.close();
}
