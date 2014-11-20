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
    E=new double*[nx];
    p = new double*[nx];


    for (int i = 0; i < nx; i++)
    {
        rho[i] = new double[nt];
        h[i] = new double[nt];
        vx[i] = new double[nt];
        Bx[i] = new double[nt];
        rhovx[i] = new double[nt];
        E[i] = new double[nt];
        p[i] = new double[nt];
    }


    for (int i=0; i<nx; i++)
    {
        h[i][0]=H[i];
        p[i][0]=2.*pow(h[i][0],gamma);
        rho[i][0]=RHO[i];
        vx[i][0]=VX[i];
        Bx[i][0]=BX[i];
        rhovx[i][0]=rho[i][0]*vx[i][0];
        E[i][0]=0.5*rho[i][0]*vx[i][0]*vx[i][0]+p[i][0]/(gamma-1);
    }


}
int sgn(double x)
{
    if (x==0)
        return 0;
    if (x>0)
        return  1 ;
    else
    {return -1;}
}

void MHD::exp()
{
    double b=dt/(2.*dx);
    double a=dt/(dx*dx);
    double P;


   for(int j=0;j<nt;j++)
    {
       if((j%50)==0)
       {
            smooth(nx,j,vx);
            smooth(nx,j,rho);
            smooth(nx,j,h);
       }

       for(int i=0;i<nx;i++)
        {
           p[i][j]=2.*pow(h[i][j],gamma); //values of pressure for current time step
        }

        for(int i=1;i<nx-1;i++) //explicit scheme
        {

            E[i][j]=0.5*rho[i][j]*vx[i][j]*vx[i][j]+p[i][j]/(gamma-1);

          //rho[i][j+1]=rho[i][j]*(1.-2.*a)+rho[i+1][j]*(a-b)+rho[i-1][j]*(a+b);

            rho[i][j+1]=rho[i][j]*(1.-b*(vx[i+1][j]-vx[i-1][j])-2.*nu*a)+rho[i+1][j]*(-b*vx[i][j]+nu*a)+rho[i-1][j]*(b*vx[i][j]+nu*a);
//          rho[i][j+1]=-b*(rho[i+1][j]*vx[i+1][j]-rho[i-1][j]*vx[i-1][j])+rho[i][j];

            P=(p[i+1][j]-p[i-1][j])*dt/(4.*dx);
            rhovx[i][j+1]=rhovx[i][j]*(1.-2.*nu*a)+rhovx[i+1][j]*(-b*vx[i+1][j]+nu*a)+rhovx[i-1][j]*(b*vx[i-1][j]+nu*a)-P;

            h[i][j+1]=h[i][j]*(1.-b*(vx[i+1][j]-vx[i-1][j])-2.*nu*a)+h[i+1][j]*(-b*vx[i][j]+nu*a)+h[i-1][j]*(b*vx[i][j]+nu*a);
         // h[i][j+1]=-b*(h[i+1][j]*vx[i+1][j]-h[i-1][j]*vx[i-1][j])+h[i][j];


            vx[i][j+1]=rhovx[i][j+1]/rho[i][j+1];

        }

        rho[0][j+1]=rho[nx-2][j+1];
        vx[0][j+1]=vx[nx-2][j+1];
        h[0][j+1]=h[nx-2][j+1];

        rho[nx-1][j+1]=rho[1][j+1];
        vx[nx-1][j+1]=vx[1][j+1];
        h[nx-1][j+1]=h[1][j+1];

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

void MHD::smooth(int nx,int j,double **B)
{

        int i=0;

        B[i][j]=(B[nx-2][j]+B[nx-1][j]+B[i][j]+B[i+1][j]+B[i+2][j])/5.;

        i=1;

        B[i][j]=(B[nx-1][j]+B[i-1][j]+B[i][j]+B[i+1][j]+B[i+2][j])/5.;

        i=nx-2;

        B[i][j]=(B[i-2][j]+B[i-1][j]+B[i][j]+B[i+1][j]+B[0][j])/5.;

        i=nx-1;

        B[i][j]=(B[i-2][j]+B[i-1][j]+B[i][j]+B[0][j]+B[1][j])/5.;

     for(i=2;i<nx-2;i++)
        {
            B[i][j]=(B[i-2][j]+B[i-1][j]+B[i][j]+B[i+1][j]+B[i+2][j])/5.;
        }

}

