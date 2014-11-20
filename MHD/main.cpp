#include <iostream>
#include <MHD.h>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

void write(double *z, double *y, int n, char *file);

int main()
{
    //plasma parcel measurements and step lengths

    int nt=800;
    int nx=500;
    double dx=.1;
    double dt=1e-4;

    //boundary conditions

    double NU=0.000;
    double *RHO, *VX,*H,*BX;
    RHO= new double [nx];
    VX= new double [nx];
    H= new double [nx];
    BX= new double [nx];
    double X=pow(-2111.3,5./3.);
    cout << X ;

    for (int i=0; i<nx;i++)
    {
        RHO[i]=1.;
        BX[i]=2.;
        H[i]=16.;
        if(i < nx/2.)
        {
            VX[i]=2.;
     //       RHO[i]=1;
        }
        else
        {
            VX[i]=-2.;
     //       RHO[i]=-1;
        }

    }
    VX[0]=-2.; //periodic conditions
    VX[nx-1]=2.;
    //RHO[0]=-1.; //periodic conditions
    //RHO[nx-1]=1.;


    //definition and setup with boundary conditions

    MHD pparcel(nt, nx, dx, dt, NU);
    pparcel.setup(RHO,H,VX,BX);
    pparcel.exp();

    double *xvalue;
    double *tvalue;
    double *yvalue;
    double *yvalue2;
    double *yvalue3;
    double *yvalue4;
    xvalue=new double[nx];
    tvalue=new double[nt];
    yvalue=new double[nx];
    yvalue2=new double[nx];
    yvalue3=new double[nx];
    yvalue4=new double[nx];



    for (int i=0;i<nx;i++)
    {
        int t =480;
        int x;
        xvalue[i]=i*dx;
        yvalue[i]=pparcel.rho[i][t];
        yvalue2[i]=pparcel.vx[i][t];
        yvalue3[i]=pparcel.h[i][t];
        yvalue4[i]=pparcel.p[i][t];
    }
   write(xvalue,yvalue,nx,"test.dat");
   write(xvalue,yvalue2,nx,"test1.dat");
   write(xvalue,yvalue3,nx,"test2.dat");
   write(xvalue,yvalue4,nx,"test3.dat");

    double *Eges;
    Eges=new double[nt];

    for(int j=0;j<nt;j++)
    {
        tvalue[j]=j*dt;

    for (int i=1;i<nx-1;i++)
    {
        Eges[j]+=pparcel.E[i][j];
    }
    }

    write(tvalue,Eges,nt,"energ.dat");

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
