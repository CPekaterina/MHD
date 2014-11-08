#ifndef MHD_H
#define MHD_H

class MHD
{
private:

    double gamma=5./3.; //ratio of specific heats for monoatomic gas

public:

    double **rho,**vx,**h,**Bx,*p,**rhovx;
    double nu, dt, dx;
    int nx, nt;
    MHD();
    MHD(int NT, int NX, double DX, double DT, double NU);


    void setup(double *RHO,double *H,double *VX,double *BX);
    void exp();
    void write(double* x, int n, char*file);


};

#endif // MHD_H
