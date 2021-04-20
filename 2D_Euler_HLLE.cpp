#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>

//  HLLE MUSCL Scheme to solve 2D compressible Euler Equations (shock-tube problem)
//
//                                      U_t + F_x + G_y = 0, for (x,y,t) in (0,1)^2x(0,0.3]
//
// The MUSCL scheme is a finite volume method that can provide highly
// accurate numerical solutions for a given system,
// even in cases where the solutions exhibit shocks, discontinuities, or large gradients.
// MUSCL stands for Monotonic Upstream-centered Scheme for Conservation Laws (van Leer, 1979),
// and the term was introduced in a seminal paper by Bram van Leer (van Leer, 1979).
//
// In this code, we integrate in time using RK-2, use the HLLE flux solver, monotonized central (MC) flux limiter (van Leer, 1977)
// Computational Region (x,y)
//
//                                                  *---------*----------*
//                                                  |         |          |
//                                                  |  reg 2  |   reg 1  |
//                                                  |         |          |
//                                                  *---------*----------*
//                                                  |         |          |
//                                                  |  reg 3  |   reg 4  |
//                                                  |         |          |
//                                                  *---------*----------*
//

using namespace std;

// Find maximum of three numbers
double maxs(double a, double b, double c){
    double k;
    k = a;
    if  (b>k){
        k=b;
    }
    if  (c>k){
        k=c;
    }
    return k;
}
// Find minimum of three numbers
double mins(double a, double b, double c){
    double k;
    k = a;
    if (b<k){
        k=b;
    }
    if(c<k){
        k=c;
    }
    return k;
}

// Sign(x) function
int sgn(double v) {
    if (v < 0) return -1;
    if (v > 0) return 1;
    return 0;
}

// Minmod Limiter
double minmod(double a, double b, double c) {
    double vs[] = {abs(a), abs(b), abs(c)};
    double v = mins(vs[0],vs[1],vs[2]);
    double mm;
    double s;
    // Using Harten's generalized definition
    // minmod: zero if opposite sign, otherwise the one of smaller magnitude.
    s = (sgn(a) + sgn(b) + sgn(c))/3;

    if (abs(s) == 1.0){
        mm = s * v;
    }

    if (abs(s) != 1){
        mm = 0;
    }
    return mm;
}

template <size_t size_t>
void HLLE(double (&flux)[size_t], double qL[4],double qR[4], double gamma, int normal_x, int normal_y){
    int nx;
    int ny;
    double rL,uL,vL,pL;
    double rR,uR,vR,pR;
    double vnL;
    double vnR;
    double vn;
    double aL,HL,HR;
    double aR,RT;
    double u,v,H,a;

    // normal vectors
    nx = normal_x;
    ny = normal_y;

    // Left State
    rL = qL[0];
    uL = qL[1]/rL;
    vL = qL[2]/rL;
    vnL = uL*nx + vL*ny;
    pL = (gamma-1)*(qL[3] - rL*(uL*uL + vL*vL)/2);
    aL= sqrt((gamma*pL)/rL);
    HL = (qL[3] + pL)/rL;

    // Right State
    rR = qR[0];
    uR = qR[1]/rR;
    vR = qR[2]/rR;
    vnR= uR*nx + vR*ny;
    pR = (gamma-1)*(qR[3] - rR*(uR*uR + vR*vR)/2);
    aR = sqrt(gamma*pR/rR);
    HR = (qR[3] + pR)/rR;

    // First compute the Roe Averages
    RT = sqrt(rR/rL);
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = (HL+RT*HR)/(1+RT);
    a = sqrt((gamma-1)*(H - (u*u + v*v)/2));
    vn = u*nx + v*ny;

    double SLm;
    double SRp;

    // Wave speed estimates
    SLm = mins((vnL-aL),(vn - a),0.0);
    SRp = maxs((vnR+aR), (vn + a), 0.0);

    if ((SRp-SLm) == 0) {
        printf("DIVIDINGS BY 0: SLm = %lf, SRp = %lf  \n", SLm, SRp);
    }

    // Left and Right fluxes
    double FL[4] = {rL*vnL, rL*vnL*uL + pL*nx, rL*vnL*vL + pL*ny, rL*vnL*HL};
    double FR[4] = {rR*vnR, rR*vnR*uR + pR*nx, rR*vnR*vR + pR*ny, rR*vnR*HR};

    // Compute HLL Flux
    for(int i=0; i<=3; i++) {
        flux[i] = ((SRp)*FL[i] - (SLm) * FR[i] + (SLm)* (SRp)*(qR[i] - qL[i])) / ((SRp) - (SLm));
    }

}

double*** FLUXSOLVER(double*** q, int nx, int ny, double dx, double dy,double gamma) {
    double dqw;
    double dqe;
    double dqc;
    double dqn;
    double dqs;
    double*** res;
    double dqdx[ny][nx][4];
    double dqdy[ny][nx][4];
    double qxL[4];
    double qxR[4];
    double qyL[4];
    double qyR[4];
    double qL[4];
    double qR[4];
    double flux[4];

    res = new double**[nx];
    for (int i=0; i<ny; i++) {
        res[i] = new double*[nx];
    }
    for (int i=0; i<ny; i++) {
        for (int j = 0; j < nx; j++) {
            res[i][j] = new double[4];
        }
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < 4; k++) {
                dqdx[i][j][k] = 0;
                dqdy[i][j][k] = 0;
            }
        }
    }

    for (int k = 0; k < 4; k++) {
        qxL[k] = 0;
        qxR[k] = 0;
        qyL[k] = 0;
        qyR[k] = 0;
        qL[k]  = 0;
        qR[k]  = 0;
    }
    // Compute and limit slopes at cells (i,j)
    for (int i = 1; i <= ny - 2; i++) {
        for (int j = 1; j <= nx - 2; j++) {
            for (int k = 0; k <= 3; k++) {
                dqw = 2*(q[i][j][k] - q[i][j - 1][k])/dx;
                dqe = 2*(q[i][j + 1][k] - q[i][j][k])/dx;
                dqc = (q[i][j + 1][k] - q[i][j - 1][k]) /(2 * dx);
                dqdx[i][j][k] = minmod(dqw, dqe, dqc);

                dqs = 2 * (q[i][j][k] - q[i - 1][j][k])/ dy;
                dqn = 2 * (q[i + 1][j][k] - q[i][j][k])/ dy;
                dqc = (q[i + 1][j][k] - q[i - 1][j][k])/(2 * dy);
                dqdy[i][j][k] = minmod(dqs, dqn, dqc);
            }
        }
    }

    // Compute residuals x-direction
    for (int i = 1; i <= ny - 2; i++) {
        for (int j = 1; j <= nx - 3; j++) {
            for (int k = 0; k <= 3; k++) {
                qxL[k] = q[i][j][k] + (dqdx[i][j][k])*dx/2;
                qxR[k] = q[i][j + 1][k] - (dqdx[i][j + 1][k])*dx/2;
            }
            HLLE(flux, qxL, qxR, gamma, 1, 0);
            for (int k = 0; k <= 3; k++) {
                res[i][j][k] = res[i][j][k] + (flux[k])/ dx;
                res[i][j + 1][k] = res[i][j + 1][k] - (flux[k])/dx;

            }
        }
    }

    // Compute residuals y-direction
    for (int i = 1; i <= ny - 3; i++) {
        for (int j = 1; j <= nx - 2; j++) {
            for (int k = 0; k <= 3; k++) {
                qyL[k] = q[i][j][k] + dqdy[i][j][k]*dy/2;
                qyR[k] = q[i + 1][j][k] - dqdy[i + 1][j][k]*dy/2;
            }
            HLLE(flux, qyL, qyR, gamma, 0, 1);
            for (int k = 0; k <= 3; k++) {
                res[i][j][k] = res[i][j][k] + (flux[k])/dy;
                res[i + 1][j][k] = res[i + 1][j][k] - (flux[k])/dy;
            }
        }
    }
    // Set Boundary Conditions
    for (int j = 1; j <= nx - 2; j++) {
        for (int k = 0; k <= 3; k++) {
            qR[k] = q[ny - 2][j][k] + dqdy[ny - 2][j][k] * dy / 2;
            qL[k] = qR[k];
        }
        HLLE(flux, qL, qR, gamma, 0, 1);
        for (int k = 0; k <= 3; k++) {
            res[ny - 2][j][k] = res[ny - 2][j][k] + (flux[k]) / dy;
        }
    }

    // Flux contribution of the MOST EAST FACE: east face of cell j=ny-2
    for (int i = 1; i <= ny - 2; i++) {
        for (int k = 0; k <= 3; k++) {
            qR[k] = q[i][nx - 2][k] + dqdx[i][nx - 2][k] * dx / 2;
            qL[k] = qR[k];
        }
        HLLE(flux, qL, qR, gamma, 1, 0);
        for (int k = 0; k <= 3; k++) {
            res[i][nx - 2][k] = res[i][nx - 2][k] + (flux[k]) / dx;
        }
    }
    // Flux contribution of the MOST SOUTH FACE: south face of cells j=1.
    for (int j = 1; j <= nx - 2; j++) {
        for (int k = 0; k <= 3; k++) {
            qR[k] = q[1][j][k] - dqdy[1][j][k]*dy/2;
            qL[k] = qR[k];
        }
        HLLE(flux, qL, qR, gamma, 0, -1);
        for (int k = 0; k <= 3; k++) {
            res[1][j][k] = res[1][j][k] + (flux[k])/ dy;
        }
    }
    // Flux contribution of the MOST WEST FACE: west face of cells j=1.
    for (int i = 1; i <= ny - 2; i++) {
        for (int k = 0; k <= 3; k++) {
            qR[k] = q[i][1][k] - dqdx[i][1][k] * dx / 2;
            qL[k] = qR[k];
        }
        HLLE(flux, qL, qR, gamma, -1, 0);
        for (int k = 0; k <= 3; k++) {
            res[i][1][k] = res[i][1][k] + (flux[k]) / dx;
        }
    }
    for (int i = 0; i <=nx-1; i++) {
        for (int j = 0; j <= nx - 1; j++) {
            for (int k = 0; k <= 3; k++) {
                res[i][j][k] = res[i][j][k];
            }
        }
    }
    return res;
}
int main(int arcv, char** argv){
    double CFL = 0.50;                  // CFL
    double tEnd = 0.3;                  // Final time
    int nx = 240;                       // # of points in x
    int ny = 240;                       // # number of points in y
    double gamma = 1.4;                 // Heat Capacity Ratio
    double Lx = 1;                      // x boundary (left)
    double Ly = 1;                      // y boundary (top)
    double dx = Lx/nx;                  // x spatial step
    double dy = Ly/ny;                  // y spatial step
    double dt = 0.0011;                 // t time step
    double xc[nx];                      // x array
    double yc[ny];                      // y array


// Initialization
    xc[0]= dx/2;
    yc[0]= dy/2;
// Fill x,y array
    for (int i=1; i<ny; i++){
        xc[i] = dx/2 + dx*i;
        yc[i] = dx/2 + dy*i;
    }
// Initial Condition
    double r0[ny][nx];
    double u0[ny][nx];
    double v0[ny][nx];
    double p0[ny][nx];

    for (int i=0; i<ny; i++){
        for (int j=0; j<nx; j++){
            // First Quadrant
            if(xc[j]>=0.5 && yc[i]>=0.5){
                p0[i][j] = 1.5;
                r0[i][j] = 1.5;
                u0[i][j] = 0.0;
                v0[i][j] = 0.0;
            }
            // Second Quadrant
            if(xc[j]<0.5 && yc[i]>=0.5){
                p0[i][j] = 0.3;
                r0[i][j] = 0.5323;
                u0[i][j] = 1.206;
                v0[i][j] = 0.0;


            }
            // Third Quadrant
            if(xc[j]<0.5 && yc[i]<0.5){
                p0[i][j] = 0.029;
                r0[i][j] = 0.138;
                u0[i][j] = 1.206;
                v0[i][j] = 1.206;
            }
            // Fourth Quadrant
            if(xc[j]>=0.5 && yc[i]<0.5){
                p0[i][j] = 0.3;
                r0[i][j] = 0.5323;
                u0[i][j] = 0.0;
                v0[i][j] = 1.206;

            }
        }
    }
// Initial Total and Internal Energy
    double E0[ny][nx];
    double c0[ny][nx];
    for (int i=0; i<ny; i++){
        for (int j=0; j<nx; j++) {
            E0[i][j] = p0[i][j]/((gamma-1)*r0[i][j]) + 0.5*(u0[i][j]*u0[i][j]+ v0[i][j]*v0[i][j]);
            c0[i][j] = sqrt(gamma*p0[i][j]/r0[i][j]);
        }
    }
// Initial Solution vector Q_0 = [rho_0, rho_0*u_0, rho_0*v_0, rho_0*E_0 ]
    double Q0[ny][nx][4];
    for (int i=0; i<ny; i++){

        for (int j=0; j<nx; j++) {
            Q0[i][j][0] = r0[i][j];
            Q0[i][j][1] = r0[i][j]*u0[i][j];
            Q0[i][j][2] = r0[i][j]*v0[i][j];
            Q0[i][j][3] = r0[i][j]*E0[i][j];

        }
    }
// Introduce Ghost points
    nx +=  2;
    ny += 2;
    double*** q;
    q = new double**[nx];
    for (int i=0; i<ny; i++) {
        q[i] = new double*[nx];
    }
    for (int i=0; i<ny; i++) {
        for (int j = 0; j < nx; j++) {
            q[i][j] = new double[4];
        }
    }

    for (int i=0; i<ny-2; i++){
        for (int j=0; j<nx-2; j++) {
            q[i+1][j+1][0] = Q0[i][j][0];
            q[i+1][j+1][1] = Q0[i][j][1];
            q[i+1][j+1][2] = Q0[i][j][2];
            q[i+1][j+1][3] = Q0[i][j][3];
        }
    }

    for (int i=0; i<ny; i++){
        for(int k =0; k<4; k++){
            q[i][0][k] = q[i][1][k];
            q[i][nx-1][k] = q[i][nx-2][k];
            q[0][i][k] = q[1][i][k];
            q[ny-1][i][k] = q[ny-2][i][k];
        }
    }
    ofstream init("rho_ic.txt");
    for (int i = 1; i < nx-1 ; i++) {
        for (int j = 1; j < nx - 1; j++) {
            init << q[i][j][0] << " ";
            // cout << q[i][j][0] << " ";
        }
        init << "\n";
    }
    init.close();

    ofstream init_p("p_ic.txt");
    for (int i = 1; i < nx-1 ; i++) {
        for (int j = 1; j < nx - 1; j++) {
            init_p << (gamma-1)*q[i][j][0]*((q[i][j][3]/q[i][j][0])-0.5*( pow(q[i][j][1]/q[i][j][0],2) + pow(q[i][j][2]/q[i][j][0],2) ) ) << " ";
            // cout << q[i][j][0] << " ";
        }
        init_p << "\n";
    }
    init_p.close();
    double*** qs;
    qs = new double**[nx];
    for (int i=0; i<ny; i++) {
        qs[i] = new double*[nx];
    }
    for (int i=0; i<ny; i++) {
        for (int j = 0; j < nx; j++) {
            qs[i][j] = new double[4];
        }
    }

    double*** res;
    res = new double**[nx];
    for (int i=0; i<ny; i++) {
        res[i] = new double*[nx];
    }
    for (int i=0; i<ny; i++) {
        for (int j = 0; j < nx; j++) {
            res[i][j] = new double[4];
        }
    }

// Integrate using RK-2
    double t = 0;
    while(t<tEnd) {
        res = FLUXSOLVER(q, nx, ny, dx, dy, gamma);
        for (int y = 0; y <=nx-1 ; y++) {
            for (int x = 0; x <= nx-1; x++) {
                for (int c = 0; c <= 3; c++) {
                    qs[y][x][c] = q[y][x][c] - dt*res[y][x][c];

                }
            }
        }
        for (int y = 0; y <= ny-1; y++) {
            for (int c = 0; c <= 3; c++) {
                qs[y][0][c] = qs[y][1][c];
                qs[y][nx - 1][c] = qs[y][nx - 2][c];
            }
        }
        for (int x = 0; x <= nx-1; x++) {
            for (int c = 0; c <= 3; c++) {
                qs[0][x][c] = qs[1][x][c];
                qs[ny - 1][x][c] = qs[ny - 2][x][c];
            }
        }

        res = FLUXSOLVER(qs, nx, ny, dx, dy, gamma);
        for (int y = 0; y <= nx-1 ; y++) {
            for (int x = 0; x <= nx-1; x++) {
                for (int c = 0; c <= 3; c++) {
                    q[y][x][c] = ((q[y][x][c] + qs[y][x][c]- dt*res[y][x][c]))/2;
                }
            }
        }

        for (int y = 0; y <= ny-1; y++) {
            for (int c = 0; c <=3; c++) {
                q[y][0][c] = q[y][1][c];
                q[y][nx - 1][c] = q[y][nx - 2][c];
            }
        }
        for (int x= 0; x <= nx-1; x++) {
            for (int c = 0; c <=3; c++) {
                q[0][x][c] = q[1][x][c];
                q[ny - 1][x][c] = q[ny - 2][x][c];
            }
        }

        if (t + dt > tEnd) {
            dt = tEnd - t;
        }
        t = t + dt;
        printf("%lf \n",t);
    }

    ofstream r("rho.txt");
    for (int i = 1; i < nx-1 ; i++) {
        for (int j = 1; j < nx - 1; j++) {
            r << q[i][j][0] << " ";
            // cout << q[i][j][0] << " ";
        }
        r << "\n";
    }
    r.close();
    ofstream p("p.txt");
    for (int i = 1; i < nx-1 ; i++) {
        for (int j = 1; j < nx - 1; j++) {
            p << (gamma-1)*q[i][j][0]*((q[i][j][3]/q[i][j][0])-0.5*( pow(q[i][j][1]/q[i][j][0],2) + pow(q[i][j][2]/q[i][j][0],2) ) ) << " ";
            // cout << q[i][j][0] << " ";
        }
        p << "\n";
    }
    p.close();
    return 0;
}

