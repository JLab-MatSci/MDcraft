#include <cmath>
#include <mdcraft/solver/potential/EAMW0644.h>
#include <algorithm>
#include <cassert>
#include <stdexcept>


namespace mdcraft {
namespace solver {
namespace potential {

EAMW0644::EAMW0644(
     const double  Rcutoff,
     const double  _aVr,
     const double  _bVr,
     const double  _rVr,
     const std::size_t  _arrsz,
     const double* _V,
     const double* _W,
     const double* _C,
     const double* _D,
     const double* _F,
     const double* _G
) : Base (Rcutoff, ::mdcraft::solver::potential::type::eam),
    aVr(_aVr),
    bVr(_bVr),
    rVr(_rVr)
{
         if (!(Rcutoff > 0) || !std::isfinite(Rcutoff) ||
        !std::isfinite(_aVr) || !std::isfinite(_bVr) || !std::isfinite(_rVr))
        throw std::invalid_argument("Invalid analytic potential parameters");
    if (_arrsz < 12 || _arrsz > arrsz)
        throw std::invalid_argument("EAMW0644: coefficient array length out of range");
    if (!(_V && _W && _C && _D && _F && _G))
        throw std::invalid_argument("EAMW0644: null coefficient array");
    for (std::size_t i = 0; i < arrsz; ++i) {
        V[i] = W[i] = 0.0; C[i] = D[i] = F[i] = G[i] = 0.0;
    }
R2cut = Rcutoff*Rcutoff;
     for (std::size_t i = 0; i < _arrsz; i++) {
          V[i] = _V[i]; W[i] = _W[i]; C[i] = _C[i];
          D[i] = _D[i]; F[i] = _F[i]; G[i] = _G[i];
     }

     // Rebuild derivatives from the energy coefficients to keep forces conservative.
     for (std::size_t i = 2; i <= 9; ++i) {
         V[i+10] = -2.0*(i+1)*V[i];
         W[i+10] = -2.0*(i-1)*W[i];
     }
     for (std::size_t i = 2; i <= 11; ++i) {
         C[i+10] = -2.0*(i+1)*C[i];
         D[i+10] = -2.0*(i-1)*D[i];
         F[i+10] = static_cast<double>(i)*F[i];
         G[i+10] = static_cast<double>(i)*G[i];
     }
}

// potential U(r)
double EAMW0644::value(const double r) {
     double u, r2 = r*r;
     value(&r2, &u, 1ul);
     return u;
}
double EAMW0644::value(const vector r) {
     double u, r2;
     r2 = r.dot(r);
     value(&r2, &u, 1ul);
     return u;
}
// virial U'(r)/r
matrix EAMW0644::virial(const vector r, const double dF) {
     double r2 = r.dot(r);
     matrix v;
     virial(&r, &r2, &dF, &v, 1ul);
     return v;
}
// force
vector EAMW0644::force(const vector r, const double dF) {
     double r2 = r.dot(r);
     vector f;
     force(&r, &r2, &dF, &f, 1ul);
     return f;
}

void EAMW0644::value(const double* r2, double* Uout, std::size_t size) {
     double X01, X02, X03, X04, X05, X06, X07, X08, X09;
     double F0, F1, G0, G1;
     double R1, R2, DV, FV, XK;
     double* FR; // out V'(r)/r

     #pragma ivdep
     #pragma vector aligned
     for (std::size_t k = 0; k < size; k++) {
          R2 = r2[k];

          X01 = R2cut - R2;
          X02 = X01*X01;
          X03 = X01*X02;
          X05 = X03*X02;
          X07 = X05*X02;
          X09 = X07*X02;

          F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01)
             + X09*(V[8] + V[9]*X01);

          F1 = X02*(V[12] + V[19]*X07)
             + X03*(V[13] + V[14]*X01)
             + X05*(V[15] + V[16]*X01)
             + X07*(V[17] + V[18]*X01);


          G0 = X01*(W[2] + W[3]*X01)
             + X03*(W[4] + W[5]*X01)
             + X05*(W[6] + W[7]*X01)
             + X07*(W[8] + W[9]*X01) + 1.0;

          G1 = (W[12] + W[13]*X01)
               + X02*(W[14] + W[19]*X05)
               + X03*(W[15] + W[16]*X01)
               + X05*(W[17] + W[18]*X01);

          FV = R2*R2*R2;                    // FV = R2**(Kpp-1)

          X02 = X02*X02;

          XK = bVr* FV*R2 + X02;          // bVr* R2**Kpp + X01**Lpp

          R1 = 1.0/R2;

          DV = R1 - rVr;

// V(r)
          Uout[k] = aVr* DV* XK* X03* X03 + F0/G0;   // aVr* X01**Mpp * DV* XK
          //Uout[k] *= 0.5;

          // FV = 2*bVr* FV* (3.0*X01 - 5.0*R2) - 14.0*X02;

// V'(r)/r
//        FR[k] = aVr2* X09* ( DV*FV - R1*R1* XK* X01)  // aVr2* X01**(Mpp-1)
//           + (F1*G0 - F0*G1)/(G0*G0);

     }
     return;
}

double EAMW0644::value(const double* r2, std::size_t size) {
     double X01, X02, X03, X04, X05, X06, X07, X08, X09;
     double F0, F1, G0, G1;
     double R1, R2, DV, FV, XK;
     double* FR; // out V'(r)/r
     double Utot = 0.0;

     #pragma ivdep
     #pragma vector aligned
     for (std::size_t k = 0; k < size; k++) {
          R2 = r2[k];

          X01 = R2cut - R2;
          X02 = X01*X01;
         X03 = X01*X02;
         X05 = X03*X02;
         X07 = X05*X02;
         X09 = X07*X02;

         F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01)
             + X09*(V[8] + V[9]*X01);

        F1 = X02*(V[12] + V[19]*X07)
             + X03*(V[13] + V[14]*X01)
             + X05*(V[15] + V[16]*X01)
             + X07*(V[17] + V[18]*X01);


        G0 = X01*(W[2] + W[3]*X01)
             + X03*(W[4] + W[5]*X01)
             + X05*(W[6] + W[7]*X01)
             + X07*(W[8] + W[9]*X01) + 1.0;

        G1 = (W[12] + W[13]*X01)
             + X02*(W[14] + W[19]*X05)
             + X03*(W[15] + W[16]*X01)
             + X05*(W[17] + W[18]*X01);

        FV = R2*R2*R2;                    // FV = R2**(Kpp-1)

        X02 = X02*X02;

        XK = bVr* FV*R2 + X02;          // bVr* R2**Kpp + X01**Lpp

        R1 = 1.0/R2;

        DV = R1 - rVr;

// V(r)
        Utot += aVr* DV* XK* X03* X03 + F0/G0;   // aVr* X01**Mpp * DV* XK

//         FV = bVr* FV*(Kpp* X01 - Mpp*R2) - (Mpp+Lpp)* X01**Lpp

//        FV = bVr2* FV* (3.0*X01 - 5.0*R2) - 14.0*X02;

// V'(r)/r
//        FR[k] = aVr2* X09* ( DV*FV - R1*R1* XK* X01)  // aVr2* X01**(Mpp-1)
//           + (F1*G0 - F0*G1)/(G0*G0);

     }
     return 0.5*Utot;
}

matrix EAMW0644::virial(
          const vector*  r,
          const double*  r2,
          const double*  dF,
          std::size_t    size
) {
     double  X01, X02, X03, X04, X05, X06, X07, X08, X09;
     double  F0, F1, G0, G1;
     double  R1, R2, DV, FV, XK, Fn;
     double  FR;
     matrix Vtot = matrix::Zero();
     #pragma ivdep
     #pragma vector aligned
     for (std::size_t k = 0; k < size; k++) {
        R2 = r2[k];

        X01 = R2cut - R2;
        X02 = X01*X01;
        X03 = X01*X02;
        X05 = X03*X02;
        X07 = X05*X02;
        X09 = X07*X02;

        F0 = X03*(C[2] + C[3]*X01)
             + X05*(C[4] + C[5]*X01)
             + X07*(C[6] + C[7]*X01)
             + X09*(C[8] + C[9]*X01 + C[10]*X02 + C[11]*X03);

        F1 = X02*(C[12] + C[21]*X09)
             + X03*(C[13] + C[14]*X01)
             + X05*(C[15] + C[16]*X01)
             + X07*(C[17] + C[18]*X01)
             + X09*(C[19] + C[20]*X01);

        G0 = X01*(D[2] + D[3]*X01)
             + X03*(D[4] + D[5]*X01)
             + X05*(D[6] + D[7]*X01)
             + X07*(D[8] + D[9]*X01)
             + X09*(D[10] + D[11]*X01) + 1.0;

        G1 = (D[12] + D[13]*X01)
             + X02*(D[14] + D[21]*X07)
             + X03*(D[15] + D[16]*X01)
             + X05*(D[17] + D[18]*X01)
             + X07*(D[19] + D[20]*X01);
// n'(r)/r
        Fn = (F1*G0 - F0*G1)/(G0*G0);

// V(r)
         F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01)
             + X09*(V[8] + V[9]*X01);

        F1 = X02*(V[12] + V[19]*X07)
             + X03*(V[13] + V[14]*X01)
             + X05*(V[15] + V[16]*X01)
             + X07*(V[17] + V[18]*X01);


        G0 = X01*(W[2] + W[3]*X01)
             + X03*(W[4] + W[5]*X01)
             + X05*(W[6] + W[7]*X01)
             + X07*(W[8] + W[9]*X01) + 1.0;

        G1 = (W[12] + W[13]*X01)
             + X02*(W[14] + W[19]*X05)
             + X03*(W[15] + W[16]*X01)
             + X05*(W[17] + W[18]*X01);


        FV = R2*R2*R2;                    // FV = R2**(Kpp-1)

        X02 = X02*X02;

        XK = bVr* FV*R2 + X02;         // bVr* R2**Kpp + X01**Lpp

        R1 = 1.0/R2;

        DV = R1 - rVr;

// V(r)/2
    //    Uout[k] = 0.5*(aVr* DV* XK* X05* X05 + F0/G0);   // aVr* X01**Mpp * DV* XK

//         FV = bVr* FV*(Kpp* X01 - Mpp*R2) - (Mpp+Lpp)* X01**Lpp

        FV = 4.0*bVr* FV* (X01 - 1.5*R2) - 10.0*X02;

// V'(r)/r
        DV = 2.0*aVr* X05* (DV*FV - R1*R1* XK* X01)  // aVr2* X01**(Mpp-1)
                + (F1*G0 - F0*G1)/(G0*G0);

// [F'(n) n'(r)+ V'(r)]/r
        FR = Fn* dF[k] + DV;

        // Vtot += matrix(FR*r[k].x*r[k].x, FR*r[k].x*r[k].y, FR*r[k].x*r[k].z,
        //                FR*r[k].y*r[k].x, FR*r[k].y*r[k].y, FR*r[k].y*r[k].z,
        //                FR*r[k].z*r[k].x, FR*r[k].z*r[k].y, FR*r[k].z*r[k].z)*0.5; // *0.5 ot not???

        Vtot(0) += FR*r[k].x()*r[k].x();
        Vtot(1) += FR*r[k].x()*r[k].y();
        Vtot(2) += FR*r[k].x()*r[k].z();

        Vtot(3) += FR*r[k].y()*r[k].x();
        Vtot(4) += FR*r[k].y()*r[k].y();
        Vtot(5) += FR*r[k].y()*r[k].z();

        Vtot(6) += FR*r[k].z()*r[k].x();
        Vtot(7) += FR*r[k].z()*r[k].y();
        Vtot(8) += FR*r[k].z()*r[k].z();

        // Apply 1/2 once after the neighbour sum.
    }
    return 0.5 * Vtot;
}

void EAMW0644::virial(
          const vector*  r,
          const double*  r2,
          const double*  dF,
          matrix*        Vout,
          std::size_t    size
) {
     double  X01, X02, X03, X04, X05, X06, X07, X08, X09;
     double  F0, F1, G0, G1;
     double  R1, R2, DV, FV, XK, Fn;
     double  FR;
     #pragma ivdep
     #pragma vector aligned
     for (std::size_t k = 0; k < size; k++) {
        R2 = r2[k];

        X01 = R2cut - R2;
        X02 = X01*X01;
        X03 = X01*X02;
        X05 = X03*X02;
        X07 = X05*X02;
        X09 = X07*X02;

        F0 = X03*(C[2] + C[3]*X01)
             + X05*(C[4] + C[5]*X01)
             + X07*(C[6] + C[7]*X01)
             + X09*(C[8] + C[9]*X01 + C[10]*X02 + C[11]*X03);

        F1 = X02*(C[12] + C[21]*X09)
             + X03*(C[13] + C[14]*X01)
             + X05*(C[15] + C[16]*X01)
             + X07*(C[17] + C[18]*X01)
             + X09*(C[19] + C[20]*X01);

        G0 = X01*(D[2] + D[3]*X01)
             + X03*(D[4] + D[5]*X01)
             + X05*(D[6] + D[7]*X01)
             + X07*(D[8] + D[9]*X01)
             + X09*(D[10] + D[11]*X01) + 1.0;

        G1 = (D[12] + D[13]*X01)
             + X02*(D[14] + D[21]*X07)
             + X03*(D[15] + D[16]*X01)
             + X05*(D[17] + D[18]*X01)
             + X07*(D[19] + D[20]*X01);
// n'(r)/r
        Fn = (F1*G0 - F0*G1)/(G0*G0);

// V(r)
         F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01)
             + X09*(V[8] + V[9]*X01);

        F1 = X02*(V[12] + V[19]*X07)
             + X03*(V[13] + V[14]*X01)
             + X05*(V[15] + V[16]*X01)
             + X07*(V[17] + V[18]*X01);


        G0 = X01*(W[2] + W[3]*X01)
             + X03*(W[4] + W[5]*X01)
             + X05*(W[6] + W[7]*X01)
             + X07*(W[8] + W[9]*X01) + 1.0;

        G1 = (W[12] + W[13]*X01)
             + X02*(W[14] + W[19]*X05)
             + X03*(W[15] + W[16]*X01)
             + X05*(W[17] + W[18]*X01);


        FV = R2*R2*R2;                    // FV = R2**(Kpp-1)

        X02 = X02*X02;

        XK = bVr* FV*R2 + X02;         // bVr* R2**Kpp + X01**Lpp

        R1 = 1.0/R2;

        DV = R1 - rVr;

// V(r)/2
    //    Uout[k] = 0.5*(aVr* DV* XK* X05* X05 + F0/G0);   // aVr* X01**Mpp * DV* XK

//         FV = bVr* FV*(Kpp* X01 - Mpp*R2) - (Mpp+Lpp)* X01**Lpp

        FV = 4.0*bVr* FV* (X01 - 1.5*R2) - 10.0*X02;

// V'(r)/r
        DV = 2.0*aVr* X05* (DV*FV - R1*R1* XK* X01)  // aVr2* X01**(Mpp-1)
                + (F1*G0 - F0*G1)/(G0*G0);

// [F'(n) n'(r)+ V'(r)]/r
        FR = Fn* dF[k] + DV;

        // Vout[k] = matrix(FR*r[k].x*r[k].x, FR*r[k].x*r[k].y, FR*r[k].x*r[k].z,
        //                  FR*r[k].y*r[k].x, FR*r[k].y*r[k].y, FR*r[k].y*r[k].z,
        //                  FR*r[k].z*r[k].x, FR*r[k].z*r[k].y, FR*r[k].z*r[k].z)*0.5; // *0.5 ot not???

        Vout[k](0) = FR*r[k].x()*r[k].x();
        Vout[k](1) = FR*r[k].x()*r[k].y();
        Vout[k](2) = FR*r[k].x()*r[k].z();

        Vout[k](3) = FR*r[k].y()*r[k].x();
        Vout[k](4) = FR*r[k].y()*r[k].y();
        Vout[k](5) = FR*r[k].y()*r[k].z();

        Vout[k](6) = FR*r[k].z()*r[k].x();
        Vout[k](7) = FR*r[k].z()*r[k].y();
        Vout[k](8) = FR*r[k].z()*r[k].z();

        Vout[k] *= 0.5; // full neighbour list: dE/dstrain convention
    }
    return;
}

void EAMW0644::virial(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
     auto& atom = *self;
     // use std::vector instead of c-array: almost no difference
     auto const nlist_size = nlist.size();

     std::vector<vector> rall    (nlist_size);
     std::vector<double> r2      (nlist_size    , 0.0);
     std::vector<double> nCharge (nlist_size + 1, 0.0);
     std::vector<double> dFemb   (nlist_size + 1, 0.0);
     std::vector<double> Femb    (nlist_size + 1, 0.0);

     #pragma ivdep
     #pragma vector aligned
     for (auto i = 0; i < nlist_size; i++) {
          auto& j = nlist[i];
          auto neib = neibs[j];

          auto const r_ji = (vector&) neib.r -
                            (vector&) atom.r;

          auto const r2ij = r_ji.dot(r_ji);

          r2[i]      = r2ij;
          rall[i]    = r_ji;
          nCharge[i] = (double&) neib.nc;
     }
     nCharge[nlist_size] = (double&) atom.nc;

     embeddingEnergy(nCharge.data(), Femb.data(), dFemb.data(), nlist_size + 1);
     if (data::owns_atom(neibs, atom)) atom.Em += Femb[nlist_size];

     std::size_t natoms = 0;
     #pragma ivdep
     #pragma vector aligned
     for (auto i = 0; i < nlist_size; i++){
          if (r2[i] > R2cut) continue;
          if (r2[i] < 1e-15) continue;

          r2[natoms] = r2[i];
          rall[natoms] = rall[i];
          dFemb[natoms] = dFemb[i] + dFemb[nlist_size]; // self

          ++natoms;
     }

     (matrix&) atom.V += virial(rall.data(), r2.data(), dFemb.data(), natoms);
     (double&) atom.Ep += value(r2.data(), natoms);

     return;
}

vector EAMW0644::force(
          const vector* r,
          const double* r2,
          const double* dF,
          std::size_t   size
) {
     double  X01, X02, X03, X04, X05, X06, X07, X08, X09;
     double  F0, F1, G0, G1, FF;
     double  R1, R2, DV, FV, XK;
     double SX = 0, SY = 0, SZ = 0; //sum of forces
     vector Ftot = vector::Zero();
     #pragma ivdep
     #pragma vector aligned
     for (std::size_t k = 0; k < size; k++) {
          // R2 = r[k].x*r[k].x + r[k].y*r[k].y + r[k].z*r[k].z;
          R2 = r2[k];

          X01 = R2cut - R2;
        X02 = X01*X01;
        X03 = X01*X02;
        X05 = X03*X02;
        X07 = X05*X02;
        X09 = X07*X02;

        F0 = X03*(C[2] + C[3]*X01)
             + X05*(C[4] + C[5]*X01)
             + X07*(C[6] + C[7]*X01)
             + X09*(C[8] + C[9]*X01 + C[10]*X02 + C[11]*X03);

        F1 = X02*(C[12] + C[21]*X09)
             + X03*(C[13] + C[14]*X01)
             + X05*(C[15] + C[16]*X01)
             + X07*(C[17] + C[18]*X01)
             + X09*(C[19] + C[20]*X01);

        G0 = X01*(D[2] + D[3]*X01)
             + X03*(D[4] + D[5]*X01)
             + X05*(D[6] + D[7]*X01)
             + X07*(D[8] + D[9]*X01)
             + X09*(D[10] + D[11]*X01) + 1.0;

        G1 = (D[12] + D[13]*X01)
             + X02*(D[14] + D[21]*X07)
             + X03*(D[15] + D[16]*X01)
             + X05*(D[17] + D[18]*X01)
             + X07*(D[19] + D[20]*X01);

// [Fn(i)+Fn(j)]*n'(r)/r
        FF = dF[k]* (F1*G0 - F0*G1)/(G0*G0);

// V(r)
         F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01)
             + X09*(V[8] + V[9]*X01);

        F1 = X02*(V[12] + V[19]*X07)
             + X03*(V[13] + V[14]*X01)
             + X05*(V[15] + V[16]*X01)
             + X07*(V[17] + V[18]*X01);


        G0 = X01*(W[2] + W[3]*X01)
             + X03*(W[4] + W[5]*X01)
             + X05*(W[6] + W[7]*X01)
             + X07*(W[8] + W[9]*X01) + 1.0;

        G1 = (W[12] + W[13]*X01)
             + X02*(W[14] + W[19]*X05)
             + X03*(W[15] + W[16]*X01)
             + X05*(W[17] + W[18]*X01);

        FV = R2*R2*R2;                    // FV = R2**(Kpp-1)

        X02 = X02*X02;

        XK = bVr* FV*R2 + X02;        // bVr* R2**Kpp + X01**Lpp

        R1 = 1.0/R2;

        DV = R1 - rVr;

//      FV = bVr* FV*(Kpp* X01 - Mpp*R2) - (Mpp+Lpp)* X01**Lpp

        FV = 4.0*bVr* FV* (X01 - 1.5*R2) - 10.0*X02;

// F'(r)/r
        FF = 2.0*aVr* X05* (DV*FV - R1*R1* XK* X01)  // aVr2* X01**(Mpp-1)
            + (F1*G0 - F0*G1)/(G0*G0) + FF;

// out: F(r)x/r, F(r)y/r, F(r)z/r
        Ftot.x() += FF* r[k].x();
        Ftot.y() += FF* r[k].y();
        Ftot.z() += FF* r[k].z();

     }

     return Ftot;
}


void EAMW0644::force(
          const vector* r,
          const double* r2,
          const double* dF,
          vector*       Fout,
          std::size_t   size
) {
     double  X01, X02, X03, X04, X05, X06, X07, X08, X09;
     double  F0, F1, G0, G1, FF;
     double  R1, R2, DV, FV, XK;
     double SX = 0, SY = 0, SZ = 0; //sum of forces
     #pragma ivdep
     #pragma vector aligned
     for (std::size_t k = 0; k < size; k++) {
          // R2 = r[k].x*r[k].x + r[k].y*r[k].y + r[k].z*r[k].z;
          R2 = r2[k];

          X01 = R2cut - R2;
        X02 = X01*X01;
        X03 = X01*X02;
        X05 = X03*X02;
        X07 = X05*X02;
        X09 = X07*X02;

        F0 = X03*(C[2] + C[3]*X01)
             + X05*(C[4] + C[5]*X01)
             + X07*(C[6] + C[7]*X01)
             + X09*(C[8] + C[9]*X01 + C[10]*X02 + C[11]*X03);

        F1 = X02*(C[12] + C[21]*X09)
             + X03*(C[13] + C[14]*X01)
             + X05*(C[15] + C[16]*X01)
             + X07*(C[17] + C[18]*X01)
             + X09*(C[19] + C[20]*X01);

        G0 = X01*(D[2] + D[3]*X01)
             + X03*(D[4] + D[5]*X01)
             + X05*(D[6] + D[7]*X01)
             + X07*(D[8] + D[9]*X01)
             + X09*(D[10] + D[11]*X01) + 1.0;

        G1 = (D[12] + D[13]*X01)
             + X02*(D[14] + D[21]*X07)
             + X03*(D[15] + D[16]*X01)
             + X05*(D[17] + D[18]*X01)
             + X07*(D[19] + D[20]*X01);

// [Fn(i)+Fn(j)]*n'(r)/r
        FF = dF[k]* (F1*G0 - F0*G1)/(G0*G0);

// V(r)
         F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01)
             + X09*(V[8] + V[9]*X01);

        F1 = X02*(V[12] + V[19]*X07)
             + X03*(V[13] + V[14]*X01)
             + X05*(V[15] + V[16]*X01)
             + X07*(V[17] + V[18]*X01);


        G0 = X01*(W[2] + W[3]*X01)
             + X03*(W[4] + W[5]*X01)
             + X05*(W[6] + W[7]*X01)
             + X07*(W[8] + W[9]*X01) + 1.0;

        G1 = (W[12] + W[13]*X01)
             + X02*(W[14] + W[19]*X05)
             + X03*(W[15] + W[16]*X01)
             + X05*(W[17] + W[18]*X01);

        FV = R2*R2*R2;                    // FV = R2**(Kpp-1)

        X02 = X02*X02;

        XK = bVr* FV*R2 + X02;        // bVr* R2**Kpp + X01**Lpp

        R1 = 1.0/R2;

        DV = R1 - rVr;

//      FV = bVr* FV*(Kpp* X01 - Mpp*R2) - (Mpp+Lpp)* X01**Lpp

        FV = 4.0*bVr* FV* (X01 - 1.5*R2) - 10.0*X02;

// F'(r)/r
        FF = 2.0*aVr* X05* (DV*FV - R1*R1* XK* X01)  // aVr2* X01**(Mpp-1)
            + (F1*G0 - F0*G1)/(G0*G0) + FF;

// out: F(r)x/r, F(r)y/r, F(r)z/r
        Fout[k].x() = FF* r[k].x(); SX += Fout[k].x();
        Fout[k].y() = FF* r[k].y(); SY += Fout[k].y();
        Fout[k].z() = FF* r[k].z(); SZ += Fout[k].z();

     }
     return;
}

void EAMW0644::force(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
     auto& atom = *self;
     // use std::vector instead of c-array: almost no difference
     auto const nlist_size = nlist.size();

     std::vector<vector> rall    (nlist_size);
     std::vector<double> r2      (nlist_size    , 0.0);
     std::vector<double> nCharge (nlist_size + 1, 0.0);
     std::vector<double> dFemb   (nlist_size + 1, 0.0);
     std::vector<double> Femb    (nlist_size + 1, 0.0);

     #pragma ivdep
     #pragma vector aligned
     for (auto i = 0; i < nlist_size; i++) {
          auto& j = nlist[i];
          auto neib = neibs[j];

          auto const r_ji = (vector&) neib.r -
                            (vector&) atom.r;

          auto const r2ij =r_ji.dot(r_ji);

          r2[i]      = r2ij;
          rall[i]    = r_ji;
          nCharge[i] = (double&) neib.nc;
     }
     nCharge[nlist_size] = (double&) atom.nc;

     embeddingEnergy(nCharge.data(), Femb.data(), dFemb.data(), nlist_size + 1);

     std::size_t natoms = 0;
     #pragma ivdep
     #pragma vector aligned
     for (auto i = 0; i < nlist_size; i++){
          if (r2[i] > R2cut) continue;
          if (r2[i] < 1e-15) continue;

          r2[natoms] = r2[i];
          rall[natoms] = rall[i];
          dFemb[natoms] = dFemb[i] + dFemb[nlist_size]; // self

          ++natoms;
     }

     (vector&) atom.f += force(rall.data(), r2.data(), dFemb.data(), natoms);
     return;
}

double EAMW0644::density(const double* r2, std::size_t size) {
     double  SD = 0.0; // sum of densities
     double  X01, X02, X03, X04, X05, X06, X07, X08, X09;
     double  R2, F0, G0;
     double n;
     #pragma ivdep
     #pragma vector aligned
     for (std::size_t k = 0; k < size; k++) {
          R2 = r2[k];

          X01 = R2cut - R2;
          X02 = X01*X01;
          X03 = X01*X02;
          X05 = X03*X02;
          X07 = X05*X02;
          X09 = X07*X02;

          F0 = X03*(C[2] + C[3]*X01)
               + X05*(C[4] + C[5]*X01)
               + X07*(C[6] + C[7]*X01)
               + X09*(C[8] + C[9]*X01 + C[10]*X02 + C[11]*X03);

          G0 = X01*(D[2] + D[3]*X01)
               + X03*(D[4] + D[5]*X01)
               + X05*(D[6] + D[7]*X01)
               + X07*(D[8] + D[9]*X01)
               + X09*(D[10] + D[11]*X01) + 1.0;

    // density n(r)
        n = F0/G0;  SD += n;
     }
     return SD;
}

void EAMW0644::density(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
     auto& atom = *self;
     std::vector<double> r2  (nlist.size());

     #pragma ivdep
     #pragma vector aligned
     for (auto i = 0; i < nlist.size(); i++) {
          auto& j = nlist[i];
          auto neib = neibs[j];

          auto const r_ji = (vector&) neib.r -
                            (vector&) atom.r;

          auto const r2ij = r_ji.dot(r_ji);

          r2[i]   = r2ij;
     }

     std::size_t natoms = 0;
     for (auto i = 0; i < nlist.size(); i++){
          if (r2[i] > R2cut) continue;
          if (r2[i] < 1e-15) continue;

          r2[natoms] = r2[i];

          ++natoms;
     }
     (double&) atom.nc += density(r2.data(), natoms);
     return;
}

// helper fucntion dRho for density derivative

void EAMW0644::dRho(const double* r2, double* dn, std::size_t size) {
     double  X01, X02, X03, X04, X05, X06, X07, X08, X09;
     double  R2, F0, F1, G0, G1;

     for (std::size_t k = 0; k < size; k++) {
          R2 = r2[k];

          X01 = R2cut - R2;
        X02 = X01*X01;
        X03 = X01*X02;
        X05 = X03*X02;
        X07 = X05*X02;
        X09 = X07*X02;

        F0 = X03*(C[2] + C[3]*X01)
             + X05*(C[4] + C[5]*X01)
             + X07*(C[6] + C[7]*X01)
             + X09*(C[8] + C[9]*X01 + C[10]*X02 + C[11]*X03);

        F1 = X02*(C[12] + C[21]*X09)
             + X03*(C[13] + C[14]*X01)
             + X05*(C[15] + C[16]*X01)
             + X07*(C[17] + C[18]*X01)
             + X09*(C[19] + C[20]*X01);

        G0 = X01*(D[2] + D[3]*X01)
             + X03*(D[4] + D[5]*X01)
             + X05*(D[6] + D[7]*X01)
             + X07*(D[8] + D[9]*X01)
             + X09*(D[10] + D[11]*X01) + 1.0;

        G1 = (D[12] + D[13]*X01)
             + X02*(D[14] + D[21]*X07)
             + X03*(D[15] + D[16]*X01)
             + X05*(D[17] + D[18]*X01)
             + X07*(D[19] + D[20]*X01);

        // n'(r)/r
        dn[k] = (F1*G0 - F0*G1)/(G0*G0);

     }
     return;
}

void EAMW0644::embeddingEnergy(
     const double*  nr,
     double*        EF,
     double*        dF,
     std::size_t    size
) {
     double  X01, X02, X03, X04, X05, X06, X07, X08, X09, X10, X11;
     double  F0, F1, G0, G1;
     #pragma ivdep
     #pragma vector always
     for (std::size_t k = 0; k < size; k++) {
        X01 = nr[k] - 1.0;

        X02 = X01*X01;  X03 = X02*X01;
        X04 = X02*X02;  X05 = X03*X02;
        X06 = X04*X02;  X07 = X04*X03;
        X08 = X04*X04;  X09 = X05*X04;
        X10 = X06*X04;  X11 = X06*X05;

        F0 = F[2]*(X02-1.0) + F[3]*(X03+1.0)
            + F[4]*(X04-1.0) + F[5]*(X05+1.0)
            + F[6]*(X06-1.0) + F[7]*(X07+1.0)
            + F[8]*(X08-1.0) + F[9]*(X09+1.0)
            + F[10]*(X10-1.0) + F[11]*(X11+1.0);

        F1 = X01*(F[12] + F[13]* X01)
            + X03*(F[14] + F[15]* X01)
            + X05*(F[16] + F[17]* X01)
            + X07*(F[18] + F[19]* X01)
            + X09*(F[20] + F[21]* X01);

        G0 = G[2]*(X02-1.0) + G[3]*(X03+1.0)
            + G[4]*(X04-1.0) + G[5]*(X05+1.0)
            + G[6]*(X06-1.0) + G[7]*(X07+1.0)
            + G[8]*(X08-1.0) + G[9]*(X09+1.0)
            + G[10]*(X10-1.0) + G[11]*(X11+1.0) + 1.0;

        G1 = X01*(G[12] + G[13]* X01)
            + X03*(G[14] + G[15]* X01)
            + X05*(G[16] + G[17]* X01)
            + X07*(G[18] + G[19]* X01)
            + X09*(G[20] + G[21]* X01);

// F(n)
        EF[k] = F0/G0;
// dF/dn
        dF[k] = (F1*G0 - F0*G1)/(G0*G0);
    }
    return;
}


} // potential
} // solver
} // mdcraft