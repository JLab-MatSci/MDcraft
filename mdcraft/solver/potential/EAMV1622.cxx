#include <cmath>
#include <mdcraft/solver/potential/EAMV1622.h>
#include <algorithm>
#include <cassert>
#include <stdexcept>


namespace mdcraft {
namespace solver {
namespace potential {

EAMV1622::EAMV1622(
     const double  Rcutoff,
     const double  _aVr,
     const double  _bVr,
     const double  _rVr,
     const std::size_t  _arrsz,
     const double* _V,
     const double* _W
) : Base (Rcutoff, ::mdcraft::solver::potential::type::pair),
    aVr(_aVr),
    bVr(_bVr),
    rVr(_rVr)
{
         if (!(Rcutoff > 0) || !std::isfinite(Rcutoff) ||
        !std::isfinite(_aVr) || !std::isfinite(_bVr) || !std::isfinite(_rVr))
        throw std::invalid_argument("Invalid analytic potential parameters");
    if (_arrsz < 10 || _arrsz > arrsz)
        throw std::invalid_argument("EAMV1622: coefficient array length out of range");
    if (!(_V && _W))
        throw std::invalid_argument("EAMV1622: null coefficient array");
    for (std::size_t i = 0; i < arrsz; ++i) {
        V[i] = W[i] = 0.0;
    }
R2cut = Rcutoff*Rcutoff;
     for (std::size_t i = 0; i < _arrsz; i++) {
          V[i] = _V[i]; W[i] = _W[i];
     }

     // Derivative tables follow from the supplied energy coefficients.
     for (std::size_t i = 2; i <= 9; ++i) {
         V[i+10] = -2.0*(i+1)*V[i];
         W[i+10] = -2.0*(i-1)*W[i];
     }
}

// potential U(r)
double EAMV1622::value(const double r) {
     double u, r2 = r*r;
     value(&r2, &u, 1ul);
     return u;
}
double EAMV1622::value(const vector r) {
     double u, r2;
     r2 = r.dot(r);
     value(&r2, &u, 1ul);
     return u;
}
// virial U'(r)/r
matrix EAMV1622::virial(const vector r, const double dF) {
     double r2 = r.dot(r);
     matrix v;
     virial(&r, &r2, &dF, &v, 1ul);
     return v;
}
// force
vector EAMV1622::force(const vector r, const double dF) {
     double r2 = r.dot(r);
     vector f;
     force(&r, &r2, &dF, &f, 1ul);
     return f;
}

void EAMV1622::value(const double* r2, double* Uout, std::size_t size) {
     double X01, X02, X03, X04, X05, X06, X07;
     double F0, F1, G0, G1;
     double R1, R2, DV, FV, XK;

     #pragma ivdep
     #pragma vector aligned

     for (std::size_t k = 0; k < size; k++) {
          R2 = r2[k];

         X01 = R2cut - R2;
         X02 = X01*X01;
         X03 = X01*X02;
         X05 = X03*X02;
         X07 = X05*X02;

        F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01 + X02*(V[8] + V[9]*X01));

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

        XK = (bVr* R2* R2 + X02)*X01;

        R1 = 1.0/R2;

        DV = R1 - rVr;

        X05 = X05*X05*X05;

// V(r)
        Uout[k] = 0.5*(aVr* DV* XK* X05 + F0/G0); // here 0.5 as usual (later will be += value)  // aVr* X01**Mpp * DV* XK

     }
}

double EAMV1622::value(const double* r2, std::size_t size) {
     double X01, X02, X03, X04, X05, X06, X07;
     double F0, F1, G0, G1;
     double R1, R2, DV, FV, XK;
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

        F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01 + X02*(V[8] + V[9]*X01));

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

        XK = (bVr* R2* R2 + X02)*X01;

        R1 = 1.0/R2;

        DV = R1 - rVr;

        X05 = X05*X05*X05;

// V(r)
        Utot += aVr* DV* XK* X05 + F0/G0; // here 0.5 as usual (later will be += value)  // aVr* X01**Mpp * DV* XK

     }
     return 0.5*Utot;
}

void EAMV1622::virial(
          const vector*  r,
          const double*  r2,
          const double*  dF,
          matrix*        Vout,
          std::size_t    size
) {
     double X01, X02, X03, X04, X05, X06, X07;
     double F0, F1, G0, G1;
     double R1, R2, DV, FV, XK;
     double FR;

     #pragma ivdep
     #pragma vector aligned

     for (std::size_t k = 0; k < size; k++) {
          R2 = r2[k];

         X01 = R2cut - R2;
         X02 = X01*X01;
         X03 = X01*X02;
         X05 = X03*X02;
         X07 = X05*X02;

        F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01 + X02*(V[8] + V[9]*X01));

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

        XK = (bVr* R2* R2 + X02)*X01;

        R1 = 1.0/R2;

        DV = R1 - rVr;

        X05 = X05*X05*X05;


        FV = 2.0*bVr* R2* (X01 - 8.0*R2) - 18.0*X02;

// V'(r)/r
        FR = 2.0*aVr* X05* ( DV*FV - R1*R1* XK)  // aVr2* X01**(Mpp-1)
               + (F1*G0 - F0*G1)/(G0*G0);

        // Vout[k] = /*0.5 **/ matrix(  FR*r[k].x()*r[k].x(), FR*r[k].x()*r[k].y(), FR*r[k].x()*r[k].z(),
        //                    FR*r[k].y()*r[k].x(), FR*r[k].y()*r[k].y(), FR*r[k].y()*r[k].z(),
        //                    FR*r[k].z()*r[k].x(), FR*r[k].z()*r[k].y(), FR*r[k].z()*r[k].z());// *0.5 ot not???

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
}

matrix EAMV1622::virial(
          const vector*  r,
          const double*  r2,
          //const double*  dF,
          std::size_t    size
) {
     double X01, X02, X03, X04, X05, X06, X07;
     double F0, F1, G0, G1;
     double R1, R2, DV, FV, XK;
     double FR;
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

        F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01 + X02*(V[8] + V[9]*X01));

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

        XK = (bVr* R2* R2 + X02)*X01;

        R1 = 1.0/R2;

        DV = R1 - rVr;

        X05 = X05*X05*X05;


        FV = 2.0*bVr* R2* (X01 - 8.0*R2) - 18.0*X02;

// V'(r)/r
        FR = 2.0*aVr* X05* ( DV*FV - R1*R1* XK)  // aVr2* X01**(Mpp-1)
               + (F1*G0 - F0*G1)/(G0*G0);

        Vtot(0) += FR*r[k].x()*r[k].x();
        Vtot(1) += FR*r[k].x()*r[k].y();
        Vtot(2) += FR*r[k].x()*r[k].z();

        Vtot(3) += FR*r[k].y()*r[k].x();
        Vtot(4) += FR*r[k].y()*r[k].y();
        Vtot(5) += FR*r[k].y()*r[k].z();

        Vtot(6) += FR*r[k].z()*r[k].x();
        Vtot(7) += FR*r[k].z()*r[k].y();
        Vtot(8) += FR*r[k].z()*r[k].z();

        // Apply half-list ownership factor once after accumulating all neighbours.

    }
    return 0.5 * Vtot;
}

void EAMV1622::virial(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
     // use std::vector instead of c-array: almost no difference
     auto& atom = *self;
     auto const nlist_size = nlist.size();

     std::vector<vector> rall    (nlist_size    , vector::Zero());
     std::vector<double> r2      (nlist_size    , 0.0);
     //std::vector<double> nCharge (nlist_size + 1, 0.0);
     //std::vector<double> dFemb   (nlist_size + 1, 0.0);
     //std::vector<double> Femb    (nlist_size + 1, 0.0);

     #pragma ivdep
     #pragma vector aligned
     for (auto i = 0; i < nlist_size; i++) {
          auto& j = nlist[i];
          auto neib = neibs[j];

          auto const r_ji = (vector&) neib.r - (vector&) atom.r;

          auto const r2ij = r_ji.dot(r_ji);

          r2[i]      = r2ij;
          rall[i]    = r_ji;
          //nCharge[i] = (double&) neib[dtype::chargeDensity];
     }
     //nCharge[nlist_size] = (double&) self[dtype::chargeDensity];

     // embeddingEnergy(nCharge.data(), Femb.data(), dFemb.data(), nlist_size + 1);
     // (double&) self[dtype::embeddingEnergy] = Femb[nlist_size];

     std::size_t natoms = 0;
     #pragma ivdep
     #pragma vector aligned
     for (auto i = 0; i < nlist_size; i++){
          if (r2[i] > R2cut) continue;
          if (r2[i] < 1e-15) continue;

          r2[natoms] = r2[i];
          rall[natoms] = rall[i];
          //dFemb[natoms] = dFemb[i] + dFemb[nlist_size]; // self

          ++natoms;
     }

     (matrix&) atom.V += virial(rall.data(), r2.data(), natoms);
     (double&) atom.Ep += value(r2.data(), natoms);

     return;
}

void EAMV1622::force(
          const vector* r,
          const double* r2,
          const double* dF,
          vector*       Fout,
          std::size_t   size
) {
     double  X01, X02, X03, X04, X05, X06, X07;
     double  F0, F1, G0, G1, FF;
     double  R1, R2, DV, FV, XK;
     double SX = 0, SY = 0, SZ = 0; //sum of forces

     #pragma ivdep
     #pragma vector aligned

     for (std::size_t k = 0; k < size; k++) {
          R2 = r2[k];

          X01 = R2cut - R2;
        X02 = X01*X01;
        X03 = X01*X02;
        X05 = X03*X02;
        X07 = X05*X02;

        F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01 + X02*(V[8] + V[9]*X01));

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

        XK = (bVr* R2* R2 + X02)*X01;

        R1 = 1.0/R2;

        DV = R1 - rVr;

        X05 = X05*X05*X05;

        FV = 2.0*bVr* R2* (X01 - 8.0*R2) - 18.0*X02;

        FV = 2.0*aVr* X05* ( DV*FV - R1*R1* XK)  // aVr2* X01**(Mpp-1)
                  + (F1*G0 - F0*G1)/(G0*G0);
// F'(r)/r
        // FF = 2.0*aVr* X05* (DV*FV - R1*R1* XK* X01)  // aVr2* X01**(Mpp-1)
        //     + (F1*G0 - F0*G1)/(G0*G0) + FF;

// out: F(r)x/r, F(r)y/r, F(r)z/r
        Fout[k].x() = FV* r[k].x(); SX += Fout[k].x();
        Fout[k].y() = FV* r[k].y(); SY += Fout[k].y();
        Fout[k].z() = FV* r[k].z(); SZ += Fout[k].z();

     }
     return;
}


vector EAMV1622::force(
          const vector* r,
          const double* r2,
          //const double* dF,
          std::size_t   size
) {
     double  X01, X02, X03, X04, X05, X06, X07;
     double  F0, F1, G0, G1, FF;
     double  R1, R2, DV, FV, XK;
     vector Ftot = vector::Zero();
     double SX = 0, SY = 0, SZ = 0; //sum of forces

     #pragma ivdep
     #pragma vector aligned

     for (std::size_t k = 0; k < size; k++) {
          R2 = r2[k];

          X01 = R2cut - R2;
        X02 = X01*X01;
        X03 = X01*X02;
        X05 = X03*X02;
        X07 = X05*X02;

        F0 = X03*(V[2] + V[3]*X01)
             + X05*(V[4] + V[5]*X01)
             + X07*(V[6] + V[7]*X01 + X02*(V[8] + V[9]*X01));

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

        XK = (bVr* R2* R2 + X02)*X01;

        R1 = 1.0/R2;

        DV = R1 - rVr;

        X05 = X05*X05*X05;

        FV = 2.0*bVr* R2* (X01 - 8.0*R2) - 18.0*X02;

        FV = 2.0*aVr* X05* ( DV*FV - R1*R1* XK)  // aVr2* X01**(Mpp-1)
                  + (F1*G0 - F0*G1)/(G0*G0);
// F'(r)/r
        // FF = 2.0*aVr* X05* (DV*FV - R1*R1* XK* X01)  // aVr2* X01**(Mpp-1)
        //     + (F1*G0 - F0*G1)/(G0*G0) + FF;

// out: F(r)x/r, F(r)y/r, F(r)z/r
        Ftot.x() += FV* r[k].x();
        Ftot.y() += FV* r[k].y();
        Ftot.z() += FV* r[k].z();

     }
     return Ftot;
}

void EAMV1622::force(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
     // use std::vector instead of c-array: almost no difference
     auto& atom = *self;
     auto const nlist_size = nlist.size();

     std::vector<vector> rall    (nlist_size);
     std::vector<double> r2      (nlist_size    , 0.0);
     // std::vector<double> nCharge (nlist_size + 1, 0.0);
     // std::vector<double> dFemb   (nlist_size + 1, 0.0);
     // std::vector<double> Femb    (nlist_size + 1, 0.0);

     #pragma ivdep
     #pragma vector aligned
     for (auto i = 0; i < nlist_size; i++) {
          auto& j = nlist[i];
          auto neib = neibs[j];

          auto const r_ji = (vector&) neib.r - (vector&) atom.r;

          auto const r2ij = r_ji.x()*r_ji.x() + r_ji.y()*r_ji.y() + r_ji.z()*r_ji.z();

          r2[i]      = r2ij;
          rall[i]    = r_ji;
          //nCharge[i] = (double&) neib[dtype::chargeDensity];
     }
     // nCharge[nlist_size] = (double&) self[dtype::chargeDensity];

     // embeddingEnergy(nCharge.data(), Femb.data(), dFemb.data(), nlist_size + 1);

     std::size_t natoms = 0;
     #pragma ivdep
     #pragma vector aligned
     for (auto i = 0; i < nlist_size; i++){
          if (r2[i] > R2cut) continue;
          if (r2[i] < 1e-15) continue;

          r2[natoms] = r2[i];
          rall[natoms] = rall[i];
          //dFemb[natoms] = dFemb[i] + dFemb[nlist_size]; // self

          ++natoms;
     }

     (vector&) atom.f += force(rall.data(), r2.data(), natoms);
     return;
}

} // potential
} // solver
} // mdcraft