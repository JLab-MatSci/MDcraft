#include <cmath>
#include <mdcraft/solver/potential/EAMMOD1044.h>
#include <algorithm>
#include <cassert>
#include <stdexcept>


namespace mdcraft {
namespace solver {
namespace potential {

EAMMOD1044::EAMMOD1044(
    const double       Rcutoff,
    const double       _aVr,
    const double       _bVr,
    const double       _rVr,
    const std::size_t  _arrsz,
    const double*      _V,
    const double*      _W,
    const double*      _C,
    const double*      _D,
    const double*      _F,
    const double*      _G
) : Base (Rcutoff, potential::type::eam),
    aVr(_aVr),
    bVr(_bVr),
    rVr(_rVr)
{
        if (!(Rcutoff > 0) || !std::isfinite(Rcutoff) ||
        !std::isfinite(_aVr) || !std::isfinite(_bVr) || !std::isfinite(_rVr))
        throw std::invalid_argument("Invalid analytic potential parameters");
    if (_arrsz < 12 || _arrsz > arrsz)
        throw std::invalid_argument("EAMMOD1044: coefficient array length out of range");
    if (!(_V && _W && _C && _D && _F && _G))
        throw std::invalid_argument("EAMMOD1044: null coefficient array");
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
double EAMMOD1044::value(const double r) {
    double u, r2 = r*r;
    value(r2, u);
    return u;
}
double EAMMOD1044::value(const vector r) {
    double u, r2;
    r2 = r.squaredNorm();
    value(r2, u);
    return u;
}
// virial U'(r)/r
matrix EAMMOD1044::virial(const vector r, const double dF) {
    double r2 = r.squaredNorm();
    matrix v;
    virial(r, r2, dF, v);
    return v;
}
// force
vector EAMMOD1044::force(const vector r, const double dF) {
    double r2 = r.squaredNorm();
    vector f;
    force(r, r2, dF, f);
    return f;
}

void EAMMOD1044::eFG(const double r2, double& F0, double& F1, double& G0, double& G1) {
    double X01, X02, X03, X04, X05, X06, X07, X08, X09;

    X01 = R2cut - r2;
    assert(X01 >= 0);
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

    return;
}

void EAMMOD1044::value(const double r2, double& Uout) {
    double X01, X02, X03, X04, X05, X06, X07, X08, X09;
    double F0, F1, G0, G1;
    double R1, DV, FV, XK;

    eFG(r2, F0, F1, G0, G1);

    X01 = R2cut - r2;
    assert(X01 >= 0);
    X02 = X01*X01;
    X03 = X01*X02;
    X05 = X03*X02;
    X07 = X05*X02;
    X09 = X07*X02;

    FV = r2*r2*r2;                    // FV = R2**(Kpp-1);

    X02 = X02*X02;

    XK = bVr* FV*r2 + X02;          // bVr* R2**Kpp + X01**Lpp

    R1 = 1.0/r2;

    DV = R1 - rVr;

    Uout = 0.5 * (aVr* DV* XK* X05* X05 + F0/G0);   // aVr* X01**Mpp * DV* XK
    return;
}

void EAMMOD1044::fill_relpos(
    Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist,
    std::vector<vector>& rall, std::vector<double>& r2
) {
    auto& atom = *self;
    if (atom.kind >= 2) throw std::out_of_range("Modified EAM requires kinds 0 and 1");
    auto const nlist_size = nlist.size();

    #pragma ivdep
    #pragma vector aligned
    for (auto i = 0; i < nlist_size; i++) {
        auto& j = nlist[i];
        auto& neib = neibs[j];

        vector const r_ji = neib.r - atom.r;

        double const r2ij = r_ji.squaredNorm();
        if (r2ij < 1e-8) continue;

        r2[i]      = r2ij;
        rall[i]    = r_ji;
    }
}

void EAMMOD1044::HydrogenCase(
    Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist,
    const std::vector<vector>& rall, const std::vector<double>& r2
) {
    auto& atom = *self;
    if (atom.kind >= 2) throw std::out_of_range("Modified EAM requires kinds 0 and 1");
    assert(1 == atom.kind);

    auto const nlist_size = nlist.size();

    double EF  = 0;

    embeddingEnergy(atom.nc_mod[0], EF);

    // add embedding energy only if the call was from real atoms, not ghosts
    if (data::owns_atom(neibs, atom))
        atom.Em += EF;
    return;
}

void EAMMOD1044::HeavyAtomCase(
    Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist,
    const std::vector<vector>& rall, const std::vector<double>& r2
) {
    auto& atom = *self;
    if (atom.kind >= 2) throw std::out_of_range("Modified EAM requires kinds 0 and 1");
    assert(0 == atom.kind);
    auto const nlist_size = nlist.size();

    double nc_  = 0;
    double EF_  = 0;
    double dFi = 0;

    // The C embedding term is PER C-H bond, including neighbours in the
    // separate MPI alien list and physical boundary-image bonds. Unlike H,
    // these terms must NOT be restricted to the local-neighbour pass.

    #pragma ivdep
    #pragma vector aligned
    for (auto i = 0; i < nlist_size; i++) {
        double const r2ij = r2[i];
        if (r2ij < 1e-8) continue;
        // embedding energy on the current heavy atom is accumulated on the fly
        if (r2ij <= R2cut) {
            density(r2ij, nc_);
            embeddingEnergy(nc_, EF_);
            atom.Em += EF_;
        }
    }
}

void EAMMOD1044::virial(
    const vector  r,
    const double  r2,
    const double  dF_sum,
    matrix&       Vout
) {
    double  X01, X02, X03, X04, X05, X06, X07, X08, X09, X10;
    double r4, r6, r8;
    double  F0, F1, G0, G1;
    double  R1, DV, FV, XK, n, dn, dF, FF;
    double  FR;

    assert(r2 != 0);

    density(r2, n);
    dRho(r2, dn);
    eFG(r2, F0, F1, G0, G1);

    embeddingEnergyDeriv(n, dF);

    r4 = r2*r2;
    r6 = r4*r2;
    r8 = r6*r2;

    X01 = R2cut - r2;
    assert(X01 >= 0);
    X02 = X01*X01;
    X03 = X01*X02;
    X04 = X01*X03;
    X05 = X03*X02;
    X07 = X05*X02;
    X09 = X07*X02;
    X10 = X01*X09;

    FV = r6;

    X02 = X02*X02;

    XK = bVr* FV*r2 + X02;

    R1 = 1.0/r2;

    DV = R1 - rVr;

    FV = 2.0*bVr* FV* (2.0*X01 - 5.0*r2) - 14.0*X02;

    double dVdr_over_r = aVr * X09 * (
          (-2.0/r4) * X01 * (bVr * r8 + X04)
        + DV * ( 8*X01*(bVr*r6 - X03) - 20 * (bVr * r8 + X04) )
    );
    FF = dVdr_over_r + dn * (dF + dF_sum) + (F1*G0 - F0*G1)/(G0*G0);

    Vout = 0.5 * FF * r * r.transpose(); // dE/dstrain convention, full list

    return;
}

void EAMMOD1044::virial(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
    auto const nlist_size = nlist.size();
    // Do not return for an empty list: a local H atom can have only alien C
    // neighbours, while its F(rho_C) must still be counted ONCE here.

    auto& atom = *self;
    if (atom.kind >= 2) throw std::out_of_range("Modified EAM requires kinds 0 and 1");
    if (!std::all_of(nlist.cbegin(), nlist.cend(), [&](auto i) { return neibs[i].kind < 2 && neibs[i].kind != atom.kind; }))
        throw std::invalid_argument("Cross modified EAM requires unlike binary kinds");

    std::vector<vector> rall    (nlist_size    , vector::Zero());
    std::vector<double> r2      (nlist_size    , 0.0);

    double Uout = 0;
    matrix Vout = matrix::Zero();

    fill_relpos(self, neibs, nlist, rall, r2);

    /* ADD EMBEDDING ENERGY */
    switch (atom.kind) {
        case 0:
            HeavyAtomCase(self, neibs, nlist, rall, r2);
            break;
        case 1:
            HydrogenCase( self, neibs, nlist, rall, r2);
            break;
        default:
            throw std::runtime_error("Only two atomic kinds are supported for this potential\n");
    }

    // On H: same summed C density for every bond.
    // On C: each neighbouring H has its own summed C density.
    double dFemb_sum = 0.0;
    if (atom.kind == 1)
        embeddingEnergyDeriv(atom.nc_mod[0], dFemb_sum);

    for (std::size_t i = 0; i < nlist_size; ++i) {
        if (r2[i] < 1e-8 || r2[i] > R2cut) continue;
        if (atom.kind == 0) {
            const auto& neib = neibs[nlist[i]];
            embeddingEnergyDeriv(neib.nc_mod[0], dFemb_sum);
        }
        value(r2[i], Uout);
        virial(rall[i], r2[i], dFemb_sum, Vout);
        atom.Ep += Uout;
        atom.V  += Vout;
    }

    return;
}

void EAMMOD1044::force(
    const vector r,
    const double r2,
    const double dF_sum,
    vector&      Fout
) {
    double X01, X02, X03, X04, X05, X06, X07, X08, X09, X10;
    double r4, r6, r8;
    double F0, F1, G0, G1;
    double R1, dn, DV, FV, XK, FF, dF, n;


    density(r2, n);
    dRho(r2, dn);

    eFG(r2, F0, F1, G0, G1);

    embeddingEnergyDeriv(n, dF);

    r4 = r2*r2;
    r6 = r4*r2;
    r8 = r6*r2;

    X01 = R2cut - r2;
    assert(X01 >= 0);
    X02 = X01*X01;
    X03 = X01*X02;
    X04 = X01*X03;
    X05 = X03*X02;
    X07 = X05*X02;
    X09 = X07*X02;
    X10 = X01*X09;

    FV = r6;

    X02 = X02*X02;

    XK = bVr* FV*r2 + X02;

    R1 = 1.0/r2;

    DV = R1 - rVr;

    FV = 2.0*bVr* FV* (2.0*X01 - 5.0*r2) - 14.0*X02;

    double dVdr_over_r = aVr * X09 * (
          (-2.0/r4) * X01 * (bVr * r8 + X04)
        + DV * ( 8*X01*(bVr*r6 - X03) - 20 * (bVr * r8 + X04) )
    );

    FF =  dVdr_over_r + (F1*G0 - F0*G1)/(G0*G0) + dn * (dF + dF_sum);

    Fout = FF * r;

    return;
}

void EAMMOD1044::force(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
    auto const nlist_size = nlist.size();
    if (nlist_size == 0) return;

    auto& atom = *self;
    if (atom.kind >= 2) throw std::out_of_range("Modified EAM requires kinds 0 and 1");
    if (!std::all_of(nlist.cbegin(), nlist.cend(), [&](auto i) { return neibs[i].kind < 2 && neibs[i].kind != atom.kind; }))
        throw std::invalid_argument("Cross modified EAM requires unlike binary kinds");

    std::vector<vector> rall    (nlist_size);
    std::vector<double> r2      (nlist_size, 0.0);

    vector Fout = vector::Zero();

    fill_relpos(self, neibs, nlist, rall, r2);

    // for each H atom we nned to precalculate the sum of n for all its C neibs

    double dFemb_sum = 0;
    double n_sum = 0;


    if (1 == atom.kind) {
        // n_sum is already calulated on each H atom by definition!
        n_sum = atom.nc_mod[0];
        embeddingEnergyDeriv(n_sum, dFemb_sum);
    }

    #pragma ivdep
    #pragma vector aligned
    for (std::size_t i = 0; i < nlist_size; i++) {

        if (0 == atom.kind) {
            // n_sum is already calulated on each H atom by definition!
            n_sum = neibs[nlist[i]].nc_mod[0];
            embeddingEnergyDeriv(n_sum, dFemb_sum);
        }

        if (r2[i] < 1e-8) continue;
        if (r2[i] > R2cut) continue;
        force(rall[i], r2[i], dFemb_sum, Fout);
        atom.f  += Fout;
    }

    return;
}

void EAMMOD1044::density(const double r2, double& n) {
    double X01, X02, X03, X04, X05, X06, X07, X08, X09;
    double F0, G0;

    X01 = R2cut - r2;
    assert(X01 >= 0);
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
    n = F0/G0;
}

void EAMMOD1044::density(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
    auto const nlist_size = nlist.size();
    if (nlist_size == 0) return;

    auto& atom = *self;
    if (atom.kind >= 2) throw std::out_of_range("Modified EAM requires kinds 0 and 1");
    if (!std::all_of(nlist.cbegin(), nlist.cend(), [&](auto i) { return neibs[i].kind < 2 && neibs[i].kind != atom.kind; }))
        throw std::invalid_argument("Cross modified EAM requires unlike binary kinds");

    if (0 == atom.kind) // if not H atom (1 is assumed to always be H)
        return;

    std::vector<double> r2(nlist_size);

    double nc_ = 0;

    #pragma ivdep
    #pragma vector aligned
    for (std::size_t i = 0; i < nlist_size; i++) {
        auto& neib = neibs[nlist[i]];

        vector const r_ji = neib.r - atom.r;

        double const r2ij = r_ji.squaredNorm();

        if (r2ij > R2cut) continue;

        r2[i]   = r2ij;

        density(r2[i], nc_);

        // charge on H atom is a sum over C neibs
        atom.nc_mod[neib.kind] += nc_;
    }

    return;
}

// dn/dr
// helper fucntion for density and dendity derivative
void EAMMOD1044::dRho(const double r2, double& dn) {
    double  X01, X02, X03, X04, X05, X06, X07, X08, X09;
    double  F0, F1, G0, G1;

    X01 = R2cut - r2;
    assert(X01 >= 0);
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

    G1 =     (D[12] + D[13]*X01)
       + X02*(D[14] + D[21]*X07)
       + X03*(D[15] + D[16]*X01)
       + X05*(D[17] + D[18]*X01)
       + X07*(D[19] + D[20]*X01);

    // r * dn/dr
    dn = (F1*G0 - F0*G1)/(G0*G0);

    return;
}

void EAMMOD1044::embGF(const double nr, double& F0, double& F1, double& G0, double& G1) {
    double  X01, X02, X03, X04, X05, X06, X07, X08, X09, X10, X11;

    X01 = nr;// - 1.0;

    X02 = X01*X01;  X03 = X02*X01;
    X04 = X02*X02;  X05 = X03*X02;
    X06 = X04*X02;  X07 = X04*X03;
    X08 = X04*X04;  X09 = X05*X04;
    X10 = X06*X04;  X11 = X06*X05;

    F0 = F[2]*X02 + F[3]*X03
       + F[4]*X04 + F[5]*X05
       + F[6]*X06 + F[7]*X07
       + F[8]*X08 + F[9]*X09
       + F[10]*X10 + F[11]*X11;

    F1 = X01*(F[12] + F[13]* X01)
       + X03*(F[14] + F[15]* X01)
       + X05*(F[16] + F[17]* X01)
       + X07*(F[18] + F[19]* X01)
       + X09*(F[20] + F[21]* X01);

    G0 = G[2]*X02 + G[3]*X03
       + G[4]*X04 + G[5]*X05
       + G[6]*X06 + G[7]*X07
       + G[8]*X08 + G[9]*X09
       + G[10]*X10 + G[11]*X11 + 1.0;

    G1 = X01*(G[12] + G[13]* X01)
       + X03*(G[14] + G[15]* X01)
       + X05*(G[16] + G[17]* X01)
       + X07*(G[18] + G[19]* X01)
       + X09*(G[20] + G[21]* X01);
}

void EAMMOD1044::embeddingEnergy(const double nr, double& EF) {
    double F0, F1, G0, G1;
    embGF(nr, F0, F1, G0, G1);
    EF = F0/G0;
    return;
}

// dG/dn
void EAMMOD1044::embeddingEnergyDeriv(const double nr, double& dF) {
    double F0, F1, G0, G1;
    embGF(nr, F0, F1, G0, G1);
    dF = (F1*G0 - F0*G1)/(G0*G0);
    return;
}

} // potential
} // solver
} // mdcraft