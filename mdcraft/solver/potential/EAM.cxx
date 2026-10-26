#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <mdcraft/solver/potential/EAM.h>

namespace mdcraft { 
namespace solver { 
namespace potential {

const double N_a       = 6.02214076e+23;
const double Qelectron = 1.602176634e-19;
const double UnMass    = 1e-3/N_a;
const double UnLeng    = 1e-9;
const double UnTime    = 1e-12;
const double UnEnrg    = UnMass*(UnLeng/UnTime)*(UnLeng/UnTime);

EAM::EAM(
	      double  Rcutoff,
	const double* _Femb,
	const double* _Rho,
	const double* _Vpar,
	      int     _Nr1,
	      double  _dr,
	      int     _Nrh1,
	      double  _drh
) : Base         (Rcutoff, potential::type::eam),
    Femb         (std::max(0, _Nrh1)),
    dFemb        (std::max(0, _Nrh1)),
    Rho          (std::max(0, _Nr1)),
    dRho         (std::max(0, _Nr1)),
    Vpar         (std::max(0, _Nr1)),
    dVpar        (std::max(0, _Nr1))
{
    if (_Nr1 < 3 || !(_dr > 0) || !(Rcutoff > 0) || !_Vpar)
        throw std::invalid_argument("Invalid potential table dimensions or grid");
    if (_Nrh1 < 3 || !(_drh > 0) || !_Femb || !_Rho)
        throw std::invalid_argument("Invalid embedding table dimensions or grid");
	Nr1  = _Nr1;
	Nrh1 = _Nrh1;
	R2cut = Rcutoff * Rcutoff;
	dr   = _dr;
	drh  = _drh;

	dr1  = 1.0/dr;           
	drh1 = 1.0/drh;

	for (std::size_t i = 0; i < Nr1; i++){
		Rho[i]  = _Rho[i]; 
		Vpar[i] = _Vpar[i]; 
	}

	for (std::size_t i = 0; i < Nrh1; i++){
		Femb[i] = _Femb[i];
	}

	init_spline(Vpar.data(),dVpar.data(),dr1 ,Nr1, 1);
	init_spline(Rho.data(),dRho.data(),dr1 ,Nr1, 1);
	init_spline(Femb.data(),dFemb.data(),drh1 ,Nrh1, 2);

}

EAM::EAM(std::string filename) : Base(0.0, potential::type::eam) {
	std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open potential file: " + filename);
	int Z; double M, Ecoh;
	std::string lattice;
	std::string line;
	// skip 4 lines
	std::getline(file, line);
	std::getline(file, line);
	std::getline(file, line);
	std::getline(file, line);
	file >> Nrh1 >> drh >> Nr1 >> dr >> R1cut;
	file >> Z >> M >> Ecoh >> lattice;
    if (!file || Nr1 < 3 || !(dr > 0) || !(R1cut > 0))
        throw std::runtime_error("Invalid potential table header: " + filename);
    if (Nrh1 < 3 || !(drh > 0))
        throw std::runtime_error("Invalid embedding table header: " + filename);
	dr /= 10;    // from A to nm
	R1cut /= 10; // from A to nm
	R2cut = R1cut * R1cut;
	dr1  = 1.0/dr;           
	drh1 = 1.0/drh;
	Femb.resize(Nrh1); dFemb.resize(Nrh1);
	Rho.resize(Nr1);   dRho.resize(Nr1);
	Vpar.resize(Nr1);  dVpar.resize(Nr1);
	for (std::size_t i = 0; i < Nrh1; ++i) {
		file >> Femb[i]; Femb[i] *= Qelectron/UnEnrg; // from eV to kJ/mole
	}
	for (std::size_t i = 0; i < Nr1; ++i) {
		file >> Rho[i];
	}
	for (std::size_t i = 0; i < Nr1; ++i) {
		file >> Vpar[i]; Vpar[i] *= Qelectron/(10*UnEnrg); // from eV*A to kJ*nm/mole
	}
    if (!file) throw std::runtime_error("Truncated potential table: " + filename);
	init_spline(Vpar.data(), dVpar.data(), dr1,   Nr1, 1);
	init_spline( Rho.data(),  dRho.data(), dr1,   Nr1, 1);
	init_spline(Femb.data(), dFemb.data(), drh1, Nrh1, 2);
}

// potential U(r)
double EAM::value(const vector r) {
	double u, r2;
	r2 = r.squaredNorm();
	value(&r2, &u, 1ul);
	return u;
}
// virial U'(r)/r
matrix EAM::virial(const vector r, const double AdF) {
	double r2 = r.squaredNorm();
	matrix v;
	virial(&r, &r2, &AdF, &v, 1ul);
	return v;
}
// force
vector EAM::force(const vector r, const double AdF) {
	double r2 = r.squaredNorm();
	vector f;
	force(&r, &r2, &AdF, &f, 1ul);
	return f;
}

// temp functions for debug
void EAM::U_arr(double* Uout, int size) {
	for (int k = 0; k < size; k++){
		Uout[k] = Vpar[k];
	}
	return;
}

void EAM::dU_arr(double* Uout, int size) {
	for (int k = 0; k < size; k++){
		Uout[k] = dVpar[k];
	}
	return;
}

void EAM::dFemb_arr(double* Uout, int size) {
	for (int k = 0; k < size; k++){
		Uout[k] = dFemb[k];
	}
	return;
}

void EAM::dRho_arr(double* Uout, int size) {
	for (int k = 0; k < size; k++){
		Uout[k] = dRho[k];
	}
	return;
}

//
void EAM::value(const double* r2, double* Uout, std::size_t size) {
	double R1, RN, VR, DV, aa, bb;
	int N1;

	#pragma ivdep
	#pragma vector aligned
	for (int k = 0; k < size; k++){
		R1 = std::sqrt(r2[k]);
		RN = R1*dr1;
		N1 = std::trunc(RN);
	
		if (N1 >= 0 && N1 + 1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;

			VR = (Vpar[N1+1] - Vpar[N1])*dr1;
			aa = (dVpar[N1+1] - VR) - (VR - dVpar[N1]);
			bb = 				aa - (VR - dVpar[N1]);
			aa = aa * DV;

			VR = Vpar[N1] + RN*(dVpar[N1] + DV*(aa - bb));
			// DV = dVpar(N1) + 2.0* DV * (1.5*aa - bb)''
		}
		else{
			VR = 0.0;
			// DV = 0.0;
		}

		Uout[k] = 0.5*VR/R1; // V(r)/2
	}

	return;
}

double EAM::value(const double* r2, std::size_t size) {
	double R1, RN, VR, DV, aa, bb, Utot = 0.0;
	int N1;

	#pragma ivdep
	#pragma vector aligned
	for (int k = 0; k < size; k++){
		R1 = std::sqrt(r2[k]);
		RN = R1*dr1;
		N1 = std::trunc(RN);
	
		if (N1 >= 0 && N1 + 1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;

			VR = (Vpar[N1+1] - Vpar[N1])*dr1;
			aa = (dVpar[N1+1] - VR) - (VR - dVpar[N1]);
			bb = 				aa - (VR - dVpar[N1]);
			aa = aa * DV;

			VR = Vpar[N1] + RN*(dVpar[N1] + DV*(aa - bb));
			// DV = dVpar(N1) + 2.0* DV * (1.5*aa - bb)''
		}
		else{
			VR = 0.0;
			// DV = 0.0;
		}

		Utot += 0.5*VR/R1; // V(r)/2
	}

	return Utot;
}

matrix EAM::virial(
	const vector* r, 
	const double* r2, 
	const double* AdF, 
	std::size_t   size
) {
	double R1, R2, RN, VR, DV, FF, aa, bb;
	int N1;
	matrix Vtot = matrix::Zero();
	#pragma ivdep
	#pragma vector aligned
	for (std::size_t k = 0; k < size; k++){
		R2 = r2[k];
		R1 = std::sqrt(R2);
		RN = R1 * dr1;
		N1 = std::trunc(RN);

		if(N1 >= 0 && N1 + 1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;
			VR = ( Rho[N1+1] - Rho[N1])*dr1;
			aa = (dRho[N1+1] - VR) - (VR - dRho[N1]);
			bb =               aa - (VR - dRho[N1]);

			FF =  dRho[N1] + 2.0* DV*(1.5* DV*aa - bb);  // n'(r)           

			VR = ( Vpar[N1+1] - Vpar[N1])*dr1;
			aa = (dVpar[N1+1] - VR) - (VR - dVpar[N1]);
			bb =                aa  - (VR - dVpar[N1]);
			aa = aa* DV;

			VR =  Vpar[N1] + RN*(dVpar[N1] + DV*(aa - bb));
			DV = dVpar[N1] + 2.0* DV*(1.5*aa - bb);
		}
		else{
			FF = 0.0;
			VR = 0.0;
			DV = 0.0;
		}
		// ! Force(r)/r = [(D(i)+D(j)*n'(r) + V'(r)]/r
		FF = ( AdF[k]* FF + (DV* R1 - VR)/R2 )/R1;
		Vtot += 0.5 * FF * r[k] * r[k].transpose();
	}
	return Vtot;
}


void EAM::virial(
	const vector* r, 
	const double* r2, 
	const double* AdF, 
	matrix*       Vout,
	std::size_t   size
) {
	double R1, R2, RN, VR, DV, FF, aa, bb;
	int N1;
	#pragma ivdep
	#pragma vector aligned
	for (std::size_t k = 0; k < size; k++){
		R2 = r2[k];
		R1 = std::sqrt(R2);
		RN = R1 * dr1;
		N1 = std::trunc(RN);

		if(N1 >= 0 && N1 + 1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;
			VR = ( Rho[N1+1] - Rho[N1])*dr1;
			aa = (dRho[N1+1] - VR) - (VR - dRho[N1]);
			bb =               aa - (VR - dRho[N1]);

			FF =  dRho[N1] + 2.0* DV*(1.5* DV*aa - bb);  // n'(r)           

			VR = ( Vpar[N1+1] - Vpar[N1])*dr1;
			aa = (dVpar[N1+1] - VR) - (VR - dVpar[N1]);
			bb =                aa  - (VR - dVpar[N1]);
			aa = aa* DV;

			VR =  Vpar[N1] + RN*(dVpar[N1] + DV*(aa - bb));
			DV = dVpar[N1] + 2.0* DV*(1.5*aa - bb);
		}
		else{
			FF = 0.0;
			VR = 0.0;
			DV = 0.0;
		}
		// ! Force(r)/r = [(D(i)+D(j)*n'(r) + V'(r)]/r
		FF = ( AdF[k]* FF + (DV* R1 - VR)/R2 )/R1;
		Vout[k] = 0.5 * FF * r[k] * r[k].transpose();
	}
	return;
}

void EAM::virial(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
	// use std::vector instead of c-array: almost no difference
	auto const nlist_size = nlist.size(); 

	auto& atom = *self;

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

		auto const r_ji = neib.r - atom.r;
		auto const r2ij = r_ji.squaredNorm();

		r2[i]      = r2ij;
		rall[i]    = r_ji;
		nCharge[i] = neib.nc;
	}
	nCharge[nlist_size] = atom.nc;

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

	atom.V  += virial(rall.data(), r2.data(), dFemb.data(), natoms);
	atom.Ep += value(r2.data(), natoms);

	return;
}

vector EAM::force(
	const vector* r, 
	const double* r2, 
	const double* AdF, 
	std::size_t   size
) {
	double R1, R2, RN, VR, DV, FF, aa, bb;
	int N1;
	vector Ftot = vector::Zero();
	#pragma ivdep
	#pragma vector aligned
	for (std::size_t k = 0; k < size; k++){
		R2 = r2[k];
		R1 = std::sqrt(R2);
		RN = R1 * dr1;
		N1 = std::trunc(RN);

		if(N1 >= 0 && N1 + 1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;
			VR = ( Rho[N1+1] - Rho[N1])*dr1;
			aa = (dRho[N1+1] - VR) - (VR - dRho[N1]);
			bb =               aa - (VR - dRho[N1]);

			FF =  dRho[N1] + 2.0* DV*(1.5* DV*aa - bb);  // n'(r)           

			VR = ( Vpar[N1+1] - Vpar[N1])*dr1;
			aa = (dVpar[N1+1] - VR) - (VR - dVpar[N1]);
			bb =                aa  - (VR - dVpar[N1]);
			aa = aa* DV;

			VR =  Vpar[N1] + RN*(dVpar[N1] + DV*(aa - bb));
			DV = dVpar[N1] + 2.0* DV*(1.5*aa - bb);
		}
		else{
			FF = 0.0;
			VR = 0.0;
			DV = 0.0;
		}
		// ! Force(r)/r = [(D(i)+D(j)*n'(r) + V'(r)]/r
		FF = ( AdF[k]* FF + (DV* R1 - VR)/R2 )/R1;
		Ftot += FF * r[k];
	}

	return Ftot;
}

void EAM::force(
	const vector* r, 
	const double* r2,
	const double* AdF, 
	vector*       Fout, 
	std::size_t   size
) {
	double R1, R2, RN, VR, DV, FF, aa, bb;
	int N1;

	#pragma ivdep
	#pragma vector aligned
	for (std::size_t k = 0; k < size; k++){
		R2 = r2[k];
		R1 = std::sqrt(R2);
		RN = R1 * dr1;
		N1 = std::trunc(RN);

		if(N1 >= 0 && N1 + 1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;
			VR = ( Rho[N1+1] - Rho[N1])*dr1;
			aa = (dRho[N1+1] - VR) - (VR - dRho[N1]);
			bb =               aa - (VR - dRho[N1]);

			FF =  dRho[N1] + 2.0* DV*(1.5* DV*aa - bb);  //! n'(r)           

			VR = ( Vpar[N1+1] - Vpar[N1])*dr1;
			aa = (dVpar[N1+1] - VR) - (VR - dVpar[N1]);
			bb =                aa  - (VR - dVpar[N1]);
			aa = aa* DV;

			VR =  Vpar[N1] + RN*(dVpar[N1] + DV*(aa - bb));
			DV = dVpar[N1] + 2.0* DV*(1.5*aa - bb);
		}
		else{
			FF = 0.0;
			VR = 0.0;
			DV = 0.0;
		}
		// ! Force(r)/r = [(D(i)+D(j)*n'(r) + V'(r)]/r
		FF = ( AdF[k]* FF + (DV* R1 - VR)/R2 )/R1;
		Fout[k] += FF * r[k];
	}

	return;
}

void EAM::force(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
	// use std::vector instead of c-array: almost no difference
	auto const nlist_size = nlist.size(); 
	if (nlist_size == 0) return;
	auto& atom = *self;

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

		vector const r_ji = neib.r - atom.r;
		double const r2ij = r_ji.squaredNorm();

		r2[i]      = r2ij;
		rall[i]    = r_ji;
		nCharge[i] = neib.nc;
	}
	nCharge[nlist_size] = atom.nc;

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

	atom.f += force(rall.data(), r2.data(), dFemb.data(), natoms);
}

double EAM::density(const double* r2, std::size_t size){
	double SD = 0.0; // sum of densities
	double nr;//Temp density
	double R1, RN, VR, DV, aa, bb;
	int N1;
	#pragma ivdep
	#pragma vector aligned
	for (std::size_t k = 0; k < size; k ++){
		R1 = std::sqrt(r2[k]);
		RN = R1 * dr1;
		N1 = std::trunc(RN);
		
		if (N1 >= 0 && N1 + 1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;  

			VR = ( Rho[N1+1] - Rho[N1])*dr1;
			aa = (dRho[N1+1] - VR) - (VR - dRho[N1]);
			bb =                aa - (VR - dRho[N1]);

			nr = Rho[N1] + RN*(dRho[N1] + DV*(DV*aa - bb));

			SD = SD + nr;
		}
		else{
			nr = 0.0;
		}

	}
	return SD;
}

void EAM::density(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
	auto const nlist_size = nlist.size(); 
	if (nlist_size == 0) return;
	std::vector<double> r2(nlist_size);
	auto& atom = *self;

	#pragma ivdep
	#pragma vector aligned
	for (auto i = 0; i < nlist.size(); i++) {
		auto& j = nlist[i];
		auto neib = neibs[j];

		vector const r_ji = neib.r - atom.r;
		double const r2ij = r_ji.squaredNorm();

		r2[i]   = r2ij;
	}

	std::size_t natoms = 0;
	for (auto i = 0; i < nlist.size(); i++){
		if (r2[i] > R2cut) continue;
		if (r2[i] < 1e-15) continue;

		r2[natoms] = r2[i];

		++natoms;
	}
	atom.nc += density(r2.data(), natoms);
}

void EAM::embeddingEnergy(
	const double* nr, 
	double*       EF, 
	double*       dF, 
	std::size_t   size
) {
	double RN, VR, DV, aa, bb;
	int N1;
	#pragma ivdep
	#pragma vector aligned
	for (std::size_t k = 0; k < size; k++){
        const double rho = nr[k];
        if (!std::isfinite(rho)) throw std::domain_error("Nonfinite EAM density");
        if (rho < 0.0) {
            EF[k] = Femb.front() + dFemb.front() * rho;
            dF[k] = dFemb.front();
            continue;
        }
        if (rho >= (Nrh1 - 1) * drh) {
            EF[k] = Femb.back() + dFemb.back() * (rho - (Nrh1 - 1) * drh);
            dF[k] = dFemb.back();
            continue;
        }
        RN = rho * drh1;
        N1 = std::min(static_cast<int>(RN), Nrh1 - 2);
        DV = RN - N1;
        const double delta = rho - N1 * drh;
        VR = (Femb[N1+1] - Femb[N1]) * drh1;
        aa = (dFemb[N1+1] - VR) - (VR - dFemb[N1]);
        bb = aa - (VR - dFemb[N1]);
        aa *= DV;
        EF[k] = Femb[N1] + delta * (dFemb[N1] + DV * (aa - bb));
        dF[k] = dFemb[N1] + 2.0 * DV * (1.5 * aa - bb);
	}
}

void EAM::init_spline(const double* y, double* z, double H1, std::size_t N, std::size_t ib){
// !
// ! the original code was taken from the book by Shikin E.V., Plis A.I. 
// !           Handbook on splines for the user (CRC,1995)(ISBN 084939404X)
// !
// ! Program constructs the interpolating cubic spline for the function
// ! y(0:N) given in a tabulated form in the domain [0,N*h] 
// !
// ! y - array containing the function values at equidistant knots
// !
// ! a,b,c,d   - work arrays (of size N)
// ! s,t,u,v,w - work arrays (of size N+1)
// !
// !  ib - type of boundary conditions
// !  ib = 1 - first derivatives are given at the end points of the segment
// !  ib = 2 - second derivatives are given at the end points
// !            (ax=0, bx=0 are natural end conditions)
// ! Output:
// ! z - array (of size N+1) containing the spline parameters
// !-----------------------------------------------------------------------
	double rp;
	std::vector<double> a(N),b(N),c(N),d(N);
	std::vector<double> s(N+1),t(N+1),u(N+1),v(N+1),w(N+1);

	std::fill(a.begin(),a.end(),2.0);
	b[0] = 1.0;
	c[0] = 0.0;
	d[0] = 3.0*H1*(y[1] - y[0]);

	b[N-1] = 0.0;

	if (ib == 1){
		c[N-1] = 0.0; //Y'(N) = 0
		d[N-1] = 0.0;
	}
	else if (ib == 2){
		c[N-1] = 1.0;
		d[N-1] = 3.0*H1*(y[N-1] - y[N-2]); // Y''(N) = 0
	}

	for (std::size_t i = 1; i < N - 1; i++){
		c[i] = 0.5;
		b[i] = 0.5;
		d[i] =  1.5*H1*(y[i+1] - y[i-1]);
	}

// !-----------------------------------------------------------------------
// !  Solving linear system of the form
// !
// !        a(1)*x(1) + b(1)*x(2)      + c(1)*x(n) = d(1)
// !        ............................................
// !        c(i)*x(i-1) + a(i)*x(i)  + b(i)*x(i+1) = d(i)
// !        .............................................
// !        b(n)*x(1)     + a(n)*x(n-1) + b(n)*x(n)= d(n)
// !
// !  with tridiagonal matrix by sweep method.
// !-----------------------------------------------------------------------

	u[0] = 0.0;
	v[0] = 0.0;
	w[0] = 1.0;

	for (std::size_t i = 0; i < N ; i++){
		rp = -1.0/(a[i] + c[i]*v[i]);
		v[i+1] =  rp * b[i];
		u[i+1] =  rp *(c[i]*u[i] - d[i]);
		w[i+1] =  rp * c[i]*w[i];
	}

	s[N-1] = 1.0;// N-1 or N
	t[N-1] = 0.0;// N-1 or N

 	for (int i = N - 2; i >= 0 ; i--){//N-2 or N-1
		s[i]= v[i+1]*s[i+1] + w[i+1];
		t[i]= v[i+1]*t[i+1] + u[i+1];
	}
	z[N-1] = (d[N-1] - b[N-1]*t[0] - c[N-1]*t[N-2])/
	         (a[N-1] + b[N-1]*s[0] + c[N-1]*s[N-2]);

	for (std::size_t i = 0; i < N-1 ; i++){
	 z[i]= s[i]*z[N-1] + t[i];
	}

}

} // potential
} // solver
} // mdcraft