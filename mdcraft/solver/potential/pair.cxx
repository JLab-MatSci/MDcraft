#include <fstream>
#include <mdcraft/solver/potential/pair.h>

namespace mdcraft { 
namespace solver { 
namespace potential {

const double N_a       = 6.02214076e+23;
const double Qelectron = 1.602176634e-19;
const double UnMass    = 1e-3/N_a;
const double UnLeng    = 1e-9;
const double UnTime    = 1e-12;
const double UnEnrg    = UnMass*(UnLeng/UnTime)*(UnLeng/UnTime);

Pair::Pair(
	      double  Rcutoff,
	const double* _Vpar,
	      int     _Nr1,
	      double  _dr
) : Base         (Rcutoff, potential::type::pair),
    Vpar         (_Nr1), 
    dVpar        (_Nr1)
{
	Nr1  = _Nr1;
	dr   = _dr;
	
	dr1  = 1.0/dr;           
	
	for (std::size_t i = 0; i < Nr1; i++){
		Vpar[i] = _Vpar[i]; 
	}

	init_spline(Vpar.data(),dVpar.data(),dr1 ,Nr1, 1);
}

Pair::Pair(std::string filename) : Base(0.0, potential::type::pair) {
	std::ifstream file(filename);
	int Z; double M, Ecoh;
	std::string lattice;
	std::string line;
	// skip 4 lines
	std::getline(file, line);
	std::getline(file, line);
	std::getline(file, line);
	std::getline(file, line);
	file >> Nr1 >> dr >> R1cut;
	dr /= 10;    // from A to nm
	R1cut /= 10; // from A to nm
	R2cut = R1cut * R1cut;
	dr1  = 1.0/dr;           
	Vpar.resize(Nr1);  dVpar.resize(Nr1);
	for (std::size_t i = 0; i < Nr1; ++i) {
		file >> Vpar[i]; Vpar[i] *= Qelectron/(10*UnEnrg); // from eV*A to kJ*nm/mole
	}
	init_spline(Vpar.data(), dVpar.data(), dr1, Nr1, 1);
}

// potential U(r)
double Pair::value(const vector r) {
	double u, r2;
	r2 = r.squaredNorm();
	value(&r2, &u, 1ul);
	return u;
}
// virial U'(r)/r
matrix Pair::virial(const vector r) {
	double r2 = r.squaredNorm();
	matrix v;
	virial(&r, &r2, &v, 1ul);
	return v;
}
// force
vector Pair::force(const vector r) {
	double r2 = r.squaredNorm();
	vector f;
	force(&r, &r2, &f, 1ul);
	return f;
}

//
void Pair::value(const double* r2, double* Uout, std::size_t size) {
	double R1, RN, VR, DV, aa, bb;
	int N1;

	#pragma ivdep
	#pragma vector aligned
	for (int k = 0; k < size; k++){
		R1 = std::sqrt(r2[k]);
		RN = R1*dr1;
		N1 = std::trunc(RN);
	
		if (N1 < Nr1){
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

double Pair::value(const double* r2, std::size_t size) {
	double R1, RN, VR, DV, aa, bb;
	double Utot = 0.0;
	int N1;

	#pragma ivdep
	#pragma vector aligned
	for (int k = 0; k < size; k++){
		R1 = std::sqrt(r2[k]);
		RN = R1*dr1;
		N1 = std::trunc(RN);
	
		if (N1 < Nr1){
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

matrix Pair::virial(
	const vector* r, 
	const double* r2, 
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

		if(N1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;

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
		FF = ((DV* R1 - VR)/R2 )/R1;
		Vtot += 0.5 * FF * r[k] * r[k].transpose();
	}
	return Vtot;
}


void Pair::virial(
	const vector* r, 
	const double* r2, 
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

		if(N1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;

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
		FF = ((DV* R1 - VR)/R2 )/R1;
		Vout[k] = 0.5 * FF * r[k] * r[k].transpose();
	}
	return;
}

void Pair::virial(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
	if (nlist.half()) {
		auto& atom = *self;

		for (std::size_t i = 0; i < nlist.size(); ++i) {
			const auto j = nlist[i];
			auto& neib = neibs[j];

			const vector r_ji = neib.r - atom.r;
			const double r2ij = r_ji.squaredNorm();

			if (r2ij > R2cut || r2ij < 1e-15) continue;

			double R2 = r2ij;
			double R1 = std::sqrt(R2);
			double RN = R1 * dr1;
			int N1 = std::trunc(RN);

			double VR = 0.0;
			double DV = 0.0;
			double FF = 0.0;
			double aa = 0.0;
			double bb = 0.0;

			if (N1 < Nr1) {
				VR = N1;
				DV = RN - VR;
				RN = R1 - VR * dr;

				VR = (Vpar[N1 + 1] - Vpar[N1]) * dr1;
				aa = (dVpar[N1 + 1] - VR) - (VR - dVpar[N1]);
				bb = aa - (VR - dVpar[N1]);
				aa = aa * DV;

				VR = Vpar[N1] + RN * (dVpar[N1] + DV * (aa - bb));
				DV = dVpar[N1] + 2.0 * DV * (1.5 * aa - bb);
			}

			FF = ((DV * R1 - VR) / R2) / R1;
			const vector F = FF * r_ji;
			const matrix V = 0.5 * F * r_ji.transpose();
			const double E = value(&R2, 1);

			atom.V += V;
			atom.Ep += E;
			neib.V += V;
			neib.Ep += E;
		}

		return;
	}

	// use std::vector instead of c-array: almost no difference
	auto const nlist_size = nlist.size();
	auto& atom = *self;

	std::vector<vector> rall    (nlist_size);
	std::vector<double> r2      (nlist_size, 0.0);

	#pragma ivdep
	#pragma vector aligned
	for (auto i = 0; i < nlist_size; i++) {
		auto& j = nlist[i];
		auto neib = neibs[j];

		auto const r_ji = neib.r - atom.r;
		auto const r2ij = r_ji.squaredNorm();

		r2[i]      = r2ij;
		rall[i]    = r_ji;
	}

	std::size_t natoms = 0;
	#pragma ivdep
	#pragma vector aligned
	for (auto i = 0; i < nlist_size; i++){
		if (r2[i] > R2cut) continue;
		if (r2[i] < 1e-15) continue;

		r2[natoms] = r2[i];
		rall[natoms] = rall[i];

		++natoms;
	}

	atom.V  += virial(rall.data(), r2.data(), natoms);
	atom.Ep += value(r2.data(), natoms);

	return;
}

vector Pair::force(
	const vector* r, 
	const double* r2, 
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

		if(N1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;

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
		FF = ((DV* R1 - VR)/R2 )/R1;
		Ftot += FF * r[k];
	}

	return Ftot;
}

void Pair::force(
	const vector* r, 
	const double* r2,
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

		if(N1 < Nr1){
			VR = N1;
			DV = RN - VR;
			RN = R1 - VR*dr;

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
		FF = ( (DV* R1 - VR)/R2 )/R1;
		Fout[k] += FF * r[k];
	}

	return;
}

void Pair::force(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist) {
	if (nlist.half()) {
		auto& atom = *self;

		for (std::size_t i = 0; i < nlist.size(); ++i) {
			const auto j = nlist[i];
			auto& neib = neibs[j];

			const vector r_ji = neib.r - atom.r;
			const double r2ij = r_ji.squaredNorm();

			if (r2ij > R2cut || r2ij < 1e-15) continue;

			double R2 = r2ij;
			double R1 = std::sqrt(R2);
			double RN = R1 * dr1;
			int N1 = std::trunc(RN);

			double VR = 0.0;
			double DV = 0.0;
			double FF = 0.0;
			double aa = 0.0;
			double bb = 0.0;

			if (N1 < Nr1) {
				VR = N1;
				DV = RN - VR;
				RN = R1 - VR * dr;

				VR = (Vpar[N1 + 1] - Vpar[N1]) * dr1;
				aa = (dVpar[N1 + 1] - VR) - (VR - dVpar[N1]);
				bb = aa - (VR - dVpar[N1]);
				aa = aa * DV;

				VR = Vpar[N1] + RN * (dVpar[N1] + DV * (aa - bb));
				DV = dVpar[N1] + 2.0 * DV * (1.5 * aa - bb);
			}

			FF = ((DV * R1 - VR) / R2) / R1;
			const vector F = FF * r_ji;

			atom.f += F;
			neib.f -= F;
		}

		return;
	}

	// use std::vector instead of c-array: almost no difference
	auto const nlist_size = nlist.size();
	auto& atom = *self;

	std::vector<vector> rall    (nlist_size);
	std::vector<double> r2      (nlist_size    , 0.0);

	#pragma ivdep
	#pragma vector aligned
	for (auto i = 0; i < nlist_size; i++) {
		auto& j = nlist[i];
		auto neib = neibs[j];

		vector const r_ji = neib.r - atom.r;
		double const r2ij = r_ji.squaredNorm();

		r2[i]      = r2ij;
		rall[i]    = r_ji;
	}

	std::size_t natoms = 0;
	#pragma ivdep
	#pragma vector aligned
	for (auto i = 0; i < nlist_size; i++){
		if (r2[i] > R2cut) continue;
		if (r2[i] < 1e-15) continue;

		r2[natoms] = r2[i];
		rall[natoms] = rall[i];

		++natoms;
	}

	atom.f += force(rall.data(), r2.data(), natoms);
}

void Pair::init_spline(const double* y, double* z, double H1, std::size_t N, std::size_t ib){
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
	 z[i]= s[i]*z[i] + t[i];
	}

}

} // potential
} // solver
} // mdcraft
