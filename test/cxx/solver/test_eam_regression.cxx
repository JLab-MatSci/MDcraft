#include <gtest/gtest.h>
#include <atomic>
#include <mdcraft/solver/multi.h>
#include <mdcraft/solver/boundary/periodic.h>
#include <mdcraft/solver/potential/EAM.h>
#include <test/cxx/support/md_test_utils.h>
#ifdef mdcraft_ENABLE_MPI
#include <mdcraft/decomp/VD3.h>
#endif
using namespace mdcraft;
namespace {
struct DensityProbe : solver::potential::Base {
    DensityProbe():Base(0.5, solver::potential::type::eam){}
    void density(data::Atoms::iterator a, data::Atoms& n, neibs::ListOne& list) override {
        for(auto j:list) if((n[j].r-a->r).squaredNorm()<0.25) a->nc+=1.0;
    }
};
struct CountThermostat : solver::thermostat::Base {
    std::atomic<int> calls{0};
    void apply_one(data::Atoms::iterator) override {++calls;}
};
}
TEST(EamRegression, PeriodicDensityIsNotDoubled) {
    data::Atoms atoms(2);
    for(int i=0;i<2;++i){atoms[i].r<< (i?0.7:-0.7),0,0;atoms[i].kind=0;atoms[i].rns=.6;}
    lattice::Domain domain(-.81,.81,-.81,.81,-.81,.81);domain.set_periodic(0,true);
    solver::boundary::Periodic boundary(domain);DensityProbe potential;
    solver::Multi multi(boundary,solver::thermostat::dummy_thermostat,1);
    multi.add_potential(0,0,potential);neibs::List list{{1},{0}};
    multi.prepare(atoms,atoms,list);multi.charge(atoms,atoms,list);
    EXPECT_DOUBLE_EQ(atoms[0].nc,1);EXPECT_DOUBLE_EQ(atoms[1].nc,1);
}
TEST(EamRegression, TabulatedEndpointsAreFinite) {
    solver::potential::EAM potential((test::source_root()/"problems/NiAl/Al.eam.alloy").string());
    double rho[]={-1,0,7,100}, energy[4], deriv[4];
    potential.embeddingEnergy(rho,energy,deriv,4);
    for(int i=0;i<4;++i){EXPECT_TRUE(std::isfinite(energy[i]));EXPECT_TRUE(std::isfinite(deriv[i]));}
    double r=potential.rcut()*(1.0-.5/10000), r2=r*r;
    EXPECT_TRUE(std::isfinite(potential.density(&r2,1)));
}
TEST(EamRegression, InvalidTableFileThrows) {
    EXPECT_THROW(solver::potential::EAM("/nonexistent/eam-file"),std::runtime_error);
}
TEST(EamRegression, InvalidKindExceptionReachesCallerWithThreads) {
    tools::Threads threads(2);DensityProbe potential;solver::Multi multi(
      solver::boundary::dummy_boundary,solver::thermostat::dummy_thermostat,1,threads);
    multi.add_potential(0,0,potential);data::Atoms atoms(8);neibs::List list(8);
    for(auto& atom:atoms) atom.kind=3;
    EXPECT_THROW(multi.charge(atoms,atoms,list),std::out_of_range);
}
#ifdef mdcraft_ENABLE_MPI
TEST(EamRegression, DistributedStepperAppliesThermostatOnce) {
    tools::Network comm;tools::Threads threads(2);
    data::Atoms atoms;
    if(comm.master()) for(int i=0;i<4;++i){data::Atom a{};a.uid=i+1;a.m=1;a.rcut=.2;a.rns=.3;a.r<< -.6+.4*i,0,0;atoms.push_back(a);}
    lattice::Domain domain(-2,2,-2,2,-2,2);
    decomp::VD3 vd(comm,atoms,domain,threads,{});vd.exchange();
    neibs::VerletList local(atoms,atoms,domain,threads),alien(atoms,vd.aliens(),domain,threads);
    local.update(-1);alien.update(-1);
    CountThermostat thermostat;solver::potential::Base potential(.2,solver::potential::type::pair);
    solver::Multi multi(solver::boundary::dummy_boundary,thermostat,1,threads);multi.add_potential(0,0,potential);
    solver::stepper::Verlet stepper(multi,threads);
    for(int i=0;i<3;++i)stepper.make_step(vd,local.get(),alien.get(),.0001);
    EXPECT_EQ(comm.sum(thermostat.calls.load()),12);
}
#endif

TEST(EamRegression, ThermostatSeesBoundaryForce) {
    struct BoundaryForce : solver::boundary::Base {
        void prepare_one(data::Atoms::iterator a, data::Atoms&, neibs::ListOne&) override { a->bc=1; }
        void force_one(data::Atoms::iterator a, data::Atoms&, neibs::ListOne&) override { a->f.x()+=2; }
    } boundary;
    struct ForceObserver : solver::thermostat::Base {
        double observed=0;
        void apply_one(data::Atoms::iterator a) override { observed=a->f.x(); }
    } thermostat;
    solver::potential::Base potential(.2,solver::potential::type::pair);
    solver::Multi multi(boundary,thermostat,1);multi.add_potential(0,0,potential);
    data::Atoms atoms(1);atoms[0].kind=0;neibs::List list(1);
    multi.prepare(atoms,atoms,list);multi.forces(atoms,atoms,list);
    EXPECT_DOUBLE_EQ(thermostat.observed,2.0);
}
