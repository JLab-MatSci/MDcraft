#include <gtest/gtest.h>
#include <iostream>
#include <cmath>
#include <mdcraft/solver/multi.h>
#include <mdcraft/solver/boundary/periodic.h>
#include <mdcraft/solver/potential/EAMWMOD1046.h>
#include <mdcraft/solver/potential/EAMWMOD0644.h>
#include <mdcraft/solver/potential/EAMMOD1044.h>
TEST(EamStageInvariants, PeriodicSplitAndOrderWithActualParallelWork) {
 using namespace mdcraft;
 tools::Threads threads(4);lattice::Domain domain(-1,1,-1,1,-1,1);domain.set_periodic(0,true);
 solver::boundary::Periodic boundary(domain);solver::Multi multi(boundary,solver::thermostat::dummy_thermostat,2,threads);
 double v[12]{},w[12]{},c[12]{},d[12]{},f[12]{},g[12]{};v[2]=.5;c[2]=3;f[2]=.2;g[2]=.01;
 solver::potential::EAMWMOD1046 cc(.6,0,0,0,12,v,w,c,d,f,g);
 solver::potential::EAMWMOD0644 hh(.6,0,0,0,12,v,w,c,d,f,g);
 solver::potential::EAMMOD1044 ch(.6,0,0,0,12,v,w,c,d,f,g);
 multi.add_potential(0,0,cc);multi.add_potential(0,1,ch);multi.add_potential(1,1,hh);
 // At least 2000 atoms are needed to create four tasks with the current
 // minimum of 500 atoms per task. Small fixtures do not exercise this race.
 data::Atoms atoms(2048);neibs::List list(2048);
 for(size_t i=0;i<atoms.size();++i){auto& a=atoms[i];a.r<<((i%2)? .94:-.94),.001*(i%31),.0007*(i/31);a.kind=(i/2)%2;a.uid=i+1;a.rns=.8;a.rcut=.6;a.m=1;a.v.setZero();for(int offset:{0,1,-1,2,-2,3,-3,511,-511,512})list[i].push_back((i+atoms.size()+offset)%atoms.size());}
 for(int step=0;step<5;++step){
  multi.prepare(atoms,atoms,list);multi.charge(atoms,atoms,list);
  auto charged=atoms;
  multi.forces(atoms,atoms,list);
  for(size_t i=0;i<atoms.size();++i)for(int k=0;k<2;++k)EXPECT_DOUBLE_EQ(atoms[i].nc_mod[k],charged[i].nc_mod[k]);
  multi.virials(atoms,atoms,list);
  for(size_t i=0;i<atoms.size();++i)for(int k=0;k<2;++k)EXPECT_DOUBLE_EQ(atoms[i].nc_mod[k],charged[i].nc_mod[k]);
 }

 auto reference=atoms;
 neibs::List first(atoms.size()),second(atoms.size());
 for(size_t i=0;i<atoms.size();++i)for(size_t k=0;k<list[i].size();++k)
   (k%2 ? first[i] : second[i]).push_back(list[i][k]);
 // Split density across owned and halo passes, then exchange FINAL density.
 multi.prepare(atoms,atoms,list);
 auto aliens=atoms;
 multi.charge(atoms,atoms,first);multi.charge(atoms,aliens,second);
 double density_error=0;
 for(size_t i=0;i<atoms.size();++i)for(int k=0;k<2;++k)
   density_error=std::max(density_error,std::abs(atoms[i].nc_mod[k]-reference[i].nc_mod[k]));
 aliens=atoms;
 multi.forces(atoms,atoms,first);multi.forces(atoms,aliens,second);
 double split_force_error=0,scale=0;
 for(size_t i=0;i<atoms.size();++i){split_force_error=std::max(split_force_error,(atoms[i].f-reference[i].f).norm());scale=std::max(scale,reference[i].f.norm());}
 for(auto& a:atoms){a.Ep=a.Em=0;a.V.setZero();}
 for(size_t i=0;i<atoms.size();++i){
   multi.virial_one(atoms.begin()+i,atoms,first[i]);boundary.virial_one(atoms.begin()+i,atoms,first[i]);
   multi.virial_one(atoms.begin()+i,aliens,second[i]);boundary.virial_one(atoms.begin()+i,aliens,second[i]);
 }
 double split_energy_error=0;
 for(size_t i=0;i<atoms.size();++i)split_energy_error=std::max(split_energy_error,std::abs(atoms[i].Ep+atoms[i].Em-reference[i].Ep-reference[i].Em));
 // Independently reverse the order of actual potential calls with frozen rho.
 solver::potential::Base* potentials[2][2]={{&cc,&ch},{&ch,&hh}};
 double order_error=0;
 for(size_t i=0;i<atoms.size();++i){
   auto neighbours=reference;
   for(auto& neighbour:neighbours)neighbour.r=atoms[i].r+domain.shortest(neighbour.r-atoms[i].r);
   atoms[i].f.setZero();
   for(int kind=1;kind>=0;--kind){
     neibs::ListOne typed;
     for(auto j:list[i])if(neighbours[j].kind==kind)typed.push_back(j);
     if(!typed.empty())potentials[atoms[i].kind][kind]->force(atoms.begin()+i,neighbours,typed);
   }
   order_error=std::max(order_error,(atoms[i].f-reference[i].f).norm());
 }
 std::cout<<"INVARIANTS density="<<density_error<<" split_force="<<split_force_error<<" split_energy="<<split_energy_error<<" reverse_order="<<order_error<<" force_scale="<<scale<<"\n";
 EXPECT_LT(density_error,1e-11);
 EXPECT_LT(split_force_error,1e-10);
 EXPECT_LT(split_energy_error,1e-10);
 EXPECT_LT(order_error,1e-10);
 std::cout<<"PASS: real 2048-atom Multi phases, 4 threads, periodic boundary, density immutable during forces/virials\n";
}
