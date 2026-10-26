#include <cmath>
#include <stdexcept>
#include <mdcraft/solver/multi.h>

namespace mdcraft {
namespace solver {

using ::mdcraft::tools::Threads;
using ::mdcraft::data::vector;
using ::mdcraft::data::matrix;

Multi::Multi(
	Boundary&   boundary,
	Thermostat& thermostat,
	std::size_t nkinds,
	Threads&    threads
) : Pair(boundary, thermostat, threads), nkinds(nkinds)
{
    if (nkinds == 0) {
        throw std::invalid_argument("Multi requires at least one atomic kind");
    }
    potentials.resize(nkinds, std::vector<Potential*>(nkinds, nullptr));
}

void Multi::add_potential(std::size_t i, std::size_t j, Potential& pot) {
    if (i >= nkinds || j >= nkinds) {
        throw std::out_of_range("Multi::add_potential: kind index out of range");
    }
	potentials[i][j] = &pot;
	potentials[j][i] = &pot;

	if (stages.empty()) {
		if (pot.type == potential::type::pair) {
			stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
				forces(atoms, neibs, nlist);
			});
		}
		else if (pot.type == potential::type::eam) {
			stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
				charge(atoms, neibs, nlist);
			});
			stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
				forces(atoms, neibs, nlist);
			});
		}
	} else if(pot.type == potential::type::eam and stages.size() == 1) {
		stages.emplace(stages.begin(), [this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
			charge(atoms, neibs, nlist);
		});
	}
}

void Multi::prepare_one(
	Atoms::iterator atomit, 
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
    // Also clear inactive atoms: diagnostics must not see stale values from
    // an earlier step after an atom has been disabled (tag == 0).
	atom.f  = vector::Zero();
	atom.V  = matrix::Zero();
	atom.Ep = 0.0;
	atom.Em = 0.0;
	atom.nc = 0.0;
    // Standard EAM accumulates into nc; modified C-H EAM uses separate
    // source-resolved densities. Both channels must be reset every step.
    for (auto& rho : atom.nc_mod) rho = 0.0;
    if (atom.tag == 0) return;
}

void Multi::charge_one(
	Atoms::iterator atomit,
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	auto const nlist_size = nlist.size();
	std::vector<NeibsListOne> multi_nlist(nkinds);

	for (std::size_t neibKind = 0; neibKind < nkinds; ++neibKind) {
		multi_nlist[neibKind].reserve(nlist_size);
	}

    const auto atomKind = static_cast<std::size_t>(atom.kind);
    if (atomKind >= nkinds) throw std::out_of_range("Multi: invalid atom kind");

    for (std::size_t i = 0; i < nlist_size; ++i) {
		auto& j = nlist[i];
		const auto& neib = neibs[j];
        const auto neibKind = static_cast<std::size_t>(neib.kind);
        if (neibKind >= nkinds) throw std::out_of_range("Multi: invalid neighbour kind");
        multi_nlist[neibKind].push_back(j);
	}

	for (std::size_t neibKind = 0; neibKind < nkinds; ++neibKind) {
        auto* pot = potentials[atomKind][neibKind];
        // Missing pair entries are intentionally non-interacting. Pair-only
        // potentials must not enter the electronic-density stage.
        if (pot != nullptr && pot->type == potential::type::eam &&
            !multi_nlist[neibKind].empty()) {
            pot->density(atomit, neibs, multi_nlist[neibKind]);
        }
	}


}

void Multi::force_one(
	Atoms::iterator atomit,
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	auto const nlist_size = nlist.size();
	std::vector<NeibsListOne> multi_nlist(nkinds);

	for (std::size_t neibKind = 0; neibKind < nkinds; ++neibKind) {
		multi_nlist[neibKind].reserve(nlist_size);
	}

    const auto atomKind = static_cast<std::size_t>(atom.kind);
    if (atomKind >= nkinds) throw std::out_of_range("Multi: invalid atom kind");

    for (std::size_t i = 0; i < nlist_size; ++i) {
		auto& j = nlist[i];
		auto& neib = neibs[j];
        const auto neibKind = static_cast<std::size_t>(neib.kind);
        if (neibKind >= nkinds) throw std::out_of_range("Multi: invalid neighbour kind");
        multi_nlist[neibKind].push_back(j);
	}

	for (std::size_t neibKind = 0; neibKind < nkinds; ++neibKind) {
        auto* pot = potentials[atomKind][neibKind];
        if (pot != nullptr && !multi_nlist[neibKind].empty()) {
            pot->force(atomit, neibs, multi_nlist[neibKind]);
        }
	}


}

void Multi::virial_one(
	Atoms::iterator atomit,
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	auto const nlist_size = nlist.size();

	std::vector<NeibsListOne> multi_nlist(nkinds);

	for (std::size_t neibKind = 0; neibKind < nkinds; ++neibKind) {
		multi_nlist[neibKind].reserve(nlist_size);
	}

    const auto atomKind = static_cast<std::size_t>(atom.kind);
    if (atomKind >= nkinds) throw std::out_of_range("Multi: invalid atom kind");

    for (std::size_t i = 0; i < nlist_size; ++i) {
		auto& j = nlist[i];
		auto& neib = neibs[j];
        const auto neibKind = static_cast<std::size_t>(neib.kind);
        if (neibKind >= nkinds) throw std::out_of_range("Multi: invalid neighbour kind");
        multi_nlist[neibKind].push_back(j);
	}

	for (std::size_t neibKind = 0; neibKind < nkinds; ++neibKind) {
        auto* pot = potentials[atomKind][neibKind];
        // EAM embedding energy F(rho) is a per-atom contribution. During the
        // local pass it must be evaluated even when the local list is empty
        // (all relevant neighbours may be in the MPI alien list).
        if (pot != nullptr &&
            (pot->type == potential::type::eam || !multi_nlist[neibKind].empty())) {
            pot->virial(atomit, neibs, multi_nlist[neibKind]);
        }
	}

}

} // solver
} // mdcraft