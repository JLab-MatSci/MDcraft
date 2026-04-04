#include <atomic>
#include <cstring>
#include <stdexcept>
#include <mdcraft/neibs/verlet-list.h>

namespace mdcraft::neibs {

Atoms dummy_atoms(1);
Atoms dummy_neibs(1);

VerletList::VerletList(
	Atoms&   atoms,
	Atoms&   neibs,
	Domain&  domain,
	Threads& pool,
	bool     half
) : grid(domain, pool),
	domain(domain),
	pool(pool),
	half(half)
{
	fetch_data(atoms, neibs);

	if (m_atoms != m_neibs)
		this->half = false;

	if (this->half && pool.active()) {
		throw std::runtime_error(
			"Half Verlet lists are currently incompatible with multithreaded pair-force evaluation. "
			"Disable threads or use a full neighbor list."
		);
	}

	nlist.resize(m_atoms->size());
	for (auto& l: nlist) {
		l.reserve(128);
		l.set_half(this->half);
	}
}

void VerletList::update(double kbuf) {
	if (m_neibs->empty()) {
		clear_lists();
		snapshot_atoms();
		return;
	}
	if (kbuf > 0.0) optimize_search_radius(kbuf);
	clear_lists();
	// put_atoms_in_cells(*m_neibs);
	// find_cell_starts();
	// count_atoms_in_cells();
	// build();
	build_linked_cells();
	snapshot_atoms();
}

void VerletList::fetch_data(
	Atoms& atoms,
	Atoms& neibs
) {
	if (atoms.data() != dummy_atoms.data()) {
		m_atoms = &atoms;
		if (neibs.data() != dummy_neibs.data()) {
			m_neibs = &neibs;
		}
		else {
			m_neibs = &atoms;
		}
	}
}

void VerletList::optimize_search_radius(double kbuf) {
	auto optimize_one = [&](Atoms::iterator p) {
		auto& atom = *p;
		auto iatom = p - m_atoms->begin();
		auto& neibs = *m_neibs;
		double hmax = 0.0;
		for (int i = 0; i < nlist[iatom].size(); ++i) {
			auto& neighbor = neibs[nlist[iatom][i]];

			auto const r_ji = neighbor.r - atom.r;

			auto const r = std::sqrt(r_ji.dot(r_ji));
			auto const h = atom.rcut + neighbor.rcut;

			hmax = std::max(h, hmax);
		}
		atom.rns = kbuf * hmax;
	};

	pool.for_each(m_atoms->begin(), m_atoms->end(), optimize_one);
}

void VerletList::clear_lists() {
	nlist.resize(m_atoms->size());
	for (auto& l: nlist) {
		l.clear();
		l.set_half(half);
	}
}

ListOne& VerletList::operator[](point_id i) {
	return i < nlist.size() ? nlist[i] : dummy_nlist; 
}

List& VerletList::get() { 
	return nlist; 
}

bool VerletList::decide() const {
	if (m_atoms == nullptr) {
		return true;
	}

	const auto& atoms = *m_atoms;

	if (atoms.size() != saved_atoms.size()) {
		return true;
	}

	if (atoms.empty()) {
		return false;
	}

	std::atomic<bool> need_rebuild = false;

	pool.parallel_for(std::size_t{0}, atoms.size(), [&](std::size_t i) {
		if (need_rebuild.load(std::memory_order_relaxed)) {
			return;
		}

		const auto& atom = atoms[i];
		const auto& saved = saved_atoms[i];

		if (atom.uid != saved.uid ||
			atom.rcut != saved.rcut ||
			atom.rns != saved.rns
		) {
			need_rebuild.store(true, std::memory_order_relaxed);
			return;
		}

		const double skin = saved.rns - saved.rcut;
		if (skin <= 0.0) {
			need_rebuild.store(true, std::memory_order_relaxed);
			return;
		}

		const auto dr = domain.shortest(atom.r - saved.r);
		const double trigger2 = 0.25 * skin * skin;

		if (dr.squaredNorm() > trigger2) {
			need_rebuild.store(true, std::memory_order_relaxed);
		}
	});

	return need_rebuild.load(std::memory_order_relaxed);
}

void VerletList::snapshot_atoms() {
	if (m_atoms == nullptr) {
		saved_atoms.clear();
		return;
	}

	const auto& atoms = *m_atoms;
	saved_atoms.resize(atoms.size());

	pool.parallel_for(std::size_t{0}, atoms.size(), [&](std::size_t i) {
		saved_atoms[i].r = atoms[i].r;
		saved_atoms[i].rcut = atoms[i].rcut;
		saved_atoms[i].rns = atoms[i].rns;
		saved_atoms[i].uid = atoms[i].uid;
	});
}

void VerletList::put_atoms_in_cells(Atoms& neibs) {
	cell_start.clear();
	cell_endin.clear();

	if (neibs.empty()) return;

	pool.for_each(neibs.begin(), neibs.end(),
		[domain=this->domain](Atom& atom) {
			domain.fit_in_period(atom.r);
		});

	grid.build(neibs);
	grid.set_indices(neibs, atoms_in_cells);

	auto cells_count = grid.cells_number();

	cell_start.resize(cells_count, 0ul);
	cell_endin.resize(cells_count, 0ul);

	pool.sort(
		atoms_in_cells.begin(),
		atoms_in_cells.end());
}

void VerletList::find_cell_starts() {
	auto find_one = [&](point_id i) {
		if (atoms_in_cells[i].cell_id != atoms_in_cells[i - 1].cell_id) {
			auto cell_id = atoms_in_cells[i].cell_id;
			cell_start[cell_id] = i;
			cell_endin[cell_id] = i;
		}
	};

	pool.parallel_for(1ul, atoms_in_cells.size(), find_one);
}

void VerletList::count_atoms_in_cells() {
	auto cell0 = atoms_in_cells[0].cell_id;
	if (pool.active()) {
		auto count_one = [&](std::size_t cell) -> void {
			auto i = cell_start[cell];
			std::size_t count = 0;
			if (i > 0 || cell == cell0) {
				while (i + count < atoms_in_cells.size() &&
					   atoms_in_cells[i + count].cell_id == cell
				) ++count;
				cell_endin[cell] += count;
			}
		};
		// cell-wise cycle for parallel version (to avoid mutex)
		pool.parallel_for(0ul, cell_start.size(), count_one);
	}
	else
		// element-wise cycle for serial version
		for (point_id i = 0; i < atoms_in_cells.size(); ++i) {
			auto cell = atoms_in_cells[i].cell_id;
			++cell_endin[cell];
		}
}

Atoms VerletList::sort(Atoms& atoms) {
	put_atoms_in_cells(atoms);

	Atoms sorted(atoms.size());

	auto func = [&](point_id i) {
        sorted[i] = atoms[atoms_in_cells[i].atom_id];
    };

    pool.parallel_for(0ul, atoms.size(), func);

	return sorted;
}

ListOne VerletList::list_for(vector r) const {
	ListOne result;
	result.reserve(32);

	auto [adj_cells, is_edge] = grid.adjacent_cells(r);

	auto& neibs = *m_neibs;

	for (auto& adjacent : adj_cells) {
		if (adjacent == empty_cell) break;

		for (auto i = head[adjacent]; i != -1; i = next[i]) {
			auto ineib = static_cast<point_id>(i);
			if (are_neibs(
					r, grid.cell_size(),
					neibs[ineib].r, 
					neibs[ineib].rns, 
					is_edge
				)
			) {
				result.push_back(ineib);
			}
		}
	}
	return result;
}

void VerletList::build() {
	auto& neibs = *m_neibs;

	auto build_one = [&](point_id iself) {
		Atom& atom = (*m_atoms)[iself];

		auto [adj_cells, is_edge] = grid.adjacent_cells(atom.r);

		for (auto adjacent: adj_cells) {
			if (adjacent == empty_cell) break;

			auto from = cell_start[adjacent];
			auto to   = cell_endin[adjacent];
			for (auto i = from; i < to; ++i) {
				auto ineib = atoms_in_cells[i].atom_id;

				if (half && ineib < iself) continue;

				if (are_neibs(atom, neibs[ineib], is_edge)) {
					nlist[iself].push_back(ineib);
				}
			}
		}
	};

	pool.parallel_for(0ul, m_atoms->size(), build_one);
}

void VerletList::build_linked_cells() {
	auto& neibs = *m_neibs;

	if (neibs.empty()) return;

	pool.for_each(neibs.begin(), neibs.end(),
		[this](Atom& atom) {
			domain.fit_in_period(atom.r);
		});

	grid.build(neibs);

	const std::size_t n_neibs = neibs.size();
	const std::size_t n_cells = grid.cells_number();

	if (n_neibs == 0 || n_cells == 0) return;

	next.assign(n_neibs, -1);
	head.assign(n_cells, -1);

	for (std::size_t i = 0; i < n_neibs; ++i) {
		auto cell_id = grid.cell_index(grid.local_index3D(neibs[i].r));
		next[i] = head[cell_id];
		head[cell_id] = i;
	}

	auto& atoms = *m_atoms;

	auto build_one = [&] (std::size_t i) {
		Atom& atom = atoms[i];

		auto [adj_cells, is_edge] = grid.adjacent_cells(atom.r);

		for (auto adjacent : adj_cells) {
			if (adjacent == empty_cell) break;

			for (auto j = head[adjacent]; j != -1; j = next[j]) {
				auto& neib = neibs[j];

				if (half && j < i) continue;

				if (are_neibs(atom, neib, is_edge)) {
					nlist[i].push_back(j);
				}
			}
		}
	};

	pool.parallel_for(0ul, m_atoms->size(), build_one);
}

inline double ns_radius(double da, double db) {
	return std::max(da, db);
}

bool VerletList::are_neibs(
	Atom& a, 
	Atom& b, 
	Flag3D is_edge
) const {
	return are_neibs(
		a.r, a.rns,
		b.r, b.rns,
		is_edge
	);
}

bool VerletList::are_neibs(
	const vector& ra,
	double radius_x, 
	const vector& rb,
	double radius_y, 
	Flag3D is_edge
) const {
	double d = ns_radius(radius_x, radius_y);
	auto d2 = d * d;

	auto xa = ra(0), ya = ra(1), za = ra(2);
	auto xb = rb(0), yb = rb(1), zb = rb(2);

	auto drx = xa - xb;
	auto dry = ya - yb;
	auto drz = za - zb;

	auto xsize = domain.xsize(), ysize = domain.ysize(), zsize = domain.zsize();

	if (is_edge[0] && domain.periodic(0) && (2*drx < -xsize)) drx += xsize;
	if (is_edge[0] && domain.periodic(0) && (2*drx >  xsize)) drx -= xsize;
	if (is_edge[1] && domain.periodic(1) && (2*dry < -ysize)) dry += ysize;
	if (is_edge[1] && domain.periodic(1) && (2*dry >  ysize)) dry -= ysize;
	if (is_edge[2] && domain.periodic(2) && (2*drz < -zsize)) drz += zsize;
	if (is_edge[2] && domain.periodic(2) && (2*drz >  zsize)) drz -= zsize;

	return drx*drx + dry*dry + drz*drz <= d2; 
}

} // namespace mdcraft::neibs
