#include <cstring>
#include <mdcraft/neibs/verlet-list.h>

namespace mdcraft::neibs {

Atoms dummy_atoms(1);
Atoms dummy_neibs(1);

VerletList::VerletList(
	Atoms&   atoms,
	Atoms&   neibs,
	Domain&  domain,
	Threads& pool
) : grid(domain, pool),
	domain(domain),
	pool(pool)
{
	fetch_data(atoms, neibs);
}

void VerletList::update(double kbuf) {
	if (m_neibs->empty()) {
		clear_lists();
		return;
	}
	if (kbuf > 0.0) optimize_search_radius(kbuf);
	clear_lists();
	put_atoms_in_cells(*m_neibs);
	find_cell_starts();
	count_atoms_in_cells();
	build();
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
		assert(atom.rcut > 0);
		auto iatom = p - m_atoms->begin();
		auto& neibs = *m_neibs;
		double hmax = 0.0;
		for (int i = 0; i < nlist[iatom].size(); ++i) {
			auto& neighbor = neibs[i];

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
	for (auto& l: nlist) l.clear();
}

ListOne& VerletList::operator[](point_id i) {
	return i < nlist.size() ? nlist[i] : dummy_nlist; 
}

List& VerletList::get() { 
	return nlist; 
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

		auto from = cell_start[adjacent];
		auto to   = cell_endin[adjacent];

		for (auto i = from; i != to; ++i) {
			auto ineib = atoms_in_cells[i].atom_id;
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
				if (are_neibs(atom, neibs[ineib], is_edge)) {
					nlist[iself].push_back(ineib);
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
	vector x, 
	double radius_x, 
	vector y, 
	double radius_y, 
	Flag3D is_edge
) const {

	double d = ns_radius(radius_x, radius_y);

	auto d2 = d * d;
	vector dist = x - y;

	vector sizes = {domain.xsize(), domain.ysize(), domain.zsize()};

	for (int i = 0; i < domain.dim(); ++i) {
		if (is_edge[i] && domain.periodic(i)) {
			dist[i] = std::fabs(dist[i] + sizes[i]) < std::fabs(dist[i]) ? dist[i] + sizes[i] : dist[i];
			dist[i] = std::fabs(dist[i] - sizes[i]) < std::fabs(dist[i]) ? dist[i] - sizes[i] : dist[i];
		}
	}
	return dist.dot(dist) <= d2;
}

} // namespace mdcraft::neibs
