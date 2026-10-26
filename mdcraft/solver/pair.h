#pragma once

#include <mdcraft/tools/threads.h>

#include <mdcraft/solver/isolver.h>
#include <mdcraft/solver/thermostat/base.h>
#include <mdcraft/solver/boundary/base.h>

namespace mdcraft::solver {

using Boundary   = boundary::Base;
using Thermostat = thermostat::Base;

using tools::Threads;

using NeibsListOne = ::mdcraft::neibs::ListOne;

/** \brief
	\~russian Базовый класс решателя для расчета сил
	          межатомного взаимодействия.
	\~english Base solver class for interatomic forces
	          evaluation.
	\~
*/
class Pair : public ISolver {
public:
	/** \brief
		\~russian Конструктор класса.
		\~english Class constructor.
		\~
		\param[in] threads
			\~russian пул тредов для параллельного счета.
			\~english thread pool for parallel calculation.
			\~
	*/
	Pair(
		Boundary&   boundary,
		Thermostat& thermostat = thermostat::dummy_thermostat,
		Threads&    threads    = tools::dummy_pool
	);
	/** \brief
		\~russian Деструктор класса.
		\~english Class destructor.
		\~
	*/
	~Pair() override = default;
#ifdef mdcraft_ENABLE_MPI
    void evaluate(Decomp& decomp, NeibsList& local, NeibsList& alien) override;
#endif

	/** \brief
		\~russian Подготовить  к расчету.
		\~english Prepare atoms to solve.
		\~
		\param[in] atoms
			\~russian данные атомов.
			\~english Atoms data.
			\~
		\details
		\~russian Функция запускает другую
			функцию, \p prepare_one, для всех элементов 
			входного атомы. 
		\~english This function executes \p prepare_one
			for all the \p atoms.
		\~
	*/
	void prepare(
		Atoms&     atoms,
		Atoms&     neibs,
		NeibsList& nlist
		) override;

#ifdef mdcraft_ENABLE_MPI
	void prepare(
		Decomp&    decomp,
		NeibsList& nlist1,
		NeibsList& nlist2) override;
#endif

	/** \brief
		\~russian Рассчитать наведенную зарядовую плотность.
		\~english Calculate charge density.
		\~
		\param[in] atoms
			\~russian данные атомов.
			\~english Atoms data.
			\~
		\details
		\~russian Функция запускает другую
			функцию, \p solve_one, для всех элементов 
			входного атомы. 
		\~english This function executes \p solve_one
			for all the \p atoms.
		\~
	*/
	void charge(
		Atoms&     atoms, 
		Atoms&     neibs, 
		NeibsList& nlist
	);

#ifdef mdcraft_ENABLE_MPI
	void charge(
		Decomp&    decomp,
		NeibsList& nlist1,
		NeibsList& nlist2
	);
#endif

	/** \brief
		\~russian Провести расчет сил.
		\~english Calculate forces actin on atoms.
		\~
		\param[in] atoms
			\~russian Данные атомов.
			\~english Atoms data.
			\~
		\details
		\~russian Функция запускает другую
			функцию, \p solve_one, для всех элементов 
			входного атомы. 
		\~english This function executes \p solve_one
			for all the \p atoms.
		\~
	*/
	void virtual forces(
		Atoms&     atoms,
		Atoms&     neibs,
		NeibsList& nlist
		);

#ifdef mdcraft_ENABLE_MPI
	void virtual forces(
		Decomp&    decomp,
		NeibsList& nlist1,
		NeibsList& nlist2
	);
#endif

	/** \brief
		\~russian Провести расчет вириалов.
		\~english Solve atoms with virial calculation.
		\~
		\param[in] atoms
			\~russian Данные атомов.
			\~english Atoms data.
			\~
		\details
		\~russian Функция запускает другую
			функцию, \p solve_one, для всех элементов 
			входного атомы. 
		\~english This function executes \p solve_one
			for all the \p atoms.
		\~
	*/
	void virials(
		Atoms&     atoms, 
		Atoms&     neibs, 
		NeibsList& nlist
	);

#ifdef mdcraft_ENABLE_MPI
	void virials(
		Decomp&    decomp,
		NeibsList& nlist1,
		NeibsList& nlist2
	);
#endif
	
	/** \brief
		\~russian Подготовить поля атомы 
			для одной частицы к решению.
			(Например, записать 0 в поля структуры, 
			где будет суммироваться решение.)
		\~english Prepare atoms fields 
			of a single particle to solve.
			Ex., write 0 to fields which a used for result 
			accomodating.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ 
				к частице.
			\~english iterator which gives access to 
				a particle.
			\~
	*/
	virtual void prepare_one(
		Atoms::iterator atom,
		Atoms&          neibs, 
		NeibsListOne&   nlist
	) = 0;

	/** \brief
		\~russian Рассчитать зарядовую плотность
			для одной частицы.
		\~english Calculate charge density 
			of a single particle.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ
				к частице.
			\~english iterator which gives access to
				a particle.
			\~
	*/
	virtual void charge_one(		
		Atoms::iterator atom, 
		Atoms&          neibs, 
		NeibsListOne&   nlist
	) = 0;


	/** \brief
		\~russian Расчет сил для атома.
		\~english Calculate forces acting on atom.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ к атому.
			\~english iterator which gives access to an atom.
			\~
		\param[in] neibs
			\~russian Возможные соседи атома \p atom.
			\~english Possible neighboring atoms with \p atom.
			\~
		\param[in] nlist
			\~russian список индексов частиц из \p neibs,
				которые являются соседними с \p atom.
			\~english a list of indices of atoms 
				from \p neibs, which are neighbors of \p atom.
			\~
	*/
	virtual void force_one(
		Atoms::iterator atom, 
		Atoms&          neibs, 
		NeibsListOne&   nlist
	) = 0;

	/** \brief
		\~russian Расчет вириала для атома.
		\~english Solve for an atom with virial calculation.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ к атому.
			\~english iterator which gives access to an atom.
			\~
		\param[in] neibs
			\~russian Возможные соседи атома \p atom.
			\~english Possible neighboring atoms with \p atom.
			\~
		\param[in] nlist
			\~russian список индексов частиц из \p neibs,
				которые являются соседними с \p atom.
			\~english a list of indices of atoms 
				from \p neibs, which are neighbors of \p atom.
			\~
	*/
	virtual void virial_one(
		Atoms::iterator atom, 
		Atoms&          neibs, 
		NeibsListOne&   nlist
	) = 0;

protected:

    void accumulate_virials(Atoms& atoms, Atoms& neibs, NeibsList& nlist, bool reset);

    // Accumulate one neighbour-list contribution. In MPI mode, thermostat
    // application is deferred until both local and alien passes complete.
    void accumulate_forces(
        Atoms& atoms, Atoms& neibs, NeibsList& nlist,
        bool apply_thermostat
    );

	Thermostat& thermostat;
	Boundary&   boundary;
	Threads&    pool;
		/**<\~russian ссылка на пул тредов для работы функции 
				в многопоточном режиме.
			\~english reference to a pool of threads for 
				multithreaded execution.
			\~
		*/

	template <typename Func>
	void process_ghosts(
		Atoms::iterator atomit,
		Atoms&           neibs,
		NeibsListOne&    nlist,
		Func&&	         func
	) {
		auto& atom = *atomit;
		auto&& ghosts = boundary.prepare_ghosts(atomit, neibs, nlist);
		std::size_t sz = ghosts.size();
		if (sz == 0) return;
		NeibsListOne nlist_ghosts(sz);
		for (std::size_t i = 0; i < sz; nlist_ghosts[i] = i, i++) {}
        // Suppress recursive ghost expansion, but retain ALL boundary flags
        // even if the callback fails. The flags need not equal 1.
        const auto saved_bc = atom.bc;
        atom.bc = 0;
        try {
            func(atomit, ghosts, nlist_ghosts);
        } catch (...) {
            atom.bc = saved_bc;
            throw;
        }
        atom.bc = saved_bc;
	}
};

} // namespace mdcraft::solver
