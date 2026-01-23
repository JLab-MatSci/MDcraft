#pragma once

#include <numeric>

#include <mdcraft/neibs/grid.h>

namespace mdcraft::neibs {

using lattice::Domain;
using lattice::dummy_domain;

using ListOne = std::vector<point_id>;
using List = std::vector<ListOne>;

extern Atoms dummy_atoms;
extern Atoms dummy_neibs;

#ifdef mdcraft_ENABLE_CUDA_THRUST
struct NeibsListOneDev
{
	point_id*    start;
	std::size_t   size;

	__host__ __device__ point_id operator[] (std::size_t i)
	{
		return start[i];
	}
};

auto inline create_nlist_device(
	List::iterator first, List::iterator last, std::size_t total
)
{
	std::array<thrust::device_vector<std::size_t>, 3> arr_dev_vecs;

	auto n_atoms = last - first;

	std::vector<std::size_t> flat, lens, offs;

	flat.reserve(total);
	lens.resize(n_atoms);
	offs.resize(n_atoms);

	for (auto it = first; it != last; ++it)
	{
		auto& v = *it;
		std::copy(/*PAR*/ v.begin(), v.end(), std::back_inserter(flat));
	}

	std::transform(/*PAR*/ first, last, lens.begin(), [] (auto&& v) { return v.size(); });

	offs[0] = 0;
	std::inclusive_scan(/*PAR*/ lens.cbegin(), lens.cend(), offs.begin() + 1);

	thrust::device_vector<std::size_t> flat_dev(flat.cbegin(), flat.cend());
	thrust::device_vector<std::size_t> lens_dev(lens.cbegin(), lens.cend());
	thrust::device_vector<std::size_t> offs_dev(offs.cbegin(), offs.cend());

	arr_dev_vecs[0] = std::move(flat_dev);
	arr_dev_vecs[1] = std::move(lens_dev);
	arr_dev_vecs[2] = std::move(offs_dev);

	return arr_dev_vecs;
}
#endif

/** 
	\brief 
		\~russian Список Верле соседних частиц. 
		\~english Verlet list of neighbor atoms.
		\~
*/
class VerletList {
public:
	/** 
		\brief 
			\~russian Конструктор класса. 
			\~english Class constructor.
			\~
		\param[in] atoms
			\~russian ссылка на хранилище с элементами, 
				для которых осуществляется поиск соседей. 
			\~english reference to an array with atoms
				for which neighbors are to be found.
			\~
			\~english indicator of looking for 
				neighbor atoms by all axes
				(separately).
			\~
		\param[in] threads
			\~russian пул тредов для работы функции 
				в многопоточном режиме.
			\~english pool of threads for 
				multithreaded execution.
			\~
	*/
	VerletList(
		Atoms&   atoms,
		Atoms&   neibs,
		Domain&  domain,
		Threads& threads = tools::dummy_pool
	);

	/** 
		\brief 
			\~russian Обновить список соседних элементов.
			\~english Update neighbors list.
			\~
		\details
			\~russian Попутно сортирует элементы.
			\~english Sorts atoms at the same time.
			\~
		\see sort().
	*/
	void update(double kBuf);

	/** 
		\brief 
			\~russian Сортировать элементы в atoms 
				по их одномерным индексам в сетке grid.
			\~english Sort atoms in atoms by their
				1D indices in grid.
			\~
	*/
	Atoms sort(Atoms& atoms);
	
	/** 
		\brief 
			\~russian Получить список соседних частиц 
				для точки x в радиусе grid.cell_size.
			\~english Get neighbors list for a point x 
				using grid.cell_size as a radius.
			\~
		\param[in] x
			\~russian координаты точки.
			\~english point coordinates.
			\~
		\see grid.cell_size.
	*/
	ListOne list_for(vector x) const;

	/** 
		\brief 
			\~russian Получить полный список соседей.
			\~english Get whole neighbors list.
			\~
		\details
			\~russian Возвращает список из списков,
				список соседей для элемента i находится 
				в позиции i общего списка. 
			\~english Returns a list of lists, 
				neighbors list for element i is 
				in position i of the outer list.
			\~
	*/
	List& get();

	/** 
		\brief 
			\~russian Получить список соседей для элемента
				i из external_data.
			\~english Get list for element i 
				of actual external_data Atoms.
			\~
		\param[in] i 
			\~russian индекс элемента
			\~english element index
	*/
	ListOne& operator[](point_id i);

private:
	void optimize_search_radius(double kBuf = -1.0);

	/** 
		\brief 
			\~russian Очистить список соседей.
			\~english Clear neighbors lists.
			\~
	*/
	void clear_lists();

	/** 
		\brief 
			\~russian Составить массив из номеров элементов, 
				упорядоченный по одномерному индексу ячеек, 
				в которых находятся элементы. Одновременно
				с этим модифицируется массив cell_endin и 
				заполняется permutation.
			\~english Make array of element indices 
				sorted by 1D containing cell index. At the 
				same time cell_endin and permutation are modified.
			\~
		\see atoms_in_cells, permutation, cell_endin

	*/    
	void count_atoms_in_cells();

	/** 
		\brief 
			\~russian Найти позиции в atoms_in_cells, 
				начиная с которых у элементов 
				меняется одномерный индекс ячейки. 
			\~english Find atoms_in_cells positions
				of atoms, which start new cell in the 
				atoms_in_cells array.
			\~
	*/
	void find_cell_starts();

	/** 
		\brief 
			\~russian Функция, обеспечивающая заполнение
				массивов grid_indices, cell_start, 
				cell_endin, atoms_in_cells, permutation
				данными об элементах internal_data.
			\~english Function dealing with filling of 
				grid_indices, cell_start, 
				cell_endin, atoms_in_cells, permutation
				with atoms data of internal_data.
			\~
	*/
	void put_atoms_in_cells(Atoms& neibs);

	/** 
		\brief 
			\~russian Функция, проверяющая, являются ли соседними 
				элементы хранилища a и b. 
			\~english The function checks whether atoms a and b 
				of a Atoms are neighbors. 
			\~
		\details
			\~russian Функция проверяет расстояние между элементами
				по полю coords и сравнивает его с neibsSearchRadius.
			\~english Функция проверяет расстояние между элементами
				по полю coords и сравнивает его с neibsSearchRadius.
			\~
		\param[in] a, b
			\~russian элементы хранилища
			\~english Atoms atoms.
			\~
		\param[in] is_edge
			\~russian является ли ячейка элемента a крайней по 
				по каждой из осей. Предполагается переменная типа 
				bool[dim].
			\~english if an atoms cell is boundary by three axes.
				Bool[dim] array type variable is expected.
			\~
		\see dim.
	*/
	bool are_neibs(
		Atom& a, 
		Atom& b, 
		Flag3D is_edge
	) const;

	/** 
		\brief 
			\~russian Функция, проверяющая, являются ли соседними 
				точка и элемент хранилища. 
			\~english The function checks whether a point with 
				neibsSearchRadius and an element are neighbors. 
			\~
		\param[in] x
			\~russian координаты точки.
			\~english point coordinates.
			\~
		\param[in] radius_x
			\~russian радиус поиска соседей точки.
			\~english neighbors search radius of the point.
			\~
		\param[in] is_edge
			\~russian является ли ячейка точки x крайней по 
				по каждой из осей. Предполагается переменная типа 
				bool[dim].
			\~english if a point's cell is boundary by three axes.
				Bool[dim] array type variable is expected.
			\~
		\see dim.
	*/
	bool are_neibs(
		vector x, 
		double radius_x, 
		vector y, 
		double radius_y, 
		Flag3D is_edge
	) const;

	/** 
		\brief 
			\~russian Обновить данные во внутренних массивах
				internal_data, internal_data_alien (если используется) по 
				данным в external_data, external_data_aliens (если используется).
			\~english Update data in inner Atoms
				internal_data, internal_data_alien (if active) using
				external_data, external_data_aliens (if active).
			\~
		\see internal_data, internal_data_alien, external_data, external_data_aliens.
	*/
	void fetch_data(
		Atoms& atoms, 
		Atoms& neibs
	); 

	/** 
		\brief 
			\~russian Составить списки соседних элементов.
			\~english Make neighbors lists.
			\~
		\see list.
	*/
	void build();
	
	std::vector<point_id> cell_start;
		/**<\~russian Вектор, i-й элемент которого показывает
				начало i-й ячейки сетки в векторе atoms_in_cells.
			\~english A vector, element i of which is equal to
				position of a beginning of i-th grid cell in 
				atoms_in_cells array. \~
			\see atoms_in_cells, cell_endin
		*/
	std::vector<point_id> cell_endin;
		/**<\~russian Вектор, i-й элемент которого показывает
				конец i-й ячейки сетки в векторе atoms_in_cells.
			\~english A vector, element i of which is equal to
				position of an end of i-th grid cell in 
				atoms_in_cells array. \~
			\see atoms_in_cells, cell_start
		*/
	
	std::vector<point_id> dummy_nlist;
		/**<\~russian Пустой массив соседей.
			\~english Empty neighbors list. \~
			\see list.
		*/
	
	std::vector<element_ids> atoms_in_cells;
		/**<\~russian Сортированный по cell_id массив индексов частиц.
			\~english Sorted indices array. \~
		*/

	Grid grid;
		/**<\~russian Сетка ячеек.
			\~english Cell grid. \~
		*/
	List nlist;
		/**<\~russian Список соседних элементов.
			\~english Neighbors list. \~
		*/
	
	Atoms* m_atoms = nullptr;
		/**<\~russian Ссылка на данные об элементах, для которых ищутся соседи
		 *            среди массива neighbors.
			\~english Reference to atoms data. \~
		*/
	
	Atoms* m_neibs = nullptr;
		/**<\~russian Ссылка на данные об элементах, среди которых ищутся соседи
		 *            для atoms.
			\~english Dummy empty storage for using in constructor. 
			\~
		*/

	Domain& domain;
	/**<\~russian ссылка на расчетную область для
	              применения периодических гран условий
		\~english reference to a domain to apply
				  periodic boundary conditions
		\~
	*/
	
	Threads& pool;
		/**<\~russian ссылка на пул тредов для работы функции 
				в многопоточном режиме.
			\~english reference to a pool of threads for 
				multithreaded execution.
			\~
		*/
};

} // namespace mdcraft::neibs
