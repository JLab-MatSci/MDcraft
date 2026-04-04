#pragma once 
#include <mdcraft/neibs/grid.h>

namespace mdcraft::neibs {

using lattice::Domain;
using lattice::dummy_domain;

class ListOne {
public:
    using value_type = std::vector<point_id>::value_type;
    using reference = std::vector<point_id>::reference;
    using const_reference = std::vector<point_id>::const_reference;
    using iterator = std::vector<point_id>::iterator;
    using const_iterator = std::vector<point_id>::const_iterator;
    using reverse_iterator = std::vector<point_id>::reverse_iterator;
    using const_reverse_iterator = std::vector<point_id>::const_reverse_iterator;
    using size_type = std::vector<point_id>::size_type;
    using difference_type = std::vector<point_id>::difference_type;

    ListOne() = default;
    explicit ListOne(bool half) : half_(half) {}
    explicit ListOne(size_type n, bool half = false) : half_(half), data_(n) {}
    ListOne(size_type n, const value_type& val, bool half = false) : half_(half), data_(n, val) {}
    template <class InputIt>
    ListOne(InputIt first, InputIt last, bool half = false) : half_(half), data_(first, last) {}
    ListOne(std::initializer_list<value_type> init, bool half = false) : half_(half), data_(init) {}
    ListOne(const ListOne& other) : half_(other.half_), data_(other.data_) {}
    ListOne(ListOne&& other) noexcept : half_(other.half_), data_(std::move(other.data_)) {}

    ListOne& operator=(const ListOne& other) {
        if (this != &other) {
            half_ = other.half_;
            data_ = other.data_;
        }
        return *this;
    }

    ListOne& operator=(ListOne&& other) noexcept {
        if (this != &other) {
            half_ = other.half_;
            data_ = std::move(other.data_);
        }
        return *this;
    }

    ListOne& operator=(std::initializer_list<value_type> init) {
        data_ = init;
        return *this;
    }

    reference operator[](size_type pos) { return data_[pos]; }
    const_reference operator[](size_type pos) const { return data_[pos]; }
    reference at(size_type pos) { return data_.at(pos); }
    const_reference at(size_type pos) const { return data_.at(pos); }
    reference front() { return data_.front(); }
    const_reference front() const { return data_.front(); }
    reference back() { return data_.back(); }
    const_reference back() const { return data_.back(); }
    value_type* data() { return data_.data(); }
    const value_type* data() const { return data_.data(); }

    iterator begin() { return data_.begin(); }
    const_iterator begin() const { return data_.begin(); }
    const_iterator cbegin() const { return data_.cbegin(); }
    iterator end() { return data_.end(); }
    const_iterator end() const { return data_.end(); }
    const_iterator cend() const { return data_.cend(); }
    reverse_iterator rbegin() { return data_.rbegin(); }
    const_reverse_iterator rbegin() const { return data_.rbegin(); }
    const_reverse_iterator crbegin() const { return data_.crbegin(); }
    reverse_iterator rend() { return data_.rend(); }
    const_reverse_iterator rend() const { return data_.rend(); }
    const_reverse_iterator crend() const { return data_.crend(); }

    bool empty() const { return data_.empty(); }
    size_type size() const { return data_.size(); }
    size_type max_size() const { return data_.max_size(); }
    void reserve(size_type new_cap) { data_.reserve(new_cap); }
    size_type capacity() const { return data_.capacity(); }
    void shrink_to_fit() { data_.shrink_to_fit(); }

    void clear() { data_.clear(); }
    iterator insert(const_iterator pos, const value_type& value) { return data_.insert(pos, value); }
    iterator insert(const_iterator pos, value_type&& value) { return data_.insert(pos, std::move(value)); }
    iterator insert(const_iterator pos, size_type count, const value_type& value) { return data_.insert(pos, count, value); }
    template <class InputIt>
    iterator insert(const_iterator pos, InputIt first, InputIt last) { return data_.insert(pos, first, last); }
    iterator insert(const_iterator pos, std::initializer_list<value_type> init) { return data_.insert(pos, init); }
    template <class... Args>
    iterator emplace(const_iterator pos, Args&&... args) { return data_.emplace(pos, std::forward<Args>(args)...); }
    iterator erase(const_iterator pos) { return data_.erase(pos); }
    iterator erase(const_iterator first, const_iterator last) { return data_.erase(first, last); }
    void push_back(const value_type& value) { data_.push_back(value); }
    void push_back(value_type&& value) { data_.push_back(std::move(value)); }
    template <class... Args>
    reference emplace_back(Args&&... args) { return data_.emplace_back(std::forward<Args>(args)...); }
    void pop_back() { data_.pop_back(); }
    void resize(size_type count) { data_.resize(count); }
    void resize(size_type count, const value_type& value) { data_.resize(count, value); }
    void swap(ListOne& other) {
        std::swap(half_, other.half_);
        data_.swap(other.data_);
    }

    bool half() const { return half_; }
    void set_half(bool half) { half_ = half; }

    bool operator==(const ListOne& other) const { return half_ == other.half_ && data_ == other.data_; }
    bool operator!=(const ListOne& other) const { return !(*this == other); }
    bool operator<(const ListOne& other) const { return data_ < other.data_; }
    bool operator<=(const ListOne& other) const { return data_ <= other.data_; }
    bool operator>(const ListOne& other) const { return data_ > other.data_; }
    bool operator>=(const ListOne& other) const { return data_ >= other.data_; }

private:
    bool half_ = false;
    std::vector<point_id> data_;
};

inline void swap(ListOne& lhs, ListOne& rhs) { lhs.swap(rhs); }

using List = std::vector<ListOne>;

extern Atoms dummy_atoms;
extern Atoms dummy_neibs;

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
		Threads& threads = tools::dummy_pool,
		bool     half = false
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

	auto& get_pool() { return pool; }

	bool decide() const;

	bool is_half() const { return half; }

	void set_half(bool b) {
		half = b;
		for (auto& l : nlist) {
			l.set_half(half);
		}
	}

private:
	void optimize_search_radius(double kBuf = -1.0);
	void snapshot_atoms();

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
		const vector& x, 
		double radius_x, 
		const vector& y, 
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

	void build_linked_cells();

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
	
	ListOne dummy_nlist;
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

	std::vector<int> head;
	std::vector<int> next;

	bool half;

	struct atom_state {
		vector r = vector::Zero();
		double rcut = 0.0;
		double rns = 0.0;
		unsigned int uid = 0;
	};

	std::vector<atom_state> saved_atoms;
};

} // namespace mdcraft::neibs
