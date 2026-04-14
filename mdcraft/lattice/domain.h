#pragma once 

#include <array>
#include <mdcraft/data/vector.h>

namespace mdcraft::lattice {

/**
	\brief
		\~russian Расчетная область моделирования.
		\~english Modeling domain.
		\~
*/
class Domain {
public:
	Domain() = default;
	/** \brief
		\~russian Конструктор класса.
		\~english Class constructor.
		\~ */
	Domain(
		double xmin, double xmax,
		double ymin, double ymax,
		double zmin, double zmax);

	void reshape(
		double xmin, double xmax,
		double ymin, double ymax,
		double zmin, double zmax);

	/**
		\brief
			\~russian Сделать границы по оси \a axis периодическими.
			\~english Make boundaries by \a axis periodic.
			\~
		\param[in] axis
			\~russian ось, по которой назначить периодичность. По умолчанию \a Ox.
			\~english axis for periodicity setting. Default is \a Ox.
			\~
		\param[in] val
			\~russian булево значение периодичности.
			\~english bool (true or false) variable for periodicity
			\~
	*/
	void set_periodic(int axis, bool val) { m_periodic[axis] = val; }

	/// @brief Dimension of the domain
	int dim() const { return m_dim; }

	double xmin() const { return m_bmin[0]; }
	double ymin() const { return m_bmin[1]; }
	double zmin() const { return m_bmin[2]; }

	double xmax() const { return m_bmax[0]; }
	double ymax() const { return m_bmax[1]; }
	double zmax() const { return m_bmax[2]; }

	double xsize() const { return m_size[0]; }
	double ysize() const { return m_size[1]; }
	double zsize() const { return m_size[2]; }

	bool xperiodic() { return m_periodic[0]; }
	bool yperiodic() { return m_periodic[1]; }
	bool zperiodic() { return m_periodic[2]; }


	double bmin(int axis) const { return m_bmin[axis]; }
	double bmax(int axis) const { return m_bmax[axis]; }
	double size(int axis) const { return m_size[axis]; }
	bool periodic(int axis) const { return m_periodic[axis]; }

	/**
		\brief
			\~russian Узнать, находится ли точка в расчетной области.
			\~english Find out if a point is in the domain.
			\~
		\param[in] p
			\~russian координаты точки.
			\~english point coordinates.
			\~
	*/
	bool belong(data::vector p) const;

	/**
		\brief
			\~russian Разместить точку в главном периоде.
			\~english Fit a point into the main period.
			\~
		\param[in] p
			\~russian координаты точки.
			\~english point coordinates.
			\~
	*/
	void fit_in_period(data::vector& p) const;

	/**
		\brief
			\~russian Вернуть кратчайший вектор с учетом периодов.
			\~english Get shortest vector using all the periods.
			\~
		\details
			\~russian Функция может быть использована для нахождения
				расстояния между точками с учетом периодичности расчетной области.
				Входной вектор может быть не кратчайшим расстоянием между точками,
				результатом же будет кратчайший.
			\~english Function may be used to find distance between points
				in periodic domain. Input vector may be any distance between points,
				resulting vector is the shortest one.
			\~
		\see periods.
		\param[in] v
			\~russian вектор.
			\~english a vector.
			\~
	*/
	data::vector shortest(data::vector v) const;

	double nearest_distance(data::vector v, int axis) const;

private:
	int m_dim;
	std::array<double, 3> m_bmin;
	std::array<double, 3> m_bmax;
	std::array<double, 3> m_size;
	std::array<bool,   3> m_periodic;
};

static Domain dummy_domain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

} // namespace mdcraft::lattice
