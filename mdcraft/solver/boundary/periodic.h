#pragma once
#include <mdcraft/solver/boundary/base.h>
#include <mdcraft/solver/pair.h>

namespace mdcraft {
namespace solver {
namespace boundary {

/** 
	\brief 
		\~russian Тип периодических граничных условий.
		\~english Periodic boundary condition type.
		\~
*/
class Periodic : public Base {
public: 
	/** 
		\brief 
			\~russian Конструктор класса.
			\~english Class constructor.
			\~
		\param[in] solver 
			\~russian ссылка на объект решателя, который будет
				использован для расчета пары частиц в периоде. 
			\~english reference to solver object, which is used to 
				solve pairs of MD particles.
			\~
	*/
	Periodic(Domain& domain);

	void prepare_one(
		Atoms::iterator atom, 
		Atoms&          neibs,
		NeibsListOne&   nlist
	) override;

	void charge_one(
		Atoms::iterator atom, 
		Atoms&          neibs,
		NeibsListOne&   nlist
	) override;

	void force_one(
		Atoms::iterator atom,
		Atoms&          neibs,
		NeibsListOne&   nlist
	) override;

	void virial_one(
		Atoms::iterator atom,
		Atoms&          neibs,
		NeibsListOne&   nlist
	) override;

	Atoms prepare_ghosts (
		Atoms::iterator atom,
		Atoms&          neibs,
		NeibsListOne&   nlist
	) override;

private:
	Atoms make_ghosts(Atoms::iterator atom, Atoms& neibs,
	                  NeibsListOne& nlist, bool with_density);

	unsigned short need_apply(
		Atoms::iterator atom
	);

	static const unsigned int type = 0x1;
};

} // boundary
} // solver
} // mdcraft
