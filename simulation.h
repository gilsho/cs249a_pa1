
#ifndef SIMULATION_H
#define SIMULATION_H

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <map>
#include "fwk/LinkedList.h"
#include "Tissue.h"

using namespace std;
using namespace boost;

#define DEFAULT_CELL_TYPE (Cell::helperCell())

class Simulation : public Fwk::NamedInterface
{

public:

  typedef Fwk::Ptr<Simulation const> PtrConst;
  typedef Fwk::Ptr<Simulation> Ptr;
	static Simulation::Ptr SimulationNew(Fwk::String _name) {
      Ptr s = new Simulation(_name);
      s->referencesDec(1);
      return s;
   }

	Cell::Ptr cellNew(Cell::Coordinates loc, Cell::CellType ctype);

	void cloneNew(Cell::Coordinates loc, CellMembrane::Side side);

	void cloneCellsNew(CellMembrane::Side side);

	void antibodyStrengthIs(Cell::Coordinates loc, CellMembrane::Side side, 
													AntibodyStrength strength);

	void infectionStart(Cell::Coordinates loc, CellMembrane::Side side, 
											AntibodyStrength strength);

	void infectedCellsDel();

	Tissue::Ptr tissue();

protected:

	static const U32 initialCytotoxicStrength = 100;
	static const U32 initialHelperStrength = 0;

	class TissueReactor : public Tissue::Notifiee
	{ 
		public:
			virtual void onCellNew( Cell::Ptr );
			virtual void onCellDel( Cell::Ptr );

			static TissueReactor *TissueReactorIs(Tissue *t) {
				return new TissueReactor(t);
			}
			Simulation *psim;
		protected:
			TissueReactor(Tissue *t) : Tissue::Notifiee() {}
	};

	Simulation(Fwk::String _name);
	~Simulation() {}
	Cell::Coordinates coordinateShifted(Cell::Coordinates loc, 
                                     CellMembrane::Side side);
	bool infectionSpreadTo(Cell::Ptr c, CellMembrane::Side side, 
                                   AntibodyStrength attack, 
                                   S32& difference,
                                   U32& attempts);
	Cell::Ptr neighbor(Cell::Ptr c, CellMembrane::Side side);
	CellMembrane::Side oppositeSide(CellMembrane::Side side);
	void stats(U32 attempts, S32 difference, U32 path);

	Tissue::Ptr t;

	U32 infectionVolume();
	U32 infectedCells(); 
	U32 cytotoxicCells_;
	U32 helperCells_;
};

#endif
