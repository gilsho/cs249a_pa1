
#ifndef SIMULATION_H
#define SIMULATION_H

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <map>
#include "fwk/LinkedList.h"
#include "Tissue.h"

using namespace std;
using namespace boost;


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

	Tissue::Ptr tissueNew(Fwk::String _name);

	Tissue::Ptr tissue(Fwk::String _name);

	Cell::Ptr cellNew(Fwk::String tissue, Cell::Coordinates loc,
							 Cell::CellType ctype);

	void infectionStart(Fwk::String tissue, Cell::Coordinates loc, 
	                    CellMembrane::Side side, AntibodyStrength strength);

	void infectedCellsDel(Fwk::String tissue);

	void cloneNew(Fwk::String tissue, Cell::Coordinates loc, 
								CellMembrane::Side side);

	void antibodyStrengthIs(Fwk::String tissue, Cell::Coordinates loc,
	                         CellMembrane::Side side, AntibodyStrength strength);

	void cloneCellsNew(Fwk::String tissue, CellMembrane::Side side);

	void commandIs(Fwk::String textLine);


protected:

	static const U32 initialCytotoxicStrength = 100;
	static const U32 initialHelperStrength = 0;

	class TissueReactor : public Tissue::Notifiee
	{ 
		public:
			virtual void onCellNew( Cell::Ptr );
			static TissueReactor *TissueReactorIs(Tissue *t) {
				return new TissueReactor(t);
			}
			Simulation *psim;
		protected:
			TissueReactor(Tissue *t) : Tissue::Notifiee() {}
	};

	Simulation(Fwk::String _name);
	~Simulation() {}
	Cell::Coordinates coordinateIs(tokenizer<char_separator<char> >::iterator token);
	CellMembrane::Side sideIs(tokenizer<char_separator<char> >::iterator token);
	Cell::Coordinates coordinateShifted(Cell::Coordinates loc, 
                                     CellMembrane::Side side);
	bool infectionSpreadTo(Cell::Ptr c, CellMembrane::Side side, 
												 AntibodyStrength strength, S32& difference);
	Cell::Ptr neighbor(Tissue::Ptr t, Cell::Ptr c, CellMembrane::Side side);
	CellMembrane::Side oppositeSide(CellMembrane::Side side);
	void stats(Fwk::String _tissue, U32 attempts, S32 difference, U32 path);

	map<Fwk::String, Tissue::Ptr> tissues_;

	U32 infectionVolume(Fwk::String _tissue);
	U32 infectedCells(Fwk::String _tissue); 
	U32 cytotoxicCells_;
	U32 helperCells_;
};

#endif
