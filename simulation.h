
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
		protected:
			TissueReactor(Tissue *t) : Tissue::Notifiee() {}
	};

	Simulation(Fwk::String _name);
	~Simulation() {}
	Cell::Coordinates coordinateIs(tokenizer<>::iterator token);
	CellMembrane::Side sideIs(tokenizer<>::iterator token);
	Cell::Coordinates coordinateShifted(Cell::Coordinates loc, 
                                     CellMembrane::Side side);
	bool infectionSpreadTo(Cell::Ptr c, CellMembrane::Side side, AntibodyStrength strength);
	Cell::Ptr neighbor(Tissue::Ptr t, Cell::Ptr c, CellMembrane::Side side);
	CellMembrane::Side oppositeSide(CellMembrane::Side side);


	map<Fwk::String, Tissue::Ptr> tissues_;

	/*
	The statistics should be in the following format: 

	a b c d e f g 

	Where a is the total number of infected cells, b is the total infection 
	attempts in that infection round, c is the total difference between disease 
	strength and antibody strength for all infection attempts in that round, d is 
	the total number of cytotoxic cells alive (infected and healthy), e is the 
	number of helper cells (infected and healthy), f is the infection spread, and 
	g is the length of the longest infection path measured from the root of 
	infection. The infection spread is defined as the volume of the smallest 
	rectangular box containing all infected cells. 
	*/

	U32 infectedCells() { return infectedCells_; }
	U32 infectedAttempts() { return infectedAttempts_; }
	U32 strengthDifference() { return strengthDifference_; }
	U32 cytotoxicCells() { return cytotoxicCells_; }
	U32 helperCells() { return helperCells_; }
	U32 infectionSpread() { return infectionSpread_; }
	U32 longestInfectionPath() { return longestInfectionPath_; }

	U32 infectedCells_;
	U32 infectedAttempts_;
	U32 strengthDifference_;
	U32 cytotoxicCells_;
	U32 helperCells_;
	U32 infectionSpread_;
	U32 longestInfectionPath_;
};

#endif
