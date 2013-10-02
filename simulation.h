
#ifndef SIMULATION_H
#define SIMULATION_H

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include "fwk/LinkedList.h"
#include "Tissue.h"
#include <map>


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

	void cytotoxicCellNew(Fwk::String tissue, Cell::Coordinates loc);

	void helperCellNew(Fwk::String tissue, Cell::Coordinates loc);

	void infectionStart(Fwk::String tissue, Cell::Coordinates loc, 
	                    CellMembrane::Side side, AntibodyStrength strength);

	void infectedCellsDel(Fwk::String tissue);

	void cloneNew(Fwk::String tissue, Cell::Coordinates loc, 
								CellMembrane::Side side);

	void setAntibodyStrength(Fwk::String tissue, Cell::Coordinates loc,
	                         CellMembrane::Side side, AntibodyStrength strength);

	void cloneCellsNew(Fwk::String tissue, CellMembrane::Side side);

	void parseCommand(Fwk::String textLine);

	int test() {return 2;}


protected:
	Simulation(Fwk::String _name);
	~Simulation() {}
	Cell::Coordinates getCoordinate(tokenizer<>::iterator token);
	CellMembrane::Side getSide(tokenizer<>::iterator token);
	Fwk::String coordToStr(Cell::Coordinates c);

	map<Fwk::String, Tissue::Ptr> tissues_;

};


#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include "Tissue.h"
#include "simulation.h"

using namespace std;
using namespace boost;

#endif
