#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include "Tissue.h"

using namespace std;
using namespace boost;


class Simulation : public Fwk::NamedInterface
{

public:
	Simulation();
	~Simulation();

	static void parseCommand(Fwk::String textLine);

	void tissueNew(Fwk::String tissue);
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

private:
	Cell::Coordinates  getCoordinate(tokenizer<>::iterator token);
	CellMembrane::Side getSide(tokenizer<>::iterator token);
	Fwk::String coordToStr(Cell::Coordinates c);


}