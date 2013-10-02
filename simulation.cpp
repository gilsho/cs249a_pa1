#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include "Tissue.h"
#include "simulation.h"

using namespace std;
using namespace boost;

Simulation::Simulation(Fwk::String _name) : Fwk::NamedInterface(_name) {}

/*
Create a new tissue named "name". Quoted names are not accepted and therefore 
names cannot contain any whitespace. NOTE: The tissue name WILL NOT be any of 
the commands, i.e. there will not be a tissue named "tissueNew".
*/
void Simulation::tissueNew(Fwk::String tissue) {
  cout << "creating new tissue named: " << tissue << "..." << endl;
}

/*
Create a new healthy CytotoxicCell in "tissue" at location "loc". In case an 
exception is thrown (ex. a cell already exists at location), your program 
should catch the exception and, if you wish, use some logging mechanism to 
report it. In any case, your simulation should continue running as if this 
entry did not exist. A new CytotoxicCell should have antibody strength of 100 
on all its membranes.
*/
void Simulation::cytotoxicCellNew(Fwk::String tissue, Cell::Coordinates loc)
{
  cout << "creating new cytoToxic cell in: " << tissue << ", at: " 
    << coordToStr(loc) << endl;
}

/*
Create a new healthy HelperCell in "tissue" at location "loc". This command is 
identical to the cytotoxicCellNew command, except that it creates a new 
HelperCell, which has antibody strength 0 on all its membranes.
*/
void Simulation::helperCellNew(Fwk::String tissue, Cell::Coordinates loc)
{
  cout << "creating new helper cell in: " << tissue << ", at: " 
    << coordToStr(loc) << endl;
}

/*
Starts an infection of strength 99 at cell at "loc" entering from the "loc"
membrane. You should proceed to the next command only when no more cells can be 
infected. At the end of the infection round, you should print statistics to 
standard out as described here.
*/
void Simulation::infectionStart(Fwk::String tissue, Cell::Coordinates loc, 
                    CellMembrane::Side side, AntibodyStrength strength)
{
  cout << "starting infection in: " << tissue << "at: " << coordToStr(loc) 
    << ", side:" << side << " with strength: " << strength << endl;
}

/*
Remove all infected cells from "tissue".
*/
void Simulation::infectedCellsDel(Fwk::String tissue)
{
  cout << "removing infected cells in: " << tissue << endl;
}

/*
Clones cell at location "loc" and places the new cell "side" of "loc" (x,y+1,z) 
in "tissue". Like the other cell creation commands, the simulation should 
continue running despite any exception that may be thrown.
*/

void Simulation::cloneNew(Fwk::String tissue, Cell::Coordinates loc, 
              CellMembrane::Side side)
{
  cout << "cloning cell in: " << tissue << "at: " << coordToStr(loc) 
    << " in direction: " << side << endl;
}

void Simulation::setAntibodyStrength(Fwk::String tissue, Cell::Coordinates loc,
                         CellMembrane::Side side, AntibodyStrength strength)
{
  cout << "setting membrane strength of cell in: " << tissue << "at: " 
    << coordToStr(loc) << " dir: " << side << " strength: " 
    << strength << endl;
}

/*
Tissue Tissue1 cloneCellsNew west — Clones all cells in "tissue" to the "loc"
direction. Equivalent to writing Cell x y z cloneNew "loc" for each cell in
the tissue. If any single cell throws an exception, you should continue the 
simulation and clone the remaining cells.
*/
void Simulation::cloneCellsNew(Fwk::String tissue, CellMembrane::Side side) 
{
  cout << "cloning cells in tissue: " << tissue << " to: " 
    << side << endl;
}


Cell::Coordinates Simulation::getCoordinate(tokenizer<>::iterator token)
{
  Cell::Coordinates loc;
  loc.x = lexical_cast<int>(*token++);
  loc.y = lexical_cast<int>(*token++);
  loc.z = lexical_cast<int>(*token);
  return loc;
}

CellMembrane::Side Simulation::getSide(tokenizer<>::iterator token)
{
  if (*token == "north")
    return CellMembrane::north_;

  if (*token == "south")
    return CellMembrane::south_;

  if (*token == "east")
    return CellMembrane::east_;

  if (*token == "west")
    return CellMembrane::west_;

  if (*token == "up")
    return CellMembrane::up_;

  if (*token == "down")
    return CellMembrane::down_;

  throw "Unrecognized membrane side";

}

void Simulation::parseCommand(Fwk::String textLine) 
{
  if (textLine == "" || textLine[0] == '#')
    return;

  tokenizer<> tokenizedLine(textLine);
  tokenizer<>::iterator token=tokenizedLine.begin(); 

  if (*token == "Tissue") {
    token++;
    if (*token == "tissueNew") {
      token++;
      tissueNew(*token);
    } else {
      Fwk::String tissue = *token;
      token++;
      if (*token == "cytotoxicCellNew") {
        token++;
        Cell::Coordinates loc = getCoordinate(token);
        cytotoxicCellNew(tissue, loc);
      } else if (*token == "helperCellNew") {
        token++;
        Cell::Coordinates loc = getCoordinate(token);
        helperCellNew(tissue, loc);
      } else if (*token == "infectionStartLocationIs") {
        token++;
        Cell::Coordinates loc = getCoordinate(token);
        token++; token++; token++;
        CellMembrane::Side side = getSide(token++);
        AntibodyStrength strength = 
          AntibodyStrength(lexical_cast<int>(*token++)); 
        infectionStart(tissue, loc, side, strength);
      } else if (*token == "infectedCellsDel") {
        infectedCellsDel(tissue);
      } else if (*token == "cloneCellsNew") {
        token++;
        CellMembrane::Side side = getSide(token++);
        cloneCellsNew(tissue, side);
      } else {
        throw "Malformed command";
      }
    }
  } else if (*token == "Cell") {
    token++;
    Fwk::String tissue = *token++;
    Cell::Coordinates loc = getCoordinate(token);
    token++; token++; token++;
    if (*token == "membrane") {
      token++;
      CellMembrane::Side side = getSide(token++);
      if (*token == "antibodyStrengthIs") {
        *token++;
        AntibodyStrength strength = 
          AntibodyStrength(lexical_cast<int>(*token++));
        setAntibodyStrength(tissue, loc, side, strength);
      } else {
        throw "Malformed command";
      }
    } else if (*token == "cloneNew") {
      token++;
      CellMembrane::Side side = getSide(token++);
      cloneNew(tissue, loc, side);
    } else {
      throw "Malformed command";
    }
  } else {
    throw "Malformed command";
  }
}

Fwk::String Simulation::coordToStr(Cell::Coordinates c) 
{
  Fwk::String s;
  s = "(" + lexical_cast<Fwk::String>(c.x) + "," + lexical_cast<Fwk::String>(c.y) + "," 
    + lexical_cast<Fwk::String>(c.z) + ")";
  return s;
}

