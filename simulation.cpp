#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include "simulation.h"
#include <queue>

using namespace std;
using namespace boost;

template <typename T>
void assertValidPtr(Fwk::Ptr<T> p) {
  if (p.ptr() == NULL) {
    throw "null pointer exception";
  }
}

void Simulation::TissueReactor::onCellNew(Cell::Ptr c)
{
  c->membraneNew(c->name() + " north", CellMembrane::north());
  c->membraneNew(c->name() + " south", CellMembrane::south());
  c->membraneNew(c->name() + " east", CellMembrane::east());
  c->membraneNew(c->name() + " west", CellMembrane::west());
  c->membraneNew(c->name() + " up", CellMembrane::up());
  c->membraneNew(c->name() + " down", CellMembrane::down());
}

Simulation::Simulation(Fwk::String _name) : Fwk::NamedInterface(_name) { }

/*
Create a new tissue named "name". Quoted names are not accepted and therefore 
names cannot contain any whitespace. NOTE: The tissue name WILL NOT be any of 
the commands, i.e. there will not be a tissue named "tissueNew".
*/
Tissue::Ptr Simulation::tissueNew(Fwk::String _name) {
  Tissue::Ptr t = Tissue::TissueNew(_name);
  TissueReactor::Ptr r = TissueReactor::TissueReactorIs(t.ptr());
  r->notifierIs(t);
  tissues_[_name] = t;
  return t;
}

Tissue::Ptr Simulation::tissue(Fwk::String _name)
{
  return tissues_[_name];
}

/*
Create a new healthy CytotoxicCell in "tissue" at location "loc". In case an 
exception is thrown (ex. a cell already exists at location), your program 
should catch the exception and, if you wish, use some logging mechanism to 
report it. In any case, your simulation should continue running as if this 
entry did not exist. A new CytotoxicCell should have antibody strength of 100 
on all its membranes.
*/
void Simulation::cytotoxicCellNew(Fwk::String _tissue, Cell::Coordinates loc)
{
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);
  Cell::Ptr c = Cell::CellNew(loc, t.ptr(), Cell::cytotoxicCell_);
  assertValidPtr(c);
  t->cellIs(c);
  //raise exception if cell already exists??
}

/*
Create a new healthy HelperCell in "tissue" at location "loc". This command is 
identical to the cytotoxicCellNew command, except that it creates a new 
HelperCell, which has antibody strength 0 on all its membranes.
*/
void Simulation::helperCellNew(Fwk::String _tissue, Cell::Coordinates loc)
{
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);
  Cell::Ptr c = Cell::CellNew(loc, t.ptr(), Cell::helperCell_);
  assertValidPtr(c);
  t->cellIs(c);
  //raise exception if cell already exists??
}

/*
Starts an infection of strength 99 at cell at "loc" entering from the "loc"
membrane. You should proceed to the next command only when no more cells can be 
infected. At the end of the infection round, you should print statistics to 
standard out as described here.
*/
void Simulation::infectionStart(Fwk::String _tissue, Cell::Coordinates loc, 
                    CellMembrane::Side side, AntibodyStrength strength)
{
  // Tissue::Ptr t = tissues_[_tissue];
  // assertValidPtr(t);
  // Cell::Ptr c = *(t->cellIter(loc));
  // CellMembrane::Ptr m = *(c->membraneIter(side));
  // if (m->antibodyStrength() >= strength)

  //   return;
}

/*
Remove all infected cells from "_tissue".
*/
void Simulation::infectedCellsDel(Fwk::String _tissue)
{
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);
  Tissue::CellIterator it = t->cellIter();
  for (int i = 0; i < t->cells(); i++, ++it) {
    Cell::Ptr c = *it;
    assertValidPtr(c);
    if (c->health() == Cell::infected()) {
      t->cellDel(c->name());
    }
  }
}

/*
Clones cell at location "loc" and places the new cell "side" of "loc" (x,y+1,z) 
in "_tissue". Like the other cell creation commands, the simulation should 
continue running despite any exception that may be thrown.
*/

void Simulation::cloneNew(Fwk::String _tissue, Cell::Coordinates loc, 
              CellMembrane::Side side)
{
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);
  Cell::Ptr c = *(t->cellIter(loc));
  assertValidPtr(c);
  Cell::Coordinates cloneLoc = loc;
    if (side == CellMembrane::north())
      cloneLoc.y++;
    else if (side == CellMembrane::south())
      cloneLoc.y--;
    else if (side == CellMembrane::east())
     cloneLoc.x++;
    else if (side == CellMembrane::west())
      cloneLoc.x--;
    else if (side == CellMembrane::up())
      cloneLoc.z++;
    else if (side == CellMembrane::down())
      cloneLoc.z--;
  Cell::Ptr clone = Cell::CellNew(cloneLoc, t.ptr(), c->cellType());
  t->cellIs(clone);
}

void Simulation::setAntibodyStrength(Fwk::String _tissue, Cell::Coordinates loc,
                         CellMembrane::Side side, AntibodyStrength strength)
{
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);
  Cell::Ptr c = *(t->cellIter(loc));
  assertValidPtr(c);
  CellMembrane::Ptr m = *(c->membraneIterConst(side));
  assertValidPtr(m);
  m->antibodyStrengthIs(strength);
}

/*
Tissue Tissue1 cloneCellsNew west — Clones all cells in "_tissue" to the "loc"
direction. Equivalent to writing Cell x y z cloneNew "loc" for each cell in
the _tissue. If any single cell throws an exception, you should continue the 
simulation and clone the remaining cells.
*/
void Simulation::cloneCellsNew(Fwk::String _tissue, CellMembrane::Side side) 
{
  cout << "cloning cells in _tissue: " << _tissue << " to: " 
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
      Fwk::String _tissue = *token;
      token++;
      if (*token == "cytotoxicCellNew") {
        token++;
        Cell::Coordinates loc = getCoordinate(token);
        cytotoxicCellNew(_tissue, loc);
      } else if (*token == "helperCellNew") {
        token++;
        Cell::Coordinates loc = getCoordinate(token);
        helperCellNew(_tissue, loc);
      } else if (*token == "infectionStartLocationIs") {
        token++;
        Cell::Coordinates loc = getCoordinate(token);
        token++; token++; token++;
        CellMembrane::Side side = getSide(token++);
        AntibodyStrength strength = 
          AntibodyStrength(lexical_cast<int>(*token++)); 
        infectionStart(_tissue, loc, side, strength);
      } else if (*token == "infectedCellsDel") {
        infectedCellsDel(_tissue);
      } else if (*token == "cloneCellsNew") {
        token++;
        CellMembrane::Side side = getSide(token++);
        cloneCellsNew(_tissue, side);
      } else {
        throw "Malformed command";
      }
    }
  } else if (*token == "Cell") {
    token++;
    Fwk::String _tissue = *token++;
    Cell::Coordinates loc = getCoordinate(token);
    token++; token++; token++;
    if (*token == "membrane") {
      token++;
      CellMembrane::Side side = getSide(token++);
      if (*token == "antibodyStrengthIs") {
        *token++;
        AntibodyStrength strength = 
          AntibodyStrength(lexical_cast<int>(*token++));
        setAntibodyStrength(_tissue, loc, side, strength);
      } else {
        throw "Malformed command";
      }
    } else if (*token == "cloneNew") {
      token++;
      CellMembrane::Side side = getSide(token++);
      cloneNew(_tissue, loc, side);
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

