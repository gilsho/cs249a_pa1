#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <queue>
#include "simulation.h"

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
  CellMembrane::Ptr m;
  AntibodyStrength strength;
  if (c->cellType() == Cell::cytotoxicCell()) {
    strength = AntibodyStrength(initialCytotoxicStrength);
    psim->cytotoxicCells_++;
  } else {
    strength = AntibodyStrength(initialHelperStrength);
    psim->helperCells_++;
  }

  m = c->membraneNew(c->name() + " north", CellMembrane::north());
  m->antibodyStrengthIs(strength);
  m = c->membraneNew(c->name() + " south", CellMembrane::south());
  m->antibodyStrengthIs(strength);
  m = c->membraneNew(c->name() + " east", CellMembrane::east());
  m->antibodyStrengthIs(strength);
  m = c->membraneNew(c->name() + " west", CellMembrane::west());
  m->antibodyStrengthIs(strength);
  m = c->membraneNew(c->name() + " up", CellMembrane::up());
  m->antibodyStrengthIs(strength);
  m = c->membraneNew(c->name() + " down", CellMembrane::down());
  m->antibodyStrengthIs(strength);
}

Simulation::Simulation(Fwk::String _name) : Fwk::NamedInterface(_name) {}

/*
Create a new tissue named "name". Quoted names are not accepted and therefore 
names cannot contain any whitespace. NOTE: The tissue name WILL NOT be any of 
the commands, i.e. there will not be a tissue named "tissueNew".
*/
Tissue::Ptr Simulation::tissueNew(Fwk::String _name) {
  Tissue::Ptr t = Tissue::TissueNew(_name);
  TissueReactor::Ptr r = TissueReactor::TissueReactorIs(t.ptr());
  r->notifierIs(t);
  ((TissueReactor *)r.ptr())->psim = this;
  tissues_[_name] = t;
  helperCells_ = 0;
  cytotoxicCells_ = 0;
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
Cell::Ptr Simulation::cellNew(Fwk::String _tissue, Cell::Coordinates loc, 
  Cell::CellType ctype)
{
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);

  // raise exception if cell already exists
  Cell::Ptr c = *(t->cellIter(loc));
  if (c)
    throw "trying to create cell in non-empty location";

  c = Cell::CellNew(loc, t.ptr(), ctype);
  assertValidPtr(c);
  t->cellIs(c);
  
  return c;
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
  U32 attempts = 0;
  S32 difference = 0;
  U32 path = 0;


  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);

  Cell::Ptr rootCell = *(t->cellIter(loc)); 
  assertValidPtr(rootCell);

  attempts++;
  if (!infectionSpreadTo(rootCell, side, strength, difference)) {
    stats(_tissue, attempts, difference, path);
    return;
  }
  
  queue<Cell::Ptr> curRound, nextRound;
  curRound.push(rootCell);

  while (!curRound.empty()) {
    Cell::Ptr c = curRound.front();
    curRound.pop();
    c->healthIs(Cell::infected());
    for (U32 rawSide = 0; rawSide < 6; rawSide++) {
      CellMembrane::Side side = CellMembrane::SideInstance(rawSide);
      Cell::Ptr nbr = neighbor(t, c, side);
      if (nbr && nbr->health() != Cell::infected()) {
        attempts++;
        S32 diff = 0;
        if (infectionSpreadTo(nbr, oppositeSide(side), strength, diff)) {
          difference += diff;
          nextRound.push(nbr);
        }
      }
    } 

    if (curRound.empty()) {
      swap(curRound,nextRound);
      path++;
    }

  }

  stats(_tissue, attempts, difference, path);
}

void Simulation::stats(Fwk::String _tissue, U32 attempts, S32 difference, 
                       U32 path)
{
  // cout << "infected|attempts|diff|cyto|helper|volume|path" << endl;
  cout << infectedCells(_tissue) << " " << attempts << " " 
    << difference << " " << cytotoxicCells_ << " " 
    << helperCells_ << " " << infectionVolume(_tissue) << " " 
    << path << endl;
}

Cell::Ptr Simulation::neighbor(Tissue::Ptr t, Cell::Ptr c, CellMembrane::Side side)
{
  Cell::Coordinates loc = c->location();

  if (side == CellMembrane::north())
    loc.y++;

  else if (side == CellMembrane::south())
    loc.y--;

  else if (side == CellMembrane::east())
    loc.x++;

  else if (side == CellMembrane::west())
    loc.x--;

  else if (side == CellMembrane::up())
    loc.z++;

  else if (side == CellMembrane::down())
    loc.z--;

  else 
    throw "invalid cell side";

  return (*t->cellIter(loc));
}

CellMembrane::Side Simulation::oppositeSide(CellMembrane::Side side) 
{
  if (side == CellMembrane::north())
    return CellMembrane::south();

  if (side == CellMembrane::south())
    return CellMembrane::north();

  if (side == CellMembrane::east())
    return CellMembrane::west();

  if (side == CellMembrane::west())
    return CellMembrane::east();

  if (side == CellMembrane::up()) 
    return CellMembrane::down();

  if (side == CellMembrane::down()) 
    return CellMembrane::up();

  throw "invalid cell side";
}

bool Simulation::infectionSpreadTo(Cell::Ptr c, CellMembrane::Side side, 
                                   AntibodyStrength attack, S32& difference)
{
  CellMembrane::Ptr m = *(c->membraneIterConst(side));
  difference = (S32)attack.value()  - (S32)m->antibodyStrength().value();
  return (attack > m->antibodyStrength());
}

/*
Remove all infected cells from "_tissue".
*/
void Simulation::infectedCellsDel(Fwk::String _tissue)
{
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);
  queue<Fwk::String> cellsQueue;
  
  for (Tissue::CellIterator it = t->cellIter(); it; ++it) {
    Cell::Ptr c = *it;
    assertValidPtr(c);
    if (c->health() == Cell::infected()) {
      cellsQueue.push(c->name());
    }
  }

  while (!cellsQueue.empty()) {
    Fwk::String s = cellsQueue.front();
    cellsQueue.pop();
    t->cellDel(s);
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
  Cell::CellType ctype = c->cellType();
  Cell::Coordinates cloneLoc = coordinateShifted(loc, side); 
  Cell::Ptr clone = cellNew(_tissue, cloneLoc, ctype);
  
  clone->healthIs(c->health());

  for (Cell::MembraneIteratorConst it = c->membraneIterConst(); it; ++it) {
    CellMembrane::Ptr m = *it;
    CellMembrane::Ptr mclone = *(clone->membraneIterConst(m->side()));
    mclone->antibodyStrengthIs(m->antibodyStrength());
  }
}

Cell::Coordinates Simulation::coordinateShifted(Cell::Coordinates loc, 
                                               CellMembrane::Side side)
{
  Cell::Coordinates shiftLoc = loc;
  if (side == CellMembrane::north())
    shiftLoc.y++;
  else if (side == CellMembrane::south())
    shiftLoc.y--;
  else if (side == CellMembrane::east())
   shiftLoc.x++;
  else if (side == CellMembrane::west())
    shiftLoc.x--;
  else if (side == CellMembrane::up())
    shiftLoc.z++;
  else if (side == CellMembrane::down())
    shiftLoc.z--;

  return shiftLoc;
}

void Simulation::antibodyStrengthIs(Fwk::String _tissue, Cell::Coordinates loc,
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
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);
  queue<Cell::Coordinates> cellsQueue;
  for (Tissue::CellIterator it = t->cellIter(); it; ++it) {
    Cell::Ptr c = *it;
    assertValidPtr(c);
    cellsQueue.push(c->location());
  }

  while (!cellsQueue.empty()) {
    Cell::Coordinates loc = cellsQueue.front();
    cellsQueue.pop();
    cloneNew(_tissue, loc, side);
  }
}


Cell::Coordinates Simulation::coordinateIs(tokenizer<>::iterator token)
{
  Cell::Coordinates loc;
  loc.x = lexical_cast<int>(*token++);
  loc.y = lexical_cast<int>(*token++);
  loc.z = lexical_cast<int>(*token);
  return loc;
}

CellMembrane::Side Simulation::sideIs(tokenizer<>::iterator token)
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

U32 Simulation::infectionVolume(Fwk::String _tissue)
{
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);

  Tissue::CellIterator it = t->cellIter();

  // Loop to find first infected cell
  for (++it; it; ++it) {
    Cell::Ptr c = *it;
    assertValidPtr(c);
    if (c->health() == Cell::infected())
      break;
  }
  if (!it)
    return 0;

  Cell::Coordinates minLoc = (*it)->location();
  Cell::Coordinates maxLoc = (*it)->location();

  for (++it; it; ++it) {
    Cell::Ptr c = *it;
    assertValidPtr(c);

    if (c->health() == Cell::healthy())
      continue;

    Cell::Coordinates loc = c->location();

    if (loc.x < minLoc.x)
      minLoc.x = loc.x;
    if (loc.y < minLoc.y)
      minLoc.y = loc.y;
    if (loc.z < minLoc.z)
      minLoc.z = loc.z;

    if (loc.x > maxLoc.x)
      maxLoc.x = loc.x;
    if (loc.y > maxLoc.y)
      maxLoc.y = loc.y;
    if (loc.z > maxLoc.z)
      maxLoc.z = loc.z;
  }

  return (maxLoc.x - minLoc.x + 1) * (maxLoc.y - minLoc.y + 1) *
          (maxLoc.z - minLoc.z + 1);
}

U32 Simulation::infectedCells(Fwk::String _tissue) 
{ 
  Tissue::Ptr t = tissues_[_tissue];
  assertValidPtr(t);
  
  U32 count = 0;
  for (Tissue::CellIterator it = t->cellIter(); it; ++it) {
    Cell::Ptr c = *it;
    assertValidPtr(c);
    if (c->health() == Cell::infected()) {
      count++;
    }
  }
  return count;
}

void Simulation::commandIs(Fwk::String textLine) 
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
        Cell::Coordinates loc = coordinateIs(token);
        cellNew(_tissue, loc, Cell::cytotoxicCell());
      } else if (*token == "helperCellNew") {
        token++;
        Cell::Coordinates loc = coordinateIs(token);
        cellNew(_tissue, loc, Cell::helperCell());
      } else if (*token == "infectionStartLocationIs") {
        token++;
        Cell::Coordinates loc = coordinateIs(token);
        token++; token++; token++;
        CellMembrane::Side side = sideIs(token++);
        AntibodyStrength strength = 
          AntibodyStrength(lexical_cast<int>(*token++)); 
        infectionStart(_tissue, loc, side, strength);
      } else if (*token == "infectedCellsDel") {
        infectedCellsDel(_tissue);
      } else if (*token == "cloneCellsNew") {
        token++;
        CellMembrane::Side side = sideIs(token++);
        cloneCellsNew(_tissue, side);
      } else {
        throw "Malformed command";
      }
    }
  } else if (*token == "Cell") {
    token++;
    Fwk::String _tissue = *token++;
    Cell::Coordinates loc = coordinateIs(token);
    token++; token++; token++;
    if (*token == "membrane") {
      token++;
      CellMembrane::Side side = sideIs(token++);
      if (*token == "antibodyStrengthIs") {
        *token++;
        AntibodyStrength strength = 
          AntibodyStrength(lexical_cast<int>(*token++));
        antibodyStrengthIs(_tissue, loc, side, strength);
      } else {
        throw "Malformed command";
      }
    } else if (*token == "cloneNew") {
      token++;
      CellMembrane::Side side = sideIs(token++);
      cloneNew(_tissue, loc, side);
    } else {
      throw "Malformed command";
    }
  } else {
    throw "Malformed command";
  }
}

