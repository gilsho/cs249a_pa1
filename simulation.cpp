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

string CoordToStr(Cell::Coordinates loc) 
{
  return "(" + lexical_cast<string>(loc.x) + ", " +
               lexical_cast<string>(loc.y) + ", " +
               lexical_cast<string>(loc.z) + ")";
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

void Simulation::TissueReactor::onCellDel(Cell::Ptr c)
{
  if (c->cellType() == Cell::cytotoxicCell()) {
    psim->cytotoxicCells_--;
  } else {
    psim->helperCells_--;
  }
}

// create new simulation object, wrapping a tissue_
Simulation::Simulation(Fwk::String _name) : 
  Fwk::NamedInterface("sim" + _name) 
{
  tissue_ = Tissue::TissueNew(_name);
  TissueReactor::Ptr r = TissueReactor::TissueReactorIs(tissue_.ptr());
  r->notifierIs(tissue_);
  ((TissueReactor *)r.ptr())->psim = this;
  helperCells_ = 0;
  cytotoxicCells_ = 0;
}

Tissue::Ptr Simulation::tissue()
{
  return tissue_;
}

/*
Create a new healthy CytotoxicCell at location "loc". In case an 
exception is thrown (ex. a cell already exists at location), your program 
should catch the exception and, if you wish, use some logging mechanism to 
report it. In any case, your simulation should continue running as if this 
entry did not exist. A new CytotoxicCell should have antibody strength of 100 
on all its membranes.
*/
Cell::Ptr Simulation::cellNew(Cell::Coordinates loc, 
  Cell::CellType ctype)
{
  // raise exception if cell already exists
  Cell::Ptr c = *(tissue_->cellIter(loc));
  if (c)
    throw "trying to create cell in non-empty location";

  c = Cell::CellNew(loc, tissue_.ptr(), ctype);
  tissue_->cellIs(c);
  
  return c;
}

void printTissue(Tissue::Ptr t)
{
  for (Tissue::CellIterator it = t->cellIter(); it; ++it) 
  {
    Cell::Ptr c = *it;
    cout << c->name() << endl;
  }
}

int printQueue(queue<Cell::Ptr> q)
{
  queue<Cell::Ptr> p;
  int count = 0;
  while (!q.empty()) {
    Cell::Ptr c = q.front();
    q.pop();
    cout << c->name() << endl;
    p.push(c);
    count++;
  }

  while(!p.empty()) {
    Cell::Ptr c = p.front(); 
    p.pop();
    q.push(c);
  }

  return count;
}

/*
Starts an infection of strength 99 at cell at "loc" entering from the "loc"
membrane. You should proceed to the next command only when no more cells can be 
infected. At the end of the infection round, you should print statistics to 
standard out as described here.
*/
void Simulation::infectionStart(Cell::Coordinates loc, 
                    CellMembrane::Side side, AntibodyStrength strength)
{

  U32 attempts = 0;
  S32 difference = 0;
  U32 path = 0;

  Cell::Ptr rootCell = *(tissue_->cellIter(loc)); 
  if (!rootCell) {
    stats(attempts, difference, path);
    return;
  }

  if (!infectionSpreadTo(rootCell, side, strength, difference, attempts)) {
    stats(attempts, difference, path);
    return;
  }
  
  queue<Cell::Ptr> curRound, nextRound;
  curRound.push(rootCell);

  while (!curRound.empty()) {
    Cell::Ptr c = curRound.front();
    curRound.pop();
      
    Cell::Ptr nbr;
    CellMembrane::Side side;

    side = CellMembrane::north(); 
    nbr = neighbor(c, side); 
    if (infectionSpreadTo(nbr, oppositeSide(side), strength, difference, attempts))
        nextRound.push(nbr);

    side = CellMembrane::east(); 
    nbr = neighbor(c, side); 
    if (infectionSpreadTo(nbr, oppositeSide(side), strength, difference, attempts))
        nextRound.push(nbr);

    side = CellMembrane::south(); 
    nbr = neighbor(c, side);
    if (infectionSpreadTo(nbr, oppositeSide(side), strength, difference, attempts))
        nextRound.push(nbr);

    side = CellMembrane::west(); 
    nbr = neighbor(c, side);
    if (infectionSpreadTo(nbr, oppositeSide(side), strength, difference, attempts))
        nextRound.push(nbr);

    side = CellMembrane::up(); 
    nbr = neighbor(c, side);
    if (infectionSpreadTo(nbr, oppositeSide(side), strength, difference, attempts))
        nextRound.push(nbr);

    side = CellMembrane::down(); 
    nbr = neighbor(c, side);
    if (infectionSpreadTo(nbr, oppositeSide(side), strength, difference, attempts))
        nextRound.push(nbr);

    if (curRound.empty()) {
      swap(curRound,nextRound);
      path++;
    }

  }

  stats(attempts, difference, path);
}

// spreads an infection to cell from specific side. updates statistics
// as well
bool Simulation::infectionSpreadTo(Cell::Ptr c, CellMembrane::Side side, 
                                   AntibodyStrength attack, 
                                   S32& difference,
                                   U32& attempts)
{
  if (!c || c->health() == Cell::infected())
    return false;

  attempts++;
  CellMembrane::Ptr m = *(c->membraneIterConst(side));
  difference += (S32)attack.value()  - (S32)m->antibodyStrength().value();
  if (attack > m->antibodyStrength()) {
    c->healthIs(Cell::infected());
    return true;
  }
  return false;
}

// print out statistics about an infection round
void Simulation::stats(U32 attempts, S32 difference, 
                       U32 path)
{
  cout << infectedCells() << " " << attempts << " " 
    << difference << " " << cytotoxicCells_ << " " 
    << helperCells_ << " " << infectionVolume() << " " 
    << path << endl;
}

//returns the neighbor of a cell in a particular direction
Cell::Ptr Simulation::neighbor(Cell::Ptr c, CellMembrane::Side side)
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
    return NULL;


  return (*tissue_->cellIter(loc));
}

// returns the inverted side. north becomes south, east becomes west, 
// and so forth
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

  return side;
}


// Remove all infected cells from "_tissue_".
void Simulation::infectedCellsDel()
{
  queue<Fwk::String> cellsQueue;
  
  for (Tissue::CellIterator it = tissue_->cellIter(); it; ++it) {
    Cell::Ptr c = *it;
    if (c->health() == Cell::infected()) {
      cellsQueue.push(c->name());
    }
  }

  while (!cellsQueue.empty()) {
    Fwk::String s = cellsQueue.front();
    cellsQueue.pop();
    tissue_->cellDel(s);
  }
}

/*
Clones cell at location "loc" and places the new cell "side" of "loc"
Like the other cell creation commands, the simulation should 
continue running despite any exception that may be thrown.
*/

void Simulation::cloneNew(Cell::Coordinates loc, 
              CellMembrane::Side side)
{

  Cell::Ptr c = *(tissue_->cellIter(loc));
  Cell::CellType ctype = c->cellType();
  Cell::Coordinates cloneLoc = coordinateShifted(loc, side); 
  // cout << CoordToStr(cloneLoc) << endl;
  Cell::Ptr clone = cellNew(cloneLoc, ctype);
  
  clone->healthIs(c->health());

  for (Cell::MembraneIteratorConst it = c->membraneIterConst(); it; ++it) {
    CellMembrane::Ptr m = *it;
    CellMembrane::Ptr mclone = *(clone->membraneIterConst(m->side()));
    mclone->antibodyStrengthIs(m->antibodyStrength());
  }
}

// returns a shifted version of given coordinate in the specied side
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

// sets the antibody stength of a cell in the simulation tissue_
void Simulation::antibodyStrengthIs(Cell::Coordinates loc,
                         CellMembrane::Side side, AntibodyStrength strength)
{
  Cell::Ptr c = *(tissue_->cellIter(loc));
  if (!c) 
    cellNew(loc, DEFAULT_CELL_TYPE);
  CellMembrane::Ptr m = *(c->membraneIterConst(side));
  assertValidPtr(m);
  m->antibodyStrengthIs(strength);
}

/*
Tissue Tissue1 cloneCellsNew west — Clones all cells to the "loc"
direction. Equivalent to writing Cell x y z cloneNew "loc" for each cell in
the tissue_. If any single cell throws an exception, you should continue the 
simulation and clone the remaining cells.
*/
void Simulation::cloneCellsNew(CellMembrane::Side side) 
{
  queue<Cell::Coordinates> cellsQueue;
  for (Tissue::CellIterator it = tissue_->cellIter(); it; ++it) {
    Cell::Ptr c = *it;
    if (c)
      cellsQueue.push(c->location());
  }

  while (!cellsQueue.empty()) {
    Cell::Coordinates loc = cellsQueue.front();
    cellsQueue.pop();
    try {
      cloneNew(loc, side);
    } catch(...) {
      //do nothing
    }
  }
}


// Calculates the bounding box around infected cells
U32 Simulation::infectionVolume()
{
  Tissue::CellIterator it = tissue_->cellIter();

  // Loop to find first infected cell
  for (; it; ++it) {
    Cell::Ptr c = *it;
    if (c->health() == Cell::infected())
      break;
  }
  if (!it)
    return 0;

  Cell::Coordinates minLoc = (*it)->location();
  Cell::Coordinates maxLoc = (*it)->location();

  for (++it; it; ++it) {
    Cell::Ptr c = *it;

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

// returns the total number of infected cells in tissue_
U32 Simulation::infectedCells() 
{ 
  U32 count = 0;
  for (Tissue::CellIterator it = tissue_->cellIter(); it; ++it) {
    Cell::Ptr c = *it;
    if (c->health() == Cell::infected()) {
      count++;
    }
  }
  return count;
}
