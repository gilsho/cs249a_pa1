#include <fstream>
#include <stdlib.h>
#include "simulation.h"
#include "Tissue.h"

using namespace std;
using namespace boost;

/*
  The main takes in one input, the file name with the rules.
  The rules are then executed and the appropriate statistics are printed
  to the console.
*/

Cell::Coordinates coordinateIs(
  tokenizer<char_separator<char> >::iterator token)
{
  Cell::Coordinates loc;
  loc.x = lexical_cast<S32>(*token++);
  loc.y = lexical_cast<S32>(*token++);
  loc.z = lexical_cast<S32>(*token);
  return loc;
}

CellMembrane::Side sideIs(
  tokenizer<char_separator<char> >::iterator token)
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

void commandIs(Fwk::String textLine, map<Fwk::String, Simulation::Ptr>& sims) 
{
  if (textLine == "" || textLine[0] == '#')
    return;

  char_separator<char> sep(" ");
  tokenizer<char_separator<char> > tokenizedLine(textLine, sep);
  tokenizer<char_separator<char> >::iterator token=tokenizedLine.begin(); 

  if (*token == "Tissue") {
    token++;
    if (*token == "tissueNew") {
      token++;
      Simulation::Ptr curSim = Simulation::SimulationNew(*token);
      sims[*token] = curSim;
    } else {
      Simulation::Ptr curSim = sims[*token];
      token++;
      if (*token == "cytotoxicCellNew") {
        token++;
        Cell::Coordinates loc = coordinateIs(token);
        curSim->cellNew(loc, Cell::cytotoxicCell());
      } else if (*token == "helperCellNew") {
        token++;
        Cell::Coordinates loc = coordinateIs(token);
        curSim->cellNew(loc, Cell::helperCell());
      } else if (*token == "infectionStartLocationIs") {
        token++;
        Cell::Coordinates loc = coordinateIs(token);
        token++; token++; token++;
        CellMembrane::Side side = sideIs(token++);
        AntibodyStrength strength = 
          AntibodyStrength(lexical_cast<int>(*token++)); 
        curSim->infectionStart(loc, side, strength);
      } else if (*token == "infectedCellsDel") {
        curSim->infectedCellsDel();
      } else if (*token == "cloneCellsNew") {
        token++;
        CellMembrane::Side side = sideIs(token++);
        curSim->cloneCellsNew(side);
      } else {
        throw "Malformed command";
      }
    }
  } else if (*token == "Cell") {
    token++;
    Simulation::Ptr curSim = sims[*token];
    token++;
    Cell::Coordinates loc = coordinateIs(token);
    token++; token++; token++;
    if (*token == "membrane") {
      token++;
      CellMembrane::Side side = sideIs(token++);
      if (*token == "antibodyStrengthIs") {
        *token++;
        AntibodyStrength strength = 
          AntibodyStrength(lexical_cast<int>(*token++));
        curSim->antibodyStrengthIs(loc, side, strength);
      } else {
        throw "Malformed command";
      }
    } else if (*token == "cloneNew") {
      token++;
      CellMembrane::Side side = sideIs(token++);
      curSim->cloneNew(loc, side);
    } else {
      throw "Malformed command";
    }
  } else {
    throw "Malformed command";
  }
}



int main(int argc, const char* argv[]) {
  map<Fwk::String, Simulation::Ptr> sims;
  ifstream infile(argv[1]);
  if(infile.fail()){
    //File error. Halt program.
    cout << "error reading file" << endl;
    return 1;
  }

  //read data in, parse it, excute commands.
  Fwk::String textLine;
  while(!infile.eof()){
    getline(infile, textLine);
    try {
      commandIs(textLine, sims);
    }
    catch (...) {
      cerr << "Excetion occurred while parseing command: [" << textLine << "]" 
        << endl;
    }
  }
  return 0;
}
