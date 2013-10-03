#include <fstream>
#include <stdlib.h>
#include "simulation.h"

using namespace std;

/*
  The main takes in one input, the file name with the rules.
  The rules are then executed and the appropriate statistics are printed
  to the console.
*/
int main(int argc, const char* argv[]) {
  Simulation::Ptr sim = Simulation::SimulationNew("main simulation");
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
    // try {
      sim->commandIs(textLine);
    // }
    // catch (...) {
    //   cerr << "Excetion occurred while parseing command: [" << textLine << "]" 
    //     << endl;
    // }
  }
  return 0;
}
