#include "gtest/gtest.h"
#include <fstream>
#include <stdlib.h>
#include <queue>
#include "simulation.h"
#include <iostream>


TEST(Simulation, tissueNew) 
{
  Simulation::Ptr sim = Simulation::SimulationNew("test");

  string str1 = "tissue1";
  string str2 = "tissue2";

  Tissue::Ptr t1 = sim->tissueNew(str1);
  Tissue::Ptr t2 = sim->tissueNew(str2);

  Tissue::Ptr t3 = sim->tissue(str1);
  Tissue::Ptr t4 = sim->tissue(str2);

  ASSERT_TRUE(t1.ptr() != NULL);
  ASSERT_TRUE(t2.ptr() != NULL);
  ASSERT_TRUE(t3.ptr() != NULL);
  ASSERT_TRUE(t4.ptr() != NULL);

  ASSERT_TRUE(t1.ptr() == t3.ptr());
  ASSERT_TRUE(t2.ptr() == t4.ptr());

}

TEST(Simulation, cytotoxicCellNew)
{
  Simulation::Ptr sim = Simulation::SimulationNew("test");
  string tissue = "tissue1";
  Tissue::Ptr t = sim->tissueNew(tissue);
  Cell::Coordinates loc;
  loc.x = 4; 
  loc.y = 5; 
  loc.z = 17;
  sim->cytotoxicCellNew(tissue, loc);
  Tissue::CellIteratorConst it = t->cellIterConst();
  ASSERT_TRUE((*it)->location().x == 4);
  ASSERT_TRUE((*it)->location().y == 5);
  ASSERT_TRUE((*it)->location().z == 17);
  ASSERT_TRUE((*it)->cellType() == Cell::cytotoxicCell_);
}

TEST(Simulation, helperCellNew)
{
  Simulation::Ptr sim = Simulation::SimulationNew("test");
  string tissue = "tissue1";
  Tissue::Ptr t = sim->tissueNew(tissue);
  Cell::Coordinates loc;
  loc.x = 14; 
  loc.y = 31; 
  loc.z = 7;
  sim->helperCellNew(tissue, loc);
  Tissue::CellIteratorConst it = t->cellIterConst();
  ASSERT_TRUE((*it)->location().x == 14);
  ASSERT_TRUE((*it)->location().y == 31);
  ASSERT_TRUE((*it)->location().z == 7);
  ASSERT_TRUE((*it)->cellType() == Cell::helperCell_);
}

TEST(Simulation, setAntibodyStrength)
{
  // Simulation::Ptr sim = Simulation::SimulationNew("test");
  // string tissue = "tissue1";
  // Tissue::Ptr t = sim->tissueNew(tissue);
  // Cell::Coordinates loc;
  // loc.x = 14; 
  // loc.y = 31; 
  // loc.z = 7;
  // sim->helperCellNew(tissue, loc);
  // Cell::CellIteratorConst it = t->cellIterConst();
  // ASSERT_TRUE(*it->location().x == 14);
  // ASSERT_TRUE(*it->location().y == 31);
  // ASSERT_TRUE(*it->location().z == 7);
}

TEST(Simulation, cloneNew)
{
  Simulation::Ptr sim = Simulation::SimulationNew("test");
  string tissue = "tissue1";
  Tissue::Ptr t = sim->tissueNew(tissue);
  Cell::Coordinates loc;
  loc.x = 1; 
  loc.y = 1; 
  loc.z = 1;
  //Cell::CellIterator it = t->cellIter();
  Cell::Ptr c = Cell::CellNew(loc, t.ptr(), Cell::helperCell_);
  t->cellIs(c);
  CellMembrane::Side side = CellMembrane::north_;
  sim->cloneNew(tissue, loc, side);
  
  ASSERT_TRUE(t->cells() == 2);

  Tissue::CellIteratorConst it = t->cellIterConst();
  bool newExists = false;
  bool oldExists = false;
  // for (int i = 0; i < 2; i++, it->advance()) {
    ASSERT_TRUE((*it)->location().x == 1);
    ASSERT_TRUE((*it)->location().y == 1);
    if ((*it)->location().z == 1)
      oldExists = true;
    if ((*it)->location().z == 2) 
      newExists = true;
  // }
  ASSERT_TRUE(newExists);
  ASSERT_TRUE(oldExists);
}


TEST(Simulation, infectionStart)
{
  ASSERT_TRUE(true);
}

TEST(Simulation, infectedCellsDel)
{
  ASSERT_TRUE(true);
}