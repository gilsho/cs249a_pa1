#include "gtest/gtest.h"
#include <fstream>
#include <stdlib.h>
#include <queue>
#include <iostream>
#include "simulation.h"


bool membraneStrength(Tissue::Ptr t, Cell::Coordinates loc,
    CellMembrane::Side side, AntibodyStrength strength)
{
  Cell::Ptr c = (*t->cellIter(loc));
  if (c.ptr() == NULL)
    return false;
  CellMembrane::Ptr m = (*c->membraneIterConst(side));
  if (m.ptr() == NULL)
    return false;
  return (m->antibodyStrength() == strength);
}

bool cellExists(Tissue::Ptr t, Cell::Coordinates loc, 
                Cell::CellType ctype, Cell::HealthId health)
{
  Tissue::CellIterator it = t->cellIter(loc);
  if (it.ptr() == NULL)
    return false;
  if (it->cellType() != ctype)
    return false;
  if (it->health() != health)
    return false;

  return true;
}

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
  Cell::Coordinates loc = { 4, 5, 17};
  loc.x = 4; 
  loc.y = 5; 
  loc.z = 17;
  sim->cytotoxicCellNew(tissue, loc);
  ASSERT_TRUE(cellExists(t, loc, Cell::cytotoxicCell(), 
    Cell::healthy()));

}

TEST(Simulation, helperCellNew)
{
  Simulation::Ptr sim = Simulation::SimulationNew("test");
  string tissue = "tissue1";
  Tissue::Ptr t = sim->tissueNew(tissue);
  Cell::Coordinates loc = {12, -13, -21};
  sim->helperCellNew(tissue, loc);
  ASSERT_TRUE(cellExists(t, loc, Cell::helperCell(), 
    Cell::healthy()));
}

TEST(Simulation, setAntibodyStrength)
{
  Simulation::Ptr sim = Simulation::SimulationNew("test");
  string tissue = "tissue1";
  Tissue::Ptr t = sim->tissueNew(tissue);
  Cell::Coordinates loc = {2, 2, 2};
  CellMembrane::Side side = CellMembrane::up();
  AntibodyStrength strength = AntibodyStrength(20);
  sim->helperCellNew(tissue, loc);
  sim->setAntibodyStrength(tissue, loc, side, strength);
  ASSERT_TRUE(membraneStrength(t, loc, side, strength));
}



TEST(Simulation, cloneNew)
{
  Simulation::Ptr sim = Simulation::SimulationNew("test");
  string tissue = "tissue1";
  Tissue::Ptr t = sim->tissueNew(tissue);
  Cell::Coordinates loc1 = {1, 1, 1};
  Cell::Coordinates loc2 = {1, 1, 2};
  Cell::Coordinates loc3 = {1, 1, 0};
  Cell::Coordinates loc4 = {0, 1, 1};
  Cell::Coordinates loc5 = {2, 1, 1};
  Cell::Coordinates loc6 = {1, 0, 1};
  Cell::Coordinates loc7 = {1, 2, 1};

  sim->helperCellNew(tissue, loc1);
  sim->cloneNew(tissue, loc1, CellMembrane::up());
  
  ASSERT_TRUE(t->cells() == 2);
  ASSERT_TRUE(cellExists(t, loc1, Cell::helperCell(), Cell::healthy()));
  ASSERT_TRUE(cellExists(t, loc2, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(tissue, loc1, CellMembrane::down());
  ASSERT_TRUE(t->cells() == 3);
  ASSERT_TRUE(cellExists(t, loc3, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(tissue, loc1, CellMembrane::west());
  ASSERT_TRUE(t->cells() == 4);
  ASSERT_TRUE(cellExists(t, loc4, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(tissue, loc1, CellMembrane::east());
  ASSERT_TRUE(t->cells() == 5);
  ASSERT_TRUE(cellExists(t, loc5, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(tissue, loc1, CellMembrane::south());
  ASSERT_TRUE(t->cells() == 6);
  ASSERT_TRUE(cellExists(t, loc6, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(tissue, loc1, CellMembrane::north());
  ASSERT_TRUE(t->cells() == 7);
  ASSERT_TRUE(cellExists(t, loc7, Cell::helperCell(), Cell::healthy()));

}


TEST(Simulation, infectionStart_FirstCellFail)
{
  ASSERT_TRUE(false);  
}

TEST(Simulation, infectedCellsDel)
{
  Simulation::Ptr sim = Simulation::SimulationNew("test");
  string tissue = "tissue1";
  Tissue::Ptr t = sim->tissueNew(tissue);
  Cell::Coordinates loc1 = {1, 2, 3};
  Cell::Coordinates loc2 = {4, 5, 6};
  sim->helperCellNew(tissue, loc1);
  sim->helperCellNew(tissue, loc2);
  Cell::Ptr c1;
  Cell::Ptr c2;

  sim->infectedCellsDel(tissue);
  c1 = *(t->cellIter(loc1));
  c2 = *(t->cellIter(loc2));
  ASSERT_TRUE(c1.ptr() != NULL);
  ASSERT_TRUE(c2.ptr() != NULL);

  c1->healthIs(Cell::infected());
  sim->infectedCellsDel(tissue);
  c1 = *(t->cellIter(loc1));
  c2 = *(t->cellIter(loc2));
  ASSERT_TRUE(c1.ptr() == NULL);
  ASSERT_TRUE(c2.ptr() != NULL);

}
