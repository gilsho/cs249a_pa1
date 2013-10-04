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
  string str1 = "tissue1";
  string str2 = "tissue2";

  Simulation::Ptr sim1 = Simulation::SimulationNew(str1);
  Simulation::Ptr sim2 = Simulation::SimulationNew(str2);

  Tissue::Ptr t1 = sim1->tissue();
  Tissue::Ptr t2 = sim2->tissue();

  ASSERT_TRUE(t1.ptr() != NULL);
  ASSERT_TRUE(t2.ptr() != NULL);

  ASSERT_TRUE(t1->name() == str1);
  ASSERT_TRUE(t2->name() == str2);  

}

TEST(Simulation, cytotoxicCellNew)
{
  string tissue = "tissue1";
  Simulation::Ptr sim = Simulation::SimulationNew(tissue);
  Cell::Coordinates loc = { 4, 5, 17};
  loc.x = 4; 
  loc.y = 5; 
  loc.z = 17;
  sim->cellNew(loc, Cell::cytotoxicCell());
  ASSERT_TRUE(cellExists(sim->tissue(), loc, Cell::cytotoxicCell(), 
    Cell::healthy()));

}

TEST(Simulation, helperCellNew)
{
  string tissue = "tissue1";
  Simulation::Ptr sim = Simulation::SimulationNew(tissue);
  Cell::Coordinates loc = {12, -13, -21};
  sim->cellNew(loc, Cell::helperCell());
  ASSERT_TRUE(cellExists(sim->tissue(), loc, Cell::helperCell(), 
    Cell::healthy()));
}

TEST(Simulation, antibodyStrength)
{
  string tissue = "tissue1";
  Simulation::Ptr sim = Simulation::SimulationNew(tissue);
  Cell::Coordinates loc = {2, 2, 2};
  CellMembrane::Side side = CellMembrane::up();
  AntibodyStrength strength = AntibodyStrength(20);
  sim->cellNew(loc, Cell::helperCell());
  sim->antibodyStrengthIs(loc, side, strength);
  ASSERT_TRUE(membraneStrength(sim->tissue(), loc, side, strength));
}

TEST(Simulation, cloneNew)
{
  string tissue = "tissue1";
  Simulation::Ptr sim = Simulation::SimulationNew(tissue);
  Tissue::Ptr t = sim->tissue();
  Cell::Coordinates loc1 = {1, 1, 1};
  Cell::Coordinates loc2 = {1, 1, 2};
  Cell::Coordinates loc3 = {1, 1, 0};
  Cell::Coordinates loc4 = {0, 1, 1};
  Cell::Coordinates loc5 = {2, 1, 1};
  Cell::Coordinates loc6 = {1, 0, 1};
  Cell::Coordinates loc7 = {1, 2, 1};

  sim->cellNew(loc1, Cell::helperCell());
  sim->cloneNew(loc1, CellMembrane::up());
  
  ASSERT_TRUE(t->cells() == 2);
  ASSERT_TRUE(cellExists(t, loc1, Cell::helperCell(), Cell::healthy()));
  ASSERT_TRUE(cellExists(t, loc2, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(loc1, CellMembrane::down());
  ASSERT_TRUE(t->cells() == 3);
  ASSERT_TRUE(cellExists(t, loc3, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(loc1, CellMembrane::west());
  ASSERT_TRUE(t->cells() == 4);
  ASSERT_TRUE(cellExists(t, loc4, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(loc1, CellMembrane::east());
  ASSERT_TRUE(t->cells() == 5);
  ASSERT_TRUE(cellExists(t, loc5, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(loc1, CellMembrane::south());
  ASSERT_TRUE(t->cells() == 6);
  ASSERT_TRUE(cellExists(t, loc6, Cell::helperCell(), Cell::healthy()));

  sim->cloneNew(loc1, CellMembrane::north());
  ASSERT_TRUE(t->cells() == 7);
  ASSERT_TRUE(cellExists(t, loc7, Cell::helperCell(), Cell::healthy()));

}


TEST(Simulation, cloneCellsNew)
{
  string tissue = "tissue1";
  Simulation::Ptr sim = Simulation::SimulationNew(tissue);
  Tissue::Ptr t = sim->tissue();

  for (int i = 0; i < 10; i++) {
    Cell::Coordinates loc = {0, 0, i};
    sim->cellNew(loc, Cell::helperCell());
  }
  Cell::Ptr c = *(t->cellIter());
  CellMembrane::Side side = CellMembrane::north();
  Cell::HealthId health = c->health();
  Cell::CellType ctype = c->cellType();

  sim->cloneCellsNew(side);

  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 2; j++) {
      Cell::Coordinates loc = {0, 0, i};
      ASSERT_TRUE(cellExists(t, loc, ctype, health));
    }
  }

}


TEST(Simulation, infectionStart1)
{
  string tissue = "tissue1";
  Simulation::Ptr sim = Simulation::SimulationNew(tissue);
  Tissue::Ptr t = sim->tissue();

  Cell::Coordinates loc = {0, 0, 0};
  sim->cellNew(loc, Cell::cytotoxicCell());
  Cell::Ptr c = *(t->cellIter(loc));
  AntibodyStrength strength = AntibodyStrength(20);
  CellMembrane::Side side = CellMembrane::up();
  CellMembrane::Ptr m = *(c->membraneIterConst(side));
  m->antibodyStrengthIs(strength);

  sim->infectionStart(loc, side, strength);
  c = *(t->cellIter(loc));
  ASSERT_TRUE(c->health() == Cell::healthy());

  strength = AntibodyStrength(15);
  sim->infectionStart(loc, side, strength);
  ASSERT_TRUE(c->health() == Cell::healthy());

  strength = AntibodyStrength(25);
  sim->infectionStart(loc, side, strength);
  ASSERT_TRUE(c->health() == Cell::infected());
}

TEST(Simulation, infectionStart2)
{
  string tissue = "tissue1";
  Simulation::Ptr sim = Simulation::SimulationNew(tissue);
  Tissue::Ptr t = sim->tissue();

  for (int i = 0; i < 10; i++) {
    Cell::Coordinates loc = {0, 0, i};
    sim->cellNew(loc, Cell::helperCell());
  }
  sim->infectedCellsDel();

  Cell::Coordinates loc = {0, 0, 0};
  CellMembrane::Side side = CellMembrane::down();
  AntibodyStrength strength = AntibodyStrength(20);
  sim->infectionStart(loc, side, strength);

  for (int i = 0; i < 10; i++) {
    Cell::Coordinates loc = {0, 0, i};
    Cell::Ptr c = *(t->cellIter(loc));
    ASSERT_TRUE(c->health() == Cell::infected());
  }

  sim->infectedCellsDel();
  ASSERT_TRUE(t->cells() == 0);
}

TEST(Simulation, infectionStart3)
{
  string tissue = "tissue1";
  Simulation::Ptr sim = Simulation::SimulationNew(tissue);
  Tissue::Ptr t = sim->tissue();

  for (int i = 0; i < 10; i++) {
    Cell::Coordinates loc = {0, 0, i};
    if (i == 5) {
      sim->cellNew(loc, Cell::cytotoxicCell());
    } else {
      sim->cellNew(loc, Cell::helperCell());
    }
  }


  Cell::Coordinates loc = {0, 0, 0};
  CellMembrane::Side side = CellMembrane::down();
  AntibodyStrength strength = AntibodyStrength(20);

  sim->infectionStart(loc, side, strength);
  
  for (int i = 0; i < 10; i++) {
    Cell::Coordinates loc = {0, 0, i};
    Cell::Ptr c = *(t->cellIter(loc));
    if (i < 5) 
      ASSERT_TRUE(c->health() == Cell::infected());
    else 
      ASSERT_TRUE(c->health() == Cell::healthy());
  }

  sim->infectedCellsDel();
  ASSERT_TRUE(t->cells() == 5);

}

TEST(Simulation, infectedCellsDel)
{
  string tissue = "tissue1";
  Simulation::Ptr sim = Simulation::SimulationNew(tissue);
  Tissue::Ptr t = sim->tissue();

  Cell::Coordinates loc1 = {1, 2, 3};
  Cell::Coordinates loc2 = {4, 5, 6};
  sim->cellNew(loc1, Cell::helperCell());
  sim->cellNew(loc2, Cell::helperCell());
  Cell::Ptr c1;
  Cell::Ptr c2;

  sim->infectedCellsDel();
  c1 = *(t->cellIter(loc1));
  c2 = *(t->cellIter(loc2));
  ASSERT_TRUE(c1.ptr() != NULL);
  ASSERT_TRUE(c2.ptr() != NULL);

  c1->healthIs(Cell::infected());
  sim->infectedCellsDel();
  c1 = *(t->cellIter(loc1));
  c2 = *(t->cellIter(loc2));
  ASSERT_TRUE(c1.ptr() == NULL);
  ASSERT_TRUE(c2.ptr() != NULL);

}
