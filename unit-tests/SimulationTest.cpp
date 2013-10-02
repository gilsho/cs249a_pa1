#include "gtest/gtest.h"
#include <fstream>
#include <stdlib.h>
#include <queue>
#include "simulation.h"
#include <iostream>


TEST(Simulation, helloWorld) {
  Simulation::Ptr sim = Simulation::SimulationNew("test simulation");
  ASSERT_TRUE(sim->test() == 2);
}