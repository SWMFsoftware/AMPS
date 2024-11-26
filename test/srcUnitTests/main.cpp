//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf



//the particle class
#include "pic.h"
#include "constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

#include <gtest/gtest.h>


int add(int a, int b) {
    return a + b;
}

TEST(MathOperationsTest, SubtractNegativeNumbers) {
    EXPECT_EQ(add(-5, -3), -8);
}



//$Id$

//TEST(AddFunction, HandlesNegativeNumbers) {
//    EXPECT_EQ(add(-2, -3), -5);  // Check if -2 + -3 equals -5
//}


void amps_init();
void amps_time_step();

void pbuffer_test_for_linker(); 
void collisions_test_for_linker();


int main(int argc,char **argv) {

  pbuffer_test_for_linker();
  collisions_test_for_linker();
  
  clock_t runtime =-clock();

PIC::Debugger::cGenericTimer t;
  
t.Start("main",__LINE__);
  amps_init();

t.SwitchTimeSegment(__LINE__,"first switch"); 


    ::testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
