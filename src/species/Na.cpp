//=======================================================================
//$Id$
//=======================================================================
//the file contains functions descpibing specific phhysical functions for sodium

#ifndef _Na_PAHYSICAL_PARAMETERS_
#define _Na_PAHYSICAL_PARAMETERS_

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

#include "constants.h"
#include "Na.h"


double SodiumGfactor__D1D2__Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance) {

  static const int SodiumGfactor__D1D2_TableLength__Combi_1997_icarus=101;
  static const double SodiumGfactor__D1D2_HeliocenticVelocity_Min__Combi_1997_icarus=-50.0E3;
  static const double SodiumGfactor__D1D2_HeliocenticVelocity_Max__Combi_1997_icarus=50.0E3;
  static const double SodiumGfactor__D1D2_dHeliocenticVelocity__Combi_1997_icarus=1.0E3;

  static const double SodiumGfactor__D1D2_Table__Combi_1997_icarus[SodiumGfactor__D1D2_TableLength__Combi_1997_icarus]={
      14.55069, 14.55069, 14.53236, 14.53236, 14.53236, 14.53236, 14.53236, 14.53236, 14.53236, 14.53236,
      14.33072, 14.20241, 14.12909, 14.05576, 13.98244, 13.85413, 13.74414, 13.63416, 13.50584, 13.34087,
      13.10257, 12.79095, 12.47933, 12.31435, 12.16770, 12.02106, 11.91108, 11.81942, 11.72777, 11.54446,
      11.30616, 11.03120, 10.75624, 10.40796, 10.09634, 9.74806, 9.38144, 8.90485, 8.50157, 8.09830,
      7.69502, 7.16343, 6.48520, 5.77030, 4.57881, 3.24067, 2.25081, 1.33428, 1.04099, 0.87601,
      0.80269, 0.85768, 1.05932, 1.46259, 2.23248, 3.38731, 4.28552, 5.45868, 6.26523, 6.90680,
      7.43839, 7.87833, 8.28160, 8.66655, 9.03316, 9.32645, 9.47310, 9.67473, 10.16966, 10.57294,
      11.03120, 11.30616, 11.61779, 11.85608, 12.09438, 12.29602, 12.46099, 12.68096, 12.88260, 12.97426,
      13.13923, 13.34087, 13.43252, 13.59750, 13.70748, 13.79914, 13.87245, 14.00077, 14.09243, 14.12909,
      14.16575, 14.25740, 14.49570, 14.49570, 14.49570, 14.49570, 14.49570, 14.49570, 14.49570, 14.49570,
      14.49570};

  double res=0.0;

  if (HeliocentricVelocity<=SodiumGfactor__D1D2_HeliocenticVelocity_Min__Combi_1997_icarus) res=SodiumGfactor__D1D2_Table__Combi_1997_icarus[0];
  else if (HeliocentricVelocity>=SodiumGfactor__D1D2_HeliocenticVelocity_Max__Combi_1997_icarus) res=SodiumGfactor__D1D2_Table__Combi_1997_icarus[SodiumGfactor__D1D2_TableLength__Combi_1997_icarus-1];
  else {
    int level;
    double x;

    x=(HeliocentricVelocity-SodiumGfactor__D1D2_HeliocenticVelocity_Min__Combi_1997_icarus)/SodiumGfactor__D1D2_dHeliocenticVelocity__Combi_1997_icarus;
    level=(int)x;
    x-=level;

    res=SodiumGfactor__D1D2_Table__Combi_1997_icarus[level]+x*(SodiumGfactor__D1D2_Table__Combi_1997_icarus[level+1]-SodiumGfactor__D1D2_Table__Combi_1997_icarus[level]);
  }

  //scale g-factor to the required heliocentric distance
  res*=pow(149598000000.0/HeliocentricDistance,2);

  return res;

}

//the solar dariation pressure as a function of heliocentric distrance and velocty
//taken from Combi-1997-icarus.pdf

//input parameters: heliocentric velociy (m/s) and heliocentric distance (m)
//return the radiation pressure (m/s^2)
//the positive velocity is in the direction out from the sun
//the acceleration is in teh direction out of the sun

double SodiumRadiationPressureAcceleration__Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance) {

  static const int sodiumRadiationPressure_TableLength__Combi_1997_icarus=1201;
  static const double sodiumRadiationPressure_HeliocenticVelocity_Min__Combi_1997_icarus=-60.0E3;
  static const double sodiumRadiationPressure_HeliocenticVelocity_Max__Combi_1997_icarus=60.0E3;
  static const double sodiumRadiationPressure_dHeliocenticVelocity__Combi_1997_icarus=0.1E3;

  static const double sodiumRadiationPressure_Table__Combi_1997_icarus[sodiumRadiationPressure_TableLength__Combi_1997_icarus]={
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,     46.4855,
    46.4064,     46.4064,     46.3362,     46.3238,     46.3238,     46.2430,     46.2430,     46.1797,     46.1656,     46.1656,     46.0883,     46.0883,     46.0391,     46.0233,     46.0233,
    45.9617,     45.9617,     45.9231,     45.9090,     45.9090,     45.8633,     45.8633,     45.8317,     45.8211,     45.8211,     45.7754,     45.7754,     45.7754,     45.7297,     45.7297,
    45.6875,     45.6717,     45.6717,     45.6031,     45.6031,     45.5469,     45.5240,     45.5240,     45.4326,     45.4326,     45.3658,     45.3377,     45.3377,     45.2374,     45.2374,
    45.1707,     45.1425,     45.1425,     45.0476,     45.0476,     44.9843,     44.9579,     44.9579,     44.8718,     44.8718,     44.8718,     44.7891,     44.7891,     44.7030,     44.7030,
    44.7030,     44.6186,     44.6186,     44.5554,     44.5272,     44.5272,     44.4375,     44.4375,     44.3672,     44.3391,     44.3391,     44.2406,     44.2406,     44.1668,     44.1369,
    44.1369,     44.0279,     44.0279,     43.9471,     43.9119,     43.9119,     43.7870,     43.7870,     43.7870,     43.6499,     43.6499,     43.5058,     43.5058,     43.5058,     43.3475,
    43.3475,     43.2316,     43.1805,     43.1805,     42.9941,     42.9941,     42.8606,     42.7972,     42.7972,     42.5827,     42.5827,     42.4352,     42.3507,     42.3507,     42.1062,
    42.1062,     41.9516,     41.8477,     41.8477,     41.5822,     41.5822,     41.5822,     41.3220,     41.3220,     41.0688,     41.0688,     41.0688,     40.8385,     40.8385,     40.6310,
    40.6310,     40.6310,     40.4553,     40.4553,     40.3147,     40.3077,     40.3077,     40.1865,     40.1865,     40.0635,     40.0793,     40.0793,     39.9757,     39.9757,     39.8703,
    39.8738,     39.8738,     39.7666,     39.7666,     39.7666,     39.6453,     39.6453,     39.5204,     39.5204,     39.5204,     39.3903,     39.3903,     39.2602,     39.2602,     39.2602,
    39.1423,     39.1423,     39.0861,     39.0333,     39.0333,     38.9419,     38.9419,     38.8962,     38.8645,     38.8645,     38.7977,     38.7977,     38.7555,     38.7432,     38.7432,
    38.6869,     38.6869,     38.6869,     38.6289,     38.6289,     38.5657,     38.5657,     38.5657,     38.4884,     38.4884,     38.3987,     38.3987,     38.3987,     38.2950,     38.2950,
    38.2107,     38.1825,     38.1825,     38.0524,     38.0524,     37.9470,     37.9083,     37.9083,     37.7571,     37.7571,     37.6376,     37.5918,     37.5918,     37.4196,     37.4196,
    37.4196,     37.2332,     37.2332,     37.0416,     37.0416,     37.0416,     36.8377,     36.8377,     36.6250,     36.6250,     36.6250,     36.4035,     36.4035,     36.1714,     36.1714,
    36.1714,     35.9323,     35.9323,     35.7531,     35.6845,     35.6845,     35.4366,     35.4366,     35.2539,     35.1834,     35.1834,     34.9320,     34.9320,     34.9320,     34.6789,
    34.6789,     34.4293,     34.4293,     34.4293,     34.1814,     34.1814,     33.9317,     33.9317,     33.9317,     33.6750,     33.6750,     33.4148,     33.4148,     33.4148,     33.1494,
    33.1494,     32.9666,     32.8733,     32.8733,     32.5902,     32.5902,     32.4005,     32.3001,     32.3001,     32.0030,     32.0030,     32.0030,     31.6971,     31.6971,     31.3876,
    31.3876,     31.3876,     31.0677,     31.0677,     30.7424,     30.7424,     30.7424,     30.4066,     30.4066,     30.0637,     30.0637,     30.0637,     29.7156,     29.7156,     29.3675,
    29.3675,     29.3675,     29.0246,     29.0246,     28.8103,     28.6817,     28.6817,     28.3459,     28.3459,     28.3459,     28.0206,     28.0206,     27.6988,     27.6988,     27.6988,
    27.3770,     27.3770,     27.0552,     27.0552,     27.0552,     26.7263,     26.7263,     26.3834,     26.3834,     26.3834,     26.0265,     26.0265,     25.6537,     25.6537,     25.6537,
    25.2651,     25.2651,     25.0191,     24.8571,     24.8571,     24.4333,     24.4333,     24.4333,     23.9938,     23.9938,     23.5366,     23.5366,     23.5366,     23.0601,     23.0601,
    22.5554,     22.5554,     22.5554,     22.0279,     22.0279,     21.4652,     21.4652,     21.4652,     20.8568,     20.8568,     20.2045,     20.2045,     20.2045,     19.4941,     19.4941,
    19.0092,     18.7275,     18.7275,     17.9047,     17.9047,     17.9047,     17.0290,     17.0290,     16.1059,     16.1059,     16.1059,     15.1529,     15.1529,     14.1770,     14.1770,
    14.1770,     13.1906,     13.1906,     12.2077,     12.2077,     12.2077,     11.2458,     11.2458,    10.30332,    10.30332,    10.30332,     9.39423,     9.39423,     8.53087,     8.53087,
    8.53087,     7.71324,     7.71324,     7.71324,     6.95540,     6.95540,     6.25737,     6.25737,     6.25737,     5.63144,     5.63144,     5.07585,     5.07585,     5.07585,     4.59587,
    4.59587,     4.18975,     4.18975,     4.18975,     3.85395,     3.85395,     3.57794,     3.57794,     3.57794,     3.35994,     3.35994,     3.18413,     3.18413,     3.18413,     3.04349,
    3.04349,     3.04349,     2.93800,     2.93800,     2.85186,     2.85186,     2.85186,     2.79034,     2.79034,     2.74464,     2.74464,     2.74464,     2.71829,     2.71829,     2.71128,
    2.71128,     2.71128,     2.72713,     2.72713,     2.76233,     2.76233,     2.76233,     2.82390,     2.82390,     2.91360,     2.91360,     2.91360,     3.03496,     3.03496,     3.18621,
    3.18621,     3.18621,     3.37615,     3.37615,     3.37615,     3.61182,     3.61182,     3.89673,     3.89673,     3.89673,     4.24319,     4.24319,     4.65471,     4.65471,     4.65471,
    5.14185,     5.14185,     5.70812,     5.70812,     5.70812,     6.35527,     6.35527,     7.07802,     7.07802,     7.07802,     7.87286,     7.87286,     8.73274,     8.73274,     8.73274,
    9.64360,     9.64360,     9.64360,     10.5914,     10.5914,     11.5655,     11.5655,     11.5655,     12.5467,     12.5467,     13.5313,     13.5313,     13.5313,     14.4966,     14.4966,
    15.4425,     15.4425,     15.4425,     16.3498,     16.3498,     17.2166,     17.2166,     17.2166,     18.0306,     18.0306,     18.7919,     18.7919,     18.7919,     19.4986,     19.4986,
    19.7117,     20.1509,     20.1509,     20.7575,     20.7575,     20.7575,     21.3183,     21.3183,     21.8423,     21.8423,     21.8423,     22.3398,     22.3398,     22.8110,     22.8110,
    22.8110,     23.2629,     23.2629,     23.6954,     23.6954,     23.6954,     24.1069,     24.1069,     24.5007,     24.5007,     24.5007,     24.8752,     24.8752,     25.0126,     25.2339,
    25.2339,     25.5786,     25.5786,     25.5786,     25.9127,     25.9127,     26.2450,     26.2450,     26.2450,     26.5738,     26.5738,     26.9009,     26.9009,     26.9009,     27.2244,
    27.2244,     27.5445,     27.5445,     27.5445,     27.8557,     27.8557,     28.1652,     28.1652,     28.1652,     28.4606,     28.4606,     28.5838,     28.7525,     28.7525,     29.0444,
    29.0444,     29.0444,     29.3293,     29.3293,     29.6142,     29.6142,     29.6142,     29.8885,     29.8885,     30.1488,     30.1488,     30.1488,     30.3809,     30.3809,     30.5779,
    30.5779,     30.5779,     30.7345,     30.7345,     30.8613,     30.8613,     30.8613,     30.9599,     30.9599,     31.0637,     31.0567,     31.0567,     31.1694,     31.1694,     31.2680,
    31.3242,     31.3242,     31.5370,     31.5370,     31.5370,     31.8218,     31.8218,     32.1734,     32.1734,     32.1734,     32.5794,     32.5794,     33.0206,     33.0206,     33.0206,
    33.4758,     33.4758,     33.9240,     33.9240,     33.9240,     34.3459,     34.3459,     34.4216,     34.7343,     34.7343,     35.0824,     35.0824,     35.1545,     35.3970,     35.3970,
    35.6853,     35.6853,     35.6853,     35.9578,     35.9578,     36.2162,     36.2162,     36.2162,     36.4693,     36.4693,     36.7225,     36.7225,     36.7225,     36.9669,     36.9669,
    37.2112,     37.2112,     37.2112,     37.4485,     37.4485,     37.5101,     37.6823,     37.6823,     37.9091,     37.9091,     37.9690,     38.1306,     38.1306,     38.3433,     38.3433,
    38.4014,     38.5525,     38.5525,     38.7547,     38.7547,     38.7547,     38.9498,     38.9498,     39.1362,     39.1362,     39.1362,     39.3173,     39.3173,     39.4966,     39.4966,
    39.4966,     39.6724,     39.6724,     39.7182,     39.8411,     39.8411,     40.0099,     40.0099,     40.0522,     40.1787,     40.1787,     40.3404,     40.3404,     40.3809,     40.5004,
    40.5004,     40.6568,     40.6568,     40.6568,     40.8062,     40.8062,     40.9521,     40.9521,     40.9521,     41.0946,     41.0946,     41.2352,     41.2352,     41.2352,     41.3794,
    41.3794,     41.4198,     41.5218,     41.5218,     41.6624,     41.6624,     41.6958,     41.8013,     41.8013,     41.9331,     41.9331,     41.9613,     42.0632,     42.0632,     42.1880,
    42.1880,     42.1880,     42.3128,     42.3128,     42.4341,     42.4341,     42.4341,     42.5624,     42.5624,     42.6942,     42.6942,     42.6942,     42.8296,     42.8296,     42.8613,
    42.9667,     42.9667,     43.1056,     43.1056,     43.1337,     43.2462,     43.2462,     43.3762,     43.3762,     43.3991,     43.4975,     43.4975,     43.6118,     43.6118,     43.6347,
    43.7225,     43.7225,     43.8297,     43.8297,     43.8297,     43.9352,     43.9352,     44.0442,     44.0442,     44.0442,     44.1549,     44.1549,     44.1796,     44.2709,     44.2709,
    44.3887,     44.3887,     44.4116,     44.4994,     44.4994,     44.6067,     44.6067,     44.6243,     44.7016,     44.7016,     44.7895,     44.7895,     44.8035,     44.8633,     44.8633,
    44.9371,     44.9371,     44.9371,     45.0004,     45.0004,     45.0636,     45.0636,     45.0636,     45.1322,     45.1322,     45.1410,     45.2007,     45.2007,     45.2710,     45.2710,
    45.2781,     45.3413,     45.3413,     45.4169,     45.4169,     45.4257,     45.4960,     45.4960,     45.5733,     45.5733,     45.5839,     45.6471,     45.6471,     45.7192,     45.7192,
    45.7280,     45.7877,     45.7877,     45.8492,     45.8492,     45.8492,     45.9037,     45.9037,     45.8931,     45.9458,     45.9458,     45.9827,     45.9827,     45.9598,     46.0090,
    46.0090,     46.0335,     46.0335,     46.0089,     46.0546,     46.0546,     46.0809,     46.0809,     46.0703,     46.1160,     46.1160,     46.1582,     46.1582,     46.1722,     46.2179,
    46.2179,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,     46.2883,
    46.2883};
  
  

  double res=0.0;
    
  if (HeliocentricVelocity<=sodiumRadiationPressure_HeliocenticVelocity_Min__Combi_1997_icarus) res=sodiumRadiationPressure_Table__Combi_1997_icarus[0];
  else if (HeliocentricVelocity>=sodiumRadiationPressure_HeliocenticVelocity_Max__Combi_1997_icarus) res=sodiumRadiationPressure_Table__Combi_1997_icarus[sodiumRadiationPressure_TableLength__Combi_1997_icarus-1];
  else {
    int level;
    double x;
      
    x=(HeliocentricVelocity-sodiumRadiationPressure_HeliocenticVelocity_Min__Combi_1997_icarus)/sodiumRadiationPressure_dHeliocenticVelocity__Combi_1997_icarus;
    level=(int)x;
    x-=level;

    res=sodiumRadiationPressure_Table__Combi_1997_icarus[level]+x*(sodiumRadiationPressure_Table__Combi_1997_icarus[level+1]-sodiumRadiationPressure_Table__Combi_1997_icarus[level]);
  }
  
  //convert the acceleration from [cm s^{-2}] to [m s^{-2}] and scale it to the required heliocentric distance 
  res*=0.01*pow(149598000000.0/HeliocentricDistance,2); 

  return res;
}





double SodiumGfactor__5891_58A__Killen_2009_AJSS(double HeliocentricVelocity,double HeliocentricDistance) {

  struct cGFactor {
    double RadiaVelocity,gFactor;
  };

  static const int nDataElements=201;
  static const double DataHeliocentricDistance=0.352*_AU_;

  static const cGFactor gFactor5891_58A__Killen_2009_AJSS__Fig_7[nDataElements] = {
  {-12.00000,42.33717}, {-11.88000,42.19121}, {-11.76000,42.00876}, {-11.64000,41.75333}, {-11.52000,41.57088}, {-11.40000,41.31545}, {-11.28000,41.09652}, {-11.16000,40.84109}, {-11.04000,40.58566}, {-10.92000,40.36672},
  {-10.80000,40.03831}, {-10.68000,39.78289}, {-10.56000,39.52746}, {-10.44000,39.23554}, {-10.32000,38.94363}, {-10.20000,38.68819}, {-10.08000,38.39628}, {-9.96000,38.10436}, {-9.84000,37.81244}, {-9.72000,37.48404},
  {-9.60000,37.15562}, {-9.48000,36.79074}, {-9.36000,36.42583}, {-9.24000,36.09743}, {-9.12000,35.73253}, {-9.00000,35.40412}, {-8.88000,35.07571}, {-8.76000,34.63784}, {-8.64000,34.23645}, {-8.52000,33.79858},
  {-8.40000,33.39719}, {-8.28000,32.92282}, {-8.16000,32.41198}, {-8.04000,31.79164}, {-7.92000,31.28079}, {-7.80000,30.73344}, {-7.68000,30.18610}, {-7.56000,29.63875}, {-7.44000,29.01843}, {-7.32000,28.47108},
  {-7.20000,27.88725}, {-7.08000,27.23043}, {-6.96000,26.53712}, {-6.84000,25.84383}, {-6.72000,25.22349}, {-6.60000,24.45721}, {-6.48000,23.72742}, {-6.36000,22.99763}, {-6.24000,22.30432}, {-6.12000,21.57453},
  {-6.00000,20.91772}, {-5.88000,20.15143}, {-5.76000,19.34866}, {-5.64000,18.58237}, {-5.52000,17.88907}, {-5.40000,17.15928}, {-5.28000,16.46597}, {-5.16000,15.77267}, {-5.04000,15.18883}, {-4.92000,14.56851},
  {-4.80000,13.91170}, {-4.68000,13.32786}, {-4.56000,12.74403}, {-4.44000,12.16019}, {-4.32000,11.57635}, {-4.20000,11.13848}, {-4.08000,10.66411}, {-3.96000,10.22624}, {-3.84000,9.78836}, {-3.72000,9.38697},
  {-3.60000,8.94910}, {-3.48000,8.58420}, {-3.36000,8.21930}, {-3.24000,7.85441}, {-3.12000,7.52600}, {-3.00000,7.16110}, {-2.88000,6.86918}, {-2.76000,6.57727}, {-2.64000,6.35833}, {-2.52000,6.13939},
  {-2.40000,5.92045}, {-2.28000,5.77449}, {-2.16000,5.62854}, {-2.04000,5.48258}, {-1.92000,5.33662}, {-1.80000,5.22715}, {-1.68000,5.08119}, {-1.56000,5.00821}, {-1.44000,4.89874}, {-1.32000,4.82576},
  {-1.20000,4.75278}, {-1.08000,4.71629}, {-0.96000,4.64331}, {-0.84000,4.57033}, {-0.72000,4.49735}, {-0.60000,4.42438}, {-0.48000,4.38789}, {-0.36000,4.31491}, {-0.24000,4.27842}, {-0.12000,4.24193},
  {0.00000,4.24193}, {0.12000,4.24193}, {0.24000,4.20544}, {0.36000,4.20544}, {0.48000,4.24193}, {0.60000,4.27842}, {0.72000,4.31491}, {0.84000,4.35140}, {0.96000,4.42438}, {1.08000,4.49735},
  {1.20000,4.57033}, {1.32000,4.67980}, {1.44000,4.78927}, {1.56000,4.86225}, {1.68000,5.00821}, {1.80000,5.11768}, {1.92000,5.26364}, {2.04000,5.40960}, {2.16000,5.62854}, {2.28000,5.88396},
  {2.40000,6.13939}, {2.52000,6.43131}, {2.64000,6.67761}, {2.76000,6.95129}, {2.88000,7.22496}, {3.00000,7.55337}, {3.12000,7.90914}, {3.24000,8.29228}, {3.36000,8.62069}, {3.48000,8.94910},
  {3.60000,9.27750}, {3.72000,9.63328}, {3.84000,10.04379}, {3.96000,10.56377}, {4.08000,11.11111}, {4.20000,11.68582}, {4.32000,12.28790}, {4.44000,12.80788}, {4.56000,13.35523}, {4.68000,14.03941},
  {4.80000,14.77833}, {4.92000,15.48987}, {5.04000,16.17406}, {5.16000,16.99507}, {5.28000,17.81609}, {5.40000,18.55501}, {5.52000,19.18446}, {5.64000,19.89601}, {5.76000,20.71702}, {5.88000,21.45594},
  {6.00000,22.05802}, {6.12000,22.85167}, {6.24000,23.59059}, {6.36000,24.32950}, {6.48000,24.93158}, {6.60000,25.64313}, {6.72000,26.29994}, {6.84000,26.87466}, {6.96000,27.47674}, {7.08000,28.07882},
  {7.20000,28.68090}, {7.32000,29.31034}, {7.44000,29.93979}, {7.56000,30.43240}, {7.68000,30.87028}, {7.80000,31.39026}, {7.92000,31.88287}, {8.04000,32.43021}, {8.16000,32.86809}, {8.28000,33.27860},
  {8.40000,33.71648}, {8.52000,34.09962}, {8.64000,34.48276}, {8.76000,34.81116}, {8.88000,35.13957}, {9.00000,35.49535}, {9.12000,35.79639}, {9.24000,36.09743}, {9.36000,36.39847}, {9.48000,36.67214},
  {9.60000,36.97318}, {9.72000,37.27422}, {9.84000,37.60263}, {9.96000,37.90367}, {10.08000,38.17734}, {10.20000,38.45101}, {10.32000,38.69732}, {10.44000,38.94363}, {10.56000,39.24466}, {10.68000,39.51834},
  {10.80000,39.81937}, {10.92000,40.06568}, {11.04000,40.31199}, {11.16000,40.50356}, {11.28000,40.77723}, {11.40000,40.99617}, {11.52000,41.21511}, {11.64000,41.43404}, {11.76000,41.65298}, {11.88000,41.87192},
  {12.00000,42.06349}};

  static const double dv=(gFactor5891_58A__Killen_2009_AJSS__Fig_7[1].RadiaVelocity-gFactor5891_58A__Killen_2009_AJSS__Fig_7[0].RadiaVelocity)*1.0E3;
  static const double vmin=gFactor5891_58A__Killen_2009_AJSS__Fig_7[0].RadiaVelocity*1.0E3;

  int Interval=(int)((HeliocentricVelocity-vmin)/dv);
  double res;

  if (Interval<0) res=gFactor5891_58A__Killen_2009_AJSS__Fig_7[0].gFactor;
  else if (Interval>=nDataElements) res=gFactor5891_58A__Killen_2009_AJSS__Fig_7[nDataElements-1].gFactor;
  else res=gFactor5891_58A__Killen_2009_AJSS__Fig_7[Interval].gFactor;

  res*=pow(DataHeliocentricDistance/HeliocentricDistance,2);

  return res;
}


double SodiumGfactor__5897_56A__Killen_2009_AJSS(double HeliocentricVelocity,double HeliocentricDistance) {

  struct cGFactor {
    double RadiaVelocity,gFactor;
  };

  static const int nDataElements=201;
  static const double DataHeliocentricDistance=0.352*_AU_;

  static const cGFactor gFactor5897_56A__Killen_2009_AJSS__Fig_7[nDataElements] = {
      {-12.00000,26.13574}, {-11.88000,25.99891}, {-11.76000,25.83470}, {-11.64000,25.67050}, {-11.52000,25.50629}, {-11.40000,25.31472}, {-11.28000,25.12315}, {-11.16000,24.98632}, {-11.04000,24.82211}, {-10.92000,24.65791},
      {-10.80000,24.52107}, {-10.68000,24.38424}, {-10.56000,24.24740}, {-10.44000,24.08320}, {-10.32000,23.89163}, {-10.20000,23.72742}, {-10.08000,23.56322}, {-9.96000,23.37165}, {-9.84000,23.20744}, {-9.72000,23.01587},
      {-9.60000,22.79693}, {-9.48000,22.57800}, {-9.36000,22.33169}, {-9.24000,22.11275}, {-9.12000,21.94855}, {-9.00000,21.75698}, {-8.88000,21.51067}, {-8.76000,21.29174}, {-8.64000,21.04543}, {-8.52000,20.77176},
      {-8.40000,20.52545}, {-8.28000,20.25178}, {-8.16000,19.97811}, {-8.04000,19.67707}, {-7.92000,19.40339}, {-7.80000,19.12972}, {-7.68000,18.91078}, {-7.56000,18.63711}, {-7.44000,18.36344}, {-7.32000,18.06240},
      {-7.20000,17.70662}, {-7.08000,17.35085}, {-6.96000,16.96771}, {-6.84000,16.61193}, {-6.72000,16.31089}, {-6.60000,16.03722}, {-6.48000,15.68145}, {-6.36000,15.29830}, {-6.24000,14.94253}, {-6.12000,14.53202},
      {-6.00000,14.14888}, {-5.88000,13.79310}, {-5.76000,13.38259}, {-5.64000,12.83525}, {-5.52000,12.31527}, {-5.40000,11.93213}, {-5.28000,11.54899}, {-5.16000,11.13848}, {-5.04000,10.81007}, {-4.92000,10.42693},
      {-4.80000,9.98905}, {-4.68000,9.46908}, {-4.56000,8.97646}, {-4.44000,8.51122}, {-4.32000,8.15545}, {-4.20000,7.82704}, {-4.08000,7.41653}, {-3.96000,7.06076}, {-3.84000,6.73235}, {-3.72000,6.29447},
      {-3.60000,5.93870}, {-3.48000,5.71976}, {-3.36000,5.44609}, {-3.24000,5.19978}, {-3.12000,4.92611}, {-3.00000,4.67980}, {-2.88000,4.48823}, {-2.76000,4.29666}, {-2.64000,4.10509}, {-2.52000,3.91352},
      {-2.40000,3.72195}, {-2.28000,3.61248}, {-2.16000,3.47564}, {-2.04000,3.36617}, {-1.92000,3.25671}, {-1.80000,3.11987}, {-1.68000,3.01040}, {-1.56000,2.92830}, {-1.44000,2.87356}, {-1.32000,2.81883},
      {-1.20000,2.76409}, {-1.08000,2.70936}, {-0.96000,2.68199}, {-0.84000,2.65463}, {-0.72000,2.62726}, {-0.60000,2.62726}, {-0.48000,2.59989}, {-0.36000,2.59989}, {-0.24000,2.62726}, {-0.12000,2.62726},
      {0.00000,2.62726}, {0.12000,2.59989}, {0.24000,2.59989}, {0.36000,2.57252}, {0.48000,2.59989}, {0.60000,2.62726}, {0.72000,2.62726}, {0.84000,2.62726}, {0.96000,2.62726}, {1.08000,2.62726},
      {1.20000,2.65463}, {1.32000,2.68199}, {1.44000,2.73673}, {1.56000,2.81883}, {1.68000,2.90093}, {1.80000,2.98303}, {1.92000,3.09250}, {2.04000,3.20197}, {2.16000,3.33881}, {2.28000,3.44828},
      {2.40000,3.55774}, {2.52000,3.69458}, {2.64000,3.85878}, {2.76000,4.05036}, {2.88000,4.26929}, {3.00000,4.48823}, {3.12000,4.67980}, {3.24000,4.89874}, {3.36000,5.11768}, {3.48000,5.39135},
      {3.60000,5.66502}, {3.72000,5.93870}, {3.84000,6.23974}, {3.96000,6.54078}, {4.08000,6.89655}, {4.20000,7.25233}, {4.32000,7.60810}, {4.44000,7.90914}, {4.56000,8.23755}, {4.68000,8.56596},
      {4.80000,8.89436}, {4.92000,9.27750}, {5.04000,9.77011}, {5.16000,10.23536}, {5.28000,10.67323}, {5.40000,11.11111}, {5.52000,11.54899}, {5.64000,11.93213}, {5.76000,12.34264}, {5.88000,12.75315},
      {6.00000,13.16366}, {6.12000,13.60153}, {6.24000,14.01204}, {6.36000,14.44992}, {6.48000,14.83306}, {6.60000,15.21620}, {6.72000,15.57198}, {6.84000,15.95512}, {6.96000,16.28353}, {7.08000,16.61193},
      {7.20000,16.96771}, {7.32000,17.26875}, {7.44000,17.59715}, {7.56000,17.92556}, {7.68000,18.22660}, {7.80000,18.55501}, {7.92000,18.85605}, {8.04000,19.10235}, {8.16000,19.34866}, {8.28000,19.62233},
      {8.40000,19.89601}, {8.52000,20.16968}, {8.64000,20.44335}, {8.76000,20.63492}, {8.88000,20.79912}, {9.00000,20.96333}, {9.12000,21.12753}, {9.24000,21.26437}, {9.36000,21.42857}, {9.48000,21.67488},
      {9.60000,21.86645}, {9.72000,22.08539}, {9.84000,22.30432}, {9.96000,22.52326}, {10.08000,22.66010}, {10.20000,22.82430}, {10.32000,22.98851}, {10.44000,23.15271}, {10.56000,23.28955}, {10.68000,23.45375},
      {10.80000,23.61795}, {10.92000,23.80952}, {11.04000,23.97373}, {11.16000,24.11056}, {11.28000,24.27477}, {11.40000,24.43897}, {11.52000,24.57581}, {11.64000,24.76738}, {11.76000,24.95895}, {11.88000,25.06842},
      {12.00000,25.12315}};


  static const double dv=(gFactor5897_56A__Killen_2009_AJSS__Fig_7[1].RadiaVelocity-gFactor5897_56A__Killen_2009_AJSS__Fig_7[0].RadiaVelocity)*1.0E3;
  static const double vmin=gFactor5897_56A__Killen_2009_AJSS__Fig_7[0].RadiaVelocity*1.0E3;

  int Interval=(int)((HeliocentricVelocity-vmin)/dv);
  double res;

  if (Interval<0) res=gFactor5897_56A__Killen_2009_AJSS__Fig_7[0].gFactor;
  else if (Interval>=nDataElements) res=gFactor5897_56A__Killen_2009_AJSS__Fig_7[nDataElements-1].gFactor;
  else res=gFactor5897_56A__Killen_2009_AJSS__Fig_7[Interval].gFactor;

  res*=pow(DataHeliocentricDistance/HeliocentricDistance,2);

  return res;
}


#endif

