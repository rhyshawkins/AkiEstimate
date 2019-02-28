
#include <stdio.h>

#include "simplified_ak135.hpp"

static std::vector<simplified_ak135entry>::const_iterator
simplified_ak135_interpolation_find(const std::vector<simplified_ak135entry>::const_iterator &start,
			 const std::vector<simplified_ak135entry>::const_iterator &end,
			 double depth,
			 double &alpha)
{
  if (depth <= start->depth) {
    alpha = 0.0;
    return start;
  } else if (depth >= end->depth) {
    alpha = 1.0;
    return end - 1;
  } else {

    if ((end - start) > 1) {
      
      std::vector<simplified_ak135entry>::const_iterator centre = start + (end - start)/2;
      
      if (depth < centre->depth) {
	return simplified_ak135_interpolation_find(start, centre, depth, alpha);
      } else {
	return simplified_ak135_interpolation_find(centre, end, depth, alpha);
      }
    } else {

      alpha = (depth - start->depth)/(end->depth - start->depth);
      return start;
      
    }
  }
}
			 
double simplified_ak135_density(double depth)
{
  double alpha;
  std::vector<simplified_ak135entry>::const_iterator i = simplified_ak135_interpolation_find(simplified_ak135_table.cbegin(),
								       simplified_ak135_table.cend() - 1,
								       depth / 1000.0,
								       alpha);

  return (i->density * (1.0 - alpha) + (i + 1)->density * alpha) * 1000.0;
}

double simplified_ak135_Vp(double depth)
{
  double alpha;
  std::vector<simplified_ak135entry>::const_iterator i = simplified_ak135_interpolation_find(simplified_ak135_table.cbegin(),
								       simplified_ak135_table.cend() - 1,
								       depth / 1000.0,
								       alpha);

  return (i->Vp * (1.0 - alpha) + (i + 1)->Vp * alpha) * 1000.0;
}

double simplified_ak135_Vs(double depth)
{
  double alpha;
  std::vector<simplified_ak135entry>::const_iterator i = simplified_ak135_interpolation_find(simplified_ak135_table.cbegin(),
								       simplified_ak135_table.cend() - 1,
								       depth / 1000.0,
								       alpha);

  return (i->Vs * (1.0 - alpha) + (i + 1)->Vs * alpha) * 1000.0;
}


const std::vector<simplified_ak135entry> simplified_ak135_table = {
  {   0.00,   1.0200,    1.4500,   1.0000},
  {  18.00,   2.9200,    6.8000,   3.9000},
  {  18.00,   3.6410,    8.0355,   4.4839},
  { 210.00,   3.3243,    8.3007,   4.5184},
  { 210.00,   3.3243,    8.3007,   4.5184},
  { 410.00,   3.5068,    9.0302,   4.8702},
  { 410.00,   3.9317,    9.3601,   5.0806},
  { 660.00,   3.9201,   10.2000,   5.6104},
  { 660.00,   4.2387,   10.7909,   5.9607},
  {2891.50,   5.7721,   13.6601,   7.2817},
  {2891.50,   9.9145,    8.0000,   0.0000},
  {5153.50,  12.1391,   10.2890,   0.0000},
  {5153.50,  12.7037,   11.0427,   3.5043},
  {6371.00,  13.0122,   11.2622,   3.6678}
};

