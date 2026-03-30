#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

#include "concaveman.h"

typedef std::array<double, 2> Point;

int orientation(Point a, Point b, Point c)
{
    double v = a[0]*(b[1]-c[1]) + b[0]*(c[1]-a[1]) + c[0]*(a[1]-b[1]);
    if (v < 0) return -1; // clockwise
    if (v > 0) return +1; // counter-clockwise
    return 0;
}

bool clockwise(Point a, Point b, Point c, bool include_collinear)
{
    int o = orientation(a, b, c);
    return o < 0 || (include_collinear && o == 0);
}

bool collinear(Point a, Point b, Point c)
{
  return orientation(a, b, c) == 0;
}

void convex_hull(std::vector<Point>& points, std::vector<int>& hull, bool include_collinear = false)
{
  Point p0 = *std::min_element(points.begin(), points.end(), [](Point a, Point b) {
    return std::make_pair(a[1], a[0]) < std::make_pair(b[1], b[0]);
  });
  
  std::sort(points.begin(), points.end(), [&p0](const Point& a, const Point& b) {
    int o = orientation(p0, a, b);
    if (o == 0)
      return (p0[0]-a[0])*(p0[0]-a[0]) + (p0[1]-a[1])*(p0[1]-a[1]) < (p0[0]-b[0])*(p0[0]-b[0]) + (p0[1]-b[1])*(p0[1]-b[1]);
    return o < 0;
  });
  
  if (include_collinear)
  {
    int i = (int)points.size()-1;
    while (i >= 0 && collinear(p0, points[i], points.back())) i--;
    std::reverse(points.begin()+i+1, points.end());
  }

  hull.clear();
  for (int i = 0; i < (int)points.size(); i++)
  {
    while (hull.size() > 1 && !clockwise(points[hull[hull.size()-2]], points[hull.back()], points[i], include_collinear))
      hull.pop_back();
    hull.push_back(i);
  }

  if (!include_collinear && hull.size() == 2 && points[hull[0]][0] == points[hull[1]][0] && points[hull[0]][1] == points[hull[1]][1])
    hull.pop_back();
}

int main(int argc, char** argv)
{
  if (argc < 2 || argc > 4)
  {
    std::cerr << "Usage: " << argv[0] << " cells.csv [boundary.csv [neighbors.csv]]\n";
    return 1;
  }
  
  std::ifstream f(argv[1]);
  if (!f.is_open())
  {
    std::cerr << "Error: Unable to open file " << argv[1] << "!\n";
    return 2;
  }
  
  std::vector<Point> points;
  std::vector<int> grows, neighbors;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  
  // skip header line
  std::string line;
  std::getline(f, line);
  
  // read cell rows
  while (std::getline(f, line))
  {
    if (!line.empty())
    {
      std::istringstream is(line);
      std::string substring;
      std::vector<std::string> substrings;
      while (std::getline(is, substring, ','))
        substrings.push_back(substring);
      points.push_back({std::stod(substrings[0]), std::stod(substrings[1])});
      grows.push_back(std::stoi(substrings[2]));
      neighbors.push_back(std::stoi(substrings[3]));
    }
  }
  
  f.close();
  
  // determine the colony boundary
  std::vector<int> hull;
  convex_hull(points, hull);
  auto boundary = concaveman<double, 32>(points, hull, 1.5, 0);
  
  long double R = nan, A = nan, C = nan, w = nan;
  
  if (boundary.size() > 2)
  {
    // compute the boundary centroid
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = xmin;
    double ymax = xmax;
    long double cx = 0, cy = 0;
    for (auto & p : boundary)
    {
      cx += p[0];
      cy += p[1];
      xmin = std::min(xmin, p[0]);
      xmax = std::max(xmax, p[0]);
      ymin = std::min(ymin, p[1]);
      ymax = std::max(ymax, p[1]);
    }
    cx /= boundary.size();
    cy /= boundary.size();
    
    // measure the colony using its boundary
    R = 0, A = 0, C = 0, w = 0;
    for (std::size_t i = boundary.size() - 1, j = 0; j < boundary.size(); i = j++)
    {
      R += std::sqrt(std::pow(boundary[i][0] - cx, 2) + std::pow(boundary[i][1] - cy, 2));
      A += boundary[i][0] * boundary[j][1] - boundary[i][1] * boundary[j][0];
      C += std::sqrt(std::pow(boundary[i][0] - boundary[j][0], 2) + std::pow(boundary[i][1] - boundary[j][1], 2));
    }
    R /= boundary.size();
    A = std::abs(A) / 2;
    for (auto & p : boundary)
    {
      const double thisR = std::sqrt(std::pow(p[0] - cx, 2) + std::pow(p[1] - cy, 2));
      w += std::pow(thisR - R, 2);
    }
    w = std::sqrt(w / boundary.size()); // mean square deviation from the radius
  }
  
  // fraction of growing cells
  std::size_t Ng = 0;
  for (std::size_t i = 0; i < grows.size(); ++i)
    Ng += grows[i];
  double g = (double)Ng / grows.size();
  
  // write boundary path to optional output file
  if (argc == 3)
  {
    std::ofstream f(argv[2]);
    if (!f.is_open())
    {
      std::cerr << "Error: Unable to write file " << argv[2] << "!\n";
      return 2;
    }
    
    f << "x,y\n";
    for (auto & p : boundary)
      f <<p[0] << "," << p[1] << "\n";
    
    f.close();
  }
  
  // write cell neighbor distribution to optional output file
  if (argc == 4)
  {
    std::ofstream f(argv[3]);
    if (!f.is_open())
    {
      std::cerr << "Error: Unable to write file " << argv[3] << "!\n";
      return 2;
    }
    
    int nmax = 0;
    for (auto& n : neighbors)
      nmax = std::max(nmax, n);
    std::vector<std::size_t> nn(nmax+1, 0);
    for (std::size_t i = 0; i < points.size(); ++i)
      ++nn[neighbors[i]];
    
    f << "n,p\n";
    for (std::size_t i = 0; i < nn.size(); ++i)
      f << i << "," << 100 * (double)nn[i] / points.size() << "\n";
    
    f.close();
  }
  
  // print measurements to stdout
  std::cout << points.size() << "," << R << "," << A << "," << C << "," << w << "," << g << "\n";

  return 0;
} 
