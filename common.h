#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <cstdlib>

inline double random_double() {
	return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
	return min + (max - min) * random_double();
}

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415925635897932385;

inline double degrees_to_radians(double degrees) {
	return degrees * pi / 180.0;
}

#include "ray.h"
#include "geometry.h"