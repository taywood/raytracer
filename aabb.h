#pragma once
#include "common.h"

class aabb {
public:
	aabb() {}
	aabb(const Point3f& mini, const Point3f& maxi) { minimum = mini; maximum = maxi; }

	Point3f min() const { return minimum; }
	Point3f max() const { return maximum; }

	bool hit(const ray& r, double t_min, double t_max) const {
		for (int a = 0; a < 3; a++) {
			auto t0 = fmin((minimum[a] - r.origin()[a]) / r.direction()[a],
				(maximum[a] - r.origin()[a]) / r.direction()[a]);
			auto t1 = fmax((minimum[a] - r.origin()[a]) / r.direction()[a],
				(maximum[a] - r.origin()[a]) / r.direction()[a]);
			t_min = fmax(t0, t_min);
			t_max = fmin(t1, t_max);
			if (t_max <= t_min)
				return false;

			}
		return true;

//		aabb surrounding_box(aabb box0, aabb box1) {
//			Point3f small(fmin(box0.min().x 1e - 3, box1.min().x - 1e-3),
//				fmin(box0.min().y - 1e-3, box1.min().y - 1e-3),
//				fmin(box0.min().z - 1e-3, box1.min().z - 1e-3));

//			Point3f big(fmax(box0.max().x + 1e-3, box1.max().x + 1e-3),
//				fmax(box0.max().y + 1e-3, box1.max().y + 1e-3),
//				fmax(box0.max().z + 1e-3, box1.max().z + 1e-3));

//			return aabb(small, big);
//		}

	}

	Point3f minimum;
	Point3f maximum;


};