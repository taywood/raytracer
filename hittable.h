#pragma once
#include "ray.h"
#include "geometry.h"

struct hit_record {
	Point3f p;
	Vec3f normal;
	double t;
	bool front_face;

	inline void set_face_normal(const ray& r, const Vec3f& outward_normal) {
		front_face = (r.direction().dotProduct(outward_normal)) < 0;
		normal = front_face ? outward_normal : outward_normal;
	}
};

class hittable {
public:
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};