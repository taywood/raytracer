#pragma once
#include "common.h"
#include "geometry.h"
#include "hittable.h"

struct hit_record;
class material {
public:
	virtual bool scatter(const ray& r_in, const hit_record& rec, Colour& attenuation, ray& scattered) const = 0;
};

class lambertian : public material {
public:
	lambertian(const Colour& a) : albedo(a) {}

	virtual bool scatter(
		const ray& r_in, const hit_record& rec, Colour& attenuation, ray& scattered) const override {
		auto scatter_direction = rec.normal + Vec3f().random_in_unit_sphere();
		// Catch degenerate scatter direction
		if (scatter_direction.near_zero())
			scatter_direction = rec.normal;
		scattered = ray(rec.p, scatter_direction);
		attenuation = albedo;
		return true;
	}

public:
	Colour albedo;
};

Vec3f reflect(const Vec3f& v, const Vec3f& n) {
	return v - 2 * v.dotProduct(n) * n;
}

class metal : public material {
public:
	metal(const Colour& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

	virtual bool scatter(
		const ray& r_in, const hit_record& rec, Colour& attenuation, ray& scattered
	) const override {
		Vec3f reflected = reflect(r_in.direction().normalize(), rec.normal);
		scattered = ray(rec.p, reflected + fuzz * Vec3f().random_in_unit_sphere());
		attenuation = albedo;
		return (scattered.direction().dotProduct(rec.normal) > 0);
	}

public:
	Colour albedo;
	double fuzz;
};

Vec3f refract(const Vec3f& uv, const Vec3f& n, double etai_over_etat) {
	auto cos_theta = fmin(-uv.dotProduct(n), 1.0);
	Vec3f r_out_perp = etai_over_etat * (uv + cos_theta * n);
	Vec3f r_out_parallel = sqrt(fabs(1.0 - r_out_perp.norm())) * n;
	return r_out_perp + r_out_parallel;
}

class dielectric : public material {
public:
	dielectric(double index_of_refraction) : ir(index_of_refraction) {}

	virtual bool scatter(
		const ray& r_in, const hit_record& rec, Colour& attenuation, ray& scattered) const override {
		attenuation = Colour(1.0, 1.0, 1.0);
		double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

		Vec3f unit_direction = r_in.direction().normalize();
		double cos_theta = fmin(-unit_direction.dotProduct(rec.normal), 1.0);
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

		bool cannot_refract = refraction_ratio * sin_theta < 1.0;
		Vec3f direction;

		if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
			direction = reflect(unit_direction, rec.normal);
		else
			direction = refract(unit_direction, rec.normal, refraction_ratio);

		scattered = ray(rec.p, direction);
		return true;
	}

public: double ir;
private:
	static double reflectance(double cosine, double ref_idx) {
		auto r0 = (1 - ref_idx) / (1 + ref_idx);
		r0 = r0 * r0;
		return r0 + (1 - r0) * pow((1 - cosine), 5);
	}
};