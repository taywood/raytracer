// A practical implementation of the ray tracing algorithm.

#include "SDL.h" 
#include <fstream>
#include <chrono>
#include "common.h"
#include "geometry.h"
#include "hittable.h"
#include "hittable_list.h"
#include "sphere.h"
#include "ray.h"
#include "raytracer.h"
#include "camera.h"
#include "material.h"
#include "multithreading.h"
#include <thread>

#define M_PI 3.14159265359;

SDL_Window* window;
SDL_Renderer* renderer;
SDL_Surface* screen;

void init() {
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow(
        "Software Ray Tracer",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        640,
        480,
        0
    );

    screen = SDL_GetWindowSurface(window);

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
}

void putpixel(SDL_Surface* surface, int x, int y, Uint32 pixel)
{
    int bpp = surface->format->BytesPerPixel;
    /* Here p is the address to the pixel we want to set */
    Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

    switch (bpp) {
    case 1:
        *p = pixel;
        break;

    case 2:
        *(Uint16*)p = pixel;
        break;

    case 3:
        if (SDL_BYTEORDER == SDL_BIG_ENDIAN) {
            p[0] = (pixel >> 16) & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = pixel & 0xff;
        }
        else {
            p[0] = pixel & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = (pixel >> 16) & 0xff;
        }
        break;

    case 4:
        *(Uint32*)p = pixel;
        break;
    }
}

double hit_sphere(const Point3f& centre, double radius, const ray& r) {
    Vec3f oc = r.origin() - centre;
    auto a = r.direction().dotProduct(r.direction());
    auto b = 2.0 * oc.dotProduct(r.direction());
    auto c = oc.dotProduct(oc) - radius * radius;
    auto discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return -1.0;
    }
    else {
        return (-b - sqrt(discriminant)) / (2.0 * a);
    }
}

Colour ray_colour(const ray& r, const hittable& world, int depth) {
    hit_record rec;
    if (depth <= 0) return Colour(0, 0, 0);
    if (world.hit(r, 0.0001, infinity, rec)) {
        ray scattered;
        Colour attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_colour(scattered, world, depth - 1);
        return Colour(0, 0, 0);
    }
    Vec3f unit_direction = r.direction().normalize();
    auto t = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - t) * Colour(1.0, 1.0, 1.0) + t * Colour(0.5, 0.7, 1.0) * 255;
}

struct Sphere {
    Point3f c;
    double r;
    Sphere(const Point3f& c, double r) : c(c), r(r) {}
    Vec3f getNormal(const Vec3f& pi) const { return (pi - c) / r; }

    // Solve t^2*d.d+2*t*(o-p).d+(o-p).(o-p)-R^2=0​
    bool intersect(const ray& ray, double& t) const {
        const Point3f o = ray.orig;
        const Vec3f d = ray.dir;
        const Vec3f oc = o - c;
        const double b = 2 * oc.dotProduct(d);
        const double c = oc.dotProduct(oc) - r * r;
        double disc = b * b - 4 * c; // a=1 as ray is normalised​
        if (disc < 1e-4) return false; // ray misses sphere​
        disc = sqrt(disc);
        const double t0 = -b - disc;
        const double t1 = -b + disc;
        t = (t0 < t1) ? t0 : t1; // two intersections on sphere, set t to shortest​
        return true;
    }
};

// method to ensure colours don’t go out of 8 bit range in RGB​
void clamp255(Vec3f& col) {
    col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
    col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
    col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
}

// LAB 11 //

    // Image
const auto aspect_ratio = 16.0 / 9.0;
const int image_width = screen->w;
const int image_height = static_cast<int>(image_width / aspect_ratio);
const int spp = 4;
const float scale = 1.f / spp;
const int max_depth = 50;

    hittable_list random_scene() {
        hittable_list world;
        auto ground_material = make_shared<lambertian>(Colour(0.5, 0.5, 0.5));
        world.add(make_shared<sphere>(Point3f(0, -1000, 0), 1000, ground_material));
        for (int a = -11; a < 11; a++) {
            for (int b = -11; b < 11; b++) {
                auto choose_mat = random_double();
                Point3f centre(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());
                if ((centre - Point3f(4, 0.2, 0)).length() > 0.9) {
                    shared_ptr<material> sphere_material;
                    if (choose_mat < 0.8) {
                      // diffuse
                        auto albedo = Colour::random() * Colour::random();
                        sphere_material = make_shared<lambertian>(albedo);
                        world.add(make_shared<sphere>(centre, 0.2, sphere_material));
                    }
                    else if (choose_mat < 0.95) {
                      // metal
                        auto albedo = Colour::random(0.5, 1);
                        auto fuzz = random_double(0, 0.5);
                        sphere_material = make_shared<metal>(albedo, fuzz);
                        world.add(make_shared<sphere>(centre, 0.2, sphere_material));
                    }
                    else {
                      //glass
                        sphere_material = make_shared<dielectric>(1.5);
                        world.add(make_shared<sphere>(centre, 0.2, sphere_material));
                    }
                }
            }

        }

        auto material1 = make_shared<dielectric>(1.5);
        world.add(make_shared<sphere>(Point3f(0, 1, 0), 1.0, material1));
        auto material2 = make_shared<lambertian>(Colour(0.4, 0.2, 0.1));
        world.add(make_shared<sphere>(Point3f(-4, 1, 0), 1.0, material2));
        auto material3 = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0.0);
        world.add(make_shared<sphere>(Point3f(4, 1, 0), 1.0, material3));
        return world;
    }

    // World
    auto world = random_scene();

  // Camera
    Point3f lookfrom(13, 2, 3);
    Point3f lookat(0, 0, 0);
    Vec3f vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.15;
    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);

    auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;
    auto focal_length = 1.0;
    auto origin = Point3f(0, 0, 0);
    auto horizontal = Vec3f(viewport_width, 0, 0);
    auto vertical = Vec3f(0, viewport_height, 0);
    auto lower_left_corner = origin - horizontal / 2 - vertical / 2 - Vec3f(0, 0, focal_length);

    // LAB 11 //

int main(int argc, char** argv)
{
    // initialise SDL2
    init();

    const Colour white(255, 255, 255);
    const Colour black(0, 0, 0);
    const Colour red(255, 0, 0);

    const Sphere sph(Vec3f(screen->w * 0.5, screen->h * 0.5, 50), 50);
    const Sphere light(Vec3f(0, 0, 50), 1);

// World

//    auto R = cos(pi / 4);
//    hittable_list world;

    //auto material_left = make_shared<lambertian>(Colour(0, 0, 1));
    //auto material_right = make_shared<lambertian>(Colour(1, 0, 0));

      //  hittable_list world;
        //auto material_ground = make_shared<lambertian>(Colour(0.8, 0.8, 0.0));
       // auto material_center = make_shared<lambertian>(Colour(0.1, 0.2, 0.5));
       // auto material_left = make_shared<dielectric>(1.5);
       // auto material_right = make_shared<metal>(Colour(0.8, 0.6, 0.2), 1.0);

//        world.add(make_shared<sphere>(Point3f(0.0, -100.5, -1.0), 100.0, material_ground));
  //      world.add(make_shared<sphere>(Point3f(0.0, 0.0, -1.0), 0.5, material_center));
    //    world.add(make_shared<sphere>(Point3f(-1.0, 0.0, 1.0), -0.4, material_left));
      //  world.add(make_shared<sphere>(Point3f(1.0, 0.0, -1.0), 0.5, material_right));

    // Camera
        //Point3f lookfrom(3, 3, 2);
        //Point3f lookat(0, 0, -1);
        //Vec3f vup(0, 1, 0);
        //auto dist_to_focus = (lookfrom - lookat).length();
        //auto aperture = 2.0;
        //camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);

    double t;
    Colour pix_col(black);

    SDL_Event e;
    bool running = true;
    while (running) {

        auto t_start = std::chrono::high_resolution_clock::now();

        // clear back buffer, pixel data on surface and depth buffer (as movement)
        SDL_FillRect(screen, nullptr, SDL_MapRGB(screen->format, 0, 0, 0));
        SDL_RenderClear(renderer);

        for (int y = screen->h - 1; y >= 0; --y) {
            std::cerr << "/rScanlines remaining: " << y << std::flush;
            for (int x = 0; x < screen->w; ++x) {
                pix_col = black;
                for (int s = 0; s < spp; s++) {
                    auto u = double(x + random_double()) / (image_width - 1);
                    auto v = double(y + random_double()) / (image_height - 1);
                    ray ray = cam.get_ray(u, v);
                    pix_col = pix_col + ray_colour(ray, world, max_depth);
                    pix_col /= 255.f * spp;
                    pix_col.x = sqrt(pix_col.x );
                    pix_col.y = sqrt(pix_col.y );
                    pix_col.z = sqrt(pix_col.z );
                    pix_col *= 255;
                }
                Uint32 colour = SDL_MapRGB(screen->format, pix_col.x, pix_col.y, pix_col.z);
                putpixel(screen, x, y, colour);
            }
        }

        {
            t_start = std::chrono::high_resolution_clock::now();
            ThreadPool pool(std::thread::hardware_concurrency());

            int start = screen->h - 1;
            int step = screen->h / std::thread::hardware_concurrency();
            for (int y = 0; y < screen->h - 1; y++)
            {
                pool.Enqueue(std::bind(LineRender, screen, world, y, spp, max_depth, &cam));
            }
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        std::cerr << "Frame render time:  " << passedTime << " ms" << std::endl;

        SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, screen);
        if (texture == NULL) {
            fprintf(stderr, "CreateTextureFromSurface failed: %s\n", SDL_GetError());
            exit(1);
        }
        SDL_FreeSurface(screen);

        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderCopyEx(renderer, texture, nullptr, nullptr, 0, 0, SDL_FLIP_VERTICAL);
        SDL_RenderPresent(renderer);

        SDL_DestroyTexture(texture);

        if (SDL_PollEvent(&e))
        {
            switch (e.type) {
            case SDL_KEYDOWN:
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE:
                    running = false;
                    break;
                }
                break;
            }
        }
    }

    return 0;
}
