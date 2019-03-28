#include "environment_light.h"

#include <algorithm>
#include <iostream>
#include <fstream>

namespace CGL { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
    	init();
}

EnvironmentLight::~EnvironmentLight() {
    delete[] pdf_envmap;
    delete[] conds_y;
    delete[] marginal_y;
}


void EnvironmentLight::init() {
	uint32_t w = envMap->w, h = envMap->h;
  pdf_envmap = new double[w * h];
	conds_y = new double[w * h];
	marginal_y = new double[h];

	std::cout << "[PathTracer] Initializing environment light...";

  // 3.3 step 1
	// Store the environment map pdf to pdf_envmap

	double sum = 0;
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
            pdf_envmap[w * j + i] = envMap->data[w * j + i].illum() * sin(PI * (j+.5) / h);
            sum += pdf_envmap[w * j + i];
		}
	}

    //normalize pdf_envmap
    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h; j++) {
        pdf_envmap[i + j * w] = pdf_envmap[i + j * w] / sum;
      }
    }
  std::cout << "Done Normalize." << std::endl;

  // TODO: 3.3 step 2
  // Store the marginal distribution for y to marginal_y. Make sure to normalize pdf_envmap.
  for (int j = 0; j < h; j++) {
    double marginal_dis = 0.0;
    for (int i = 0; i < w; i++) {
      for (int j_prime = 0; j_prime < j + 1; j_prime++) {
        marginal_dis += pdf_envmap[i + j_prime * w];
      }
    }
    marginal_y[j] = marginal_dis;
  }

  std::cout << "Done Calculating marginal_y." << std::endl;

  // TODO: 3.3 step 3
  // Store the conditional distribution for x given y to conds_y
  //p(j) = F(j) - F(j - 1)
  for (int j = 0; j < h; j++) {
    double p_j = j == 0 ? marginal_y[j] : marginal_y[j] - marginal_y[j - 1];
    for (int i = 0; i < w; i++) {
      double sum = 0;
      for (int i_prime = 0; i_prime < i; i_prime++) {
        sum += pdf_envmap[i_prime + j * w];
      }
      conds_y[i + j * w] = sum / p_j;
    }
  }

	if (true) {
    std::cout << "Saving out probability_debug image for debug." << std::endl;
    save_probability_debug();
  }

	std::cout << "done." << std::endl;
}

// Helper functions

void EnvironmentLight::save_probability_debug() {
	uint32_t w = envMap->w, h = envMap->h;
	uint8_t* img = new uint8_t[4*w*h];

	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
			img[4 * (j * w + i) + 3] = 255;
			img[4 * (j * w + i) + 0] = 255 * marginal_y[j];
			img[4 * (j * w + i) + 1] = 255 * conds_y[j * w + i];
			img[4 * (j * w + i) + 2] = 0;
		}
	}

    lodepng::encode("probability_debug.png", img, w, h);
    delete[] img;
}

Vector2D EnvironmentLight::theta_phi_to_xy(const Vector2D &theta_phi) const {
    uint32_t w = envMap->w, h = envMap->h;
    double x = theta_phi.y / 2. / PI * w;
    double y = theta_phi.x / PI * h;
    return Vector2D(x, y);
}

Vector2D EnvironmentLight::xy_to_theta_phi(const Vector2D &xy) const {
    uint32_t w = envMap->w, h = envMap->h;
    double x = xy.x;
    double y = xy.y;
    double phi = x / w * 2.0 * PI;
    double theta = y / h * PI;
    return Vector2D(theta, phi);
}

Vector2D EnvironmentLight::dir_to_theta_phi(const Vector3D &dir) const {
    Vector3D unit_dir = dir.unit();
    double theta = acos(unit_dir.y);
    double phi = atan2(-unit_dir.z, unit_dir.x) + PI;
    return Vector2D(theta, phi);
}

Vector3D EnvironmentLight::theta_phi_to_dir(const Vector2D& theta_phi) const {
    double theta = theta_phi.x;
    double phi = theta_phi.y;

    double y = cos(theta);
    double x = cos(phi - PI) * sin(theta);
    double z = -sin(phi - PI) * sin(theta);

    return Vector3D(x, y, z);
}

// Credits to Luowen Qian from Spring 2018 for this more robust bilerp
Spectrum EnvironmentLight::bilerp(const Vector2D& xy) const {
    long right = lround(xy.x), left, v = lround(xy.y);
    double u1 = right - xy.x + .5, v1;
    if (right == 0 || right == envMap->w) {
        left = envMap->w - 1;
        right = 0;
    } else left = right - 1;
    if (v == 0) v1 = v = 1; else if (v == envMap->h) {
        v = envMap->h - 1;
        v1 = 0;
    } else v1 = v - xy.y + .5;
    auto bottom = envMap->w * v, top = bottom - envMap->w;
    auto u0 = 1 - u1;
    return (envMap->data[top + left] * u1 + envMap->data[top + right] * u0) * v1 +
        (envMap->data[bottom + left] * u1 + envMap->data[bottom + right] * u0) * (1 - v1);
}


Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  // TODO: 3.2
	// First implement uniform sphere sampling for the environment light
//	*wi = sampler_uniform_sphere.get_sample();
//	*distToLight = INF_D;
//	*pdf = 1.0 / (4 * PI);
//
//    Vector2D x_y = theta_phi_to_xy(dir_to_theta_phi(*wi));
//    return bilerp(x_y);

  // TODO: 3.3
	// Later implement full importance sampling
    uint32_t w = envMap->w, h = envMap->h;
	Vector2D sample = sampler_uniform2d.get_sample();
	double r1 = sample[0];
	double r2 = sample[1];
	int y = std::upper_bound(marginal_y, marginal_y + h, r1) - marginal_y;
	int x = std::upper_bound(conds_y + y * w, conds_y + (y + 1) * w, r2) - (conds_y + y * w);

	*wi = theta_phi_to_dir(xy_to_theta_phi(Vector2D(x, y)));
	*distToLight = INF_D;
	*pdf = pdf_envmap[x + y * w] * (w * h) / (2 * pow(PI, 2) * sin_theta(*wi));
  return bilerp(Vector2D(x, y));
}

Spectrum EnvironmentLight::sample_dir(const Ray& r) const {
  // TODO: 3.1
	// Use the helper functions to convert r.d into (x,y)
	// then bilerp the return value
    Vector2D x_y = theta_phi_to_xy(dir_to_theta_phi(r.d));
	return bilerp(x_y);
}

} // namespace StaticScene
} // namespace CGL
