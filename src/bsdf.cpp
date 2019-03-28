#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>
#include <math.h>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
  else h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 1.2
  // Using BSDF::reflect(), implement sample_f for a mirror surface
  *pdf = 1.0;
  reflect(wo, wi);
  return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: 2.2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
  //double alpha = 0.5;
  double theta = acos(cos_theta(h));
  double tan_theta_square = pow(tan(theta), 2);
  double numerator = exp(-(tan_theta_square / pow(alpha, 2)));
  double deno = PI * pow(alpha, 2) * pow(cos(theta), 4);

  return numerator / deno;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: 2.3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.
  double cos_theta_i = cos_theta(wi);
  Spectrum eta_square = eta * eta;
  Spectrum k_square = k * k;

  Spectrum R_s_nom = (eta_square + k_square) - 2.0 * eta * cos_theta_i + (float) pow(cos_theta_i, 2);
  Spectrum R_s_denom = (eta_square + k_square) + 2.0 * eta * cos_theta_i + (float) pow(cos_theta_i, 2);
  Spectrum R_s = R_s_nom / R_s_denom;

  Spectrum R_p_nom = (eta_square + k_square) * pow(cos_theta_i, 2) - 2.0 * eta * cos_theta_i + (float) 1.0;
  Spectrum R_p_denom = (eta_square + k_square) * pow(cos_theta_i, 2) + 2.0 * eta * cos_theta_i + (float) 1.0;
  Spectrum R_p = R_p_nom / R_p_denom;
  return (R_s + R_p) / 2.0;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: 2.1
  // Implement microfacet model here

  if (wo.z < 0 || wi.z < 0) {
    return Spectrum(); //check if wo and wi are valid.
  }
  Vector3D h = (wi + wo).unit();
  return (F(wi) * G(wo, wi) * D(h)) / (4 * cos_theta(wo) * cos_theta(wi));
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 2.4
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
  //*wi = cosineHemisphereSampler.get_sample(pdf);
  //return MicrofacetBSDF::f(wo, *wi);
  Vector2D random_num = sampler.get_sample();
  double r1 = random_num[0];
  double r2 = random_num[1];
  double theta_h = atan(sqrt(-pow(alpha, 2) * log(1.0 - r1)));
  double phi_h = 2 * PI * r2;

  Vector3D h = Vector3D(cos(phi_h) * sin(theta_h), sin(theta_h) * sin(phi_h), cos(theta_h)).unit();
  *wi =  2 * dot(wo, h) * h - wo;

  if (wi->z < 0) {
    //std::cout << *pdf << std::endl;
    *pdf = 0.0;
    return Spectrum();
  }

  double p_theta = ((2.0 * sin(theta_h)) / (pow(alpha, 2.0) * pow(cos(theta_h), 3.0))) /
      exp(-(pow(tan(theta_h), 2.0) / pow(alpha, 2.0)));
  double p_phi = 1.0 / (2 * PI);
  double p_h = (p_theta * p_phi) / sin(theta_h);
  *pdf = p_h / (4.0 * dot(*wi, h));

  return MicrofacetBSDF::f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 1.4
  // Compute Fresnel coefficient and either reflect or refract based on it.
  if (!refract(wo, wi, ior)) {//Total Internal reflection
    reflect(wo, wi);
    *pdf = 1.0;
    return reflectance / abs_cos_theta(*wi);
  } else {
    double R0 = pow((1.0 - ior) / (1.0 + ior), 2);
    double R = R0 + (1.0 - R0) * pow((1 - abs_cos_theta(wo)), 5);

    if (coin_flip(R)) {
      reflect(wo, wi);
      *pdf = R;
      return R * reflectance / abs_cos_theta(*wi);
    } else {
      *pdf = 1.0 - R;
      double eta = wo.z >= 0 ? 1.0 / ior : ior;
      Spectrum result =  (1.0 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
      return result;
    }
  }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO: 1.1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 1.3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  double factor = wo.z >= 0 ? 1.0 / ior : ior;

  if (1 - factor * (1 - pow((wo.z),2)) < 0) { //Total internal reflection
    return false;
  } else {
    *wi = Vector3D(-factor * wo.x, -factor * wo.y, sqrt(1 - factor * (1 - pow((wo.z),2))));
    wi->z = wo.z >= 0 ? -wi->z : wi->z;
  }
  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
