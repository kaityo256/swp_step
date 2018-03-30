#include <cstdio>
#include <x86intrin.h>
#include <random>
#include "conf.hpp"
//----------------------------------------------------------------------
double q[N][D];
double p[N][D] = {};
int particle_number = 0;
int number_of_pairs = 0;
int number_of_partners[N];
int i_particles[MAX_PAIRS];
int j_particles[MAX_PAIRS];
int pointer[N];
int sorted_list[MAX_PAIRS];
//----------------------------------------------------------------------
class AoSDataManager : public DataManager {
private:
  double (*p)[D];
  double (*q)[D];
public:
  AoSDataManager(double _p[N][D], double _q[N][D]) {
    p = _p;
    q = _q;
  }
  void add_particle(double x, double y, double z, int &particle_number) {
    static std::mt19937 mt(2);
    std::uniform_real_distribution<double> ud(0.0, 0.1);
    q[particle_number][X] = x + ud(mt);
    q[particle_number][Y] = y + ud(mt);
    q[particle_number][Z] = z + ud(mt);
    particle_number++;
  }
  double calc_distance(int i, int j) {
    const double dx = q[i][X] - q[j][X];
    const double dy = q[i][Y] - q[j][Y];
    const double dz = q[i][Z] - q[j][Z];
    return dx * dx + dy * dy + dz * dz;
  }
  void print_results(int particle_number) {
    for (int i = 0; i < 5; i++) {
      printf("%.10f %.10f %.10f\n", p[i][X], p[i][Y], p[i][Z]);
    }
    for (int i = particle_number - 5; i < particle_number; i++) {
      printf("%.10f %.10f %.10f\n", p[i][X], p[i][Y], p[i][Z]);
    }
  }
};
//----------------------------------------------------------------------
void
force_pair(void) {
  for (int k = 0; k < number_of_pairs; k++) {
    const int i = i_particles[k];
    const int j = j_particles[k];
    double dx = q[j][X] - q[i][X];
    double dy = q[j][Y] - q[i][Y];
    double dz = q[j][Z] - q[i][Z];
    double r2 = (dx * dx + dy * dy + dz * dz);
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    if (r2 > CL2) df = 0.0;
    p[i][X] += df * dx;
    p[i][Y] += df * dy;
    p[i][Z] += df * dz;
    p[j][X] -= df * dx;
    p[j][Y] -= df * dy;
    p[j][Z] -= df * dz;
  }
}
//----------------------------------------------------------------------
void
force_sorted(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - qix;
      double dy = q[j][Y] - qiy;
      double dz = q[j][Z] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//------------------------------------------------------------------------
void
force_sorted_swp1(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    if (np == 0) continue; // 相互作用粒子ゼロなら飛ばす
    for (int k = 0; k < np; k++) {
      // ------- 8< ---------
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - qix;
      double dy = q[j][Y] - qiy;
      double dz = q[j][Z] - qiz;
      // ------- 8< ---------
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//------------------------------------------------------------------------
void
force_sorted_swp2(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    if (np == 0) continue; // 相互作用粒子ゼロなら飛ばす
    int j = sorted_list[kp];
    double dx = q[j][X] - qix;
    double dy = q[j][Y] - qiy;
    double dz = q[j][Z] - qiz;
    for (int k = 1; k < np; k++) {
      // ------- 8< ---------
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
      // ------- 8< ---------
      j = sorted_list[kp + k];
      dx = q[j][X] - qix;
      dy = q[j][Y] - qiy;
      dz = q[j][Z] - qiz;
    }
    double r2 = (dx * dx + dy * dy + dz * dz);
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    if (r2 > CL2) df = 0.0;
    pfx += df * dx;
    pfy += df * dy;
    pfz += df * dz;
    p[j][X] -= df * dx;
    p[j][Y] -= df * dy;
    p[j][Z] -= df * dz;

    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//------------------------------------------------------------------------
void
force_sorted_swp3(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    if (np == 0) continue; // 相互作用粒子ゼロなら飛ばす
    int j = sorted_list[kp];
    double dx = q[j][X] - qix;
    double dy = q[j][Y] - qiy;
    double dz = q[j][Z] - qiz;
    for (int k = 1; k < np; k++) {
      int j2 = sorted_list[kp + k];
      // ------- 8< ---------
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
      // ------- 8< ---------
      //j = sorted_list[kp + k];
      j = j2;
      dx = q[j][X] - qix;
      dy = q[j][Y] - qiy;
      dz = q[j][Z] - qiz;
    }
    double r2 = (dx * dx + dy * dy + dz * dz);
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    if (r2 > CL2) df = 0.0;
    pfx += df * dx;
    pfy += df * dy;
    pfz += df * dz;
    p[j][X] -= df * dx;
    p[j][Y] -= df * dy;
    p[j][Z] -= df * dz;

    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//------------------------------------------------------------------------
void
force_sorted_swp4(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    if (np == 0) continue; // 相互作用粒子ゼロなら飛ばす
    int j = sorted_list[kp];
    double dx = q[j][X] - qix;
    double dy = q[j][Y] - qiy;
    double dz = q[j][Z] - qiz;
    for (int k = 1; k < np; k++) {
      int j2 = sorted_list[kp + k];
      double dx2 = q[j2][X] - qix;
      double dy2 = q[j2][Y] - qiy;
      double dz2 = q[j2][Z] - qiz;
      // ------- 8< ---------
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
      // ------- 8< ---------
      //j = sorted_list[kp + k];
      //dx = q[j][X] - qix;
      //dy = q[j][Y] - qiy;
      //dz = q[j][Z] - qiz;
      j = j2;
      dx = dx2;
      dy = dy2;
      dz = dz2;
    }
    double r2 = (dx * dx + dy * dy + dz * dz);
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    if (r2 > CL2) df = 0.0;
    pfx += df * dx;
    pfy += df * dy;
    pfz += df * dz;
    p[j][X] -= df * dx;
    p[j][Y] -= df * dy;
    p[j][Z] -= df * dz;

    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//------------------------------------------------------------------------
int
main(void) {
  AoSDataManager aosdm(p, q);
  init(&aosdm, particle_number);
  check_pairlist(particle_number, number_of_pairs, number_of_partners, i_particles, j_particles, &aosdm);
  sortpair(particle_number, number_of_pairs, number_of_partners, i_particles, j_particles, pointer, sorted_list);
  //measure(&force_sorted, "sorted", particle_number);
  //measure(&force_sorted_swp1, "sorted_swp1", particle_number);
  //measure(&force_sorted_swp2, "sorted_swp2", particle_number);
  //measure(&force_sorted_swp3, "sorted_swp3", particle_number);
  measure(&force_sorted_swp4, "sorted_swp4", particle_number);
  aosdm.print_results(particle_number);
}
//----------------------------------------------------------------------
