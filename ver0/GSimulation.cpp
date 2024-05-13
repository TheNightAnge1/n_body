/*
    This file is part of the example codes which have been used
    for the "Code Optmization Workshop".

    Copyright (C) 2016  Fabio Baruffa <fbaru-dev@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "GSimulation.hpp"
#include "cpu_time.hpp"

GSimulation ::GSimulation()
{
  std::cout << "===============================" << std::endl;
  std::cout << " Initialize Gravity Simulation" << std::endl;
  set_npart(10000);
  set_nsteps(10);
  set_tstep(0.1);
  set_sfreq(1);
}

void GSimulation ::set_number_of_particles(int N)
{
  set_npart(N);
}

void GSimulation ::set_number_of_steps(int N)
{
  set_nsteps(N);
}

void GSimulation ::init_pos()
{
  std::random_device rd; // random number generator
  std::mt19937 gen(42);
  std::uniform_real_distribution<real_type> unif_d(0, 1.0);

  for (int i = 0; i < get_npart(); ++i)
  {
    particles[i].pos[0] = unif_d(gen);
    particles[i].pos[1] = unif_d(gen);
    particles[i].pos[2] = unif_d(gen);
  }
}

void GSimulation ::init_vel()
{
  std::random_device rd; // random number generator
  std::mt19937 gen(42);
  std::uniform_real_distribution<real_type> unif_d(-1.0, 1.0);

  for (int i = 0; i < get_npart(); ++i)
  {
    particles[i].vel[0] = unif_d(gen) * 1.0e-3f;
    particles[i].vel[1] = unif_d(gen) * 1.0e-3f;
    particles[i].vel[2] = unif_d(gen) * 1.0e-3f;
  }
}

void GSimulation ::init_acc()
{
  for (int i = 0; i < get_npart(); ++i)
  {
    particles[i].acc[0] = 0.f;
    particles[i].acc[1] = 0.f;
    particles[i].acc[2] = 0.f;
  }
}

void GSimulation ::init_mass()
{
  real_type n = static_cast<real_type>(get_npart());
  std::random_device rd; // random number generator
  std::mt19937 gen(42);
  std::uniform_real_distribution<real_type> unif_d(0.0, 1.0);

  for (int i = 0; i < get_npart(); ++i)
  {
    particles[i].mass = n * unif_d(gen);
  }
}

// prevents explosion in the case the particles are really close to each other
static const double softeningSquared = 1e-6; // если стоит 1e-9 - энергия системы резко подскакивает
static const double G = 6.67259e-11;

void count_impulse(Particle *particles, int num_parts, double impulse[])
{
  for (int i = 0; i < num_parts; ++i)
  {
    // Have to add
  }
}

double count_k_energy(Particle *particles, int num_parts)
{
  // Have to fix
  double energy = 0.;
  for (int i = 0; i < num_parts; ++i)
  {
    energy += particles[i].mass * (particles[i].vel[0] * particles[i].vel[0] +
                                   particles[i].vel[1] * particles[i].vel[1] +
                                   particles[i].vel[2] * particles[i].vel[2]); // 7flops
  }
  return energy;
}

double count_p_energy(Particle *particles, int num_parts)
{
  // Have to fix
  double energy = 0.;
  for (int i = 0; i < num_parts; ++i)
  {
    for (int j = 0; j < num_parts; ++j)
    {
      if (i == j)
        continue;
      double dx = particles[j].pos[0] - particles[i].pos[0]; // 1flop
      double dy = particles[j].pos[1] - particles[i].pos[1]; // 1flop
      double dz = particles[j].pos[2] - particles[i].pos[2]; // 1flop

      double distanceSqr = dx * dx + dy * dy + dz * dz + softeningSquared; // 6flops
      double distanceInv = 1.0 / sqrt(distanceSqr);                        // 1div+1sqrt
      energy += -G * particles[j].mass * distanceInv * particles[i].mass;
    }
  }
  return energy;
}

void computeAcc(Particle *src, Particle *temp, real_type step, int part_num) // part_num - число частиц
{
  for (int i = 0; i < part_num; i++) // update acceleration
  {
    src[i].acc[0] = 0.;
    src[i].acc[1] = 0.;
    src[i].acc[2] = 0.;

    for (int j = 0; j < part_num; j++)
    {
      if (i != j)
      {
        real_type dx, dy, dz;
        real_type distanceSqr = 0.0;
        real_type distanceInv = 0.0;

        dx = src[j].pos[0] - src[i].pos[0]; // 1flop
        dy = src[j].pos[1] - src[i].pos[1]; // 1flop
        dz = src[j].pos[2] - src[i].pos[2]; // 1flop

        distanceSqr = dx * dx + dy * dy + dz * dz + softeningSquared; // 6flops
        distanceInv = 1.0 / sqrt(distanceSqr);                        // 1div+1sqrt

        temp[i].acc[0] = src[i].acc[0] + dx * G * src[j].mass * distanceInv * distanceInv * distanceInv; // 6flops
        temp[i].acc[1] = src[i].acc[1] + dy * G * src[j].mass * distanceInv * distanceInv * distanceInv; // 6flops
        temp[i].acc[2] = src[i].acc[2] + dz * G * src[j].mass * distanceInv * distanceInv * distanceInv; // 6flops
      }
    }
  }
}

// void update_pos(Particle *dst, const Particle *src_1, const Particle *src_2, double coef, int part_num)
// {
//   for (int i = 0; i < part_num; ++i)
//   {
//     dst[i].pos[0] = src_1[i].pos[0] + src_2[i].vel[0] * coef; // 2 flop
//     dst[i].pos[1] = src_1[i].pos[1] + src_2[i].vel[1] * coef;
//     dst[i].pos[2] = src_1[i].pos[2] + src_2[i].vel[2] * coef;

//     dst[i].vel[0] = src_1[i].vel[0] + src_2[i].acc[0] * coef;
//     dst[i].vel[1] = src_1[i].vel[1] + src_2[i].acc[1] * coef;
//     dst[i].vel[2] = src_1[i].vel[2] + src_2[i].acc[2] * coef;
//   }
// }

void updatePos(Particle *dist, Particle *source, Particle *newparts, real_type step, int part_num)
{

  for (int i = 0; i < part_num; ++i) // update position and velocity
  {

    newparts[i].vel[0] += newparts[i].acc[0] * step;
    newparts[i].vel[1] += newparts[i].acc[1] * step;
    newparts[i].vel[2] += newparts[i].acc[2] * step;

    dist[i].pos[0] = source[i].pos[0] + newparts[i].vel[0] * step;
    dist[i].pos[1] = source[i].pos[1] + newparts[i].vel[1] * step;
    dist[i].pos[2] = source[i].pos[2] + newparts[i].vel[2] * step;

    dist[i].vel[0] = source[i].vel[0] + newparts[i].acc[0] * step;
    dist[i].vel[1] = source[i].vel[1] + newparts[i].acc[1] * step;
    dist[i].vel[2] = source[i].vel[2] + newparts[i].acc[2] * step;
  }
}

void GSimulation ::start()
{
  real_type energy_k, energy_p;
  real_type dt = get_tstep();
  int n = get_npart(); // n - число частиц
  int i, j;

  // allocate particles
  particles = new Particle[n];

  init_pos();
  init_vel();
  init_acc();
  init_mass();

  std::cout << "Initial system energy k: " << count_k_energy(particles, n) << " p:" << count_p_energy(particles, n) << " Sum: " << count_k_energy(particles, n) + count_p_energy(particles, n) << std::endl;

  print_header();

  _totTime = 0.;

  CPUTime time;
  double ts0 = 0;
  double ts1 = 0;
  double nd = double(n);
  double gflops = 1e-9 * ((11. + 18.) * nd * nd + nd * 19.);
  double av = 0.0, dev = 0.0;
  int nf = 0;

  const double t0 = time.start();
  for (int s = 1; s <= get_nsteps(); ++s) // s - номер итерации
  {
    ts0 += time.start();

    // Advanced Euler Method
    Particle *tmp = new Particle[n];
    Particle *k1 = new Particle[n];
    computeAcc(particles, k1, dt / 2, n);
    updatePos(tmp, particles, k1, dt / 2, n);

    Particle *k2 = new Particle[n];
    computeAcc(particles, k2, dt / 2, n);
    updatePos(tmp, particles, k2, dt / 2, n);

    updatePos(particles, particles, k1, dt / 6, n);
    updatePos(particles, particles, k2, dt / 3, n);

    delete[] k1;
    delete[] k2;
    delete[] tmp;

    energy_k = count_k_energy(particles, n);
    energy_p = count_p_energy(particles, n);

    double impulse[3] = {};
    count_impulse(particles, n, impulse);

    double impulse_len = sqrt(impulse[0] * impulse[0] + impulse[1] * impulse[1] + impulse[2] * impulse[2]);

    ts1 += time.stop();
    if (!(s % get_sfreq())) // если остаток от деления s на sfreq == 0
    {
      nf += 1;
      std::cout << " "
                << std::left << std::setw(8) << s
                << std::left << std::setprecision(5) << std::setw(8) << s * get_tstep()
                << std::left << std::setprecision(9) << std::setw(12) << energy_p + energy_k
                << std::left << std::setprecision(9) << std::setw(12) << impulse_len
                << std::left << std::setprecision(5) << std::setw(12) << (ts1 - ts0)
                << std::endl;
      if (nf > 2)
      {
        av += gflops * get_sfreq() / (ts1 - ts0);
        dev += gflops * get_sfreq() * gflops * get_sfreq() / ((ts1 - ts0) * (ts1 - ts0));
      }
      ts0 = 0;
      ts1 = 0;
    }
  } // end of the time step loop

  const double t1 = time.stop();
  _totTime = (t1 - t0);
  _totFlops = gflops * get_nsteps();

  av /= (double)(nf - 2);
  dev = sqrt(dev / (double)(nf - 2) - av * av);

  int nthreads = 1;

  std::cout << std::endl;
  std::cout << "# Number Threads     : " << nthreads << std::endl;
  std::cout << "# Total Time (s)     : " << _totTime << std::endl;
  std::cout << "# Average Perfomance : " << av << " +- " << dev << std::endl;
  std::cout << "===============================" << std::endl;
}

void GSimulation ::print_header()
{

  std::cout << " nPart = " << get_npart() << "; "
            << "nSteps = " << get_nsteps() << "; "
            << "dt = " << get_tstep() << std::endl;

  std::cout << "------------------------------------------------" << std::endl;
  std::cout << " "
            << std::left << std::setw(8) << "s"
            << std::left << std::setw(8) << "dt"
            << std::left << std::setw(12) << "s_energy"
            << std::left << std::setw(12) << "impulse"
            << std::left << std::setw(12) << "time (s)"
            << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
}

GSimulation ::~GSimulation()
{
  delete particles;
}
