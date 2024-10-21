#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <vector>
#include <iostream>

#include "vec3.hpp"
#include "zmorton.hpp"

#include "params.hpp"
#include "state.hpp"
#include "interact.hpp"
#include "binhash.hpp"

/* Define this to use the bucketing version of the code */
#define USE_BUCKETING
// #define USE_2H
#define USE_OMP

extern int particle_neighbour_map[3380][27];
extern omp_lock_t locks[3380][3380];
extern bool particle_updated[3380][3380];

/*@T
 * \subsection{Density computations}
 *
 * The formula for density is
 * \[
 *   \rho_i = \sum_j m_j W_{p6}(r_i-r_j,h)
 *          = \frac{315 m}{64 \pi h^9} \sum_{j \in N_i} (h^2 - r^2)^3.
 * \]
 * We search for neighbors of node $i$ by checking every particle,
 * which is not very efficient.  We do at least take advange of
 * the symmetry of the update ($i$ contributes to $j$ in the same
 * way that $j$ contributes to $i$).
 *@c*/

inline void update_density(particle_t *pi, particle_t *pj, float h2, float C)
{
    float r2 = vec3_dist2(pi->x, pj->x);
    float z = h2 - r2;
    if (z > 0)
    {
        float rho_ij = C * z * z * z;
        pi->rho += rho_ij;
        pj->rho += rho_ij;
    }
}

inline void update_density_single(particle_t *pi, particle_t *pj, float h2, float C)
{
    if (pi == pj)
    {
        return;
    }
    float r2 = vec3_dist2(pi->x, pj->x);
    float z = h2 - r2;
    if (z > 0)
    {
        float rho_ij = C * z * z * z;
#ifdef USE_OMP
#pragma omp atomic
#endif
        pi->rho += rho_ij;
        // pj->rho += rho_ij;
    }
}

void compute_density(sim_state_t *s, sim_param_t *params)
{
    int n = s->n;
    particle_t *p = s->part;
    particle_t **hash = s->hash;

    float h = params->h;
    float h2 = h * h;
    float h3 = h2 * h;
    float h9 = h3 * h3 * h3;
    float C = (315.0 / 64.0 / M_PI) * s->mass / h9;

    // Clear densities
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; ++i)
        p[i].rho = 0;

        // Accumulate density info
#ifdef USE_BUCKETING
    /* BEGIN TASK */
    // std::cout<<"in com"<<std::endl;
    // #ifdef USE_OMP
    // omp_set_num_threads(16);
    // #pragma omp parallel for schedule(guided)
    // #endif

    unsigned buckets[27];
#pragma omp parallel for private(buckets) schedule(guided, 4)
    for (int i = 0; i < n; ++i)
    {
        // std::cout<<"in for"<< i <<std::endl;
        particle_t *pi = p + i;
        pi->rho += (315.0 / 64.0 / M_PI) * s->mass / h3;

#ifdef USE_2H
        unsigned num = particle_neighborhood(buckets, pi, 2 * h);
#else
        unsigned num = particle_neighborhood(buckets, pi, h, i);
#endif

        // float local_rho = 0.0f;
        // #pragma omp parallel for reduction(+:local_rho)
        for (size_t j = 0; j < num; j++)
        {
            particle_t *p_n = hash[particle_neighbour_map[i][j]];
            while (p_n)
            {
                if (particle_updated[i][p_n->index] == false)
                {
                    if (pi != p_n)
                    {
                        float r2 = vec3_dist2(pi->x, p_n->x);
                        float z = h2 - r2;
                        if (z > 0)
                        {
                            float rho_ij = C * z * z * z;
                            int li = (i > p_n->index) ? p_n->index : i;
                            int lj = (i > p_n->index) ? i : p_n->index;
                            omp_set_lock(&locks[li][lj]);
                            if (particle_updated[i][p_n->index] == false)
                            {
                                pi->rho += rho_ij;
                                p_n->rho += rho_ij;
                                particle_updated[p_n->index][i] = true;
                            }
                            omp_unset_lock(&locks[li][lj]);
                        }
                    }
                }
                p_n = p_n->next;
            }
        }
    }
    /* END TASK */
#else
    for (int i = 0; i < n; ++i)
    {
        particle_t *pi = s->part + i;
        pi->rho += (315.0 / 64.0 / M_PI) * s->mass / h3;
        for (int j = i + 1; j < n; ++j)
        {
            particle_t *pj = s->part + j;
            update_density(pi, pj, h2, C);
        }
    }
#endif
}

/*@T
 * \subsection{Computing forces}
 *
 * The acceleration is computed by the rule
 * \[
 *   \bfa_i = \frac{1}{\rho_i} \sum_{j \in N_i}
 *     \bff_{ij}^{\mathrm{interact}} + \bfg,
 * \]
 * where the pair interaction formula is as previously described.
 * Like [[compute_density]], the [[compute_accel]] routine takes
 * advantage of the symmetry of the interaction forces
 * ($\bff_{ij}^{\mathrm{interact}} = -\bff_{ji}^{\mathrm{interact}}$)
 * but it does a very expensive brute force search for neighbors.
 *@c*/

inline void update_forces(particle_t *pi, particle_t *pj, float h2,
                          float rho0, float C0, float Cp, float Cv)
{
    float dx[3];
    vec3_diff(dx, pi->x, pj->x);
    float r2 = vec3_len2(dx);
    if (r2 < h2)
    {
        const float rhoi = pi->rho;
        const float rhoj = pj->rho;
        float q = sqrt(r2 / h2);
        float u = 1 - q;
        float w0 = C0 * u / rhoi / rhoj;
        float wp = w0 * Cp * (rhoi + rhoj - 2 * rho0) * u / q;
        float wv = w0 * Cv;
        float dv[3];
        vec3_diff(dv, pi->v, pj->v);

        // Equal and opposite pressure forces
        vec3_saxpy(pi->a, wp, dx);
        vec3_saxpy(pj->a, -wp, dx);

        // Equal and opposite viscosity forces
        vec3_saxpy(pi->a, wv, dv);
        vec3_saxpy(pj->a, -wv, dv);
    }
}

inline void update_forces_single(particle_t *pi, particle_t *pj, float h2,
                                 float rho0, float C0, float Cp, float Cv)
{
    if (pi == pj)
    {
        return;
    }
    float dx[3];
    vec3_diff(dx, pi->x, pj->x);
    float r2 = vec3_len2(dx);
    if (r2 < h2)
    {
        const float rhoi = pi->rho;
        const float rhoj = pj->rho;
        float q = sqrt(r2 / h2);
        float u = 1 - q;
        float w0 = C0 * u / rhoi / rhoj;
        float wp = w0 * Cp * (rhoi + rhoj - 2 * rho0) * u / q;
        float wv = w0 * Cv;
        float dv[3];
        vec3_diff(dv, pi->v, pj->v);

        // Equal and opposite pressure forces
        vec3_saxpy(pi->a, wp, dx);
        // vec3_saxpy(pj->a, -wp, dx);

        // Equal and opposite viscosity forces
        vec3_saxpy(pi->a, wv, dv);
        // vec3_saxpy(pj->a, -wv, dv);
    }
}

void compute_accel(sim_state_t *state, sim_param_t *params)
{
    // Unpack basic parameters
    const float h = params->h;
    const float rho0 = params->rho0;
    const float k = params->k;
    const float mu = params->mu;
    const float g = params->g;
    const float mass = state->mass;
    const float h2 = h * h;

    // Unpack system state
    particle_t *p = state->part;
    particle_t **hash = state->hash;
    int n = state->n;

    // Rehash the particles
#ifdef USE_2H
    hash_particles(state, 2 * h);
#else
    hash_particles(state, h);
#endif

// Compute density and color
#pragma omp parallel for
    for (int i = 0; i < 3380; ++i)
    {
        memset(particle_updated[i], 0, sizeof(bool) * 3380);
    }

    compute_density(state, params);
// #pragma omp parallel for
//     for (int i = 0; i < 3380; ++i)
//     {
//         memset(particle_updated[i], 0, sizeof(bool) * 3380);
//     }
// memset(particle_updated, 0, sizeof(particle_updated));
// Start with gravity and surface forces
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; ++i)
        vec3_set(p[i].a, 0, -g, 0);

    // Constants for interaction term
    float C0 = 45 * mass / M_PI / ((h2) * (h2)*h);
    float Cp = k / 2;
    float Cv = -mu;

    // Accumulate forces
#ifdef USE_BUCKETING
    /* BEGIN TASK */
#ifdef USE_OMP
    // omp_set_num_threads(16);
#pragma omp parallel for schedule(guided, 4)
#endif
    for (int i = 0; i < n; ++i)
    {
        particle_t *pi = p + i;

        for (size_t j = 0; j < pi->num; j++)
        {

            particle_t *p_n = hash[particle_neighbour_map[i][j]];
            while (p_n)
            {

                update_forces_single(pi, p_n, h2, rho0, C0, Cp, Cv);

                p_n = p_n->next;
            }
        }
        // delete[] buckets;
    }
    /* END TASK */
#else
    for (int i = 0; i < n; ++i)
    {
        particle_t *pi = p + i;
        for (int j = i + 1; j < n; ++j)
        {
            particle_t *pj = p + j;
            update_forces(pi, pj, h2, rho0, C0, Cp, Cv);
        }
    }
#endif
}
