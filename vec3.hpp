#ifndef VEC3_HPP
#define VEC3_HPP

// #define USE_OMP

/*@T
 * \subsection{Vector operations}
 * 
 * A few inline functions save us from typing the same annoying bits
 * of code repeatedly.
 *
 *@c*/

#ifndef USE_OMP
#pragma omp declare simd
#endif
inline void vec3_set(float* result, float x, float y, float z)
{
    result[0] = x;
    result[1] = y;
    result[2] = z;
}

#ifndef USE_OMP
#pragma omp declare simd
#endif
inline void vec3_copy(float* result, float* v)
{
    vec3_set(result, v[0], v[1], v[2]);
}

#ifndef USE_OMP
#pragma omp declare simd
#endif
inline void vec3_diff(float* result, float* a, float* b)
{
    result[0] = a[0]-b[0];
    result[1] = a[1]-b[1];
    result[2] = a[2]-b[2];
}

#ifndef USE_OMP
#pragma omp declare simd
#endif
inline void vec3_scale(float* result, float alpha, float* v)
{
    result[0] = alpha*v[0];
    result[1] = alpha*v[1];
    result[2] = alpha*v[2];
}

#ifndef USE_OMP
#pragma omp declare simd
#endif
inline float vec3_dist2(float* a, float* b)
{
    float dx = a[0]-b[0];
    float dy = a[1]-b[1];
    float dz = a[2]-b[2];
    return dx*dx + dy*dy + dz*dz;
}

#ifndef USE_OMP
#pragma omp declare simd
#endif
inline float vec3_len2(float* a)
{
    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
}

#ifndef USE_OMP
#pragma omp declare simd
#endif
inline void vec3_saxpy(float* result, float alpha, float* v)
{
    result[0] += alpha*v[0];
    result[1] += alpha*v[1];
    result[2] += alpha*v[2];
}

#ifndef USE_OMP
#pragma omp declare simd
#endif
inline void vec3_scalev(float* result, float alpha)
{
    result[0] *= alpha;
    result[1] *= alpha;
    result[2] *= alpha;
}

/*@q*/
#endif /* VEC3_HPP */
