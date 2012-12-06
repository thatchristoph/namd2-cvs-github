/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

//#include "InfoStream.h"
#include "Settle.h"
#include <string.h>
#include <math.h>
//#include <charm++.h> // for CkPrintf

#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE)
#include <emmintrin.h>  // SSE2
#if defined(__INTEL_COMPILER)
#define __align(X) __declspec(align(X) )
#elif defined(__PGI)
#define __align(X)  __attribute__((aligned(X) ))
#define MISSING_mm_cvtsd_f64
#elif defined(__GNUC__)
#define __align(X)  __attribute__((aligned(X) ))
#if (__GNUC__ < 4)
#define MISSING_mm_cvtsd_f64
#endif
#else
#define __align(X) __declspec(align(X) )
#endif
#endif

//
// XXX static and global variables are unsafe for shared memory builds.
// The global and static vars should be eliminated.
// Unfortunately, the routines that use these below are actually
// in use in NAMD.
//

// Initialize various properties of the waters
// settle1() assumes all waters are identical, 
// and will generate bad results if they are not.
void settle1init(BigReal pmO, BigReal pmH, BigReal hhdist, BigReal ohdist,
                 BigReal &mOrmT, BigReal &mHrmT, BigReal &ra,
                 BigReal &rb, BigReal &rc, BigReal &rra) {

    BigReal rmT = 1.0 / (pmO+pmH+pmH);
    mOrmT = pmO * rmT;
    mHrmT = pmH * rmT;
    BigReal t1 = 0.5*pmO/pmH;
    rc = 0.5*hhdist;
    ra = sqrt(ohdist*ohdist-rc*rc)/(1.0+t1);
    rb = t1*ra;
    rra = 1.0 / ra;
}


int settle1(const Vector *ref, Vector *pos, Vector *vel, BigReal invdt,
                 BigReal mOrmT, BigReal mHrmT, BigReal ra,
                 BigReal rb, BigReal rc, BigReal rra) {
#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE)
  // SSE acceleration of some of the costly parts of settle using
  // the Intel C/C++ classes.  This implementation uses the SSE units
  // less efficiency than is potentially possible, but in order to do
  // better, the settle algorithm will have to be vectorized and operate
  // on multiple waters at a time.  Doing so could give us the ability to
  // do two (double precison) or four (single precision) waters at a time.
  // This code achieves a modest speedup without the need to reorganize
  // the NAMD structure.  Once we have water molecules sorted in a single
  // block we can do far better.

  // vectors in the plane of the original positions
  Vector b0, c0;

  __m128d REF0xy = _mm_loadu_pd((double *) &ref[0].x);  // ref0.y and ref0.x
  __m128d REF1xy = _mm_loadu_pd((double *) &ref[1].x);  // ref1.y and ref1.x

  __m128d B0xy = _mm_sub_pd(REF1xy, REF0xy);
  _mm_storeu_pd((double *) &b0.x, B0xy);
  b0.z = ref[1].z - ref[0].z;

  __m128d REF2xy = _mm_loadu_pd((double *) &ref[2].x);  // ref2.y and ref2.x

  __m128d C0xy = _mm_sub_pd(REF2xy, REF0xy);
  _mm_storeu_pd((double *) &c0.x, C0xy);
  c0.z = ref[2].z - ref[0].z;

  // new center of mass
  // Vector d0 = pos[0] * mOrmT + ((pos[1] + pos[2]) * mHrmT);
  __align(16) Vector a1;
  __align(16) Vector b1;
  __align(16) Vector c1;
  __align(16) Vector d0;

  __m128d POS1xy = _mm_loadu_pd((double *) &pos[1].x);
  __m128d POS2xy = _mm_loadu_pd((double *) &pos[2].x);
  __m128d PMHrmTxy = _mm_mul_pd(_mm_add_pd(POS1xy, POS2xy), _mm_set1_pd(mHrmT));

  __m128d POS0xy = _mm_loadu_pd((double *) &pos[0].x);
  __m128d PMOrmTxy = _mm_mul_pd(POS0xy, _mm_set1_pd(mOrmT));
  __m128d D0xy = _mm_add_pd(PMOrmTxy, PMHrmTxy);

  d0.z = pos[0].z * mOrmT + ((pos[1].z + pos[2].z) * mHrmT);
  a1.z = pos[0].z - d0.z;
  b1.z = pos[1].z - d0.z;
  c1.z = pos[2].z - d0.z;

  __m128d A1xy = _mm_sub_pd(POS0xy, D0xy);
  _mm_store_pd((double *) &a1.x, A1xy); // must be aligned

  __m128d B1xy = _mm_sub_pd(POS1xy, D0xy);
  _mm_store_pd((double *) &b1.x, B1xy); // must be aligned

  __m128d C1xy = _mm_sub_pd(POS2xy, D0xy);
  _mm_store_pd((double *) &c1.x, C1xy); // must be aligned

  _mm_store_pd((double *) &d0.x, D0xy); // must be aligned
  
  // Vectors describing transformation from original coordinate system to
  // the 'primed' coordinate system as in the diagram.  
  Vector n0 = cross(b0, c0);
  Vector n1 = cross(a1, n0); 
  Vector n2 = cross(n0, n1); 
#else
  // vectors in the plane of the original positions
  Vector b0 = ref[1]-ref[0];
  Vector c0 = ref[2]-ref[0];
  
  // new center of mass
  Vector d0 = pos[0]*mOrmT + ((pos[1] + pos[2])*mHrmT);
 
  Vector a1 = pos[0] - d0;
  Vector b1 = pos[1] - d0;
  Vector c1 = pos[2] - d0;
  
  // Vectors describing transformation from original coordinate system to
  // the 'primed' coordinate system as in the diagram.  
  Vector n0 = cross(b0, c0);
  Vector n1 = cross(a1, n0); 
  Vector n2 = cross(n0, n1); 
#endif

#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE) && ! defined(MISSING_mm_cvtsd_f64)
  __m128d l1 = _mm_set_pd(n0.x, n0.y);
  l1 = _mm_mul_pd(l1, l1);
  // n0.x^2 + n0.y^2
  double l1xy0 = _mm_cvtsd_f64(_mm_add_sd(l1, _mm_shuffle_pd(l1, l1, 1)));

  __m128d l3 = _mm_set_pd(n1.y, n1.z);
  l3 = _mm_mul_pd(l3, l3);
  // n1.y^2 + n1.z^2
  double l3yz1 = _mm_cvtsd_f64(_mm_add_sd(l3, _mm_shuffle_pd(l3, l3, 1)));

  __m128d l2 = _mm_set_pd(n1.x, n0.z);
  // len(n1)^2 and len(n0)^2 
  __m128d ts01 = _mm_add_pd(_mm_set_pd(l3yz1, l1xy0), _mm_mul_pd(l2, l2));

  __m128d l4 = _mm_set_pd(n2.x, n2.y);
  l4 = _mm_mul_pd(l4, l4);
  // n2.x^2 + n2.y^2
  double l4xy2 = _mm_cvtsd_f64(_mm_add_sd(l4, _mm_shuffle_pd(l4, l4, 1)));
  double ts2 = l4xy2 + (n2.z * n2.z);              // len(n2)^2

  double  invlens[4];
  // since rsqrt_nr() doesn't work with current compiler
  // this is the next best option 
  static const __m128d fvecd1p0 = _mm_set1_pd(1.0);

  // 1/len(n1) and 1/len(n0)
  __m128d invlen12 = _mm_div_pd(fvecd1p0, _mm_sqrt_pd(ts01));

  // invlens[0]=1/len(n0), invlens[1]=1/len(n1)
  _mm_storeu_pd(invlens, invlen12);

  n0 = n0 * invlens[0];

  // shuffle the order of operations around from the normal algorithm so
  // that we can double pump sqrt() with n2 and cosphi at the same time
  // these components are usually computed down in the canonical water block
  BigReal A1Z = n0 * a1;
  BigReal sinphi = A1Z * rra;
  BigReal tmp = 1.0-sinphi*sinphi;

  __m128d n2cosphi = _mm_sqrt_pd(_mm_set_pd(tmp, ts2));
  // invlens[2] = 1/len(n2), invlens[3] = cosphi
  _mm_storeu_pd(invlens+2, n2cosphi);

  n1 = n1 * invlens[1];
  n2 = n2 * (1.0 / invlens[2]);
  BigReal cosphi = invlens[3];

  b0 = Vector(n1*b0, n2*b0, n0*b0); // note: b0.z is never referenced again
  c0 = Vector(n1*c0, n2*c0, n0*c0); // note: c0.z is never referenced again
 
  b1 = Vector(n1*b1, n2*b1, n0*b1);
  c1 = Vector(n1*c1, n2*c1, n0*c1);

  // now we can compute positions of canonical water 
  BigReal sinpsi = (b1.z - c1.z)/(2.0*rc*cosphi);
  tmp = 1.0-sinpsi*sinpsi;
  BigReal cospsi = sqrt(tmp);
#else
  n0 = n0.unit();
  n1 = n1.unit();
  n2 = n2.unit();

  b0 = Vector(n1*b0, n2*b0, n0*b0); // note: b0.z is never referenced again
  c0 = Vector(n1*c0, n2*c0, n0*c0); // note: c0.z is never referenced again
 
  BigReal A1Z = n0 * a1;
  b1 = Vector(n1*b1, n2*b1, n0*b1);
  c1 = Vector(n1*c1, n2*c1, n0*c1);

  // now we can compute positions of canonical water 
  BigReal sinphi = A1Z * rra;
  BigReal tmp = 1.0-sinphi*sinphi;
  BigReal cosphi = sqrt(tmp);
  BigReal sinpsi = (b1.z - c1.z)/(2.0*rc*cosphi);
  tmp = 1.0-sinpsi*sinpsi;
  BigReal cospsi = sqrt(tmp);
#endif

  BigReal rbphi = -rb*cosphi;
  BigReal tmp1 = rc*sinpsi*sinphi;
  BigReal tmp2 = rc*sinpsi*cosphi;
 
  Vector a2(0, ra*cosphi, ra*sinphi);
  Vector b2(-rc*cospsi, rbphi - tmp1, -rb*sinphi + tmp2);
  Vector c2( rc*cosphi, rbphi + tmp1, -rb*sinphi - tmp2);

  // there are no a0 terms because we've already subtracted the term off 
  // when we first defined b0 and c0.
  BigReal alpha = b2.x*(b0.x - c0.x) + b0.y*b2.y + c0.y*c2.y;
  BigReal beta  = b2.x*(c0.y - b0.y) + b0.x*b2.y + c0.x*c2.y;
  BigReal gama  = b0.x*b1.y - b1.x*b0.y + c0.x*c1.y - c1.x*c0.y;
 
  BigReal a2b2 = alpha*alpha + beta*beta;
  BigReal sintheta = (alpha*gama - beta*sqrt(a2b2 - gama*gama))/a2b2;
  BigReal costheta = sqrt(1.0 - sintheta*sintheta);
  
#if 0
  Vector a3( -a2.y*sintheta, 
              a2.y*costheta,
              a2.z);
  Vector b3(b2.x*costheta - b2.y*sintheta,
              b2.x*sintheta + b2.y*costheta,
              b2.z);
  Vector c3(c2.x*costheta - c2.y*sintheta,
              c2.x*sintheta + c2.y*costheta,
              c2.z);
  
#else
  Vector a3( -a2.y*sintheta, 
              a2.y*costheta,
              A1Z);
  Vector b3(b2.x*costheta - b2.y*sintheta,
              b2.x*sintheta + b2.y*costheta,
              b1.z);
  Vector c3(-b2.x*costheta - c2.y*sintheta,
            -b2.x*sintheta + c2.y*costheta,
              c1.z);

#endif

  // undo the transformation; generate new normal vectors from the transpose.
  Vector m1(n1.x, n2.x, n0.x);
  Vector m2(n1.y, n2.y, n0.y);
  Vector m0(n1.z, n2.z, n0.z);

  pos[0] = Vector(a3*m1, a3*m2, a3*m0) + d0;
  pos[1] = Vector(b3*m1, b3*m2, b3*m0) + d0;
  pos[2] = Vector(c3*m1, c3*m2, c3*m0) + d0;

  // dt can be negative during startup!
  if (invdt != 0) {
    vel[0] = (pos[0]-ref[0])*invdt;
    vel[1] = (pos[1]-ref[1])*invdt;
    vel[2] = (pos[2]-ref[2])*invdt;
  }

  return 0;
}



static int settlev(const Vector *pos, BigReal ma, BigReal mb, Vector *vel,
				   BigReal dt, Tensor *virial) {
  
  Vector rAB = pos[1]-pos[0];
  Vector rBC = pos[2]-pos[1];
  Vector rCA = pos[0]-pos[2];
 
  Vector AB = rAB.unit();
  Vector BC = rBC.unit();
  Vector CA = rCA.unit();
  
  BigReal cosA = -AB * CA;
  BigReal cosB = -BC * AB;
  BigReal cosC = -CA * BC;

  BigReal vab = (vel[1]-vel[0])*AB;
  BigReal vbc = (vel[2]-vel[1])*BC;
  BigReal vca = (vel[0]-vel[2])*CA;

  BigReal mab = ma+mb;
  
  BigReal d = (2*mab*mab + 2*ma*mb*cosA*cosB*cosC - 2*mb*mb*cosA*cosA
               - ma*mab*(cosB*cosB + cosC*cosC))*0.5/mb;

  BigReal tab = (vab*(2*mab - ma*cosC*cosC) +
                vbc*(mb*cosC*cosA - mab*cosB) +
                vca*(ma*cosB*cosC - 2*mb*cosA))*ma/d;
            
  BigReal tbc = (vbc*(mab*mab - mb*mb*cosA*cosA) +
                vca*ma*(mb*cosA*cosB - mab*cosC) +
                vab*ma*(mb*cosC*cosA - mab*cosB))/d;
  
  BigReal tca = (vca*(2*mab - ma*cosB*cosB) +
                vab*(ma*cosB*cosC - 2*mb*cosA) +
                vbc*(mb*cosA*cosB - mab*cosC))*ma/d;
 
  Vector ga = tab*AB - tca*CA;
  Vector gb = tbc*BC - tab*AB;
  Vector gc = tca*CA - tbc*BC;
#if 0
  if (virial) {
    *virial += 0.5*outer(tab, rAB)/dt;
    *virial += 0.5*outer(tbc, rBC)/dt;
    *virial += 0.5*outer(tca, rCA)/dt;
  }
#endif
  vel[0] += (0.5/ma)*ga;
  vel[1] += (0.5/mb)*gb;
  vel[2] += (0.5/mb)*gc;

  return 0;
}


int settle2(BigReal mO, BigReal mH, const Vector *pos,
                  Vector *vel, BigReal dt, Tensor *virial) {

  settlev(pos, mO, mH, vel, dt, virial);
  return 0;
}

