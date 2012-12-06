
#ifdef NAMD_CUDA

#define NAME(X) SLOWNAME( X )

#undef SLOW
#undef SLOWNAME
#ifdef DO_SLOW
#define SLOW(X) X
#define SLOWNAME(X) ENERGYNAME( X ## _slow )
#else
#define SLOW(X)
#define SLOWNAME(X) ENERGYNAME( X )
#endif

#undef ENERGY
#undef ENERGYNAME
#ifdef DO_ENERGY
#define ENERGY(X) X
#define ENERGYNAME(X) PAIRLISTNAME( X ## _energy )
#else
#define ENERGY(X)
#define ENERGYNAME(X) PAIRLISTNAME( X )
#endif

#undef GENPAIRLIST
#undef USEPAIRLIST
#undef PAIRLISTNAME
#ifdef MAKE_PAIRLIST
#define GENPAIRLIST(X) X
#define USEPAIRLIST(X)
#define PAIRLISTNAME(X) LAST( X ## _pairlist )
#else
#define GENPAIRLIST(X)
#define USEPAIRLIST(X) X
#define PAIRLISTNAME(X) LAST( X )
#endif

#define LAST(X) X

__device__ __forceinline__ static void NAME(dev_sum_forces)(
        const int force_list_index,
	const atom *atoms,
	const force_list *force_lists,
	const float4 *force_buffers,
	const float *virial_buffers,
	float4 *forces, float *virials);

__global__ static void NAME(dev_nonbonded)(
	const patch_pair *patch_pairs,
	const atom *atoms,
	const atom_param *atom_params,
	float4 *force_buffers,
	float4 *slow_force_buffers,
	unsigned int *block_flags,
	float *virial_buffers,
	float *slow_virial_buffers,
        const unsigned int *overflow_exclusions,
        unsigned int *force_list_counters,
        const force_list *force_lists,
        float4 *forces, float *virials,
        float4 *slow_forces, float *slow_virials,
        int lj_table_size,
        float3 lata, float3 latb, float3 latc,
	float cutoff2, float plcutoff2, int doSlow) {
// call with one block per patch_pair
// call with BLOCK_SIZE threads per block
// call with no shared memory

#ifdef __DEVICE_EMULATION__
  #define myPatchPair (*(patch_pair*)(&pp.i))
#else
  #define myPatchPair pp.pp
#endif
  __shared__ union {
#ifndef __DEVICE_EMULATION__
    patch_pair pp;
#endif
    unsigned int i[PATCH_PAIR_SIZE];
  } pp;

 { // start of nonbonded calc

  #define pl plu.c
  __shared__ union {
    unsigned int i[BLOCK_SIZE];
    char c[4*BLOCK_SIZE];
  } plu;

  __shared__ volatile union {
    float a2d[32][3];
    float a1d[32*3];
  } sumf;

  __shared__ volatile union {
    float a2d[32][3];
    float a1d[32*3];
  } sumf_slow;

#ifdef __DEVICE_EMULATION__
  #define jpqs ((atom*)(jpqu.i))
#else
  #define jpqs jpqu.d
#endif
  __shared__ union {
#ifndef __DEVICE_EMULATION__
    atom d[SHARED_SIZE];
#endif
    unsigned int i[4*SHARED_SIZE];
    float f[4*SHARED_SIZE];
  } jpqu;

#ifdef __DEVICE_EMULATION__
  #define japs ((atom_param*)(japu.i))
#else
  #define japs japu.d
#endif
  __shared__ union {
#ifndef __DEVICE_EMULATION__
    atom_param d[SHARED_SIZE];
#endif
    unsigned int i[4*SHARED_SIZE];
  } japu;

  if ( threadIdx.x < PATCH_PAIR_USED ) {
    unsigned int tmp = ((unsigned int*)patch_pairs)[
			PATCH_PAIR_SIZE*blockIdx.x+threadIdx.x];
    pp.i[threadIdx.x] = tmp;
  }

  if ( threadIdx.x < 96 ) { // initialize net force in shared memory
    sumf.a1d[threadIdx.x] = 0.f;
    sumf_slow.a1d[threadIdx.x] = 0.f;
  }

  __syncthreads();

  // convert scaled offset with current lattice
  if ( threadIdx.x == 0 ) {
    float offx = myPatchPair.offset.x * lata.x
               + myPatchPair.offset.y * latb.x
               + myPatchPair.offset.z * latc.x;
    float offy = myPatchPair.offset.x * lata.y
               + myPatchPair.offset.y * latb.y
               + myPatchPair.offset.z * latc.y;
    float offz = myPatchPair.offset.x * lata.z
               + myPatchPair.offset.y * latb.z
               + myPatchPair.offset.z * latc.z;
    myPatchPair.offset.x = offx;
    myPatchPair.offset.y = offy;
    myPatchPair.offset.z = offz;
  }

  __syncthreads();

  ENERGY(
  float totalev = 0.f;
  float totalee = 0.f;
  SLOW( float totales = 0.f; )
  )

  for ( int blocki = 0;
        blocki < myPatchPair.patch1_force_size;
        blocki += BLOCK_SIZE ) {

  atom ipq;
  struct {
    int vdw_type;
    int index; } iap;

  // load patch 1
  if ( blocki + threadIdx.x < myPatchPair.patch1_force_size ) {
    int i = myPatchPair.patch1_atom_start + blocki + threadIdx.x;
    float4 tmpa = ((float4*)atoms)[i];

    ipq.position.x = tmpa.x + myPatchPair.offset.x;
    ipq.position.y = tmpa.y + myPatchPair.offset.y;
    ipq.position.z = tmpa.z + myPatchPair.offset.z;
    ipq.charge = tmpa.w;

    uint4 tmpap = ((uint4*)atom_params)[i];

    iap.vdw_type = tmpap.x;
    iap.index = tmpap.y;
  }

  // avoid syncs by having all warps load pairlist
  USEPAIRLIST( {
    int i_pl = (blocki >> 2) + myPatchPair.block_flags_start;
    plu.i[threadIdx.x] = block_flags[i_pl + (threadIdx.x & 31)];
  } )
  GENPAIRLIST(
    plu.i[threadIdx.x] = 0;
  )
  int pli = 4 * ( threadIdx.x & 96 );

  float4 ife, ife_slow;
  ife.x = 0.f;
  ife.y = 0.f;
  ife.z = 0.f;
  ife.w = 0.f;
  ife_slow.x = 0.f;
  ife_slow.y = 0.f;
  ife_slow.z = 0.f;
  ife_slow.w = 0.f;

  for ( int blockj = 0;
        blockj < myPatchPair.patch2_size;
        blockj += SHARED_SIZE, ++pli ) {

#ifdef __DEVICE_EMULATION__
  USEPAIRLIST( if ( threadIdx.x == 0 ) printf("%d %d %d %d %d %d %d\n", blockIdx.x, blocki, blockj, pli, pl[pli], (pli+128)&255, pl[(pli+128)&255]); )
#endif
  USEPAIRLIST( if ( pl[pli] == 0 ) continue; )

  int shared_size = myPatchPair.patch2_size - blockj;
  if ( shared_size > SHARED_SIZE ) shared_size = SHARED_SIZE;

  // load patch 2
  __syncthreads();

  if ( threadIdx.x < 4 * shared_size ) {
    int j = myPatchPair.patch2_atom_start + blockj;
    jpqu.i[threadIdx.x] = ((unsigned int *)(atoms + j))[threadIdx.x];
    int aptmp = ((unsigned int *)(atom_params + j))[threadIdx.x];
    // scale vdw_type field, which is first in struct
    if ( (threadIdx.x & 3) == 0 ) aptmp *= lj_table_size;
    japu.i[threadIdx.x] = aptmp;
  }
  __syncthreads();

  USEPAIRLIST( if ( (pl[pli] & (1 << (threadIdx.x >> 5))) == 0 ) continue; )

  // calc forces on patch 1
  if ( blocki + threadIdx.x < myPatchPair.patch1_force_size ) {

    GENPAIRLIST( bool plpli = 0; )

    for ( int j = 0; j < shared_size; ++j ) {
      /* actually calculate force */
      float tmpx = jpqs[j].position.x - ipq.position.x;
      float tmpy = jpqs[j].position.y - ipq.position.y;
      float tmpz = jpqs[j].position.z - ipq.position.z;
      float r2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
      GENPAIRLIST( if(r2<plcutoff2) ) { GENPAIRLIST( plpli=1; )
      if ( r2 < cutoff2 ) {
        ENERGY( float rsqrtfr2; )
        float4 fi = tex1D(force_table, ENERGY(rsqrtfr2 =) rsqrtf(r2));
        ENERGY( float4 ei = tex1D(energy_table, rsqrtfr2); )
        float2 ljab = tex1Dfetch(lj_table,
                /* lj_table_size * */ japs[j].vdw_type + iap.vdw_type);
        bool excluded = false;
        int indexdiff = (int)(iap.index) - (int)(japs[j].index);
        if ( abs(indexdiff) <= (int) japs[j].excl_maxdiff ) {
          indexdiff += japs[j].excl_index;
          int indexword = indexdiff >> 5;
          if ( indexword < MAX_CONST_EXCLUSIONS )
               indexword = const_exclusions[indexword];
          else indexword = overflow_exclusions[indexword];
          excluded = ((indexword & (1<<(indexdiff&31))) != 0);
        }
        float f_slow = ipq.charge * jpqs[j].charge;
        float f = ljab.x * fi.z + ljab.y * fi.y + f_slow * fi.x;
        ENERGY(
        float ev = ljab.x * ei.z + ljab.y * ei.y;
        float ee = f_slow * ei.x;
        SLOW( float es = f_slow * ei.w; )
        )
        SLOW( f_slow *= fi.w; )
        if ( ! excluded ) { \
          ENERGY(
          totalev += ev;
          totalee += ee;
          SLOW( totales += es; )
          if ( blockj + j >= myPatchPair.patch2_force_size ) {
            /* add fixed atoms twice */
            totalev += ev;
            totalee += ee;
            SLOW( totales += es; )
          }
          )
          /* ife.w += r2 * f; */
          ife.x += tmpx * f;
          ife.y += tmpy * f;
          ife.z += tmpz * f;
          SLOW(
          /* ife_slow.w += r2 * f_slow; */
          ife_slow.x += tmpx * f_slow;
          ife_slow.y += tmpy * f_slow;
          ife_slow.z += tmpz * f_slow;
          )
        } else ife.w += 1.f;
      } }  /* cutoff */
    }

    GENPAIRLIST ( if ( plpli ) pl[pli] = 1; )

/*
    if ( plcutoff2 == 0 ) {  // use pairlist
      if ( doSlow ) {
        FORCE_INNER_LOOP(ipq,iap,1,{)
      } else {
        FORCE_INNER_LOOP(ipq,iap,0,{)
      }
    } else {  // create pairlist
      bool plpli = 0;
      if ( doSlow ) {
        FORCE_INNER_LOOP(ipq,iap,1,if(r2<plcutoff2){plpli=1;)
      } else {
        FORCE_INNER_LOOP(ipq,iap,0,if(r2<plcutoff2){plpli=1;)
      }
      if ( plpli ) pl[pli] = 1;
    }
*/

  } // if
  } // blockj loop

  if ( blocki + threadIdx.x < myPatchPair.patch1_force_size ) {
    int i_out = myPatchPair.patch1_force_start + blocki + threadIdx.x;
    force_buffers[i_out] = ife;
    if ( doSlow ) {
      slow_force_buffers[i_out] = ife_slow;
    }
    // accumulate net force to shared memory, warp-synchronous
    const int subwarp = threadIdx.x >> 2;  // 32 entries in table
    const int thread = threadIdx.x & 3;  // 4 threads share each entry
    for ( int g = 0; g < 4; ++g ) {
      if ( thread == g ) {
        sumf.a2d[subwarp][0] += ife.x;
        sumf.a2d[subwarp][1] += ife.y;
        sumf.a2d[subwarp][2] += ife.z;
        if ( doSlow ) {
          sumf_slow.a2d[subwarp][0] += ife_slow.x;
          sumf_slow.a2d[subwarp][1] += ife_slow.y;
          sumf_slow.a2d[subwarp][2] += ife_slow.z;
        }
      }
    }
  }
  if ( plcutoff2 != 0 ) {
    __syncthreads();  // all shared pairlist writes complete
    unsigned int pltmp;
    if ( threadIdx.x < 32 ) {
      pltmp = plu.i[threadIdx.x];
      pltmp |= plu.i[threadIdx.x+32] << 1;
      pltmp |= plu.i[threadIdx.x+64] << 2;
      pltmp |= plu.i[threadIdx.x+96] << 3;
    }
    __syncthreads();  // all shared pairlist reads complete
    if ( threadIdx.x < 32 ) {
      int i_pl = (blocki >> 2) + myPatchPair.block_flags_start;
      block_flags[i_pl + threadIdx.x] = pltmp;
    }
  }

  } // blocki loop

  __syncthreads();
  if ( threadIdx.x < 24 ) { // reduce forces, warp-synchronous
                            // 3 components, 8 threads per component
    const int i_out = myPatchPair.virial_start + threadIdx.x;
    {
      float f;
      f = sumf.a1d[threadIdx.x] + sumf.a1d[threadIdx.x + 24] + 
          sumf.a1d[threadIdx.x + 48] + sumf.a1d[threadIdx.x + 72];
      sumf.a1d[threadIdx.x] = f;
      f += sumf.a1d[threadIdx.x + 12];
      sumf.a1d[threadIdx.x] = f;
      f += sumf.a1d[threadIdx.x + 6];
      sumf.a1d[threadIdx.x] = f;
      f += sumf.a1d[threadIdx.x + 3];
      f *= 0.5f;  // compensate for double-counting
      // calculate virial contribution on first 3 threads
      sumf.a2d[threadIdx.x][0] = f * myPatchPair.offset.x;
      sumf.a2d[threadIdx.x][1] = f * myPatchPair.offset.y;
      sumf.a2d[threadIdx.x][2] = f * myPatchPair.offset.z;
      if ( threadIdx.x < 9 ) {  // write out output buffer
        virial_buffers[i_out] = sumf.a1d[threadIdx.x];
      }
    }
    if ( doSlow ) { // repeat above for slow forces
      float fs;
      fs = sumf_slow.a1d[threadIdx.x] + sumf_slow.a1d[threadIdx.x + 24] + 
           sumf_slow.a1d[threadIdx.x + 48] + sumf_slow.a1d[threadIdx.x + 72];
      sumf_slow.a1d[threadIdx.x] = fs;
      fs += sumf_slow.a1d[threadIdx.x + 12];
      sumf_slow.a1d[threadIdx.x] = fs;
      fs += sumf_slow.a1d[threadIdx.x + 6];
      sumf_slow.a1d[threadIdx.x] = fs;
      fs += sumf_slow.a1d[threadIdx.x + 3];
      fs *= 0.5f;
      sumf_slow.a2d[threadIdx.x][0] = fs * myPatchPair.offset.x;
      sumf_slow.a2d[threadIdx.x][1] = fs * myPatchPair.offset.y;
      sumf_slow.a2d[threadIdx.x][2] = fs * myPatchPair.offset.z;
      if ( threadIdx.x < 9 ) {
        slow_virial_buffers[i_out] = sumf_slow.a1d[threadIdx.x];
      }
    }
  }
#if ENERGY(1 +) 0
  if ( threadIdx.x < 512 )  // workaround for compiler bug
  {
    __syncthreads();
    // accumulate energies to shared memory, warp-synchronous
    const int subwarp = threadIdx.x >> 2;  // 32 entries in table
    const int thread = threadIdx.x & 3;  // 4 threads share each entry
    if ( thread == 0 ) {
      sumf.a2d[subwarp][0] = totalev;
      sumf.a2d[subwarp][1] = totalee;
      sumf.a2d[subwarp][2] = 0.f SLOW( + totales ) ;
    }
    for ( int g = 1; g < 4; ++g ) {
      if ( thread == g ) {
        sumf.a2d[subwarp][0] += totalev;
        sumf.a2d[subwarp][1] += totalee;
        SLOW( sumf.a2d[subwarp][2] += totales; )
      }
    }
    __syncthreads();
    if ( threadIdx.x < 24 ) { // reduce energies, warp-synchronous
                             // 3 components, 8 threads per component
      const int i_out = myPatchPair.virial_start + threadIdx.x;
      float f;
      f = sumf.a1d[threadIdx.x] + sumf.a1d[threadIdx.x + 24] + 
          sumf.a1d[threadIdx.x + 48] + sumf.a1d[threadIdx.x + 72];
      sumf.a1d[threadIdx.x] = f;
      f += sumf.a1d[threadIdx.x + 12];
      sumf.a1d[threadIdx.x] = f;
      f += sumf.a1d[threadIdx.x + 6];
      sumf.a1d[threadIdx.x] = f;
      f += sumf.a1d[threadIdx.x + 3];
      f *= 0.5f;  // compensate for double-counting
      if ( threadIdx.x < 3 ) {  // write out output buffer
        virial_buffers[i_out+9] = f;
      }
    }
  }
#endif

 } // end of nonbonded calc

 { // start of force sum

  // make sure forces are visible in global memory
  __threadfence();
  __syncthreads();

  __shared__ bool sumForces;

  if (threadIdx.x == 0) {
    int fli = myPatchPair.patch1_force_list_index;
    int fls = myPatchPair.patch1_force_list_size;
    int old = atomicInc(force_list_counters+fli,fls-1);
    sumForces = ( old == fls - 1 );
  }

  __syncthreads();

  if ( sumForces ) {
    NAME(dev_sum_forces)(myPatchPair.patch1_force_list_index,
       atoms,force_lists,force_buffers,
       virial_buffers,forces,virials);

    if ( doSlow ) {
      NAME(dev_sum_forces)(myPatchPair.patch1_force_list_index,
         atoms,force_lists,slow_force_buffers,
         slow_virial_buffers,slow_forces,slow_virials);
    }
  }

 } // end of force sum
}


__device__ __forceinline__ static void NAME(dev_sum_forces)(
        const int force_list_index,
	const atom *atoms,
	const force_list *force_lists,
	const float4 *force_buffers,
	const float *virial_buffers,
	float4 *forces, float *virials) {
// call with one block per patch
// call BLOCK_SIZE threads per block
// call with no shared memory

  #define myForceList fl.fl
  __shared__ union {
    force_list fl;
    unsigned int i[FORCE_LIST_SIZE];
  } fl;

  if ( threadIdx.x < FORCE_LIST_USED ) {
    unsigned int tmp = ((unsigned int*)force_lists)[
                        FORCE_LIST_SIZE*force_list_index+threadIdx.x];
    fl.i[threadIdx.x] = tmp;
  }

  __shared__ volatile union {
    float a3d[32][3 ENERGY(+1)][3];
    float a2d[32][9 ENERGY(+3)];
    float a1d[32*(9 ENERGY(+3))];
  } virial;

  for ( int i = threadIdx.x; i < 32*(9 ENERGY(+3)); i += BLOCK_SIZE ) {
    virial.a1d[i] = 0.f;
  }

  __syncthreads();

  float vxx = 0.f;
  float vxy = 0.f;
  float vxz = 0.f;
  float vyx = 0.f;
  float vyy = 0.f;
  float vyz = 0.f;
  float vzx = 0.f;
  float vzy = 0.f;
  float vzz = 0.f;

  for ( int j = threadIdx.x; j < myForceList.patch_size; j += BLOCK_SIZE ) {

    const float4 *fbuf = force_buffers + myForceList.force_list_start + j;
    float4 fout;
    fout.x = 0.f;
    fout.y = 0.f;
    fout.z = 0.f;
    fout.w = 0.f;
    for ( int i=0; i < myForceList.force_list_size; ++i ) {
      float4 f = *fbuf;
      fout.x += f.x;
      fout.y += f.y;
      fout.z += f.z;
      fout.w += f.w;
      fbuf += myForceList.patch_stride;
    }

    // compiler will use st.global.f32 instead of st.global.v4.f32
    // if forcedest is directly substituted in the assignment
    const int forcedest = myForceList.force_output_start + j;
    forces[forcedest] = fout;

    float4 pos = ((float4*)atoms)[myForceList.atom_start + j];

    // accumulate per-atom virials to registers
    vxx += fout.x * pos.x;
    vxy += fout.x * pos.y;
    vxz += fout.x * pos.z;
    vyx += fout.y * pos.x;
    vyy += fout.y * pos.y;
    vyz += fout.y * pos.z;
    vzx += fout.z * pos.x;
    vzy += fout.z * pos.y;
    vzz += fout.z * pos.z;

  }

  { // accumulate per-atom virials to shared memory, warp-synchronous
    const int subwarp = threadIdx.x >> 2;  // 32 entries in table
    const int thread = threadIdx.x & 3;  // 4 threads share each entry
    for ( int g = 0; g < 4; ++g ) {
      if ( thread == g ) {
        virial.a3d[subwarp][0][0] += vxx;
        virial.a3d[subwarp][0][1] += vxy;
        virial.a3d[subwarp][0][2] += vxz;
        virial.a3d[subwarp][1][0] += vyx;
        virial.a3d[subwarp][1][1] += vyy;
        virial.a3d[subwarp][1][2] += vyz;
        virial.a3d[subwarp][2][0] += vzx;
        virial.a3d[subwarp][2][1] += vzy;
        virial.a3d[subwarp][2][2] += vzz;
      }
    }
  }
  __syncthreads();
  { // accumulate per-compute virials to shared memory, data-parallel
    const int halfwarp = threadIdx.x >> 4;  // 8 half-warps
    const int thread = threadIdx.x & 15;
    if ( thread < (9 ENERGY(+3)) ) {
      for ( int i = halfwarp; i < myForceList.force_list_size; i += 8 ) {
        virial.a2d[halfwarp][thread] +=
          virial_buffers[myForceList.virial_list_start + 16*i + thread];
      }
    }
  }
  __syncthreads();
  { // reduce virials in shared memory, warp-synchronous
    const int subwarp = threadIdx.x >> 3;  // 16 quarter-warps
    const int thread = threadIdx.x & 7;  // 8 threads per component
    if ( subwarp < (9 ENERGY(+3)) ) {  // 9 components
      float v;
      v = virial.a2d[thread][subwarp] + virial.a2d[thread+8][subwarp] +
          virial.a2d[thread+16][subwarp] + virial.a2d[thread+24][subwarp];
      virial.a2d[thread][subwarp] = v;
      v += virial.a2d[thread+4][subwarp];
      virial.a2d[thread][subwarp] = v;
      v += virial.a2d[thread+2][subwarp];
      virial.a2d[thread][subwarp] = v;
      v += virial.a2d[thread+1][subwarp];
      virial.a2d[thread][subwarp] = v;
    }
  }
  __syncthreads();
  if ( threadIdx.x < (9 ENERGY(+3)) ) {  // 9 components
    virials[myForceList.virial_output_start + threadIdx.x] =
                                              virial.a2d[0][threadIdx.x];
  }

}

#endif // NAMD_CUDA

