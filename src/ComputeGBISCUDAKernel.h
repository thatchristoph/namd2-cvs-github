#include <stdio.h>
#define GBIS_CUDA
#include "ComputeGBIS.inl"

#define DORED 1

/*******************************************************************************
  These GBIS kernel was patterned after Jim's nonbonded kernel
*******************************************************************************/

/*
  TODO

  using the namd force indices to sum psi and dEda require
  every atom to have a force, i.e., fixed atoms won't work
  consider using numFree/numFixed atoms to correct this
*/

//00000000000000000000000000000000000000000000000000000000000
//
// GBIS Psi / dEda Reduction CUDA Kernal
//
//00000000000000000000000000000000000000000000000000000000000

// Reduce Energy
__device__ __forceinline__ static void GBIS_Energy_Reduction_Kernel (
  const int list_index,
  const force_list *lists,
  const float *sumBuffers,
  float *sum
) {
//one thread block (128 threads) is doing this single patch
  int tx = threadIdx.x;

  #define myList fl.fl
  __shared__ union {
    force_list fl;
    unsigned int i[FORCE_LIST_SIZE];
  } fl;

  if ( tx < FORCE_LIST_USED ) {
    unsigned int tmp = ((unsigned int*)lists)[
                        FORCE_LIST_SIZE*list_index+tx];
    fl.i[tx] = tmp;
  }
  __syncthreads();

  __shared__ float energySh[BLOCK_SIZE];

  if (tx < myList.force_list_size) {
    energySh[tx] = sumBuffers[myList.virial_list_start/16 + tx];
  }
  __syncthreads();

  float nf = 1.f * myList.force_list_size; // exact number to sum
  float p2f = log2(nf);
  int p2i = ceil(p2f);
  int n2i = 1<<p2i; // next power of 2 over number to sum
  int b = n2i / 2; // first midpoint
  if (tx + b < myList.force_list_size) { // handle non power of 2
    energySh[tx] += energySh[tx+b];
  }
  __syncthreads();
  b /= 2;

  for ( ; b > 0; b /= 2) { // handle remaining powers of 2
    if (tx < b) {
      energySh[tx] += energySh[tx+b];
    }
    __syncthreads();
  }

    if (threadIdx.x == 0) {
      const int dest = myList.virial_output_start/16;
      sum[dest] = energySh[0];
    }
} // GBIS_Energy_Reduction

// Reduce Forces
__device__ __forceinline__ static void GBIS_Force_Reduction_Kernel (
  const int list_index,
  const force_list *lists,
  const float4 *sumBuffers,
  float4 *sum
) {
  #define myList fl.fl
  __shared__ union {
    force_list fl;
    unsigned int i[FORCE_LIST_SIZE];
  } fl;

  if ( threadIdx.x < FORCE_LIST_USED ) {
    unsigned int tmp = ((unsigned int*)lists)[
                        FORCE_LIST_SIZE*list_index+threadIdx.x];
    fl.i[threadIdx.x] = tmp;
  }
  __syncthreads();
  for ( int j = threadIdx.x; j < myList.patch_size; j += BLOCK_SIZE ) {
    const float4 *sumBuf = sumBuffers + myList.force_list_start + j;
    float4 sumOut; sumOut.x = 0.f; sumOut.y = 0.f; sumOut.z = 0.f;
    for ( int i = 0; i < myList.force_list_size; ++i ) {
      float4 sumVal = *sumBuf;
      sumOut.x += sumVal.x;
      sumOut.y += sumVal.y;
      sumOut.z += sumVal.z;
      sumBuf += myList.patch_stride;
    } // i
    const int dest = myList.force_output_start + j;
    sum[dest].x += sumOut.x; // adding onto NAMD nonbonded, sync?
    sum[dest].y += sumOut.y;
    sum[dest].z += sumOut.z;
  } // j
} // GBIS_Force_Reduction

// Reduce psi / dEda
__device__ __forceinline__ static void GBIS_Atom_Reduction_Kernel (
  const int list_index,
  const force_list *lists,
  const GBReal *sumBuffers,
  GBReal *sum
) {
  #define myList fl.fl
  __shared__ union {
    force_list fl;
    unsigned int i[FORCE_LIST_SIZE];
  } fl;

  if ( threadIdx.x < FORCE_LIST_USED ) {
    unsigned int tmp = ((unsigned int*)lists)[
                        FORCE_LIST_SIZE*list_index+threadIdx.x];
    fl.i[threadIdx.x] = tmp;
  }
  __syncthreads();
  for ( int j = threadIdx.x; j < myList.patch_size; j += BLOCK_SIZE ) {
    const GBReal *sumBuf = sumBuffers + myList.force_list_start + j;
    GBReal sumOut = 0.f;
    for ( int i = 0; i < myList.force_list_size; ++i ) {
      GBReal sumVal = *sumBuf;
      sumOut += sumVal;
      sumBuf += myList.patch_stride;
    } // i
    const int dest = myList.force_output_start + j;
    sum[dest] = sumOut;
  } // j
} // GBIS_Atom_Reduction

//1111111111111111111111111111111111111111111111111111111111
//
// GBIS Phase 1 CUDA Kernal
//
//1111111111111111111111111111111111111111111111111111111111
__global__ static void GBIS_P1_Kernel (
  const patch_pair *patch_pairs, // atoms pointers and such
  const atom *atoms,             // position & charge
  const atom_param *atom_params,
  const float *intRad0,          // read in intrinsic radius
  const float *intRadS,          // read in intrinsic radius
  GBReal *psiSumBuffers,         // write out psi summation
  GBReal *psiSum,
  const float a_cut,             // P1 interaction cutoff
  const float rho_0,             // H(i,j) parameter
  float3 lata,
  float3 latb,
  float3 latc,
  const force_list *force_lists,
  unsigned int *P1_counters
) {
  int tx = threadIdx.x;
  int bx = blockIdx.x;
#ifdef __DEVICE_EMULATION__
//  printf("GBIS_P1: Blk(%3i) Thd(%3i)\n",bx,tx);
#endif
// one block per patch_pair
// BLOCK_SIZE 128 threads per block
// SHARED_SIZE 32
// PATCH_PAIR_SIZE 16
// PATCH_PAIR_USED 15

//macro for patch pair
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
  if ( tx < PATCH_PAIR_USED ) {
    unsigned int tmp = ((unsigned int*)patch_pairs)[PATCH_PAIR_SIZE*bx+tx];
    pp.i[tx] = tmp;
  } // ?


//macro for list of j atoms
#ifdef __DEVICE_EMULATION__
  #define jatoms ((atom*)(jpqu.i))
#else
  #define jatoms jpqu.d
#endif
  __shared__ union {
#ifndef __DEVICE_EMULATION__
    atom d[SHARED_SIZE];
#endif
    int i[4*SHARED_SIZE];
    float f[4*SHARED_SIZE];
  } jpqu;

  // convert scaled offset with current lattice
  if ( tx == 0 ) {
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

  int patch1size = myPatchPair.patch1_size;
  int patch2size = myPatchPair.patch2_size;

  //iterate over chunks of atoms within Patch 1
  for ( int blocki = 0; blocki < patch1size; blocki += BLOCK_SIZE ) {

    //this thread calculates only the force on atomi; iterating over js
    atom atomi;
    //int iindex;
    int i;

    // load BLOCK of Patch i atoms
    if ( blocki + tx < patch1size ) {
      i = myPatchPair.patch1_atom_start + blocki + tx;
      float4 tmpa = ((float4*)atoms)[i];
      atomi.position.x = tmpa.x + myPatchPair.offset.x;
      atomi.position.y = tmpa.y + myPatchPair.offset.y;
      atomi.position.z = tmpa.z + myPatchPair.offset.z;
      atomi.charge = intRad0[i]; // overwrite charge with radius
      //iindex =  ((uint4 *)(atom_params))[i].y; // store index
    } // load patch 1

    //init intermediate variables
    GBReal psiSumI = 0.f; // each thread accumulating single psi

    //iterate over chunks of atoms within Patch 2
    for ( int blockj = 0; blockj < patch2size; blockj += SHARED_SIZE ) {

      int shared_size = patch2size - blockj;
      if ( shared_size > SHARED_SIZE )
        shared_size = SHARED_SIZE;
      __syncthreads();

      //load smaller chunk of j atoms into shared memory: coordinates

      int j = myPatchPair.patch2_atom_start + blockj;
      if ( tx < 4 * shared_size ) {
        jpqu.i[tx] = ((unsigned int *)(atoms       + j))[tx];
      }
      __syncthreads();

      //load smaller chunk of j atoms into shared memory: radii
      if ( tx < shared_size ) { // every 4th thread
        jatoms[tx].charge = intRadS[j+tx]; // load radius into charge float
        //printf("jatom2 = %7.3f\n", jatoms[tx].charge);
      }
      __syncthreads();

      if ( blocki + tx < patch1size ) {
        //each thread loop over shared atoms
        for ( int j = 0; j < shared_size; ++j ) {
          float dx = atomi.position.x - jatoms[j].position.x;
          float dy = atomi.position.y - jatoms[j].position.y;
          float dz = atomi.position.z - jatoms[j].position.z;
          float r2 = dx*dx + dy*dy + dz*dz;

          // within cutoff        different atoms
          if (r2 < (a_cut+FS_MAX)*(a_cut+FS_MAX) && r2 > 0.01f) {
            //int jj = myPatchPair.patch2_atom_start + blockj + j;
            //int jindex =  ((uint4 *)(atom_params))[jj].y; // store index

            // calculate H(i,j) [and not H(j,i)]
            float r_i = 1.f / sqrt(r2);
            float r  = r2 * r_i;
            float hij;
            int dij;
            CalcH(r,r2,r_i,a_cut,atomi.charge,jatoms[j].charge,hij,dij);

          psiSumI += hij;

          } // cutoff
        } // for j
      } // if thread in bounds
    } // for block j
    //psiSumI now contains contributions from all j in 2nd patch

    // write psiSum to global memory buffer; to be accumulated later
    if ( blocki + tx < myPatchPair.patch1_force_size) {
      int i_out = myPatchPair.patch1_force_start + blocki + tx;
      psiSumBuffers[i_out] = psiSumI;
    }
  } // for block i

 { // start of force sum
  // make sure psiSums are visible in global memory
  __threadfence();
  __syncthreads();
  __shared__ bool doSum;

  if (threadIdx.x == 0) {
    int fli = myPatchPair.patch1_force_list_index;
    int fls = myPatchPair.patch1_force_list_size;
    int old = atomicInc(P1_counters+fli,fls-1);
    doSum = ( old == fls - 1 );
  }

  __syncthreads();

  if ( doSum ) {
    GBIS_Atom_Reduction_Kernel(myPatchPair.patch1_force_list_index,
       force_lists, psiSumBuffers, psiSum);
  }
 } // end of force sum
} //GBIS_P1

//2222222222222222222222222222222222222222222222222222222222
//
// GBIS Phase 2 CUDA Kernal
//
//2222222222222222222222222222222222222222222222222222222222
__global__ static void GBIS_P2_Kernel (
  const patch_pair *patch_pairs,// atoms pointers and such
  const atom *atoms,            // position & charge
  const atom_param *atom_params,
  const float *bornRad,         // read in Born radius
  GBReal *dEdaSumBuffers,       // write out dEda summation
  GBReal *dEdaSum,
  const float a_cut,            // P1 interaction cutoff
  const float r_cut,            // P1 interaction cutoff
  const float scaling,          // scale nonbonded
  const float kappa,
  const float smoothDist,       // use interaction cutoff smoothing?
  const float epsilon_p,        // protein dielectric
  const float epsilon_s,        // solvent dielectric
  float3 lata,
  float3 latb,
  float3 latc,
  const int doEnergy,           // calculate energy too?
  const int doFullElec,         // calc dEdaSum for P3 full electrostatics
  const force_list *force_lists,
  float4 *forceBuffers,
  float4 *forces,
  float *energyBuffers,
  float *energy,
  unsigned int *P2_counters
) {
  int tx = threadIdx.x;
  int bx = blockIdx.x;
#ifdef __DEVICE_EMULATION__
//  printf("GBIS_P2: Blk(%3i) Thd(%3i)\n",bx,tx);
#endif
// one block per patch_pair
// BLOCK_SIZE 128 threads per block
// SHARED_SIZE 32
// PATCH_PAIR_SIZE 16
// PATCH_PAIR_USED 15

//macro for patch pair
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
  if ( tx < PATCH_PAIR_USED ) {
    unsigned int tmp = ((unsigned int*)patch_pairs)[PATCH_PAIR_SIZE*bx+tx];
    pp.i[tx] = tmp;
  } // ?


//macro for list of j atoms
#ifdef __DEVICE_EMULATION__
  #define jatoms ((atom*)(jpqu.i))
#else
  #define jatoms jpqu.d
#endif
  __shared__ union {
#ifndef __DEVICE_EMULATION__
    atom d[SHARED_SIZE];
#endif
    int i[4*SHARED_SIZE];
    float f[4*SHARED_SIZE];
  } jpqu;

  __shared__ float jBornRad[SHARED_SIZE];
  float energyT = 0.f; // total energy for this thread; to be reduced

  // convert scaled offset with current lattice
  if ( tx == 0 ) {
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

  int patch1size = myPatchPair.patch1_size;
  int patch2size = myPatchPair.patch2_size;

  //values used in loop
  float r_cut2 = r_cut*r_cut;
  float r_cut_2 = 1.f / r_cut2;
  float r_cut_4 = 4.f*r_cut_2*r_cut_2;
  float epsilon_s_i = 1.f / epsilon_s;
  float epsilon_p_i = 1.f / epsilon_p;

  //iterate over chunks of atoms within Patch 1
  for ( int blocki = 0; blocki < patch1size; blocki += BLOCK_SIZE ) {

    //this thread calculates only the force on atomi; iterating over js
    atom atomi;
    float bornRadI;
    //int iindex;
    int i;

    // load BLOCK of Patch i atoms
    if ( blocki + tx < patch1size ) {
      i = myPatchPair.patch1_atom_start + blocki + tx;
      float4 tmpa = ((float4*)atoms)[i];
      atomi.position.x = tmpa.x + myPatchPair.offset.x;
      atomi.position.y = tmpa.y + myPatchPair.offset.y;
      atomi.position.z = tmpa.z + myPatchPair.offset.z;
      atomi.charge = - tmpa.w * scaling;
      bornRadI = bornRad[i];
      //iindex =  ((uint4 *)(atom_params))[i].y; // store index
    } // load patch 1

    //init intermediate variables
    GBReal dEdaSumI = 0.f; // each thread accumulating single psi
    float4 forceI; forceI.x = forceI.y = forceI.z = forceI.w = 0.f;

    //iterate over chunks of atoms within Patch 2
    for ( int blockj = 0; blockj < patch2size; blockj += SHARED_SIZE ) {

      int shared_size = patch2size - blockj;
      if ( shared_size > SHARED_SIZE )
        shared_size = SHARED_SIZE;
      __syncthreads();

      //load smaller chunk of j atoms into shared memory: coordinates
      int j = myPatchPair.patch2_atom_start + blockj;
      if ( tx < 4 * shared_size ) {
        jpqu.i[tx] = ((unsigned int *)(atoms + j))[tx];
      }
      __syncthreads();

      //load smaller chunk of j atoms into shared memory: radii
      if ( tx < shared_size ) {
        jBornRad[tx] = bornRad[j+tx]; // load Born radius into shared
        //printf("BornRadius = %7.3f\n", jBornRadii[tx]);
      }
      __syncthreads();

      if ( blocki + tx < patch1size ) {
        //each thread loop over shared atoms
        for ( int j = 0; j < shared_size; ++j ) {
          float dx = atomi.position.x - jatoms[j].position.x;
          float dy = atomi.position.y - jatoms[j].position.y;
          float dz = atomi.position.z - jatoms[j].position.z;
          float r2 = dx*dx + dy*dy + dz*dz;

          // within cutoff different atoms
          if (r2 < r_cut*r_cut && r2 > 0.01f) {
            //int jj = myPatchPair.patch2_atom_start + blockj + j;
            //int jindex =  ((uint4 *)(atom_params))[jj].y; // store index

            float r_i = 1.f / sqrt(r2);
            float r  = r2 * r_i;
            float bornRadJ = jBornRad[j];

            //calculate GB energy
            float qiqj = atomi.charge*jatoms[j].charge;
            float aiaj = bornRadI*bornRadJ;
            float aiaj4 = 4*aiaj;
            float expr2aiaj4 = exp(-r2/aiaj4);
            float fij = sqrt(r2+aiaj*expr2aiaj4);
            float f_i = 1/fij;
            float expkappa = exp(-kappa*fij);
            float Dij = epsilon_p_i - expkappa*epsilon_s_i;
            float gbEij = qiqj*Dij*f_i;

            //calculate energy derivatives
            float ddrfij = r*f_i*(1.f - 0.25f*expr2aiaj4);
            float ddrf_i = -ddrfij*f_i*f_i;
            float ddrDij = kappa*expkappa*ddrfij*epsilon_s_i;
            float ddrGbEij = qiqj*(ddrDij*f_i+Dij*ddrf_i);

            //NAMD smoothing function
            float scale = 1.f;
            float ddrScale = 0.f;
            float forcedEdr;
            if (smoothDist > 0.f) {
              scale = r2 * r_cut_2 - 1.f;
              scale *= scale;
              ddrScale = r*(r2-r_cut2)*r_cut_4;
              energyT += gbEij * scale * 0.5f;
              //printf("PAIRENERGY = %15.7e\n", gbEij*scale*0.5);
              forcedEdr = -(ddrGbEij)*scale-(gbEij)*ddrScale;
            } else {
              energyT += gbEij * 0.5f;
              //printf("PAIRENERGY = %15.7e\n", gbEij*0.5);
              forcedEdr = -ddrGbEij;
            }

            //add dEda
            float dEdai = 0.f;
            if (doFullElec) {
              //gbisParams->dEdaSum[0][i] += dEdai*scale;
              dEdai = 0.5f*qiqj*f_i*f_i
                      *(kappa*epsilon_s_i*expkappa-Dij*f_i)
                      *(aiaj+0.25f*r2)*expr2aiaj4/bornRadI*scale;//0
              dEdaSumI += dEdai;
            }

            forcedEdr *= r_i;
            forceI.x += dx*forcedEdr;
            forceI.y += dy*forcedEdr;
            forceI.z += dz*forcedEdr;


          } // within cutoff
          if (r2 < 0.01f) {
            // GB Self Energy
            if (doEnergy) {
              float fij = bornRadI;//inf
              float expkappa = exp(-kappa*fij);//0
              float Dij = epsilon_p_i - expkappa*epsilon_s_i;
              float gbEij = atomi.charge*(atomi.charge / (-scaling) )*Dij/fij;
              energyT += 0.5f*gbEij;
            }
          } //same atom or within cutoff
        } // for j
      } // if thread in bounds
    } // for block j
    //psiSumI now contains contributions from all j in 2nd patch

    // write psiSum to global memory buffer; to be accumulated later
    if ( blocki + tx < myPatchPair.patch1_force_size) {
      int i_out = myPatchPair.patch1_force_start + blocki + tx;
      dEdaSumBuffers[i_out] = dEdaSumI;
      forceBuffers[i_out] = forceI;
    }
  } // for block i
  //printf("energyT[%03i,%03i] = %15.7e\n", blockIdx.x, tx, energyT);

//Energy Reduction
  if (doEnergy) {
    const int i_out = myPatchPair.virial_start / 16;
    __shared__ float energySh[BLOCK_SIZE];
    energySh[tx] = energyT;
    __syncthreads();

    for ( unsigned int b = blockDim.x/2; b > 0 ; b >>= 1 ) {
      if ( tx < b ) {
        //printf("RED eSh[%i] += eSh[%i]\n", tx, tx+b);
        energySh[tx] += energySh[tx+b];
      }
      __syncthreads();
    }
    //const int i_out = myPatchPair.virial_start + tx;
    if (tx == 0) {
      energyBuffers[i_out] = energySh[0]; // TODO
      //printf("eB[%i] = %15.7e\n", i_out, energySh[0]);
    }
  } //end Energy Reduction

 { // start of reduction
  // make sure psiSums are visible in global memory
  __threadfence();
  __syncthreads();

  __shared__ bool doSum;

  if (threadIdx.x == 0) {
    int fli = myPatchPair.patch1_force_list_index;
    int fls = myPatchPair.patch1_force_list_size;
    int old = atomicInc(P2_counters+fli,fls-1);
    doSum = ( old == fls - 1 );
  }

  __syncthreads();

  if ( doSum ) {
    GBIS_Atom_Reduction_Kernel(myPatchPair.patch1_force_list_index,
       force_lists, dEdaSumBuffers, dEdaSum);
    GBIS_Force_Reduction_Kernel(myPatchPair.patch1_force_list_index,
       force_lists, forceBuffers, forces);
    if ( doEnergy ) {
      GBIS_Energy_Reduction_Kernel(myPatchPair.patch1_force_list_index,
        force_lists, energyBuffers, energy);
    }
  }
 } // end of sum
} //GBIS_P2


//3333333333333333333333333333333333333333333333333333333333
//
// GBIS Phase 3 CUDA Kernal
//
//3333333333333333333333333333333333333333333333333333333333
__global__ static void GBIS_P3_Kernel (
  const patch_pair *patch_pairs,  // atoms pointers and such
  const atom *atoms,              // position & charge
  const atom_param *atom_params,
  const float *intRad0,           // read in intrinsic radius
  const float *intRadS,           // read in intrinsic radius
  const float *dHdrPrefix,        // read in prefix
  const float a_cut,              // P1 interaction cutoff
  const float rho_0,              // H(i,j) parameter
  const float scaling,            // scale nonbonded
  float3 lata,
  float3 latb,
  float3 latc,
  const force_list *force_lists,
  float4 *forceBuffers,
  float4 *forces,
  unsigned int *P3_counters
) {
  int tx = threadIdx.x;
  int bx = blockIdx.x;
#ifdef __DEVICE_EMULATION__
//  printf("GBIS_P3: Blk(%3i) Thd(%3i)\n",bx,tx);
#endif

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
  if ( tx < PATCH_PAIR_USED ) {
    unsigned int tmp = ((unsigned int*)patch_pairs)[PATCH_PAIR_SIZE*bx+tx];
    pp.i[tx] = tmp;
  } // all threads load single shared patch pair

//macro for list of j atoms
#ifdef __DEVICE_EMULATION__
  #define jatoms ((atom*)(jpqu.i))
#else
  #define jatoms jpqu.d
#endif
  __shared__ union {
#ifndef __DEVICE_EMULATION__
    atom d[SHARED_SIZE];
#endif
    int i[4*SHARED_SIZE];
    float f[4*SHARED_SIZE];
  } jpqu;

  __shared__ float jDHdrPrefix[SHARED_SIZE];
  __shared__ float intRadJ0[SHARED_SIZE];

  // convert scaled offset with current lattice
  if ( tx == 0 ) {
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

  //iterate over chunks of atoms within Patch 1
  for ( int blocki = 0; blocki < myPatchPair.patch1_size; blocki += BLOCK_SIZE ) {
    //this thread calculates only the force on atomi; iterating over js
    atom atomi;
    float intRadIS;
    //int iindex;
    int i;
    float dHdrPrefixI;

    // load BLOCK of Patch i atoms
    if ( blocki + tx < myPatchPair.patch1_size ) {
      i = myPatchPair.patch1_atom_start + blocki + tx;
      float4 tmpa = ((float4*)atoms)[i];
      atomi.position.x = tmpa.x + myPatchPair.offset.x;
      atomi.position.y = tmpa.y + myPatchPair.offset.y;
      atomi.position.z = tmpa.z + myPatchPair.offset.z;
      atomi.charge = intRad0[i]; // overwrite charge with radius
      intRadIS = intRadS[i];
      //iindex =  ((uint4 *)(atom_params))[i].y; // store index
      dHdrPrefixI = dHdrPrefix[i];
    } // load patch 1

    //init intermediate variables
    float4 forceI; forceI.x = forceI.y = forceI.z = forceI.w = 0.f;

    //iterate over chunks of atoms within Patch 2
    for ( int blockj = 0; blockj < myPatchPair.patch2_size; blockj += SHARED_SIZE ) {

      int shared_size = myPatchPair.patch2_size - blockj;
      if ( shared_size > SHARED_SIZE )
        shared_size = SHARED_SIZE;
      __syncthreads();

      //load smaller chunk of j atoms into shared memory: coordinates

      int j = myPatchPair.patch2_atom_start + blockj;
      if ( tx < 4 * shared_size ) {
        jpqu.i[tx] = ((unsigned int *)(atoms       + j))[tx];
      }
      __syncthreads();

      //load smaller chunk of j atoms into shared memory: radii
      if ( tx < shared_size ) { // every 4th thread
        jatoms[tx].charge = intRadS[j+tx]; // load radius into charge float
        jDHdrPrefix[tx] = dHdrPrefix[j+tx]; // load dHdrPrefix into shared
        intRadJ0[tx] = intRad0[j+tx];
      }
      __syncthreads();

      if ( blocki + tx < myPatchPair.patch1_size ) {
        //each thread loop over shared atoms
        for ( int j = 0; j < shared_size; ++j ) {
          float dx = atomi.position.x - jatoms[j].position.x;
          float dy = atomi.position.y - jatoms[j].position.y;
          float dz = atomi.position.z - jatoms[j].position.z;
          float r2 = dx*dx + dy*dy + dz*dz;

          // within cutoff        different atoms
          if (r2 < (a_cut+FS_MAX)*(a_cut+FS_MAX) && r2 > 0.01f) {
            //int jj = myPatchPair.patch2_atom_start + blockj + j;
            //int jindex =  ((uint4 *)(atom_params))[jj].y; // store index

            float r_i = 1.f / sqrt(r2);
            float r  = r2 * r_i;
            float dHdrPrefixJ = jDHdrPrefix[j];
            float dhij, dhji;
            int dij, dji;
            CalcDH(r,r2,r_i,a_cut,atomi.charge,
                    jatoms[j].charge,dhij,dij);
            CalcDH(r,r2,r_i,a_cut,intRadJ0[j],
                    intRadIS,dhji,dji);

            float forceAlpha = -r_i*(dHdrPrefixI*dhij+dHdrPrefixJ*dhji);
            forceI.x += dx * forceAlpha;
            forceI.y += dy * forceAlpha;
            forceI.z += dz * forceAlpha;
          } // cutoff
        } // for j
      } // if thread in bounds
    } // for block j

    // write psiSum to global memory buffer; to be accumulated later
    if ( blocki + tx < myPatchPair.patch1_force_size) {
      int i_out = myPatchPair.patch1_force_start + blocki + tx;
      forceBuffers[i_out] = forceI;
    }
  } // for block i

 { // start of force sum
  // make sure forces are visible in global memory
  __threadfence();
  __syncthreads();
  __shared__ bool doSum;

  if (threadIdx.x == 0) {
    int fli = myPatchPair.patch1_force_list_index;
    int fls = myPatchPair.patch1_force_list_size;
    int old = atomicInc(P3_counters+fli,fls-1);
    doSum = ( old == fls - 1 );
  }
  __syncthreads();

  if ( doSum ) {
    GBIS_Force_Reduction_Kernel(myPatchPair.patch1_force_list_index,
       force_lists, forceBuffers, forces);
  }
 } // end of force sum
} //GBIS_P3
