
#include "ComputeNonbondedCUDAKernel.h"
#include <stdio.h>

#ifdef NAMD_CUDA

#ifdef WIN32
//  not supported by nvcc on Windows
// #define __thread __declspec(thread)
#define __thread
#endif

__constant__ unsigned int const_exclusions[MAX_CONST_EXCLUSIONS];

static __thread unsigned int *overflow_exclusions;

#define SET_EXCL(EXCL,BASE,DIFF) \
         (EXCL)[((BASE)+(DIFF))>>5] |= (1<<(((BASE)+(DIFF))&31))

void cuda_bind_exclusions(const unsigned int *t, int n) {

  cudaMalloc((void**) &overflow_exclusions, n*sizeof(unsigned int));
  cuda_errcheck("malloc overflow_exclusions");
  cudaMemcpy(overflow_exclusions, t,
		n*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to overflow_exclusions");
  int nconst = ( n < MAX_CONST_EXCLUSIONS ? n : MAX_CONST_EXCLUSIONS );
  cudaMemcpyToSymbol(const_exclusions, t, nconst*sizeof(unsigned int), 0);
  cuda_errcheck("memcpy to const_exclusions");
}


texture<float2, 1, cudaReadModeElementType> lj_table;
static __thread int lj_table_size;

void cuda_bind_lj_table(const float2 *t, int _lj_table_size) {
    static __thread float2 *ct;
    static __thread int lj_table_alloc;
    lj_table_size = _lj_table_size;
    if ( ct && lj_table_alloc < lj_table_size ) {
      cudaFree(ct);
      cuda_errcheck("freeing lj table");
      ct = 0;
    }
    if ( ! ct ) {
      lj_table_alloc = lj_table_size;
      cudaMalloc((void**) &ct, lj_table_size*lj_table_size*sizeof(float2));
      cuda_errcheck("allocating lj table");
    }
    cudaMemcpy(ct, t, lj_table_size*lj_table_size*sizeof(float2),
                                            cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to lj table");

    lj_table.normalized = false;
    lj_table.addressMode[0] = cudaAddressModeClamp;
    lj_table.filterMode = cudaFilterModePoint;

    cudaBindTexture((size_t*)0, lj_table, ct,
        lj_table_size*lj_table_size*sizeof(float2));
    cuda_errcheck("binding lj table to texture");
}


texture<float4, 1, cudaReadModeElementType> force_table;
texture<float4, 1, cudaReadModeElementType> energy_table;

void cuda_bind_force_table(const float4 *t, const float4 *et) {
    static __thread cudaArray *ct;
    static __thread cudaArray *ect;
    if ( ! ct ) {
      cudaMallocArray(&ct, &force_table.channelDesc, FORCE_TABLE_SIZE, 1);
      cuda_errcheck("allocating force table");
    }
    if ( ! ect ) {
      cudaMallocArray(&ect, &energy_table.channelDesc, FORCE_TABLE_SIZE, 1);
      cuda_errcheck("allocating energy table");
    }
    cudaMemcpyToArray(ct, 0, 0, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    // cudaMemcpy(ct, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to force table");
    cudaMemcpyToArray(ect, 0, 0, et, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to energy table");

    force_table.normalized = true;
    force_table.addressMode[0] = cudaAddressModeClamp;
    force_table.addressMode[1] = cudaAddressModeClamp;
    force_table.filterMode = cudaFilterModeLinear;

    energy_table.normalized = true;
    energy_table.addressMode[0] = cudaAddressModeClamp;
    energy_table.addressMode[1] = cudaAddressModeClamp;
    energy_table.filterMode = cudaFilterModeLinear;

    cudaBindTextureToArray(force_table, ct);
    cuda_errcheck("binding force table to texture");

    cudaBindTextureToArray(energy_table, ect);
    cuda_errcheck("binding energy table to texture");
}

static __thread int patch_pairs_size;
static __thread patch_pair *patch_pairs;
static __thread float *virial_buffers;  // one per patch pair
static __thread float *slow_virial_buffers;  // one per patch pair

static __thread int block_flags_size;
static __thread unsigned int *block_flags;

static __thread int force_lists_size;
static __thread force_list *force_lists;
static __thread unsigned int *force_list_counters;
static __thread unsigned int *GBIS_P1_counters;
static __thread unsigned int *GBIS_P2_counters;
static __thread unsigned int *GBIS_P3_counters;

static __thread int force_buffers_size;
static __thread float4 *force_buffers;
static __thread float4 *slow_force_buffers;

static __thread int atoms_size;
static __thread atom *atoms;
static __thread atom_param *atom_params;
static __thread float4 *forces;
static __thread float4 *slow_forces;
static __thread float *virials;  // one per patch
static __thread float *slow_virials;  // one per patch
static __thread float *energy_gbis;  // one per patch
static __thread float *energy_gbis_buffers;  // one per pair

//GBIS arrays
static __thread float  *intRad0D;      // one per patch
static __thread float  *intRadSD;      // one per patch
static __thread GBReal *psiSumD;     // one per patch
static __thread GBReal *psiSumD_buffers; // one per patch
static __thread float  *bornRadD;     // one per patch
static __thread GBReal *dEdaSumD;    // one per patch
static __thread GBReal *dEdaSumD_buffers; // one per patch
static __thread float  *dHdrPrefixD;  // one per patch

static __thread int patch_pairs_alloc;
static __thread int block_flags_alloc;
static __thread int force_buffers_alloc;
static __thread int force_lists_alloc;
static __thread int atoms_alloc;

static __thread int max_atoms_per_patch;

__thread cudaStream_t stream;
__thread cudaStream_t stream2;
 
void cuda_init() {
  forces = 0;
  slow_forces = 0;
  virials = 0;
  energy_gbis = 0;
  slow_virials = 0;
  atom_params = 0;
  atoms = 0;
  force_buffers = 0;
  slow_force_buffers = 0;
  force_lists = 0;
  force_list_counters = 0;
  GBIS_P1_counters = 0;
  GBIS_P2_counters = 0;
  GBIS_P3_counters = 0;
  patch_pairs = 0;
  virial_buffers = 0;
  energy_gbis_buffers = 0;
  slow_virial_buffers = 0;
  block_flags = 0;

  intRad0D = 0;
  intRadSD = 0;
  psiSumD = 0;
  psiSumD_buffers = 0;
  bornRadD = 0;
  dEdaSumD = 0;
  dEdaSumD_buffers = 0;
  dHdrPrefixD = 0;

  patch_pairs_alloc = 0;
  block_flags_alloc = 0;
  force_buffers_alloc = 0;
  force_lists_alloc = 0;
  atoms_alloc = 0;

  cudaStreamCreate(&stream);
  cudaStreamCreate(&stream2);
  cuda_errcheck("cudaStreamCreate");
}

void cuda_bind_patch_pairs(const patch_pair *pp, int npp,
                        const force_list *fl, int nfl,
                        int atoms_size_p, int force_buffers_size_p,
                        int block_flags_size_p, int max_atoms_per_patch_p) {

  patch_pairs_size = npp;
  force_buffers_size = force_buffers_size_p;
  force_lists_size = nfl;
  atoms_size = atoms_size_p;
  block_flags_size = block_flags_size_p;
  max_atoms_per_patch = max_atoms_per_patch_p;

#if 0
 printf("%d %d %d %d %d %d %d %d\n",
      patch_pairs_size , patch_pairs_alloc ,
      force_buffers_size , force_buffers_alloc ,
      force_lists_size , force_lists_alloc ,
      atoms_size , atoms_alloc );
#endif

 if ( patch_pairs_size > patch_pairs_alloc ||
      block_flags_size > block_flags_alloc ||
      force_buffers_size > force_buffers_alloc ||
      force_lists_size > force_lists_alloc ||
      atoms_size > atoms_alloc ) {

  block_flags_alloc = (int) (1.2 * block_flags_size);
  patch_pairs_alloc = (int) (1.2 * patch_pairs_size);
  force_buffers_alloc = (int) (1.2 * force_buffers_size);
  force_lists_alloc = (int) (1.2 * force_lists_size);
  atoms_alloc = (int) (1.2 * atoms_size);

  // if ( forces ) cudaFree(forces);
  // if ( slow_forces ) cudaFree(slow_forces);
  forces = slow_forces = 0;
  if ( atom_params ) cudaFree(atom_params);
  if ( atoms ) cudaFree(atoms);
  if ( force_buffers ) cudaFree(force_buffers);
  if ( slow_force_buffers ) cudaFree(slow_force_buffers);
  if ( force_lists ) cudaFree(force_lists);
  if ( force_list_counters ) cudaFree(force_list_counters);
  if ( GBIS_P1_counters ) cudaFree(GBIS_P1_counters);
  if ( GBIS_P2_counters ) cudaFree(GBIS_P2_counters);
  if ( GBIS_P3_counters ) cudaFree(GBIS_P3_counters);
  // if ( virials ) cudaFree(virials);
  virials = slow_virials = 0;
  energy_gbis = 0;
  if ( patch_pairs ) cudaFree(patch_pairs);
  if ( virial_buffers ) cudaFree(virial_buffers);
  if ( energy_gbis_buffers ) cudaFree(energy_gbis_buffers);
  if ( slow_virial_buffers ) cudaFree(slow_virial_buffers);
  if ( block_flags ) cudaFree(block_flags);
  if ( intRad0D ) cudaFree(intRad0D); // GBIS memory
  if ( intRadSD ) cudaFree(intRadSD);
  //if ( psiSumD ) cudaFree(psiSumD);
  if ( psiSumD_buffers ) cudaFree(psiSumD_buffers);
  if ( bornRadD ) cudaFree(bornRadD);
  //if ( dEdaSumD ) cudaFree(dEdaSumD);
  if ( dEdaSumD_buffers ) cudaFree(dEdaSumD_buffers);
  if ( dHdrPrefixD ) cudaFree(dHdrPrefixD);
  cuda_errcheck("free everything");

#if 0
  int totalmem = patch_pairs_alloc * sizeof(patch_pair) +
		force_lists_alloc * sizeof(force_list) +
		2 * force_buffers_alloc * sizeof(float4) +
		atoms_alloc * sizeof(atom) +
		atoms_alloc * sizeof(atom_param) +
		2 * atoms_alloc * sizeof(float4);
  // printf("allocating %d MB of memory on GPU\n", totalmem >> 20);
  printf("allocating %d MB of memory for block flags\n",
				(block_flags_alloc * 4) >> 20);
#endif

  cudaMalloc((void**) &block_flags, block_flags_alloc * 4);
  cudaMalloc((void**) &energy_gbis_buffers, patch_pairs_alloc * sizeof(float));
  cudaMalloc((void**) &virial_buffers, patch_pairs_alloc * 16*sizeof(float));
  cudaMalloc((void**) &slow_virial_buffers, patch_pairs_alloc * 16*sizeof(float));
  cudaMalloc((void**) &patch_pairs, patch_pairs_alloc * sizeof(patch_pair));
  // cudaMalloc((void**) &virials, 2 * force_lists_alloc * 16*sizeof(float));
  // slow_virials = virials + force_lists_size * 16;
  cudaMalloc((void**) &force_lists, force_lists_alloc * sizeof(force_list));
  cudaMalloc((void**) &force_list_counters, force_lists_alloc * sizeof(unsigned int));
  cudaMalloc((void**) &GBIS_P1_counters, force_lists_alloc * sizeof(unsigned int));
  cudaMalloc((void**) &GBIS_P2_counters, force_lists_alloc * sizeof(unsigned int));
  cudaMalloc((void**) &GBIS_P3_counters, force_lists_alloc * sizeof(unsigned int));
  cudaMalloc((void**) &force_buffers, force_buffers_alloc * sizeof(float4));
  cudaMalloc((void**) &slow_force_buffers, force_buffers_alloc * sizeof(float4));
  cudaMalloc((void**) &atoms, atoms_alloc * sizeof(atom));
  cudaMalloc((void**) &atom_params, atoms_alloc * sizeof(atom_param));
  // cudaMalloc((void**) &forces, atoms_alloc * sizeof(float4));
  // cudaMalloc((void**) &slow_forces, atoms_alloc * sizeof(float4));
  cudaMalloc((void**) &intRad0D, atoms_alloc * sizeof(float));
  cudaMalloc((void**) &intRadSD, atoms_alloc * sizeof(float));
  //cudaMalloc((void**) &psiSumD, atoms_alloc * sizeof(GBReal));
  cudaMalloc((void**) &psiSumD_buffers, force_buffers_alloc * sizeof(GBReal));
  cudaMalloc((void**) &bornRadD, atoms_alloc * sizeof(float));
  //cudaMalloc((void**) &dEdaSumD, atoms_alloc * sizeof(GBReal));
  cudaMalloc((void**) &dEdaSumD_buffers, force_buffers_alloc * sizeof(GBReal));
  cudaMalloc((void**) &dHdrPrefixD, atoms_alloc * sizeof(float));
  cuda_errcheck("malloc everything");

 } //if sizes grew

  cudaMemcpy(patch_pairs, pp, npp * sizeof(patch_pair),
				cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to patch_pairs");

  cudaMemcpy(force_lists, fl, nfl * sizeof(force_list),
				cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to force_lists");

  cudaMemset(force_list_counters, 0, nfl * sizeof(unsigned int));
  cudaMemset(GBIS_P1_counters, 0, nfl * sizeof(unsigned int));
  cudaMemset(GBIS_P2_counters, 0, nfl * sizeof(unsigned int));
  cudaMemset(GBIS_P3_counters, 0, nfl * sizeof(unsigned int));
  cuda_errcheck("memset force_list_counters");
} // bind patch pairs

void cuda_bind_atom_params(const atom_param *t) {
  cudaMemcpyAsync(atom_params, t, atoms_size * sizeof(atom_param),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to atom_params");
}

void cuda_bind_atoms(const atom *a) {
  cuda_errcheck("before memcpy to atoms");
  cudaMemcpyAsync(atoms, a, atoms_size * sizeof(atom),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to atoms");
}

void cuda_bind_forces(float4 *f, float4 *f_slow) {
  cudaHostGetDevicePointer(&forces, f, 0);
  cuda_errcheck("cudaHostGetDevicePointer forces");
  cudaHostGetDevicePointer(&slow_forces, f_slow, 0);
  cuda_errcheck("cudaHostGetDevicePointer slow_forces");
}

void cuda_bind_virials(float *v) {
  cudaHostGetDevicePointer(&virials, v, 0);
  cuda_errcheck("cudaHostGetDevicePointer virials");
  slow_virials = virials + force_lists_size*16;
}

//GBIS bindings
void cuda_bind_GBIS_energy(float *e) {
  cudaHostGetDevicePointer(&energy_gbis, e, 0);
  cuda_errcheck("cudaHostGetDevicePointer energy_gbis");
}
void cuda_bind_GBIS_intRad(float *intRad0H, float *intRadSH) {
  cudaMemcpyAsync(intRad0D, intRad0H, atoms_size * sizeof(float),
				cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(intRadSD, intRadSH, atoms_size * sizeof(float),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to intRad");
}

void cuda_bind_GBIS_psiSum(GBReal *psiSumH) {
  cudaHostGetDevicePointer(&psiSumD, psiSumH, 0);
  cuda_errcheck("cudaHostGetDevicePointer psiSum");
}

void cuda_bind_GBIS_bornRad(float *bornRadH) {
  cudaMemcpyAsync(bornRadD, bornRadH, atoms_size * sizeof(float),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to bornRad");
}

void cuda_bind_GBIS_dEdaSum(GBReal *dEdaSumH) {
  cudaHostGetDevicePointer(&dEdaSumD, dEdaSumH, 0);
  cuda_errcheck("cudaHostGetDevicePointer dEdaSum");
}

void cuda_bind_GBIS_dHdrPrefix(float *dHdrPrefixH) {
  cudaMemcpyAsync(dHdrPrefixD, dHdrPrefixH, atoms_size * sizeof(float),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to dHdrPrefix");
}
// end GBIS methods

#if 0
void cuda_load_forces(float4 *f, float4 *f_slow, int begin, int count) {
  // printf("load forces %d %d %d\n",begin,count,atoms_size);
  cudaMemcpyAsync(f+begin, forces+begin, count * sizeof(float4),
				cudaMemcpyDeviceToHost, stream);
  if ( f_slow ) {
    cudaMemcpyAsync(f_slow+begin, slow_forces+begin, count * sizeof(float4),
				cudaMemcpyDeviceToHost, stream);
  }
  cuda_errcheck("memcpy from forces");
}

void cuda_load_virials(float *v, int doSlow) {
  int count = force_lists_size;
  if ( doSlow ) count *= 2;
  cudaMemcpyAsync(v, virials, count * 16*sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cuda_errcheck("memcpy from virials");
}
#endif

#if 0
__host__ __device__ static int3 patch_coords_from_id(
        dim3 PATCH_GRID, int id) {

  return make_int3( id % PATCH_GRID.x,
                ( id / PATCH_GRID.x ) % PATCH_GRID.y,
                id / ( PATCH_GRID.x * PATCH_GRID.y ) );
}

__host__ __device__ static int patch_id_from_coords(
        dim3 PATCH_GRID, int3 coords) {

  // handles periodic boundaries
  int x = (coords.x + 4 * PATCH_GRID.x) % PATCH_GRID.x;
  int y = (coords.y + 4 * PATCH_GRID.y) % PATCH_GRID.y;
  int z = (coords.z + 4 * PATCH_GRID.z) % PATCH_GRID.z;

  return ( z * PATCH_GRID.y + y ) * PATCH_GRID.x + x;
}

__host__ __device__ static int3 patch_offset_from_neighbor(int neighbor) {

  // int3 coords = patch_coords_from_id(make_uint3(3,3,3), 13 + neighbor);
  int3 coords = patch_coords_from_id(make_uint3(3,3,3), neighbor);
  return make_int3(coords.x - 1, coords.y - 1, coords.z - 1);

}
#endif
 
#define BLOCK_SIZE 128
#define SHARED_SIZE 32


#define MAKE_PAIRLIST
#define DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef MAKE_PAIRLIST
#define DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"


void cuda_nonbonded_forces(float3 lata, float3 latb, float3 latc,
		float cutoff2, float plcutoff2,
		int cbegin, int ccount, int pbegin, int pcount,
		int doSlow, int doEnergy, int usePairlists, int savePairlists,
		cudaStream_t &strm) {

 if ( ccount ) {
   if ( usePairlists ) {
     if ( ! savePairlists ) plcutoff2 = 0.;
   } else {
     plcutoff2 = cutoff2;
   }
   int grid_dim = 65535;  // maximum allowed
   for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
     if ( grid_dim > ccount - cstart ) grid_dim = ccount - cstart;
     // printf("%d %d %d\n",cbegin+cstart,grid_dim,patch_pairs_size);

#define CALL(X) X<<< grid_dim, BLOCK_SIZE, 0, strm \
	>>>(patch_pairs+cbegin+cstart,atoms,atom_params,force_buffers, \
	     (doSlow?slow_force_buffers:0), block_flags, \
             virial_buffers, (doSlow?slow_virial_buffers:0), \
             overflow_exclusions, force_list_counters, force_lists, \
             forces, virials, \
             (doSlow?slow_forces:0), (doSlow?slow_virials:0), \
             lj_table_size, \
	     lata, latb, latc, cutoff2, plcutoff2, doSlow)
//end definition

     if ( doEnergy ) {
       if ( doSlow ) {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_slow_energy_pairlist);
         else CALL(dev_nonbonded_slow_energy);
       } else {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_energy_pairlist);
         else CALL(dev_nonbonded_energy);
       }
     } else {
       if ( doSlow ) {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_slow_pairlist);
         else CALL(dev_nonbonded_slow);
       } else {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_pairlist);
         else CALL(dev_nonbonded);
       }
     }

     cuda_errcheck("dev_nonbonded");
   }
 }

#if 0
 if ( pcount ) {
  // printf("%d %d %d\n",pbegin,pcount,force_lists_size);
  dev_sum_forces<<< pcount, BLOCK_SIZE, 0, stream
	>>>(atoms,force_lists+pbegin,force_buffers,
                virial_buffers,forces,virials);
  if ( doSlow ) {
    dev_sum_forces<<< pcount, BLOCK_SIZE, 0, stream
	>>>(atoms,force_lists+pbegin,slow_force_buffers,
                slow_virial_buffers,slow_forces,slow_virials);
  }
  cuda_errcheck("dev_sum_forces");
 }
#endif

}

//import GBIS Kernel definitions
#include "ComputeGBISCUDAKernel.h"

//////////////////////////////////////////
//  GBIS P1
//////////////////////////////////////////
void cuda_GBIS_P1(
	int cbegin,
  int ccount,
  int pbegin,
  int pcount,
  float a_cut,
  float rho_0,
  float3 lata,
  float3 latb,
  float3 latc,
  cudaStream_t &strm
) {

  int grid_dim = 65535;  // maximum allowed
  for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
    if (grid_dim > ccount - cstart) {
      grid_dim = ccount - cstart;
    }

    GBIS_P1_Kernel<<<grid_dim, BLOCK_SIZE, 0, strm>>>(
      patch_pairs+cbegin+cstart,
      atoms,
      atom_params,
      intRad0D,
      intRadSD,
      psiSumD_buffers,
      psiSumD,
      a_cut,
      rho_0,
      lata,
      latb,
      latc,
      force_lists,
      GBIS_P1_counters 
      );
    cuda_errcheck("dev_GBIS_P1");
  } // end for
} // end GBIS P1

//////////////////////////////////////////
//  GBIS P2
//////////////////////////////////////////
void cuda_GBIS_P2(
	int cbegin,
  int ccount,
  int pbegin,
  int pcount,
  float a_cut,
  float r_cut,
  float scaling,
  float kappa,
  float smoothDist,
  float epsilon_p,
  float epsilon_s,
  float3 lata,
  float3 latb,
  float3 latc,
  int doEnergy,
  int doFullElec,
  cudaStream_t &strm
) {
  int grid_dim = 65535;  // maximum allowed
  for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
    if (grid_dim > ccount - cstart)
      grid_dim = ccount - cstart;

    GBIS_P2_Kernel<<<grid_dim, BLOCK_SIZE, 0, strm>>>(
      patch_pairs+cbegin+cstart,
      atoms,
      atom_params,
      bornRadD,
      dEdaSumD_buffers,
      dEdaSumD,
      a_cut,
      r_cut,
      scaling,
      kappa,
      smoothDist,
      epsilon_p,
      epsilon_s,
      lata,
      latb,
      latc,
      doEnergy,
      doFullElec,
      force_lists,
      force_buffers,
      forces,
      energy_gbis_buffers,
      energy_gbis,
      GBIS_P2_counters 
      );
    cuda_errcheck("dev_GBIS_P2");
  } // end for
} // end P2

//////////////////////////////////////////
//  GBIS P3
//////////////////////////////////////////
void cuda_GBIS_P3(
	int cbegin,
  int ccount,
  int pbegin,
  int pcount,
  float a_cut,
  float rho_0,
  float scaling,
  float3 lata,
  float3 latb,
  float3 latc,
  cudaStream_t &strm
) {
  int grid_dim = 65535;  // maximum allowed
  for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
    if (grid_dim > ccount - cstart)
      grid_dim = ccount - cstart;

    GBIS_P3_Kernel<<<grid_dim, BLOCK_SIZE, 0, strm>>>(
      patch_pairs+cbegin+cstart,
      atoms,
      atom_params,
      intRad0D,
      intRadSD,
      dHdrPrefixD,
      a_cut,
      rho_0,
      scaling,
      lata,
      latb,
      latc,
      force_lists,
      slow_force_buffers,
      slow_forces,
      GBIS_P3_counters 
      );
    cuda_errcheck("dev_GBIS_P3");
  }
}

#if 0
int cuda_stream_finished() {
  return ( cudaStreamQuery(stream) == cudaSuccess );
}
#endif


#else  // NAMD_CUDA

// for make depends
#include "ComputeNonbondedCUDAKernelBase.h"

#endif  // NAMD_CUDA

