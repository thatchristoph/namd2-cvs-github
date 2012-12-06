#ifdef NAMD_CUDA
//this type defined in multiple files
typedef float GBReal;

void cuda_errcheck(const char *msg);

#ifndef __CUDACC__
#undef __align__(X)
#define __align__(X)
#endif

#define PATCH_PAIR_SIZE 16
#define PATCH_PAIR_USED 15

struct __align__(16) patch_pair {  // must be multiple of 16!
  float4 offset;
  unsigned int patch1_size;
  unsigned int patch2_size;
  unsigned int patch1_force_size;  // non-fixed atoms at start of list
  unsigned int patch1_atom_start;
  unsigned int patch2_atom_start;
  unsigned int patch1_force_start;
  unsigned int block_flags_start;
  unsigned int virial_start;  // virial output location padded to 16
  unsigned int patch1_force_list_index;
  unsigned int patch1_force_list_size;
  unsigned int patch2_force_size;  // used for fixed-atom energy
  unsigned int pad2;
};

#define FORCE_LIST_SIZE 8
#define FORCE_LIST_USED 8

struct __align__(16) force_list {  // must be multiple of 16!
  unsigned int force_list_start;  // beginning of compute output
  unsigned int force_list_size;  // number of computes for this patch
  unsigned int patch_size;  // real number of atoms in patch
  unsigned int patch_stride;  // padded number of atoms in patch
  unsigned int force_output_start;  // output array
  unsigned int atom_start;  // atom positions
  unsigned int virial_list_start;  // beginning of compute virial output
  unsigned int virial_output_start;  // virial output location padded to 16
};

struct __align__(16) atom {  // must be multiple of 16!
  float3 position;
  float charge;
};

struct __align__(16) atom_param {  // must be multiple of 16!
  int vdw_type;
  int index;
  int excl_index;
  int excl_maxdiff;  // maxdiff == 0 -> only excluded from self
};

#define COPY_ATOM( DEST, SOURCE ) { \
  DEST.position.x = SOURCE.position.x; \
  DEST.position.y = SOURCE.position.y; \
  DEST.position.z = SOURCE.position.z; \
  DEST.charge = SOURCE.charge; \
  }

#define COPY_PARAM( DEST, SOURCE ) { \
  DEST.sqrt_epsilon = SOURCE.sqrt_epsilon; \
  DEST.half_sigma = SOURCE.half_sigma; \
  DEST.index = SOURCE.index; \
  DEST.excl_index = SOURCE.excl_index; \
  DEST.excl_maxdiff = SOURCE.excl_maxdiff; \
  }

#define COPY_ATOM_TO_SHARED( ATOM, PARAM, SHARED ) { \
    COPY_ATOM( SHARED, ATOM ) \
    COPY_PARAM( SHARED, PARAM ) \
  }

#define COPY_ATOM_FROM_SHARED( ATOM, PARAM, SHARED ) { \
    COPY_ATOM( ATOM, SHARED ) \
    COPY_PARAM( PARAM, SHARED ) \
  }

// 2^11 ints * 2^5 bits = 2^16 bits = range of unsigned short excl_index
// 2^26 ints * 2^5 bits = 2^32 bits = range of int excl_index
#define MAX_EXCLUSIONS (1<<26)
#define MAX_CONST_EXCLUSIONS 2048  // cache size is 8k

void cuda_bind_exclusions(const unsigned int *t, int n);

void cuda_bind_lj_table(const float2 *t, int _lj_table_size);

// #define FORCE_TABLE_SIZE 512
// maximum size of CUDA array 1D texture reference is 2^13 = 8192
// #define FORCE_TABLE_SIZE 8192
// CUDA docs lie, older devices can only handle 4096
#define FORCE_TABLE_SIZE 4096

void cuda_bind_force_table(const float4 *t, const float4 *et);

void cuda_init();

void cuda_bind_patch_pairs(const patch_pair *pp, int npp,
			const force_list *fl, int nfl,
			int atoms_size_p, int force_buffers_size_p,
			int block_flags_size_p, int max_atoms_per_patch_p);

void cuda_bind_atom_params(const atom_param *t);

void cuda_bind_atoms(const atom *a);

void cuda_bind_forces(float4 *f, float4 *f_slow);

void cuda_bind_virials(float *v);

void cuda_nonbonded_forces(float3 lata, float3 latb, float3 latc,
                float cutoff2, float plcutoff2,
                int cbegin, int ccount, int pbegin, int pcount,
                int doSlow, int doEnergy, int usePairlists, int savePairlists,
		cudaStream_t &strm);

//GBIS methods
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
  );
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
  );
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
  );

void cuda_bind_GBIS_intRad(float *intRad0H, float *intRadSH);
void cuda_bind_GBIS_energy(float *energy_gbis);
void cuda_bind_GBIS_psiSum(GBReal *psiSumH);
void cuda_bind_GBIS_bornRad(float *bornRadH);
void cuda_bind_GBIS_dEdaSum(GBReal *dEdaSumH);
void cuda_bind_GBIS_dHdrPrefix(float *dHdrPrefixH);

//end GBIS methods

int cuda_stream_finished();

#endif  // NAMD_CUDA

