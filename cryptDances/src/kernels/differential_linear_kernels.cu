#include "differential_linear_kernels.cuh"
#include <inttypes.h>
#include "algorithms/chacha.cuh"
#include "types.cuh"
#include "arx_cryptanalysis.cuh"
#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>

#define U32C(v) (v##U)

#define U32V(v) ((uint32_t)(v) &U32C(0xFFFFFFFF))

//#define ROTATE(v, c) ((v<<c)^(v>>(32-c)))
#define ROT(v, c) ((v<<c)^(v>>(32-c)))
#define XOR(v, w) (((v) ^ (w)))
#define XOR1(v) (((v) ^ (33554464)))
#define XOR2(v) (((v) ^ (268439552)))
#define XOR3(v) (((v) ^ (35651584)))
#define XOR4(v) (((v) ^ (2105344)))
#define UNO 1
#define PLUS(v, w) (((v) + (w)))
#define ADD(v, w) (((v) + (w)))
#define MINUS(v, w) (((v) - (w)))
#define PLUSONE(v,c) ((v<<c)^(v>>(32-c)))


__global__ void setup_kernel(curandState *state, int my_rank)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(my_rank, id, 0, &state[id]);
}

#define DOUBLE_ROUND(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15) \
    /* Diagonal round */                                                              \
      x2 = ADD( x2,  x7);  x3 = ADD( x3,  x4); \
  x15 = XOR(x15,  x0);  x13 = XOR(x13,  x2); x14 = XOR(x14,  x3); \
  x15 = ROT(x15,  16);  x13 = ROT(x13,  16); x14 = ROT(x14,  16); \
  x10 = ADD(x10, x15); x11 = ADD(x11, x12);  x8 = ADD( x8, x13);  x9 = ADD( x9, x14); \
   x5 = XOR( x5, x10);  x6 = XOR( x6, x11);  x7 = XOR( x7,  x8);  x4 = XOR( x4,  x9); \
   x5 = ROT( x5,  12);  x6 = ROT( x6,  12);  x7 = ROT( x7,  12);  x4 = ROT( x4,  12); \
   x0 = ADD( x0,  x5);  x1 = ADD( x1,  x6);  x2 = ADD( x2,  x7);  x3 = ADD( x3,  x4); \
  x15 = XOR(x15,  x0); x12 = XOR(x12,  x1); x13 = XOR(x13,  x2); x14 = XOR(x14,  x3); \
  x15 = ROT(x15,   8); x12 = ROT(x12,   8); x13 = ROT(x13,   8); x14 = ROT(x14,   8); \
  x10 = ADD(x10, x15); x11 = ADD(x11, x12);  x8 = ADD( x8, x13);  x9 = ADD( x9, x14); \
   x5 = XOR( x5, x10);  x6 = XOR( x6, x11);  x7 = XOR( x7,  x8);  x4 = XOR( x4,  x9); \
   x5 = ROT( x5,   7);  x6 = ROT( x6,   7);  x7 = ROT( x7,   7);  x4 = ROT( x4,   7); \
   x0 = ADD( x0,  x4);  x1 = ADD( x1,  x5);  x2 = ADD( x2,  x6);  x3 = ADD( x3,  x7); \
  x12 = XOR(x12,  x0); x13 = XOR(x13,  x1); x14 = XOR(x14,  x2); x15 = XOR(x15,  x3); \
  x12 = ROT(x12,  16); x13 = ROT(x13,  16); x14 = ROT(x14,  16); x15 = ROT(x15,  16); \
   x8 = ADD( x8, x12);  x9 = ADD( x9, x13); x10 = ADD(x10, x14); x11 = ADD(x11, x15); \
   x4 = XOR( x4,  x8);  x5 = XOR( x5,  x9);  x6 = XOR( x6, x10);  x7 = XOR( x7, x11); \
   x4 = ROT( x4,  12);  x5 = ROT( x5,  12);  x6 = ROT( x6,  12);  x7 = ROT( x7,  12); \
   x0 = ADD( x0,  x4);  x1 = ADD( x1,  x5);  x2 = ADD( x2,  x6);  x3 = ADD( x3,  x7); \
  x12 = XOR(x12,  x0); x13 = XOR(x13,  x1);  x15 = XOR(x15,  x3); \
  x12 = ROT(x12,   8); x13 = ROT(x13,   8);  x15 = ROT(x15,   8); \
   x8 = ADD( x8, x12);    x11 = ADD(x11, x15); \
         x7 = XOR( x7, x11); \
       x7 = ROT( x7,   7); \
  /* x0 = ADD( x0,  x5);  x1 = ADD( x1,  x6);  x2 = ADD( x2,  x7);  x3 = ADD( x3,  x4); \
  x15 = XOR(x15,  x0); x12 = XOR(x12,  x1); x13 = XOR(x13,  x2); x14 = XOR(x14,  x3); \
  x15 = ROT(x15,  16); x12 = ROT(x12,  16); x13 = ROT(x13,  16); x14 = ROT(x14,  16); \
  x10 = ADD(x10, x15); x11 = ADD(x11, x12);  x8 = ADD( x8, x13);  x9 = ADD( x9, x14); \
   x5 = XOR( x5, x10);  x6 = XOR( x6, x11);  x7 = XOR( x7,  x8);  x4 = XOR( x4,  x9); \
   x5 = ROT( x5,  12);  x6 = ROT( x6,  12);  x7 = ROT( x7,  12);  x4 = ROT( x4,  12); \*/

#define DOUBLE_ROUND_ALT(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15) \
    /* Diagonal round */ \
    x3 = XOR1(x3); x7 = XOR2(x7); x11 = XOR3(x11); x15 = XOR4(x15);\
    x2 = ADD( x2,  x7);  x3 = ADD( x3,  x4); \
  x15 = XOR(x15,  x0);  x13 = XOR(x13,  x2); x14 = XOR(x14,  x3); \
  x15 = ROT(x15,  16);  x13 = ROT(x13,  16); x14 = ROT(x14,  16); \
  x10 = ADD(x10, x15); x11 = ADD(x11, x12);  x8 = ADD( x8, x13);  x9 = ADD( x9, x14); \
   x5 = XOR( x5, x10);  x6 = XOR( x6, x11);  x7 = XOR( x7,  x8);  x4 = XOR( x4,  x9); \
   x5 = ROT( x5,  12);  x6 = ROT( x6,  12);  x7 = ROT( x7,  12);  x4 = ROT( x4,  12); \
   x0 = ADD( x0,  x5);  x1 = ADD( x1,  x6);  x2 = ADD( x2,  x7);  x3 = ADD( x3,  x4); \
  x15 = XOR(x15,  x0); x12 = XOR(x12,  x1); x13 = XOR(x13,  x2); x14 = XOR(x14,  x3); \
  x15 = ROT(x15,   8); x12 = ROT(x12,   8); x13 = ROT(x13,   8); x14 = ROT(x14,   8); \
  x10 = ADD(x10, x15); x11 = ADD(x11, x12);  x8 = ADD( x8, x13);  x9 = ADD( x9, x14); \
   x5 = XOR( x5, x10);  x6 = XOR( x6, x11);  x7 = XOR( x7,  x8);  x4 = XOR( x4,  x9); \
   x5 = ROT( x5,   7);  x6 = ROT( x6,   7);  x7 = ROT( x7,   7);  x4 = ROT( x4,   7); \
   x0 = ADD( x0,  x4);  x1 = ADD( x1,  x5);  x2 = ADD( x2,  x6);  x3 = ADD( x3,  x7); \
  x12 = XOR(x12,  x0); x13 = XOR(x13,  x1); x14 = XOR(x14,  x2); x15 = XOR(x15,  x3); \
  x12 = ROT(x12,  16); x13 = ROT(x13,  16); x14 = ROT(x14,  16); x15 = ROT(x15,  16); \
   x8 = ADD( x8, x12);  x9 = ADD( x9, x13); x10 = ADD(x10, x14); x11 = ADD(x11, x15); \
   x4 = XOR( x4,  x8);  x5 = XOR( x5,  x9);  x6 = XOR( x6, x10);  x7 = XOR( x7, x11); \
   x4 = ROT( x4,  12);  x5 = ROT( x5,  12);  x6 = ROT( x6,  12);  x7 = ROT( x7,  12); \
   x0 = ADD( x0,  x4);  x1 = ADD( x1,  x5);  x2 = ADD( x2,  x6);  x3 = ADD( x3,  x7); \
  x12 = XOR(x12,  x0); x13 = XOR(x13,  x1);  x15 = XOR(x15,  x3); \
  x12 = ROT(x12,   8); x13 = ROT(x13,   8);  x15 = ROT(x15,   8); \
   x8 = ADD( x8, x12);    x11 = ADD(x11, x15); \
       x7 = XOR( x7, x11); \
       x7 = ROT( x7,   7); \
   /* x0 = ADD( x0,  x5);  x1 = ADD( x1,  x6);  x2 = ADD( x2,  x7);  x3 = ADD( x3,  x4); \
 x15 = XOR(x15,  x0); x12 = XOR(x12,  x1); x13 = XOR(x13,  x2); x14 = XOR(x14,  x3); \
  x15 = ROT(x15,  16); x12 = ROT(x12,  16); x13 = ROT(x13,  16); x14 = ROT(x14,  16); \
  x10 = ADD(x10, x15); x11 = ADD(x11, x12);  x8 = ADD( x8, x13);  x9 = ADD( x9, x14); \
   x5 = XOR( x5, x10);  x6 = XOR( x6, x11);  x7 = XOR( x7,  x8);  x4 = XOR( x4,  x9); \
   x5 = ROT( x5,  12);  x6 = ROT( x6,  12);  x7 = ROT( x7,  12);  x4 = ROT( x4,  12); \*/

#define DOUBLE_ROUND_128(v0,v1,v2,v3) \
v0 += v1; v3 ^= v0; v3 <<= (16, 16, 16, 16); \
v2 += v3; v1 ^= v2; v1 <<= (12, 12, 12, 12); \
v0 += v1; v3 ^= v0; v3 <<= ( 8, 8, 8, 8); \
v2 += v3; v1 ^= v2; v1 <<= ( 7, 7, 7, 7); \
v1 >>= 32; v2 >>= 64; v3 >>= 96; \
v0 += v1; v3 ^= v0; v3 <<= (16, 16, 16, 16); \
v2 += v3; v1 ^= v2; v1 <<= (12, 12, 12, 12); \
v0 += v1; v3 ^= v0; v3 <<= ( 8, 8, 8, 8); \
v2 += v3; v1 ^= v2; v1 <<= ( 7, 7, 7, 7); \
v1 <<= 32; v2 <<= 64; v3 <<= 96; \

void get_differential_from_position(int position, differential_t *diff)
{
    uint32_t iv_positions[4];

    diff->input.words[0] = position/(32*512);
    position -= diff->input.words[0] * 32 * 512;

    diff->input.bits[0] = position/512;
    position -= diff->input.bits[0] * 512;

    diff->output.words[0] = position/32;
    position -= diff->output.words[0] * 32;

    diff->output.bits[0] = position;

    get_alg_iv_positions(iv_positions, diff->alg_type);
    diff->input.words[0] = iv_positions[diff->input.words[0]];

    diff->input.number_of_bits = 1;
    diff->output.number_of_bits = 1;

    lob_compute_mask_from_list_of_bits(&(diff->input));
    lob_compute_mask_from_list_of_bits(&(diff->output));
}

/*
    Computes the differential correlation for every output bit for every 128 possible single input differentials.
    The test executes I iterations for n keys.

    The organization is the following:
        - blockIdx.y - indicate which word for the input differential, from 0 to 3.
        - blockIdx.z - indicate which bit for the input differential, from 0 to 31.
        - blockIdx.x * thread.x * ntest_for_each_key - number of tests

 * . For increased performance we tried to minimize all unnecessary operations,
 * including some that were not very straightforward. 
 * . For instance, the initialization of the state is bypassed and a completely
 * random state is considered. This is useful to avoid copying memory from one place to another.
 * . Additionally, the CUDA PRNG is called only once to initialize the first state.
 * After that, the previous resulting state is considered as a starting point since it is itself random.
 * . These choices do not seem to affect the results, as the replication of previous papers shows.
 * . Therefore, this kernel is not useful for researchers that are trying to analyze some particular
 * behaviours of the constants of the algorithms.
*/
__global__ void differential_correlation_exhaustion_kernel(unsigned long long seed,
int subrounds, int last_subround, int n_test_for_each_thread, unsigned long long *d_result, int alg_type)
{
    uint32_t id[STATE_SIZE] = { 0 };
    algorithm alg;
    uint32_t observed_od[STATE_SIZE];
    uint8_t observed_od_bits[STATE_SIZE_IN_BITS];
    curandState_t rng;
    int sum[STATE_SIZE_IN_BITS] = { 0 };
    uint32_t state[STATE_SIZE] = { 0 }, alt_state[STATE_SIZE] = { 0 };

    int word = blockIdx.y, bit = blockIdx.z;
    const unsigned long long blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
    const unsigned long long tid = blockId * blockDim.x + threadIdx.x;

    define_alg(&alg, alg_type);
    curand_init(seed, tid, 0, &rng);

    //comput id - each block may test a different id
    id[alg.iv_positions[word]] = 1 << bit;

    GENERATE_RANDOM_STATE(state);
    for (int t = 0; t < n_test_for_each_thread; t++)
    {
        xor_array(alt_state, state, id, STATE_SIZE);
        alg.subrounds(state, subrounds,last_subround);
        alg.subrounds(alt_state, subrounds,last_subround);
        xor_array(observed_od, state, alt_state, STATE_SIZE);
        transform_state_to_bits(observed_od, observed_od_bits);
        update_result(sum, observed_od_bits);
    }

    for (int i = 0; i < 512; i++)
        atomicAdd(&d_result[32 * 512 * word + 512 * bit + i], sum[i]);
}


void compute_all_single_bit_differential_correlation(int alg_type, int subrounds, int last_subround,
uint64_t number_of_trials, const char *out_file_name)
{
    unsigned long long int *d_results;
    uint64_t numblocks = NUMBER_OF_CUDA_BLOCKS/(4*32);
    uint64_t iterations = number_of_trials / NUMBER_OF_CUDA_THREADS / numblocks / NUMBER_OF_TESTS_PER_THREAD / (num_procs-1);
    unsigned long long int h_results[NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS] = { 0 }, seed;
    uint64_t acc_results[NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS] = {0}, result[NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS] = {0};
    int L = NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS * sizeof(unsigned long long int);

    srand_by_rank(); //initialize prng with different internal state for each MPI process

    memset(result,0x00,NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS*sizeof(uint64_t));
    memset(acc_results,0x00,NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS*sizeof(uint64_t));

    if (my_rank != 0)
    {
        cudaSetDevice((my_rank-1)%NUMBER_OF_DEVICES_PER_MACHINE);

        cudaMalloc((void **)&d_results, L);
        cudaMemcpy(d_results, h_results, L, cudaMemcpyHostToDevice);

        for (int i = 0; i < iterations; i++)
        {
            memset(h_results, 0, L);
            cudaMemcpy(d_results, h_results, L, cudaMemcpyHostToDevice);
            seed = seed_by_rank();

            differential_correlation_exhaustion_kernel <<<dim3(numblocks, 4, 32), NUMBER_OF_CUDA_THREADS>>> (seed,
                subrounds, last_subround, NUMBER_OF_TESTS_PER_THREAD, d_results, alg_type);

            cudaMemcpy(h_results, d_results, L, cudaMemcpyDeviceToHost);

            for(int j=0;j<NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS; j++)
                acc_results[j]+= (uint64_t) h_results[j];
        }

        cudaFree(d_results);
    }

    MPI_Allreduce(&acc_results, &result, NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

    if(my_rank == 0)
    {
        differential_t *diff = NULL;
        diff = (differential_t * ) malloc(sizeof(differential_t) * NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS);
        if(diff == NULL)
        {
            printf("Not enough memory\n");
            return;
        }
        for(int position = 0; position < NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS; position++)
        {
            memset(&diff[position],0x00, sizeof(differential_t));
            diff[position].alg_type = alg_type;
            get_differential_from_position(position, &diff[position]);
            diff[position].correlation.number_of_trials = number_of_trials;
            diff[position].correlation.correlation_count = number_of_trials - result[position];
            ct_compute_and_test_correlation(&(diff[position].correlation));
        }

        update_single_bit_differentials_from_file(out_file_name, diff);

        free(diff);
    }
}

/*
Computes the linear bias given the number of rounds, ID and OD mask.
*/
__global__ void linear_correlation_kernel(unsigned long long seed, int subrounds,
                                          int last_subround, uint32_t *id_mask, uint32_t *od_mask, int n_test_for_each_thread,
                                          unsigned long long int *d_result, int alg_type)
{
    algorithm alg;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t input_state[STATE_SIZE] = { 0 }, output_state[STATE_SIZE] = { 0 };
    curandState_t rng;
    unsigned long long int sum_parity = 0;

    define_alg(&alg, alg_type);
    curand_init(seed, tid, 0, &rng);

    GENERATE_RANDOM_STATE(output_state);
    for (int t = 0; t < n_test_for_each_thread; t++)
    {
        for(int i=0;i<STATE_SIZE;i++)
            input_state[i] = output_state[i];
        alg.subrounds(output_state, subrounds, last_subround);
        sum_parity += check_parity_of_linear_relation(id_mask, input_state, od_mask, output_state);
    }

    atomicAdd(d_result, sum_parity);
}
/**
 * Computes the differential bias given the number of subrounds, ID and OD mask. Some important observations:
 * . For increased performance we tried to minimize all unnecessary operations,
 * including some that were not very straightforward. 
 * . For instance, the initialization of the state is bypassed and a completely
 * random state is considered. This is useful to avoid copying memory from one place to another.
 * . Additionally, the CUDA PRNG is called only once to initialize the first state.
 * After that, the previous resulting state is considered as a starting point since it is itself random.
 * . These choices do not seem to affect the results, as the replication of previous papers shows.
 * . Therefore, this kernel is not useful for researchers that are trying to analyze some particular
 * behaviours of the constants of the algorithms.
 * */
#include <thrust/random.h>
#include <curand.h>

__global__ void differential_correlation_kernel(uint32_t *id,
    uint32_t *od, int n_test_for_each_thread, unsigned long long int *d_result, int alg_type, curandState *global_state, uint32_t *d_od_non_zero_words, uint32_t *d_od_non_zero_bits, int size_non_zero_bits)
{
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t state[STATE_SIZE] = { 0 }, alt_state[STATE_SIZE] = { 0 }, observed_od[STATE_SIZE] = {0};
    curandState_t  rng;
    unsigned long long int sum_parity = 0;

    rng = global_state[tid];

    for (int t = 0; t < n_test_for_each_thread; t++)
    {
        for (int i = 0; i < 16; i++) {
            state[i] = curand(&rng);
            alt_state[i] = state[i];
        }

        state[0] = ADD( state[0],  state[5]);
        alt_state[0] = state[0];
        state[1] = ADD( state[1],  state[6]);
        alt_state[1] = state[1];

        state[12] = XOR(state[12],  state[1]);
        state[12] = ROT(state[12],  16);
        alt_state[12] = state[12];

        DOUBLE_ROUND(state[0],state[1],state[2],state[3],state[4],state[5],state[6],state[7],state[8],state[9],state[10],state[11],state[12],state[13],state[14],state[15]);
        DOUBLE_ROUND_ALT(alt_state[0],alt_state[1],alt_state[2],alt_state[3],alt_state[4],alt_state[5],alt_state[6],alt_state[7],alt_state[8],alt_state[9],alt_state[10],alt_state[11],alt_state[12],alt_state[13],alt_state[14],alt_state[15]);
        xor_array_dist2(observed_od, state, alt_state, STATE_SIZE, d_od_non_zero_words, d_od_non_zero_bits, size_non_zero_bits);
        sum_parity += check_parity_of_equation_dist2(observed_od, od, d_od_non_zero_words, d_od_non_zero_bits, size_non_zero_bits);

    }
    global_state[tid] = rng;
    atomicAdd(&d_result[blockIdx.x], sum_parity);
}



void compute_differential_or_linear_correlation(diff_lin_t *diff_lin, int type)
{
    uint64_t result = 0, seed, local_sum=0;
    uint64_t iterations;
    unsigned long long int * d_sum_parity;
    uint32_t *d_id, *d_od;
    uint32_t *d_od_non_zero_words, *d_od_non_zero_bits;
    unsigned long long int * local_sum_parity = ( unsigned long long int *)malloc(NUMBER_OF_CUDA_BLOCKS * sizeof( unsigned long long int));
    /*for (int i=0;i<NUMBER_OF_CUDA_BLOCKS;i++)
        local_sum_parity[i] = 0;*/

    int subrounds = diff_lin->output.subround - diff_lin->input.subround;
    int last_subround = diff_lin->input.subround;

    srand_by_rank(); //initialize prng with different internal state for each MPI process
    iterations = diff_lin->correlation.number_of_trials / NUMBER_OF_TESTS_PER_THREAD / NUMBER_OF_CUDA_THREADS / NUMBER_OF_CUDA_BLOCKS / (num_procs);
    cudaSetDevice((my_rank)%NUMBER_OF_DEVICES_PER_MACHINE);
    cudaMalloc((void **)&d_sum_parity, NUMBER_OF_CUDA_BLOCKS * sizeof(unsigned long long int));
    cudaMalloc(&d_id, STATE_SIZE * sizeof(uint32_t));
    cudaMalloc(&d_od, STATE_SIZE * sizeof(uint32_t));
    cudaMalloc(&d_od_non_zero_words, diff_lin->output.number_of_bits * sizeof(uint32_t));
    cudaMalloc(&d_od_non_zero_bits, diff_lin->output.number_of_bits * sizeof(uint32_t));
    cudaMemcpy(d_id, diff_lin->input.mask, STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_od, diff_lin->output.mask, STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_od_non_zero_words, diff_lin->output.words, diff_lin->output.number_of_bits * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_od_non_zero_bits, diff_lin->output.bits, diff_lin->output.number_of_bits * sizeof(uint32_t), cudaMemcpyHostToDevice);

    curandState_t *devStates;

    cudaMalloc((void **)&devStates, NUMBER_OF_CUDA_BLOCKS * NUMBER_OF_CUDA_THREADS *
                                    sizeof(curandState_t));

    setup_kernel<<<NUMBER_OF_CUDA_BLOCKS, NUMBER_OF_CUDA_THREADS>>>(devStates, my_rank);
    //cudaDeviceSynchronize();


    for (int i = 0; i < iterations; i++)
    {
        memset(local_sum_parity, 0, NUMBER_OF_CUDA_BLOCKS *
                                    sizeof(unsigned long long int));
        cudaMemset(d_sum_parity, 0, NUMBER_OF_CUDA_BLOCKS *
        sizeof(unsigned long long int));
        if(type == TYPE_DIFFERENTIAL) {
            differential_correlation_kernel <<< NUMBER_OF_CUDA_BLOCKS, NUMBER_OF_CUDA_THREADS >>>(d_id, d_od,
                                                                                                  NUMBER_OF_TESTS_PER_THREAD,
                                                                                                  d_sum_parity,
                                                                                                  diff_lin->alg_type,
                                                                                                  devStates, d_od_non_zero_words, d_od_non_zero_bits,
                                                                                                  diff_lin->output.number_of_bits);


        }else
            linear_correlation_kernel <<< NUMBER_OF_CUDA_BLOCKS, NUMBER_OF_CUDA_THREADS >>> ((unsigned long long)seed,subrounds, last_subround, d_id, d_od, NUMBER_OF_TESTS_PER_THREAD, d_sum_parity, diff_lin->alg_type);

        cudaMemcpyAsync(local_sum_parity,d_sum_parity,
                   NUMBER_OF_CUDA_BLOCKS *sizeof(unsigned long long int), cudaMemcpyDeviceToHost);

        for (int j = 0;j<NUMBER_OF_CUDA_BLOCKS;j++) {
            local_sum += (uint64_t) local_sum_parity[j];
        }
    }

    cudaFree(d_sum_parity);
    cudaFree(d_id);
    cudaFree(d_od);
    cudaFree(d_od_non_zero_bits);
    cudaFree(d_od_non_zero_words);
    cudaFree(devStates);
    free(local_sum_parity);


    MPI_Allreduce(&local_sum,&result,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);


    diff_lin->correlation.correlation_count = diff_lin->correlation.number_of_trials-result;
    ct_compute_and_test_correlation(&(diff_lin->correlation));
}


__global__ void differential_correlation_kernel_original(unsigned long long seed, int subrounds, int last_subround, uint32_t *id,
                                                uint32_t *od, int n_test_for_each_thread, unsigned long long int *d_result, int alg_type)
{
    algorithm alg;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t state[STATE_SIZE] = { 0 }, alt_state[STATE_SIZE] = { 0 }, observed_od[STATE_SIZE];
    curandState_t rng;
    unsigned long long int sum_parity = 0;

    define_alg(&alg, alg_type);
    curand_init(seed, tid, 0, &rng);

    GENERATE_RANDOM_STATE(state);
    for (int t = 0; t < n_test_for_each_thread; t++)
    {
        xor_array(alt_state, state, id, STATE_SIZE);
        alg.subrounds(state, subrounds,last_subround);
        alg.subrounds(alt_state, subrounds,last_subround);
        xor_array(observed_od, state, alt_state, STATE_SIZE);
        sum_parity += check_parity_of_equation(observed_od, od);
    }

    atomicAdd(d_result, sum_parity);
}


void compute_differential_or_linear_correlation_original(diff_lin_t *diff_lin, int type)
{
    uint64_t result = 0, seed, local_sum=0;
    uint64_t iterations;
    unsigned long long int *d_sum_parity;
    uint32_t *d_id, *d_od;
    unsigned long long int local_sum_parity = 0;

    int subrounds = diff_lin->output.subround - diff_lin->input.subround;
    int last_subround = diff_lin->input.subround;

    srand_by_rank(); //initialize prng with different internal state for each MPI process
    iterations = diff_lin->correlation.number_of_trials / TOTAL_EXECUTIONS_PER_KERNEL / (num_procs);

        cudaSetDevice((my_rank)%NUMBER_OF_DEVICES_PER_MACHINE);

        cudaMalloc(&d_sum_parity, sizeof(unsigned long long int));
        cudaMalloc(&d_id, STATE_SIZE * sizeof(uint32_t));
        cudaMalloc(&d_od, STATE_SIZE * sizeof(uint32_t));

        cudaMemcpy(d_id, diff_lin->input.mask, STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_od, diff_lin->output.mask, STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);

        for (int i = 0; i < iterations; i++)
        {
            seed = seed_by_rank();
            local_sum_parity = 0;
            cudaMemcpy(d_sum_parity, &local_sum_parity,
                       sizeof(unsigned long long int), cudaMemcpyHostToDevice);

            if(type == TYPE_DIFFERENTIAL)
                differential_correlation_kernel_original <<< NUMBER_OF_CUDA_BLOCKS, NUMBER_OF_CUDA_THREADS >>> ((unsigned long long)seed,subrounds, last_subround, d_id, d_od, NUMBER_OF_TESTS_PER_THREAD, d_sum_parity, diff_lin->alg_type);
            else
                linear_correlation_kernel <<< NUMBER_OF_CUDA_BLOCKS, NUMBER_OF_CUDA_THREADS >>> ((unsigned long long)seed,
                                                                                                 subrounds, last_subround, d_id, d_od, NUMBER_OF_TESTS_PER_THREAD, d_sum_parity, diff_lin->alg_type);
            cudaMemcpy(&local_sum_parity, d_sum_parity,
                       sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
            local_sum += (uint64_t) local_sum_parity;
        }

        cudaFree(d_sum_parity);
        cudaFree(d_id);
        cudaFree(d_od);


    MPI_Allreduce(&local_sum,&result,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);

    diff_lin->correlation.correlation_count = diff_lin->correlation.number_of_trials-result;
    ct_compute_and_test_correlation(&(diff_lin->correlation));
}

void search_until_find_correlation(diff_lin_t *diff_lin, int type)
{
    uint64_t count = 0, total = 0;
    int level = 34;
    double correlation = 0;

    diff_lin->correlation.number_of_trials = 1;
    diff_lin->correlation.number_of_trials <<= level;

    while (1)
    {
        compute_differential_or_linear_correlation(diff_lin, type);
        count += diff_lin->correlation.correlation_count;
        total += diff_lin->correlation.number_of_trials;

        correlation = ((double)count) / total;
        correlation = 2 * correlation - 1;

        if(test_significance_of_correlation(correlation, total))
            break;
        if(level>MAX_LEVEL) {
            correlation = 0;
            break;
        }
        diff_lin->correlation.number_of_trials = 1;
        diff_lin->correlation.number_of_trials <<= level; //first repeat previous level since 2^level + 2^level = 2^{level+1}
        level++;
    }

    diff_lin->correlation.observed = correlation;
    diff_lin->correlation.number_of_trials = total;
    diff_lin->correlation.correlation_count = count;
    diff_lin->correlation.is_significant = TRUE;
}