#define NUMBER_OF_DEVICES 8
#define NUMBER_OF_THREADS (1<<7)
#define NUMBER_OF_BLOCKS (1<<7)
#define NUMBER_OF_TEST_PER_THREAD (1<<15)


#define five_key_recovery 0

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <inttypes.h>
#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>

void random_uint32(uint32_t *x)
{
    *x = 0;
    *x |= (rand() & 0xFF);
    *x |= (rand() & 0xFF) << 8;
    *x |= (rand() & 0xFF) << 16;
    *x |= (rand() & 0xFF) << 24;
}

void random_uint64(uint64_t *x)
{
    *x = 0;
    *x |= (rand() & 0xFFFF);
    *x |= (((uint64_t)rand() & 0xFFFF)) << 16;
    *x |= ((uint64_t)(rand() & 0xFFFF)) << 32;
    *x |= ((uint64_t)(rand() & 0xFFFF)) << 48;
}

void random_uint32_array(uint32_t *x, uint32_t size)
{
    for (uint32_t i = 0; i < size; i++)
        random_uint32(&x[i]);
}

void transform_state_to_bits(uint32_t state[16], uint8_t bits[512])
{
    int count = 0;
    for (int i = 0; i < 16; i++)
    {
        for (int b = 0; b < 32; b++)
        {
            bits[count] = (state[i] >> b) & 1;
            count++;
        }
    }
}

typedef struct chacha_ctx chacha_ctx;

#define U32C(v) (v##U)

#define U32V(v) ((uint32_t)(v) &U32C(0xFFFFFFFF))

#define ROTATE(v, c) ((v<<c)^(v>>(32-c)))
#define XOR(v, w) ((v) ^ (w))
#define PLUS(v, w) (((v) + (w)))
#define MINUS(v, w) (((v) - (w)))
#define PLUSONE(v) (PLUS((v), 1))

#define QUARTERROUND(a, b, c, d) \
a = PLUS(a, b);              \
d = ROTATE(XOR(d, a), 16);   \
c = PLUS(c, d);              \
b = ROTATE(XOR(b, c), 12);   \
a = PLUS(a, b);              \
d = ROTATE(XOR(d, a), 8);    \
c = PLUS(c, d);              \
b = ROTATE(XOR(b, c), 7);

#define HALF_1_QUARTERROUND(a, b, c, d) \
a = PLUS(a, b);              \
d = ROTATE(XOR(d, a), 16);   \
c = PLUS(c, d);              \
b = ROTATE(XOR(b, c), 12);   \


#define HALF_2_QUARTERROUND(a, b, c, d) \
a = PLUS(a, b);              \
d = ROTATE(XOR(d, a), 8);    \
c = PLUS(c, d);              \
b = ROTATE(XOR(b, c), 7);

#define INVERT_QUARTERROUND(a,b,c,d)\
b = XOR(ROTATE(b,25), c); \
c = MINUS(c, d);              \
d = XOR(ROTATE(d,24), a); \
a = MINUS(a, b);              \
b = XOR(ROTATE(b,20), c); \
c = MINUS(c, d);              \
d = XOR(ROTATE(d,16), a); \
a = MINUS(a, b);

#define HALF_HALF_1_QUARTERROUND(a, b, c, d) \
a = PLUS(a, b);              \
d = ROTATE(XOR(d, a), 16);   \

#define INVERT_HALF_HALF_QUARTERROUND(a,b,c,d)\
d = XOR(ROTATE(d,16), a); \
a = MINUS(a, b);

#define LOAD32_LE(v) (*((uint32_t *) (v)))
#define STORE32_LE(c,x) (memcpy(c,&x,4))

__host__ __device__ void chacha_init(uint32_t state[16], uint32_t k[8], uint32_t nonce[2], uint32_t ctr[2])
{
    state[0] = U32C(0x61707865);
    state[1] = U32C(0x3320646e);
    state[2] = U32C(0x79622d32);
    state[3] = U32C(0x6b206574);
    state[4] = k[0];
    state[5] = k[1];
    state[6] = k[2];
    state[7] = k[3];
    state[8] = k[4];
    state[9] = k[5];
    state[10] = k[6];
    state[11] = k[7];
    state[12] = ctr[0];
    state[13] = ctr[1];
    state[14] = nonce[0];
    state[15] = nonce[1];
}

__host__ __device__ void chacha_odd_round(uint32_t x[16])
{
    QUARTERROUND(x[0], x[4], x[8], x[12])
    QUARTERROUND(x[1], x[5], x[9], x[13])
    QUARTERROUND(x[2], x[6], x[10], x[14])
    QUARTERROUND(x[3], x[7], x[11], x[15])
}

__host__ __device__ void pretty_print_state(uint32_t * state, int state_len) {
    for (int i=0;i<state_len;i++) {
        printf("0x%08x ", state[i]);
        if ((i+1)%4==0)
            printf("\n");
    }
}

__host__ __device__ void chacha_odd_half_round(uint32_t x[16], int half)
{
    if (half == 1)
    {
        HALF_1_QUARTERROUND(x[0], x[4], x[8], x[12])
        HALF_1_QUARTERROUND(x[1], x[5], x[9], x[13])
        HALF_1_QUARTERROUND(x[2], x[6], x[10], x[14])
        HALF_1_QUARTERROUND(x[3], x[7], x[11], x[15])
    }
    else
    {
        HALF_2_QUARTERROUND(x[0], x[4], x[8], x[12])
        HALF_2_QUARTERROUND(x[1], x[5], x[9], x[13])
        HALF_2_QUARTERROUND(x[2], x[6], x[10], x[14])
        HALF_2_QUARTERROUND(x[3], x[7], x[11], x[15])
    }
}

__host__ __device__ void chacha_even_round(uint32_t x[16])
{
    QUARTERROUND(x[0], x[5], x[10], x[15])
    QUARTERROUND(x[1], x[6], x[11], x[12])
    QUARTERROUND(x[2], x[7], x[8], x[13])
    QUARTERROUND(x[3], x[4], x[9], x[14])
}

__host__ __device__ void chacha_even_half_round(uint32_t x[16], int half)
{
    if (half == 1)
    {
        HALF_1_QUARTERROUND(x[0], x[5], x[10], x[15])
        HALF_1_QUARTERROUND(x[1], x[6], x[11], x[12])
        HALF_1_QUARTERROUND(x[2], x[7], x[8], x[13])
        HALF_1_QUARTERROUND(x[3], x[4], x[9], x[14])
    }
    else
    {
        HALF_2_QUARTERROUND(x[0], x[5], x[10], x[15])
        HALF_2_QUARTERROUND(x[1], x[6], x[11], x[12])
        HALF_2_QUARTERROUND(x[2], x[7], x[8], x[13])
        HALF_2_QUARTERROUND(x[3], x[4], x[9], x[14])
    }
}

__host__ __device__ void chacha_invert_odd_round(uint32_t x[16])
{
    INVERT_QUARTERROUND(x[3], x[7], x[11], x[15])
    INVERT_QUARTERROUND(x[2], x[6], x[10], x[14])
    INVERT_QUARTERROUND(x[1], x[5], x[9], x[13])
    INVERT_QUARTERROUND(x[0], x[4], x[8], x[12])
}


__host__ __device__ void chacha_invert_even_half_half_round(uint32_t x[16])
{
    INVERT_HALF_HALF_QUARTERROUND(x[3], x[4], x[9], x[14])
    INVERT_HALF_HALF_QUARTERROUND(x[2], x[7], x[8], x[13])
    INVERT_HALF_HALF_QUARTERROUND(x[1], x[6], x[11], x[12])
    INVERT_HALF_HALF_QUARTERROUND(x[0], x[5], x[10], x[15])
}

__host__ __device__ void chacha_invert_even_round(uint32_t x[16])
{
    INVERT_QUARTERROUND(x[3], x[4], x[9], x[14])
    INVERT_QUARTERROUND(x[2], x[7], x[8], x[13])
    INVERT_QUARTERROUND(x[1], x[6], x[11], x[12])
    INVERT_QUARTERROUND(x[0], x[5], x[10], x[15])
}

__host__ __device__ void chacha_rounds(uint32_t state[16], uint32_t rounds, uint32_t lastRound)
{
    uint32_t i;

    for (i = 1; i <= rounds; i++) {
        if ((i + lastRound) % 2)
            chacha_odd_round(state);
        else
        if (i==8 && five_key_recovery){
            HALF_HALF_1_QUARTERROUND(state[0], state[5], state[10], state[15])
            HALF_HALF_1_QUARTERROUND(state[1], state[6], state[11], state[12])
            HALF_HALF_1_QUARTERROUND(state[2], state[7], state[8],  state[13])
            HALF_HALF_1_QUARTERROUND(state[3], state[4], state[9],  state[14])
        } else {
            chacha_even_round(state);
        }
    }
}

__host__ __device__ void chacha_half_rounds(uint32_t state[16], uint32_t rounds, uint32_t lastRound, int half)
{
    uint32_t i;

    for (i = 1; i <= rounds; i++) {
        if ((i + lastRound) % 2)
            chacha_odd_half_round(state, half);
        else
            chacha_even_half_round(state, half);
    }
}

__host__ __device__ void chacha_invert_rounds(uint32_t state[16], uint32_t rounds, uint32_t lastRound)
{
    uint32_t i;

    lastRound = lastRound % 2;

    if (lastRound)
    {
        for (i = 1; i <= rounds; i++) {
            if (i % 2)
                chacha_invert_odd_round(state);
            else
                chacha_invert_even_round(state);
        }
    }
    else
    {
        for (i = 1; i <= rounds; i++) {
            if (i % 2)
                if (i==1 && five_key_recovery){
                    chacha_invert_even_half_half_round(state);
                } else {
                    chacha_invert_even_round(state);
                }
            else
                chacha_invert_odd_round(state);
        }
    }
}

__host__ __device__ void chacha_encrypt(uint32_t output[16], uint32_t input[16], uint32_t rounds)
{
    uint32_t x[16];
    uint32_t i;

    for (i = 0; i < 16; ++i) x[i] = input[i];
    chacha_rounds(x, rounds, 0);
    for (i = 0; i < 16; ++i) x[i] = PLUS(x[i], input[i]);

    memcpy(output, x, 64);
}

__host__ __device__ void chacha_invert(uint32_t output[16], uint32_t input[16], uint32_t intermediate[16], uint32_t rounds, uint32_t lastRound)
{
    for (int i = 0; i < 16; ++i) intermediate[i] = MINUS(output[i], input[i]);
    chacha_invert_rounds(intermediate, rounds, lastRound);
}

#define ALG_TYPE_SALSA 0
#define ALG_TYPE_CHACHA 1

typedef struct {
    uint32_t algType;
    uint32_t key_positions[8];
    uint32_t iv_positions[4];
    void(*init)(uint32_t *, uint32_t *, uint32_t *, uint32_t *);
    void(*encrypt)(uint32_t *, uint32_t *, uint32_t);
    void(*rounds)(uint32_t *, uint32_t, uint32_t);
    void(*halfrounds)(uint32_t *, uint32_t, uint32_t, int);
    void(*invert)(uint32_t *, uint32_t *, uint32_t *, uint32_t, uint32_t);
    char name[20];
} ALGORITHM;

__host__ __device__ void define_alg(ALGORITHM *alg, uint32_t type)
{
    uint32_t chacha_iv_positions[4] = { 12,13,14,15 };
    uint32_t chacha_key_positions[8] = { 4,5,6,7,8,9,10,11 };

    switch (type)
    {
        case ALG_TYPE_CHACHA:
            memcpy(alg->key_positions, chacha_key_positions, 8 * sizeof(uint32_t));
            memcpy(alg->iv_positions, chacha_iv_positions, 4 * sizeof(uint32_t));
            alg->algType = ALG_TYPE_CHACHA;
            alg->init = &chacha_init;
            alg->encrypt = &chacha_encrypt;
            alg->invert = &chacha_invert;
            alg->rounds = &chacha_rounds;
            alg->halfrounds = &chacha_half_rounds;
            alg->name[0] = 'C'; alg->name[1] = 'h'; alg->name[2] = 'a'; alg->name[3] = 'c'; alg->name[4] = 'h'; alg->name[5] = 'a'; alg->name[6] = 0;
            break;

        default:
            break;
    }
}

__host__ __device__ void xor_array(uint32_t *z, uint32_t *x, uint32_t *y, int size)
{
    for (int i = 0; i < size; i++)
        z[i] = x[i] ^ y[i];
}

__host__ __device__ void sub_array(uint32_t *z, uint32_t *x, uint32_t *y, int size)
{
    for (int i = 0; i < size; i++)
        z[i] = x[i] - y[i];
}

__host__ __device__ uint8_t get_bit_in_position(uint32_t state[16], uint32_t pos)
{
    int w = pos / 32;
    int bit = pos % 32;

    return((state[w] >> bit) & 1);
}

__host__ __device__ uint8_t get_bit_from_word_and_bit(uint32_t state[16], uint32_t w, uint32_t bit)
{
    return((state[w] >> bit) & 1);
}

__host__ __device__ void set_bit(uint32_t state[16], uint32_t w, uint32_t bit)
{
    state[w] ^= (1 << bit);
}

__host__ __device__ void set_list_of_bits(uint32_t state[16], uint32_t *w, uint32_t *bit, uint32_t numberOfBits)
{
    for (uint32_t i = 0; i < numberOfBits; i++)
        set_bit(state, w[i], bit[i]);
}

__host__ __device__ void and_array(uint32_t *z, uint32_t *x, uint32_t *y, int size)
{
    for (int i = 0; i < size; i++)
        z[i] = x[i] & y[i];
}

__host__ __device__ uint8_t xor_bits_of_state(uint32_t state[16])
{
    uint32_t x = state[0];
    for (int i = 1; i < 16; i++)
        x ^= state[i];

    x = x ^ (x >> 16);
    x = x ^ (x >> 8);
    x = x ^ (x >> 4);
    x = x ^ (x >> 2);
    return ((x ^ (x >> 1)) & 1);
}

__host__ __device__ uint8_t check_parity_of_equation(uint32_t state[16], uint32_t ODmask[16])
{
    uint32_t aux[16];

    and_array(aux, state, ODmask, 16);
    return(xor_bits_of_state(aux));
}

__device__ uint8_t check_parity_of_linear_relation_cuda(uint32_t inputMask[16], uint32_t inputState[16], uint32_t outputMask[16], uint32_t outputState[16])
{
    uint32_t aux[16], aux2[16];

    and_array(aux, inputState, inputMask, 16);
    and_array(aux2, outputState, outputMask, 16);
    xor_array(aux, aux, aux2, 16);

    return(xor_bits_of_state(aux));
}


__device__ uint8_t chacha_test1_out[64] =
        {
                0x76,0xb8,0xe0,0xad,0xa0,0xf1,0x3d,0x90,0x40,0x5d,0x6a,0xe5,0x53,0x86,0xbd,0x28,0xbd,0xd2,0x19,0xb8,0xa0,0x8d,0xed,0x1a,0xa8,0x36,0xef,0xcc,0x8b,0x77,0x0d,0xc7,0xda,0x41,0x59,0x7c,0x51,0x57,0x48,0x8d,0x77,0x24,0xe0,0x3f,0xb8,0xd8,0x4a,0x37,0x6a,0x43,0xb8,0xf4,0x15,0x18,0xa1,0x1c,0xc3,0x87,0xb6,0x69,0xb2,0xee,0x65,0x86
        };

//---------------------------------------------------------------------------------------
//----------------------  Kernels
//---------------------------------------------------------------------------------------
__device__ int cudaCmp(uint8_t *v1, uint8_t *v2, int len)
{
    for (int i = 0; i < len; i++)
        if (v1[i] != v2[i])
            return 1;

    return 0;
}

/*
This function executes a test vector for chacha to check gpu implementation
*/
__global__ void testChachaKernel(int *rv)
{
    uint32_t k[8] = { 0 }, ctr[2] = { 0 }, nonce[2] = { 0 }, state[16] = { 0 }, state_final[16], inverted[16];
    ALGORITHM alg;
    int tx = threadIdx.x;

    //printf("Teste chacha %d\n", tx);
    define_alg(&alg, ALG_TYPE_CHACHA);

    alg.init(state, k, nonce, ctr);
    alg.encrypt(state_final, state, 20);

    if (cudaCmp((uint8_t *)state_final, chacha_test1_out, 64))
        rv[tx] = 1;
    else
        rv[tx] = 0;

    alg.invert(state_final, state, inverted, 20, 20);
    if (cudaCmp((uint8_t *)inverted, (uint8_t *)state, 64))
        rv[tx] = 2;
}

/*
Computes the differential bias given the number of rounds, ID and OD mask.
*/
__global__ void differential_bias_kernel(unsigned long long seed, int rounds, uint32_t *ID,
                                         uint32_t *ODmask, int ntestForEachThread, unsigned long long int *d_result, int algType, int addExtraHalfRound, int startFromSecondRound)
{
    ALGORITHM alg;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t K[8] = { 0 }, state[16] = { 0 }, alt_state[16] = { 0 };
    uint32_t nonce[2] = { 0 }, ctr[2] = { 0 };
    curandState_t rng;
    unsigned long long int sumParity = 0;

    //printf("Ola sou a td %d\n", tid);
    define_alg(&alg, algType);
    curand_init(seed, tid, 0, &rng);

    for (int t = 0; t < ntestForEachThread; t++)
    {
        if (startFromSecondRound)
        {
            for (int i = 0; i < 16; i++)
                state[i] = curand(&rng);
        }
        else
        {
            for (int i = 0; i < 8; i++)
                K[i] = curand(&rng);

            nonce[0] = curand(&rng); nonce[1] = curand(&rng);
            ctr[0] = curand(&rng); ctr[1] = curand(&rng);

            alg.init(state, K, nonce, ctr);
        }
        xor_array(alt_state, state, ID, 16);

        alg.rounds(state, rounds, startFromSecondRound);
        alg.rounds(alt_state, rounds, startFromSecondRound);
        if (addExtraHalfRound)
        {
            alg.halfrounds(state, 1, rounds+1, 1);
            alg.halfrounds(alt_state, 1, rounds+1, 1);
        }

        xor_array(state, state, alt_state, 16);
        sumParity += check_parity_of_equation(state, ODmask);
    }

    atomicAdd(d_result, sumParity);
}

/*
Computes the linear bias given the number of rounds, ID and OD mask.
*/
__global__ void linear_bias_kernel(unsigned long long seed, int outputRound, int inputRound, uint32_t *IDmask, uint32_t *ODmask, int ntestForEachThread, int *d_result, int algType, int special_linear_relation)
{
    ALGORITHM alg;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t inputState[16] = { 0 }, outputState[16] = { 0 };
    curandState_t rng;
    uint32_t sumParity = 0;

    //printf("Ola sou a td %d\n", tid);
    define_alg(&alg, algType);
    curand_init(seed, tid, 0, &rng);

    for (int t = 0; t < ntestForEachThread; t++)
    {
        for (int i = 0; i < 16; i++)
            inputState[i] = curand(&rng);

        for (int i = 0; i < 16; i++)
            outputState[i] = inputState[i];
        // This "if (special_linear_relation == 5)" was added since the original linear_bias_kernel does not consider half_rounds
        if (special_linear_relation == 5) {
            //chacha_half_rounds(outputState, 3, 1, 2);
            // (uint32_t state[16], uint32_t rounds, uint32_t lastRound, int half)
            alg.halfrounds(outputState, 1, 2, 2);
            // (uint32_t state[16], uint32_t rounds, uint32_t lastRound)
            alg.rounds(outputState, 1, 3);
        } else {
            alg.rounds(outputState, outputRound - inputRound, inputRound);
        }

        sumParity += check_parity_of_linear_relation_cuda(IDmask, inputState, ODmask, outputState);
    }

    atomicAdd(d_result, (int)sumParity);
}

/*This function computes \varepsilon_a from a PNB attack as presented in aumasson 2008*/
__global__ void compute_bias_of_g_for_random_key_kernel(
        unsigned long long seed, uint32_t enc_rounds, uint32_t dec_rounds,
        uint32_t *IDmask, uint32_t *ODmask,
        uint32_t *pnb, uint32_t number_of_pnb, int ntestForEachThread,
        int *d_result, int algType
)
{
    ALGORITHM alg;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t K_with_zeros[8] = { 0 }, state[16] = { 0 }, alt_state[16] = { 0 };
    uint32_t final_state[16] = { 0 }, alt_final_state[16] = { 0 }, aux[16];
    uint32_t intermediate_state[16] = { 0 }, alt_intermediate_state[16] = { 0 };
    uint32_t nonce[2] = { 0 }, ctr[2] = { 0 };
    curandState_t rng;
    uint32_t f_parity, g_parity;
    uint32_t sumParity = 0;
    uint32_t mask;

    uint32_t Krand[8];

    //printf("Ola sou a td %d\n", tid);
    define_alg(&alg, algType);
    curand_init(seed, tid, 0, &rng);

    for (int i = 0; i < 8; i++)
        Krand[i] = curand(&rng);

    for (int i = 0; i < 8; i++)
        K_with_zeros[i] = Krand[i];

    for (uint32_t j = 0; j < number_of_pnb; j++)
    {
        mask = ~(1 << (pnb[j] % 32));
        K_with_zeros[pnb[j] / 32] = K_with_zeros[pnb[j] / 32] & mask;
    }

    for (int t = 0; t < ntestForEachThread; t++)
    {
        nonce[0] = curand(&rng); nonce[1] = curand(&rng);
        ctr[0] = curand(&rng); ctr[1] = curand(&rng);

        //compute for f
        alg.init(state, Krand, nonce, ctr);
        xor_array(alt_state, state, IDmask, 16);

        alg.encrypt(final_state, state, enc_rounds);
        alg.encrypt(alt_final_state, alt_state, enc_rounds);

        alg.invert(final_state, state, intermediate_state, dec_rounds, enc_rounds);
        alg.invert(alt_final_state, alt_state, alt_intermediate_state, dec_rounds, enc_rounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, 16);
        f_parity = check_parity_of_equation(aux, ODmask);

        //compute for g
        alg.init(state, K_with_zeros, nonce, ctr);
        xor_array(alt_state, state, IDmask, 16);

        //use the same final state
        alg.invert(final_state, state, intermediate_state, dec_rounds, enc_rounds);
        alg.invert(alt_final_state, alt_state, alt_intermediate_state, dec_rounds, enc_rounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, 16);
        g_parity = check_parity_of_equation(aux, ODmask);

        if (f_parity == g_parity)
            sumParity++;
    }

    atomicAdd(d_result, (int)sumParity);
}

int testChachaOnGPU()
{
    int *d_rvs;
    const int numThreads = 1024;
    int results[numThreads];

    cudaMalloc((void **)&d_rvs, sizeof(int) * numThreads);
    testChachaKernel <<< 1, numThreads >>> (d_rvs);
    cudaDeviceSynchronize(); //make it flush

    cudaMemcpy(results, d_rvs, sizeof(int)*numThreads, cudaMemcpyDeviceToHost);

    for (int i = 0; i < numThreads; i++)
        if (results[i])
            return 1;

    cudaFree(d_rvs);

    return 0;
}


/*
Compute the differential bias for chacha given the number of rounds, the input differential and the output differential.

addExtraHalfRound - use 1 to add half round
startFromSecondRound - use 1 to start the computation from the second round. This is used when applying a technique from bierle et. al. crypto 2020
*/
double compute_differential_bias(
        int rounds,
        int algType,
        uint32_t ID[16],
        uint32_t ODmask[16],
        uint64_t N, //a power of 2
        int addExtraHalfRound,
        int startFromSecondRound
)
{
    int nTestsForEachThread = NUMBER_OF_TEST_PER_THREAD, nThreads = NUMBER_OF_THREADS, nBlocks = NUMBER_OF_BLOCKS;
    int executionsPerKernel = nTestsForEachThread * nThreads*nBlocks;
    uint64_t iterations;
    unsigned long long int *dSumParity;
    uint32_t *dID, *dODmask;
    unsigned long long int localSumParity = 0;
    double prob = 0;

    uint64_t seed = rand();
    random_uint64(&seed);

    iterations = N / (executionsPerKernel);

    cudaMalloc(&dSumParity, sizeof(unsigned long long int));
    cudaMalloc(&dID, 16 * sizeof(uint32_t));
    cudaMalloc(&dODmask, 16 * sizeof(uint32_t));

    cudaMemcpy(dID, ID, 16 * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(dODmask, ODmask, 16 * sizeof(uint32_t), cudaMemcpyHostToDevice);

    for (int i = 0; i < iterations; i++)
    {
        random_uint64(&seed);
        localSumParity = 0;
        cudaMemcpy(dSumParity, &localSumParity, sizeof(unsigned long long int), cudaMemcpyHostToDevice);

        differential_bias_kernel <<< nBlocks, nThreads >>> ((unsigned long long)seed, rounds, dID, dODmask, nTestsForEachThread, dSumParity, algType, addExtraHalfRound, startFromSecondRound);
        cudaDeviceSynchronize(); //make it flush

        cudaMemcpy(&localSumParity, dSumParity, sizeof(unsigned long long int), cudaMemcpyDeviceToHost);

        prob += ((double)(executionsPerKernel - localSumParity)) / executionsPerKernel;
    }
    prob /= iterations;

    cudaFree(dSumParity);
    cudaFree(dID);
    cudaFree(dODmask);

    return(2 * prob - 1);
}

/*
Compute the differential bias for chacha given the number of rounds, the input differential and the output differential.

addExtraHalfRound - use 1 to add half round
startFromSecondRound - use 1 to start the computation from the second round. This is used when applying a technique from bierle et. al. crypto 2020
*/
double compute_differential_bias_multiple_devices(
        int rounds,
        int algType,
        uint32_t ID[16],
        uint32_t ODmask[16],
        uint64_t N, //a power of 2
        int addExtraHalfRound,
        int startFromSecondRound
)
{
    //uint64_t nTestsForEachThread = (1 << 18), nThreads = (1 << 9), nBlocks = (1 << 8);
    uint64_t nTestsForEachThread = NUMBER_OF_TEST_PER_THREAD, nThreads = NUMBER_OF_THREADS, nBlocks = NUMBER_OF_BLOCKS;
    uint64_t executionsPerKernel = nTestsForEachThread * nThreads*nBlocks;
    uint64_t iterations;
    unsigned long long int  *dSumParityVec[8];
    uint32_t *dIDVec[8], *dODmaskVec[8];
    unsigned long long int localSumParityVec[8] = { 0 };
    double prob = 0;

    uint64_t seed = rand();
    random_uint64(&seed);

    iterations = N / (executionsPerKernel * NUMBER_OF_DEVICES);
    //printf("size of long long int = %d\n",sizeof(unsigned long long int));

    for (int d = 0; d < NUMBER_OF_DEVICES; d++)
    {
        cudaSetDevice(d);
        cudaMalloc(&(dSumParityVec[d]), sizeof(unsigned long long int));
        cudaMalloc(&(dIDVec[d]), 16 * sizeof(uint32_t));
        cudaMalloc(&(dODmaskVec[d]), 16 * sizeof(uint32_t));

        cudaMemcpy(dIDVec[d], ID, 16 * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(dODmaskVec[d], ODmask, 16 * sizeof(uint32_t), cudaMemcpyHostToDevice);
    }

    //printf("iterations %d = \n", iterations);
    for (int i = 0; i < iterations; i++)
    {
        for (int d = 0; d < NUMBER_OF_DEVICES; d++)
        {
            cudaSetDevice(d);
            random_uint64(&seed);
            localSumParityVec[d] = 0;
            cudaMemcpy(dSumParityVec[d], &(localSumParityVec[d]), sizeof(unsigned long long int), cudaMemcpyHostToDevice);

            differential_bias_kernel <<< nBlocks, nThreads >>> ((unsigned long long)seed, rounds,
                                                                dIDVec[d], dODmaskVec[d], nTestsForEachThread, dSumParityVec[d], algType,
                                                                addExtraHalfRound, startFromSecondRound);
        }
        cudaDeviceSynchronize(); //make it flush

        for (int d = 0; d < NUMBER_OF_DEVICES; d++)
        {
            cudaSetDevice(d);
            cudaMemcpy(&(localSumParityVec[d]), dSumParityVec[d],
                       sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
            //printf("%" PRIu64 "\n", executionsPerKernel);
            prob += ((double)(executionsPerKernel - localSumParityVec[d])) / executionsPerKernel;
        }
    }
    prob /= (iterations*NUMBER_OF_DEVICES);

    for (int d = 0; d < NUMBER_OF_DEVICES; d++)
    {
        cudaFree(dSumParityVec[d]);
        cudaFree(dIDVec[d]);
        cudaFree(dODmaskVec[d]);
    }

    return(2 * prob - 1);
}

//compute linear correlation using cuda
double compute_linear_bias_cuda(
        int inputRound,
        int outputRound,
        int algType,
        uint32_t IDmask[16],
        uint32_t ODmask[16],
        uint64_t N,
        int special_linear_approx
)
{
    //int nTestsForEachThread = (1 << 6), nThreads = (1 << 6), nBlocks = (1 << 6);
    int nTestsForEachThread = NUMBER_OF_TEST_PER_THREAD, nThreads = NUMBER_OF_THREADS, nBlocks = NUMBER_OF_BLOCKS;
    int executionsPerKernel = nTestsForEachThread * nThreads*nBlocks;
    uint64_t iterations;
    int *dSumParity;
    uint32_t *dID, *dODmask;
    int localSumParity = 0;
    double prob = 0;

    uint64_t seed = rand();
    random_uint64(&seed);

    iterations = N / (executionsPerKernel);

    cudaMalloc(&dSumParity, sizeof(int));
    cudaMalloc(&dID, 16 * sizeof(uint32_t));
    cudaMalloc(&dODmask, 16 * sizeof(uint32_t));

    cudaMemcpy(dID, IDmask, 16 * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(dODmask, ODmask, 16 * sizeof(uint32_t), cudaMemcpyHostToDevice);

    for (int i = 0; i < iterations; i++)
    {
        random_uint64(&seed);
        localSumParity = 0;
        cudaMemcpy(dSumParity, &localSumParity, sizeof(int), cudaMemcpyHostToDevice);

        linear_bias_kernel <<< nBlocks, nThreads >>> ((unsigned long long)seed, outputRound, inputRound, dID, dODmask, nTestsForEachThread, dSumParity, algType, special_linear_approx);

        //cudaDeviceSynchronize(); //make it flush

        cudaMemcpy(&localSumParity, dSumParity, sizeof(uint32_t), cudaMemcpyDeviceToHost);

        prob += ((double)(executionsPerKernel - localSumParity)) / executionsPerKernel;
    }
    prob /= iterations;

    cudaFree(dSumParity);
    cudaFree(dID);
    cudaFree(dODmask);

    return(2 * prob - 1);
}

/*
compute \varepsion_a in a PNB attack as in aumasson 2008
*/
double compute_mean_bias_of_g_cuda(
        uint64_t N,
        uint32_t ID[16],
        uint32_t ODmask[16],
        uint32_t enc_rounds,
        uint32_t dec_rounds,
        uint32_t *pnb,
        uint32_t number_of_pnb,
        ALGORITHM alg
)
{
    //int nTestsForEachThread = (1 << 7), nThreads = (1 << 8), nBlocks = (1 << 8);
    int nTestsForEachThread = NUMBER_OF_TEST_PER_THREAD, nThreads = NUMBER_OF_THREADS, nBlocks = NUMBER_OF_BLOCKS;
    int executionsPerKernel = nTestsForEachThread * nThreads*nBlocks;
    uint64_t iterations;
    int *dSumParity;
    uint32_t *dID, *dODmask, *dPNB;
    int localSumParity = 0;
    double prob = 0;

    uint64_t seed = rand();
    random_uint64(&seed);

    iterations = N / (executionsPerKernel);

    cudaMalloc(&dSumParity, sizeof(int));
    cudaMalloc(&dID, 16 * sizeof(uint32_t));
    cudaMalloc(&dODmask, 16 * sizeof(uint32_t));
    cudaMalloc(&dPNB, number_of_pnb * sizeof(uint32_t));

    cudaMemcpy(dID, ID, 16 * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(dODmask, ODmask, 16 * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(dPNB, pnb, number_of_pnb * sizeof(uint32_t), cudaMemcpyHostToDevice);

    for (int i = 0; i < iterations; i++)
    {
        random_uint64(&seed);
        localSumParity = 0;
        cudaMemcpy(dSumParity, &localSumParity, sizeof(int), cudaMemcpyHostToDevice);

        compute_bias_of_g_for_random_key_kernel <<< nBlocks, nThreads >>> ((unsigned long long)seed,
                                                                           enc_rounds, dec_rounds, dID, dODmask, dPNB, number_of_pnb, nTestsForEachThread,
                                                                           dSumParity, alg.algType);
        //cudaDeviceSynchronize(); //make it flush

        cudaMemcpy(&localSumParity, dSumParity, sizeof(uint32_t), cudaMemcpyDeviceToHost);

        prob += ((double)(localSumParity)) / executionsPerKernel;
    }
    prob /= iterations;

    cudaFree(dSumParity);
    cudaFree(dID);
    cudaFree(dODmask);

    return(2 * prob - 1);
}

double get_max(double x, double y)
{
    if (x>y)
        return x;
    else
        return y;
}


//compute complexity of pnb attack
void compute_complexity_of_the_attack(double *data_complexity, double *time_complexity, double bias_of_g, int number_of_pnb)
{
    int alpha;
    int m = 256 - number_of_pnb;
    double N, tc, minN = 256, minTC = 256;

    for (alpha = 0; alpha < 256; alpha++)
    {
        N = (sqrt(alpha*log(4)) + 3 * sqrt(1 - bias_of_g * bias_of_g)) / bias_of_g;
        N = N * N;
        tc = get_max(256 - alpha, m + log2(N));

        if (tc < minTC)
        {
            minTC = tc;
            minN = N;
            printf("alpha=%d\n", alpha);
        }
    }

    *data_complexity = log2(minN);
    *time_complexity = minTC;
}


/*
Alg of sec 3.3 of New features of latin dances.
neutrality_measure - is the neutrality measure of this neutral bit.
N - number of tests to be executed
input_differential - the input differential (ID) under analysis
enc_rounds - in the paper is R
dec_rounds - in the paper is r-R
od_word - identify the word of the output differential (0 to 16)
od_bit_of_word - identify the bit of the word for the output differential (0 to 31)
neutral_word - identify the word in which is located the neutral bit being tested (0 to 16)
neutral_bit - identify the bit of the word for the neutral bit being tested (0 to 31)
init, encrypt and invert are pointer to salsa or chacha or any other initialization, encryption, or inversion functions
*/
void compute_neutrality_of_bit(
        double *neutrality_measure,
        uint64_t N, uint32_t ID[16], uint32_t ODmask[16],
        uint32_t enc_rounds, uint32_t dec_rounds,
        uint32_t neutral_word, uint32_t neutral_bit,
        ALGORITHM alg)
{
    uint32_t K[8] = { 0 }, state[16] = { 0 }, alt_state[16] = { 0 };
    uint32_t final_state[16] = { 0 }, alt_final_state[16] = { 0 };
    uint32_t intermediate_state[16] = { 0 }, alt_intermediate_state[16] = { 0 }, aux[16];
    uint32_t nonce[2] = { 0 }, ctr[2] = { 0 };
    uint32_t diff, new_diff, neutral_diff;
    uint64_t count = 0;

    neutral_diff = 1 << neutral_bit;

    for (uint64_t i = 0; i < N; i++)
    {
        random_uint32_array(K, 8);
        random_uint32_array(nonce, 2);
        random_uint32_array(ctr, 2);

        alg.init(state, K, nonce, ctr);
        xor_array(alt_state, state, ID, 16);

        alg.encrypt(final_state, state, enc_rounds);
        alg.encrypt(alt_final_state, alt_state, enc_rounds);

        alg.invert(final_state, state, intermediate_state, dec_rounds, enc_rounds);
        alg.invert(alt_final_state, alt_state, alt_intermediate_state, dec_rounds, enc_rounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, 16);
        diff = check_parity_of_equation(aux, ODmask);

        state[neutral_word] ^= neutral_diff;
        alt_state[neutral_word] ^= neutral_diff;

        alg.invert(final_state, state, intermediate_state, dec_rounds, enc_rounds);
        alg.invert(alt_final_state, alt_state, alt_intermediate_state, dec_rounds, enc_rounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, 16);
        new_diff = check_parity_of_equation(aux, ODmask);

        if (diff == new_diff)
            count++;
    }

    *neutrality_measure = 2 * ((double)count) / N - 1;
}

void compute_neutrality_for_every_key_bit(
        double neutrality_list[256],
        uint64_t number_of_tests_for_each_key_bit,
        uint32_t ID[16], uint32_t ODmask[16],
        uint32_t enc_rounds, uint32_t dec_rounds,
        ALGORITHM alg)
{
    int count = 0;
    double neutrality_measure = 0;

    for (uint32_t neutral_word = 0; neutral_word < 8; neutral_word++)
    {
        for (uint32_t neutral_bit = 0; neutral_bit < 32; neutral_bit++)
        {
            //printf("%d %d\n", neutral_word, neutral_bit);
            compute_neutrality_of_bit(&neutrality_measure, number_of_tests_for_each_key_bit, ID, ODmask,
                                      enc_rounds, dec_rounds, alg.key_positions[neutral_word], neutral_bit, alg);

            neutrality_list[count] = neutrality_measure;
            count++;
        }
    }
}

void get_pnb_list(uint32_t pnb[256], uint32_t *number_of_pnb,
                  double neutrality_measure_threshold, double neutrality_list[256])
{
    *number_of_pnb = 0;
    for (uint32_t i = 0; i < 256; i++)
    {
        if (neutrality_list[i] > neutrality_measure_threshold)
        {
            pnb[*number_of_pnb] = i;
            (*number_of_pnb)++;
        }
    }
}
// key-recovery attack using Differential-Linear Distinguisher 1
void first_key_recovery_attack()
{
    printf("First key-recovery attack\n");
    uint32_t ODmask[16] = { 0 };
    uint32_t ID[16] = { 0 };
    uint64_t N = (1 << 20);
    ALGORITHM alg;
    double neutrality_list[256];
    uint32_t pnb[256], number_of_pnb;
    double threshold = 0.24, varepsilon_a, varepsilon_d; //24

    define_alg(&alg, ALG_TYPE_CHACHA);

    //Go 4, back 2
    memset(ODmask, 0x00, sizeof(uint32_t) * 16);
    N = 1;
    N <<= 17;
    ID[14] = (1 << 21);
    ID[14] = (1 << 9);
    uint32_t finalListOfWords[2] = { 11, 12 };
    uint32_t finalListOfBits[2] = { 0, 0 };
    set_list_of_bits(ODmask, finalListOfWords, finalListOfBits, 2);

    compute_neutrality_for_every_key_bit(neutrality_list, N, ID, ODmask, 7, 3, alg);
    get_pnb_list(pnb, &number_of_pnb, threshold, neutrality_list);
    for (int i = 0; i < number_of_pnb; i++)
        printf("%d, ", pnb[i]);

    printf("\nNumber of pnb: %d.\n", number_of_pnb);

    N = 1;
    N<<=46;
    varepsilon_a = compute_mean_bias_of_g_cuda(N, ID, ODmask,
                                               7, 3, pnb, number_of_pnb, alg);
    double comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);
    printf("comp=%f\n", comp);
    if (fabs(varepsilon_a)>comp) {
        printf("Found\n");
    }

    varepsilon_d = 0.00000227;
    double time_complexity, data_complexity;
    double e = varepsilon_a * varepsilon_d;

    compute_complexity_of_the_attack(&data_complexity, &time_complexity, e, number_of_pnb);
    //printf("limiar= %f, varepsilon_a = %f \\varepsilon = %f, data_complexity %f, time_complexity %f.\n", threshold, varepsilon_a, e, data_complexity, time_complexity);
    printf("comp=%.17g, limiar= %f, varepsilon_a = %.17g \\varepsilon = %f, data_complexity %f, time_complexity %f.\n", comp, threshold, varepsilon_a, e, data_complexity, time_complexity);
}


// key-recovery attack using Differential-Linear Distinguisher 4.
void second_key_recovery_attack()
{
    printf("Second key-recovery attack\n");
    uint32_t ODmask[16] = {
            0x00000001,0x00000000,0x00000000,0x00000000,
            0x00000000,0x00000000,0x00000000,0x00000000,
            0x00000000,0x00000000,0x00000000,0x00000000,
            0x00000000,0x00000000,0x00000000,0x00000000
    };
    uint32_t ID[16] = { 0 };
    uint64_t N = (1 << 20);
    ALGORITHM alg;
    double neutrality_list[256];
    uint32_t pnb[256], number_of_pnb;
    double threshold = 0.27, varepsilon_a, varepsilon_d;

    define_alg(&alg, ALG_TYPE_CHACHA);

    //Go 4, back 3

    N = 1;
    N <<= 14;
    ID[14] = (1 << 6);

    compute_neutrality_for_every_key_bit(neutrality_list, N, ID, ODmask, 7, 3, alg);
    get_pnb_list(pnb, &number_of_pnb, threshold, neutrality_list);
    for (int i = 0; i < number_of_pnb; i++)
        printf("%d, ", pnb[i]);

    printf("\nNumber of pnb: %d.\n", number_of_pnb);

    N = 1;
    N<<=42;
    varepsilon_a = compute_mean_bias_of_g_cuda(N, ID, ODmask,
                                               7, 3, pnb, number_of_pnb, alg);

    varepsilon_d = pow(2,-39.88);
    double time_complexity, data_complexity;
    double e = varepsilon_a * varepsilon_d;
    double comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

    compute_complexity_of_the_attack(&data_complexity, &time_complexity, e, number_of_pnb);
    printf("comp=%.17g, threshold = %f, varepsilon_a = %.17g \\varepsilon = %f, data_complexity %f, time_complexity %f.\n", comp, threshold, varepsilon_a, e, data_complexity, time_complexity);
}

typedef struct {
    uint32_t inputRound;
    uint32_t inputMask[16];
    uint32_t outputRound;
    uint32_t outputMask[16];
    uint64_t parityCount;
    uint64_t totalCount;
    double bias;
} LINEAR_RELATION;

void print_latex_linear_relation(LINEAR_RELATION *L)
{
    double prob;

    printf("$$ \\begin{array}{cl}\n");
    if (fabs(L->bias) > 0)
    {
        printf("%f = (1+%f)/2 = &\\\\ \\Pr(", (1 + L->bias) / 2, L->bias);
    }
    else
    {
        prob = ((double)L->totalCount - L->parityCount) / L->totalCount;
        printf("%f = (1+%f)/2 = \\Pr(", prob, 2 * prob - 1);
    }
    for (int w = 0; w < 16; w++)
    {
        for (int b = 0; b < 32; b++)
        {
            if (get_bit_from_word_and_bit(L->inputMask, w, b))
                printf("x^{(%d)}_{%d,%d} \\oplus ", L->inputRound, w, b);
        }
    }
    printf(" = & ");

    int count = 0;
    for (int w = 0; w < 16; w++)
    {
        for (int b = 0; b < 32; b++)
        {
            if (get_bit_from_word_and_bit(L->outputMask, w, b))
            {
                printf("x^{(%d)}_{%d,%d} \\oplus ", L->outputRound, w, b);
                count++;
                if (count == 8)
                {
                    count = 0;
                    printf("\\\\ &");
                }
            }
        }
    }
    printf(") \\end{array} $$ \n\n");
}

uint8_t check_parity_of_linear_relation(LINEAR_RELATION *L, uint32_t inputState[16], uint32_t outputState[16])
{
    uint32_t aux[16], aux2[16];

    and_array(aux, inputState, L->inputMask, 16);
    and_array(aux2, outputState, L->outputMask, 16);
    xor_array(aux, aux, aux2, 16);

    return(xor_bits_of_state(aux));
}

int test_linear_relation(LINEAR_RELATION *L, uint32_t N, ALGORITHM alg)
{
    uint32_t state[16], finalState[16];
    L->totalCount = N;
    for (int i = 0; i < N; i++)
    {
        random_uint32_array(state, 16);
        memcpy(finalState, state, 16 * sizeof(uint32_t));

        alg.rounds(finalState, L->outputRound - L->inputRound, L->inputRound);

        L->parityCount += check_parity_of_linear_relation(L, state, finalState);
    }

    double prob;
    prob = ((double)L->totalCount - L->parityCount) / L->totalCount;

    if (fabs(prob - 0.5) > 4 * sqrt(0.25 / N))
    {
        print_latex_linear_relation(L);
        return 1;
    }
    return 0;
}

int test_linear_relation_cuda(LINEAR_RELATION *L, uint64_t N, ALGORITHM alg, int special_linear_relation)
{
    L->bias = compute_linear_bias_cuda(L->inputRound, L->outputRound,
                                       alg.algType, L->inputMask, L->outputMask, N, special_linear_relation);

    double comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

    printf("%.17g\n", fabs(L->bias));
    if (fabs(L->bias) > comp)
    {
        print_latex_linear_relation(L);
        return 1;
    }
    return 0;
}

/*
Parameters:
w = 0 for set (0,4,8,12)
w = 1 for set (1,5,9,13)
w = 2 for set (2,6,10,14)
w = 3 for set (3,7,11,15)

letter 0,1,2,3 = A,B,C,D respectivelly
*/
#define LetterA 0
#define LetterB 1
#define LetterC 2
#define LetterD 3
#define typeLinear 0
#define typeNonLinear 1
#define typeSpecialC 2
#define typeSpecialA 3
void add_elements_to_linear_equation(LINEAR_RELATION *L, int w, int bit, int letter, unsigned char type, int round)
{
    int size;
    uint32_t abcd[4] = { 0,4,8,12 };

    uint32_t shift[4] = { 0 };
    if ((round % 2) == 0)
    {
        shift[1] = 1; shift[2] = 2; shift[3] = 3;
    }

    for (int i = 0; i < 4; i++)
        abcd[i] = abcd[i] + (shift[i] + w) % 4;

    uint32_t listOfWords[9] = { 0 }, listOfBits[9] = { 0 };
    uint32_t listOfWordsA[9] = { abcd[0], abcd[1],abcd[1], abcd[2], abcd[3], abcd[1], abcd[2], abcd[3], abcd[1] };
    uint32_t listOfBitsA[9] = { 0, 7,19,12,0,  18,11,31, 6 };

    uint32_t listOfWordsB[9] = { abcd[1],abcd[2],abcd[3],abcd[2],abcd[3] };
    uint32_t listOfBitsB[9] = { 19,12,0,0,31 };

    uint32_t listOfWordsC[9] = { abcd[3],abcd[2],abcd[3],abcd[0],abcd[0],abcd[3],abcd[3] };
    uint32_t listOfBitsC[9] = { 0,0,8,0,31,7,31 };

    uint32_t listOfWordsD[9] = { abcd[3],abcd[0],abcd[0],abcd[2],abcd[1],abcd[2],abcd[1] };
    uint32_t listOfBitsD[9] = { 24,16,0,0,7,31,6 };

    for (int j = 0; j < 9; j++)
    {
        switch (letter)
        {
            case 0:
                listOfWords[j] = listOfWordsA[j];
                listOfBits[j] = listOfBitsA[j];
                break;
            case 1:
                listOfWords[j] = listOfWordsB[j];
                listOfBits[j] = listOfBitsB[j];
                break;
            case 2:
                listOfWords[j] = listOfWordsC[j];
                listOfBits[j] = listOfBitsC[j];
                break;
            case 3:
                listOfWords[j] = listOfWordsD[j];
                listOfBits[j] = listOfBitsD[j];
                break;
        }
    }

    switch (letter)
    {
        case 0:
            if (type == 0)
                size = 5;
            else if (type == 3)
            {
                size = 9;
                listOfWords[8] = abcd[2];
                listOfBits[8] = (uint32_t)-1;
            }
            else
                size = 9;
            break;
        case 1:
            if (type == 0)
                size = 4;
            else
                size = 5;
            break;
        case 2:
            if (type == 0)
                size = 4;
            else if (type == 1)
                size = 7;
            else
                size = 6;
            break;
        case 3:
            if (type == 0)
                size = 5;
            else
                size = 7;
            break;
    }

    set_bit(L->inputMask, letter * 4 + (w + shift[letter]) % 4, bit);
    for (int j = 0; j < size; j++)
    {
        listOfWords[j] = listOfWords[j];
        listOfBits[j] = (listOfBits[j] + 32 + bit) % 32;
    }
    set_list_of_bits(L->outputMask, listOfWords, listOfBits, size);
}



void lemma_1()
{
    LINEAR_RELATION L;
    ALGORITHM alg;
    uint64_t N = 1;

#ifdef REDUCE_NUMBER_OF_ITERATIONS
    N <<= 33; //use for faster test
#else
    N <<= 33; //used for paper result
#endif

    printf("\n::Lemma 1:\n");
    // This linear trail starts at round 3.5 but remember that the transition to round 4 occurs with probability 1. So, we are starting from round 4.

    define_alg(&alg, ALG_TYPE_CHACHA);
    memset(&L, 0x00, sizeof(LINEAR_RELATION));

    {
        //------------------------------
        L.inputRound = 4;
        L.outputRound = 6;
        set_bit(L.inputMask, 12, 0);
        set_bit(L.inputMask, 11, 0);
        uint32_t listOfWords[43] ={ 0, 0, 1, 1, 1, 2, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 14, 14, 15, 15, 15};
        uint32_t listOfBits[43] ={ 8, 24, 0, 8, 24, 0, 0, 7, 19, 26, 2, 3, 14, 15, 19, 22, 23, 30, 31, 0, 6, 7, 12, 19, 0, 7, 8, 12, 27, 28, 0, 23, 24, 0, 8, 16, 0, 8, 0, 7, 14, 16, 24};
        set_list_of_bits(L.outputMask, listOfWords, listOfBits, 43);

        test_linear_relation_cuda(&L, N, alg, -1);
        //------------------------------
    }
}

void lemma_3()
{
    LINEAR_RELATION L;
    ALGORITHM alg;
    uint64_t N = 1;

#ifdef REDUCE_NUMBER_OF_ITERATIONS
    N <<= 33; //use for faster test
#else
    N <<= 33; //used for paper result
#endif

    printf("\n::Lemma 3:\n");
    // This linear trail starts at round 3.5 but remember that the transition to round 4 occurs with probability 1. So, we are starting from round 4.

    define_alg(&alg, ALG_TYPE_CHACHA);
    memset(&L, 0x00, sizeof(LINEAR_RELATION));

    {
        //------------------------------
        L.inputRound = 4;
        L.outputRound = 6;
        set_bit(L.inputMask, 2, 0);
        uint32_t listOfWords[25] ={ 0, 0, 2, 3, 3, 4, 6, 6, 7, 7, 8, 9, 10, 11, 11, 11, 11, 12, 12, 13, 14, 15, 15, 15, 15,  };
        uint32_t listOfBits[25] ={ 11, 12, 0, 0, 16, 7, 6, 26, 7, 19, 12, 0, 12, 6, 7, 18, 31, 7, 19, 0, 24, 11, 12, 19, 20,  };
        set_list_of_bits(L.outputMask, listOfWords, listOfBits, 25);

        test_linear_relation_cuda(&L, N, alg, -1);
        //------------------------------
    }
}

void toy_example()
{
    LINEAR_RELATION L;
    ALGORITHM alg;
    uint64_t N = 1;

#ifdef REDUCE_NUMBER_OF_ITERATIONS
    N <<= 49; //use for faster test
#else
    N <<= 49; //used for paper result
#endif

    printf("\n::Toy example:\n");
    // This linear trail starts at round 6

    define_alg(&alg, ALG_TYPE_CHACHA);
    memset(&L, 0x00, sizeof(LINEAR_RELATION));

    {
        //------------------------------
        L.inputRound = 6;
        L.outputRound = 7;
        set_bit(L.inputMask, 0, 11);
        set_bit(L.inputMask, 0, 12);
        set_bit(L.inputMask, 2, 0);
        set_bit(L.inputMask, 3, 0);
        set_bit(L.inputMask, 3, 16);
        set_bit(L.inputMask, 4, 7);
        set_bit(L.inputMask, 6, 6);
        set_bit(L.inputMask, 6, 26);
        set_bit(L.inputMask, 7, 7);
        set_bit(L.inputMask, 7, 19);
        set_bit(L.inputMask, 8, 12);
        set_bit(L.inputMask, 9, 0);
        set_bit(L.inputMask, 10, 12);
        set_bit(L.inputMask, 11, 7);
        set_bit(L.inputMask, 11, 6);
        set_bit(L.inputMask, 11, 31);
        set_bit(L.inputMask, 11, 18);
        set_bit(L.inputMask, 12, 7);
        set_bit(L.inputMask, 12, 19);
        set_bit(L.inputMask, 13, 0);
        set_bit(L.inputMask, 14, 24);
        set_bit(L.inputMask, 15, 11);
        set_bit(L.inputMask, 15, 12);
        set_bit(L.inputMask, 15, 19);
        set_bit(L.inputMask, 15, 20);
        uint32_t listOfWords[99]={0,0,0,0,0,0,0,1,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,6,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,8,8,10,10,10,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,12,12,13,13,13,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};
        uint32_t listOfBits[99]={3,7,10,11,12,19,23,16,0,8,10,11,24,0,3,4,6,7,12,16,17,18,19,20,28,30,31,13,14,17,18,25,29,30,7,7,13,19,25,30,31,2,3,6,7,18,22,23,27,12,18,22,23,5,11,12,18,23,24,25,26,11,19,20,27,28,7,9,10,12,19,20,31,0,8,24,0,6,10,11,16,18,19,26,0,4,6,7,11,12,15,16,17,18,19,25,26,30,31};
        set_list_of_bits(L.outputMask, listOfWords, listOfBits, 99);

        test_linear_relation_cuda(&L, N, alg, -1);
        //------------------------------
    }
}

void toy_example_groupA()
{
    LINEAR_RELATION L;
    ALGORITHM alg;
    uint64_t N = 1;

#ifdef REDUCE_NUMBER_OF_ITERATIONS
    N <<= 34; //use for faster test
#else
    N <<= 34; //used for paper result
#endif

    printf("\n::Toy example group A:\n");
    // This linear trail starts at round 6

    define_alg(&alg, ALG_TYPE_CHACHA);
    memset(&L, 0x00, sizeof(LINEAR_RELATION));

    {
        //------------------------------
        L.inputRound = 6;
        L.outputRound = 7;
        set_bit(L.inputMask, 2, 0);
        set_bit(L.inputMask, 3, 0);
        set_bit(L.inputMask, 9, 0);
        set_bit(L.inputMask, 13, 0);
        uint32_t listOfWords[15]={1,2,3,5,6,6,7,7,10,11,13,13,13,14,15};
        uint32_t listOfBits[15]={16,0,0,7,7,19,7,19,12,12,0,8,24,0,0};

        set_list_of_bits(L.outputMask, listOfWords, listOfBits, 15);

        test_linear_relation_cuda(&L, N, alg, -1);
        //------------------------------
    }
}


void toy_example_groupB()
{
    LINEAR_RELATION L;
    ALGORITHM alg;
    uint64_t N = 1;

#ifdef REDUCE_NUMBER_OF_ITERATIONS
    N <<= 34; //use for faster test
#else
    N <<= 34; //used for paper result
#endif

    printf("\n::Toy example group B:\n");
    // This linear trail starts at round 6

    define_alg(&alg, ALG_TYPE_CHACHA);
    memset(&L, 0x00, sizeof(LINEAR_RELATION));

    {
        //------------------------------
        L.inputRound = 6;
        L.outputRound = 7;
        set_bit(L.inputMask, 11, 7);
        set_bit(L.inputMask, 11, 6);
        set_bit(L.inputMask, 15, 19);
        set_bit(L.inputMask, 15, 20);
        set_bit(L.inputMask, 0, 11);
        set_bit(L.inputMask, 0, 12);
        set_bit(L.inputMask, 15, 11);
        set_bit(L.inputMask, 15, 12);
        uint32_t listOfWords[38]={0,0,3,3,3,3,3,3,3,3,3,3,4,4,4,4,7,7,7,7,8,8,11,11,11,11,11,11,12,12,15,15,15,15,15,15,15,15};
        uint32_t listOfBits[38]={11,12,3,4,6,7,11,12,19,20,27,28,18,19,30,31,18,19,26,27,23,24,6,7,11,12,19,20,11,12,3,4,6,7,11,12,14,15};

        set_list_of_bits(L.outputMask, listOfWords, listOfBits, 38);

        test_linear_relation_cuda(&L, N, alg, -1);
        //------------------------------
    }
}


void toy_example_groupA_and_groupB()
{
    LINEAR_RELATION L;
    ALGORITHM alg;
    uint64_t N = 1;

#ifdef REDUCE_NUMBER_OF_ITERATIONS
    N <<= 34; //use for faster test
#else
    N <<= 34; //used for paper result
#endif

    printf("\n::Toy example group A and group B:\n");
    // This linear trail starts at round 6

    define_alg(&alg, ALG_TYPE_CHACHA);
    memset(&L, 0x00, sizeof(LINEAR_RELATION));

    {
        //------------------------------
        L.inputRound = 6;
        L.outputRound = 7;
        set_bit(L.inputMask, 11, 7);
        set_bit(L.inputMask, 11, 6);
        set_bit(L.inputMask, 15, 19);
        set_bit(L.inputMask, 15, 20);
        set_bit(L.inputMask, 0, 11);
        set_bit(L.inputMask, 0, 12);
        set_bit(L.inputMask, 15, 11);
        set_bit(L.inputMask, 15, 12);
        set_bit(L.inputMask, 2, 0);
        set_bit(L.inputMask, 3, 0);
        set_bit(L.inputMask, 9, 0);
        set_bit(L.inputMask, 13, 0);

        //uint32_t listOfWords[15]={1,2,3,5,6,6,7,7,10,11,13,13,13,14,15};
        //uint32_t listOfBits[15]={16,0,0,7,7,19,7,19,12,12,0,8,24,0,0};
        uint32_t listOfWords[38+15]={0,0,3,3,3,3,3,3,3,3,3,3,4,4,4,4,7,7,7,7,8,8,11,11,11,11,11,11,12,12,15,15,15,15,15,15,15,15, 1,2,3,5,6,6,7,7,10,11,13,13,13,14,15};
        uint32_t listOfBits[38+15]={11,12,3,4,6,7,11,12,19,20,27,28,18,19,30,31,18,19,26,27,23,24,6,7,11,12,19,20,11,12,3,4,6,7,11,12,14,15, 16,0,0,7,7,19,7,19,12,12,0,8,24,0,0};

        set_list_of_bits(L.outputMask, listOfWords, listOfBits, 38+15);

        test_linear_relation_cuda(&L, N, alg, -1);
        //------------------------------
    }
}


void linear_trail_4()
{
    LINEAR_RELATION L;
    ALGORITHM alg;
    uint64_t N = 1;

#ifdef REDUCE_NUMBER_OF_ITERATIONS
    N <<= 48; //use for faster test
#else
    N <<= 48; //used for paper result
#endif

    printf("\nLinear Trail 4:\n");
    // This linear trail starts at round 3.5 but remember that the transition to round 4 occurs with probability 1. So, we are starting from round 4.

    define_alg(&alg, ALG_TYPE_CHACHA);
    memset(&L, 0x00, sizeof(LINEAR_RELATION));

    {
        //------------------------------
        L.inputRound = 4;
        L.outputRound = 7;
        set_bit(L.inputMask, 2, 0);
        uint32_t listOfWords[83] = {0,0,0,0,0,0,0,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,5,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,10,10,10,10,10,10,11,11,11,11,12,12,12,12,12,12,13,13,13,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15};
        uint32_t listOfBits[83] = {23,22,19,11,7,3,2,16,24,23,12,11,8,0,31,30,28,20,19,16,12,11,7,6,4,0,31,30,19,14,7,31,25,19,13,7,27,26,23,7,6,3,24,23,12,6,26,25,24,18,11,5,28,20,18,15,31,30,20,11,10,7,24,8,0,26,20,19,16,12,6,0,31,30,19,16,15,14,12,7,6,4,0};
        set_list_of_bits(L.outputMask, listOfWords, listOfBits, 83);

        test_linear_relation_cuda(&L, N, alg, -1);
        //------------------------------
    }
}

void linear_trail_5()
{
    LINEAR_RELATION L;
    ALGORITHM alg;
    uint64_t N = 1;

#ifdef REDUCE_NUMBER_OF_ITERATIONS
    N <<= 38; //use for faster test
#else
    N <<= 38; //used for paper result
#endif

    printf("\nLinear Trail 5:\n");

    define_alg(&alg, ALG_TYPE_CHACHA);
    memset(&L, 0x00, sizeof(LINEAR_RELATION));

    {
        //------------------------------
        L.inputRound = 2;
        L.outputRound = 4;
        uint32_t input_temp[16] = {
//                    0x00000010,0x00200020,0x00001000,0x10000000,
//                    0x00000010,0x20302020,0x00001000,0x10000000,
//                    0x00000000,0x30002000,0x00100000,0x00000000,
//                    0x00000000,0x00200030,0x00001800,0x18000000
                0x00000011,0x00220020,0x00001800,0x18000000,
                0x00000011,0x22222020,0x00001000,0x10000000,
                0x00000000,0x22002000,0x00100000,0x00000000,
                0x00000000,0x00330030,0x00001000,0x10000000
        };
        for (int i=0; i<16; i++)
            L.inputMask[i] = input_temp[i];
        uint32_t listOfWords[1] = {0};
        uint32_t listOfBits[1] = {0};
        set_list_of_bits(L.outputMask, listOfWords, listOfBits, 1);

        test_linear_relation_cuda(&L, N, alg, 5);
        //------------------------------
    }
}

void linear_trail_6()
{
    LINEAR_RELATION L;
    ALGORITHM alg;
    uint64_t N = 1;

#ifdef REDUCE_NUMBER_OF_ITERATIONS
    N <<= 38; //use for faster test
#else
    N <<= 38; //used for paper result
#endif

    printf("\nLinear Trail 6:\n");

    define_alg(&alg, ALG_TYPE_CHACHA);
    memset(&L, 0x00, sizeof(LINEAR_RELATION));

    {
        //------------------------------
        L.inputRound = 2;
        L.outputRound = 4;
        uint32_t input_temp[16] = {
                0x00220020,0x00001000,0x10000000,0x00000011,
                0x22222020,0x00001800,0x18000000,0x00000011,
                0x22002000,0x00100000,0x00000000,0x00000000,
                0x00330030,0x00001000,0x10000000,0x00000000
        };
        for (int i=0; i<16; i++)
            L.inputMask[i] = input_temp[i];
        uint32_t listOfWords[1] = {3};
        uint32_t listOfBits[1] = {0};
        set_list_of_bits(L.outputMask, listOfWords, listOfBits, 1);

        test_linear_relation_cuda(&L, N, alg, 5);
        //------------------------------
    }
}
void print_latex_eq_and_bias(int rounds, double bias, uint32_t ID[16], uint32_t ODmask[16])
{
    printf("Bias = %f, for ", bias);
    printf("$$\\Delta = (");
    for (int w = 0; w < 16; w++)
    {
        for (int b = 0; b< 32; b++)
        {
            if (get_bit_from_word_and_bit(ODmask, w, b))
                printf("\\Delta X^{(%d)}_{%d,%d} \\oplus ", rounds, w, b);
        }
    }
    printf(" | ");
    for (int w = 0; w < 16; w++)
    {
        for (int b = 0; b< 32; b++)
        {
            if (get_bit_from_word_and_bit(ID, w, b))
                printf("\\Delta X^{(%d)}_{%d,%d}, ", 0, w, b);
        }
    }
    printf(")$$\n");
}



/*Test to find the bias used by Maitra, the differential bias should be 
approximately 0.0272 as in page 16 of the paper.*/
void search_until_find_significant_bias_maitra()
{
    uint32_t ID[16] = { 0 };
    uint32_t ODmask[16] = { 0 };
    uint64_t N = 1;
    uint32_t listOfWords[5] = { 0 };
    uint32_t listOfBits[5] = { 0 };

    //Init odmask
    memset(ODmask, 0x00, sizeof(uint32_t) * 16);
    listOfWords[0] =11; listOfBits[0] =0;
    set_list_of_bits(ODmask, listOfWords, listOfBits, 1);

    printf("w=%d\n", listOfWords[0]);
    memset(ID, 0x00, sizeof(uint32_t) * 16);
    ID[13] = 0x00002000;

    int level = 38, rounds = 3, halfround = 0;
    N = 1;
    N <<= level;
    double comp;
    int flagFirst = 1;
    double bias = 0;
    while (1)
    {
        bias += compute_differential_bias_multiple_devices(rounds,
                                                           ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 0);

        if (flagFirst)
            flagFirst = 0;
        else
            bias /= 2;

        printf("Level %d\n", level);
        N = 1;
        N <<= level;
        comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

        level++;
        if ((fabs(bias) < comp))
            continue;

        //Check again to make sure
        double newbias = compute_differential_bias_multiple_devices(rounds,
                                                                    ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 0);
        if ((fabs(newbias) < comp))
            continue;

        printf("Bias = (%.15f e %.15f)\n", bias, newbias);
        return;
    }
}


// Experimental verification of Differential-Linear Distinguisher 1
void differential_linear_distinguisher_1(int od_word, int od_bit)
{
    uint32_t ID[16] = { 0 };
    uint32_t ODmask[16] = { 0 };
    uint64_t N = 1;
    uint32_t listOfWords[1] = { 0 };
    uint32_t listOfBits[1] = { 0 };

    //Init odmask
    memset(ODmask, 0x00, sizeof(uint32_t) * 16);
    printf("Testing differential bias for OD = X_{%d,%d}\n", od_word,od_bit);
    listOfWords[0] =od_word; listOfBits[0] =od_bit;
    set_list_of_bits(ODmask, listOfWords, listOfBits, 1);

    printf("w=%d\n", listOfWords[0]);
    memset(ID, 0x00, sizeof(uint32_t) * 16);
    ID[2] = 1<<17;
    ID[2] ^= 1<<5;
    ID[6] = 1<<24;
    ID[6] ^= 1<<8;
    ID[10] = 1<<5;
    ID[10] ^= 1<<1;
    ID[14] = 1<<25;
    ID[14] ^= 1<<1;

    int level = 46, rounds = 3, halfround = 1;
    N = 1;
    N <<= level;
    double comp;
    int flagFirst = 1;
    double bias = 0;
    while (1)
    {
        bias += compute_differential_bias_multiple_devices(rounds - 1,
                                                           ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if (flagFirst)
            flagFirst = 0;
        else
            bias /= 2;

        printf("Level %d\n", level);
        N = 1;
        N <<= level;
        comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

        level++;
        if ((fabs(bias) < comp))
            continue;

        //Check again to make sure
        double newbias = compute_differential_bias_multiple_devices(rounds - 1,
                                                                    ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if ((fabs(newbias) < comp))
            continue;

        printf("Bias = (%.15f e %.15f)\n", bias, newbias);
        return;
    }
}


void differential_linear_distinguisher_2()
{
    printf("Differential Distinguisher 2 with no groups\n");
    uint32_t ID[16] = { 0 };
    uint32_t ODmask[16] = { 0 };
    uint64_t N = 1;

    //Init odmask
    memset(ID, 0x00, sizeof(uint32_t) * 16);
    ID[3] = 1 << 25;
    ID[3] = 1 << 5;
    ID[7] = 1 << 28;
    ID[7] = 1 << 12;
    ID[11] = 1 << 25;
    ID[11] = 1 << 21;
    ID[15] = 1 << 21;
    ID[15] = 1 << 13;

    memset(ODmask, 0x00, sizeof(uint32_t) * 16);
    ODmask[13] = 1<<4;
    ODmask[8] = 1<<20;
    ODmask[8] = 1<<19;
    ODmask[2] = 1<<3;
    ODmask[2] = 1<<0;
    ODmask[2] = 1<<4;
    ODmask[7] = 1<<0;
    ODmask[7] = 1<<4;
    ODmask[7] = 1<<20;

    int level = 52, rounds = 3, halfround = 0;
    N = 1;
    N <<= level;
    double comp;
    int flagFirst = 1;
    double bias = 0;
    while (1)
    {
        bias += compute_differential_bias_multiple_devices(rounds - 1,
                                                           ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if (flagFirst)
            flagFirst = 0;
        else
            bias /= 2;

        printf("Level %d\n", level);
        N = 1;
        N <<= level;
        comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

        level++;
        if ((fabs(bias) < comp)) {
            printf("Noo\n");
            continue;
        }

        //Check again to make sure
        double newbias = compute_differential_bias_multiple_devices(rounds - 1,
                                                                    ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if ((fabs(newbias) < comp))
            continue;

        printf("Bias = (%.17g e %.17g)\n", bias, newbias);
        return;
    }
}

void differential_linear_distinguisher_2_group_1()
{
    uint32_t ID[16] = { 0 };
    uint32_t ODmask[16] = { 0 };
    uint64_t N = 1;

    //Init odmask
    memset(ID, 0x00, sizeof(uint32_t) * 16);
    ID[3] = 1 << 25;
    ID[3] = 1 << 5;
    ID[7] = 1 << 28;
    ID[7] = 1 << 12;
    ID[11] = 1 << 25;
    ID[11] = 1 << 21;
    ID[15] = 1 << 21;
    ID[15] = 1 << 13;

    memset(ODmask, 0x00, sizeof(uint32_t) * 16);
    ODmask[2] = 1<<3;
    ODmask[8] = 1<<19;
    ODmask[2] = 1<<0;
    ODmask[13] = 1<<4;
    ODmask[2] = 1<<4;

    int level = 43, rounds = 3, halfround = 0;
    N = 1;
    N <<= level;
    double comp;
    int flagFirst = 1;
    double bias = 0;
    while (1)
    {
        bias += compute_differential_bias_multiple_devices(rounds - 1,
                                                           ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if (flagFirst)
            flagFirst = 0;
        else
            bias /= 2;

        printf("Level %d\n", level);
        N = 1;
        N <<= level;
        comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

        level++;
        if ((fabs(bias) < comp))
            continue;

        //Check again to make sure
        double newbias = compute_differential_bias_multiple_devices(rounds - 1,
                                                                    ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if ((fabs(newbias) < comp))
            continue;

        printf("Bias = (%.15f e %.15f)\n", bias, newbias);
        return;
    }
}

void differential_linear_distinguisher_2_group_2()
{
    uint32_t ID[16] = { 0 };
    uint32_t ODmask[16] = { 0 };
    uint64_t N = 1;

    //Init odmask
    memset(ID, 0x00, sizeof(uint32_t) * 16);
    ID[3] = 1 << 25;
    ID[3] = 1 << 5;
    ID[7] = 1 << 28;
    ID[7] = 1 << 12;
    ID[11] = 1 << 25;
    ID[11] = 1 << 21;
    ID[15] = 1 << 21;
    ID[15] = 1 << 13;

    memset(ODmask, 0x00, sizeof(uint32_t) * 16);
    ODmask[7] = 1<<0;

    int level = 42, rounds = 3, halfround = 0;
    N = 1;
    N <<= level;
    double comp;
    int flagFirst = 1;
    double bias = 0;
    while (1)
    {
        bias += compute_differential_bias_multiple_devices(rounds - 1,
                                                           ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if (flagFirst)
            flagFirst = 0;
        else
            bias /= 2;

        printf("Level %d\n", level);
        N = 1;
        N <<= level;
        comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

        level++;
        if ((fabs(bias) < comp))
            continue;
        //Check again to make sure
        double newbias = compute_differential_bias_multiple_devices(rounds - 1,
                                                                    ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if ((fabs(newbias) < comp))
            continue;

        printf("Bias = (%.15f e %.15f)\n", bias, newbias);
        return;
    }
}

void differential_linear_distinguisher_2_group_3()
{
    uint32_t ID[16] = { 0 };
    uint32_t ODmask[16] = { 0 };
    uint64_t N = 1;

    //Init odmask
    memset(ID, 0x00, sizeof(uint32_t) * 16);
    ID[3] = 1 << 25;
    ID[3] = 1 << 5;
    ID[7] = 1 << 28;
    ID[7] = 1 << 12;
    ID[11] = 1 << 25;
    ID[11] = 1 << 21;
    ID[15] = 1 << 21;
    ID[15] = 1 << 13;

    memset(ODmask, 0x00, sizeof(uint32_t) * 16);
    ODmask[7] = 1<<20;
    ODmask[8] = 1<<20;
    ODmask[7] = 1<<4;

    int level = 46, rounds = 3, halfround = 0;
    N = 1;
    N <<= level;
    double comp;
    int flagFirst = 1;
    double bias = 0;
    while (1)
    {
        bias += compute_differential_bias_multiple_devices(rounds - 1,
                                                           ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if (flagFirst)
            flagFirst = 0;
        else
            bias /= 2;

        printf("Level %d\n", level);
        N = 1;
        N <<= level;
        comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

        level++;
        if ((fabs(bias) < comp))
            continue;

        //Check again to make sure
        double newbias = compute_differential_bias_multiple_devices(rounds - 1,
                                                                    ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 1);
        if ((fabs(newbias) < comp))
            continue;

        printf("Bias = (%.15f e %.15f)\n", bias, newbias);
        return;
    }
}



void differential_linear_distinguisher_3()
{
    uint32_t ID[16] = { 0 };
    uint32_t ODmask[16] = {
            0x00000011,0x00220020,0x00001800,0x18000000,
            0x00000011,0x22222020,0x00001000,0x10000000,
            0x00000000,0x22002000,0x00100000,0x00000000,
            0x00000000,0x00330030,0x00001000,0x10000000
    };
    uint64_t N = 1;
    memset(ID, 0x00, sizeof(uint32_t) * 16);
    ID[2] = 0x00000004;
    ID[6] = 0x20020220;
    ID[10] = 0x40400400;
    ID[14] = 0x40000400;

    int level = 40, rounds = 2, halfround = 1;
    N = 1;
    N <<= level;
    double comp;
    int flagFirst = 1;
    double bias = 0;
    while (1)
    {
        bias += compute_differential_bias_multiple_devices(rounds - 1,
                                                           ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 0);
        if (flagFirst)
            flagFirst = 0;
        else
            bias /= 2;

        printf("Level %d\n", level);
        N = 1;
        N <<= level;
        comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

        level++;
        if ((fabs(bias) < comp))
            continue;

        //Check again to make sure
        double newbias = compute_differential_bias_multiple_devices(rounds - 1,
                                                                    ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 0);
        if ((fabs(newbias) < comp))
            continue;

        printf("Bias = (%.15f e %.15f)\n", bias, newbias);
        return;
    }
}

void differential_linear_distinguisher_5()
{
    uint32_t ID[16] = { 0 };
    uint32_t ODmask[16] = {
            0x00220020,0x00001000,0x10000000,0x00000011,
            0x22222020,0x00001800,0x18000000,0x00000011,
            0x22002000,0x00100000,0x00000000,0x00000000,
            0x00330030,0x00001000,0x10000000,0x00000000
    };
    uint64_t N = 1;
    memset(ID, 0x00, sizeof(uint32_t) * 16);
    uint32_t i = 22;
    ID[3] = 1<< ((uint32_t)(i+28)%32);
    ID[7] = (1<< (uint32_t)(i+31))^(1<<(uint32_t)(i+23))^(1<<(uint32_t)(i+11))^(1<<(uint32_t)(i+3));
    ID[11] = (1<< (uint32_t)(i + 24)) ^ (1<<(uint32_t)(i + 16)) ^ (1<<(uint32_t)(i + 4));
    ID[15] = (1<< (uint32_t)(i + 24)) ^ (1<<(uint32_t)(i + 4));

    int level = 40, rounds = 2, halfround = 1;
    N = 1;
    N <<= level;
    double comp;
    int flagFirst = 1;
    double bias = 0;
    while (1)
    {
        bias += compute_differential_bias_multiple_devices(rounds - 1,
                                                           ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 0);
        if (flagFirst)
            flagFirst = 0;
        else
            bias /= 2;

        printf("Level %d\n", level);
        N = 1;
        N <<= level;
        comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

        level++;
        if ((fabs(bias) < comp))
            continue;

        //Check again to make sure
        double newbias = compute_differential_bias_multiple_devices(rounds - 1,
                                                                    ALG_TYPE_CHACHA, ID, ODmask, N, halfround, 0);
        if ((fabs(newbias) < comp))
            continue;

        printf("Bias = (%.15f e %.15f)\n", bias, newbias);
        return;
    }
}

// key-recovery attack using Differential-Linear Distinguisher 2. 
//Third key-recovery attack
//3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 23, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 51, 52, 53, 54, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 78, 79, 80, 81, 82, 83, 84, 85, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 102, 103, 104, 105, 1
//06, 107, 108, 109, 110, 111, 112, 115, 116, 117, 118, 119, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 135, 136, 140, 141, 142, 143, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 159, 160, 161, 162, 163, 164, 167, 168, 169, 176, 179, 180, 184, 188, 189, 190
//, 191, 192, 193, 194, 195, 196, 197, 198, 204, 205, 206, 207, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 231, 232, 233, 234, 235, 236, 237, 238, 243, 244, 245, 246, 247, 250, 251, 252, 255,
//Number of pnb: 163.
//alpha=76
//comp=1.9073486328125e-06, threshold = 0.340000, varepsilon_a = 5.324321591615444e-06 \varepsilon = 0.000000, data_complexity 87.096922, time_complexity 180.096922.
//alpha=80
//comp=1.9073486328125e-06, threshold = 0.340000, varepsilon_a = 5.4466314395540394e-05 \varepsilon = 0.000000, data_complexity 80.164962, time_complexity 176.164962.
void third_key_recovery_attack()
{
    printf("Third key-recovery attack\n");
    uint32_t ODmask[16] = { 0 };
    uint32_t ID[16] = { 0 };
    uint64_t N = (1 << 20);
    ALGORITHM alg;
    double neutrality_list[256];
    uint32_t pnb[256], number_of_pnb;
    double threshold = 0.34, varepsilon_a, varepsilon_d;

    define_alg(&alg, ALG_TYPE_CHACHA);

    //Go 4, back 2
    memset(ODmask, 0x00, sizeof(uint32_t) * 16);
    N = 1;
    N <<= 16;
    ID[15] = (1 << 29);
    ID[15] = (1 << 9);
    uint32_t finalListOfWords[5] = { 2, 6, 6, 10, 14 };
    uint32_t finalListOfBits[5] = { 0, 7, 19, 12, 0 };
    //uint32_t finalListOfWords[1] = { 2};
    //uint32_t finalListOfBits[1] = { 0};
    set_list_of_bits(ODmask, finalListOfWords, finalListOfBits, 5);

    compute_neutrality_for_every_key_bit(neutrality_list, N, ID, ODmask, 7, 2, alg);
    get_pnb_list(pnb, &number_of_pnb, threshold, neutrality_list);
    for (int i = 0; i < number_of_pnb; i++)
        printf("%d, ", pnb[i]);

    printf("\nNumber of pnb: %d.\n", number_of_pnb);

    N = 1;
    N<<=42;
    varepsilon_a = compute_mean_bias_of_g_cuda(N, ID, ODmask,
                                               7, 2, pnb, number_of_pnb, alg);

    varepsilon_d = pow(2,-34.15);
    double comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);
    double time_complexity, data_complexity;
    double e = varepsilon_a * varepsilon_d;

    compute_complexity_of_the_attack(&data_complexity, &time_complexity, e, number_of_pnb);
    printf("comp=%.17g, threshold = %f, varepsilon_a = %.17g \\varepsilon = %f, data_complexity %f, time_complexity %f.\n", comp, threshold, varepsilon_a, e, data_complexity, time_complexity);
}

void fourth_key_recovery_attack()
{
    printf("Fourth key-recovery attack\n");
    uint32_t ODmask[16] = {
            0x00000000,0x00000000,0x00000000,0x00000001,
            0x00000000,0x00000000,0x00000000,0x00000000,
            0x00000000,0x00000000,0x00000000,0x00000000,
            0x00000000,0x00000000,0x00000000,0x00000000
    };
    uint32_t ID[16] = { 0 };
    uint64_t N = (1 << 21);
    ALGORITHM alg;
    double neutrality_list[256];
    uint32_t pnb[256], number_of_pnb;
    double threshold = 0.265, varepsilon_a, varepsilon_d;

    define_alg(&alg, ALG_TYPE_CHACHA);

    //Go 4, back 3

    N = 1;
    N <<= 14;
    ID[15] = (1 << 22);

    compute_neutrality_for_every_key_bit(neutrality_list, N, ID, ODmask, 7, 3, alg);
    get_pnb_list(pnb, &number_of_pnb, threshold, neutrality_list);
    for (int i = 0; i < number_of_pnb; i++)
        printf("%d, ", pnb[i]);

    printf("\nNumber of pnb: %d.\n", number_of_pnb);

    N = 1;
    N<<=43;
    varepsilon_a = compute_mean_bias_of_g_cuda(N, ID, ODmask,
                                               7, 3, pnb, number_of_pnb, alg);

    varepsilon_d = pow(2,-37.98);
    double time_complexity, data_complexity;
    double e = varepsilon_a * varepsilon_d;
    double comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);

    compute_complexity_of_the_attack(&data_complexity, &time_complexity, e, number_of_pnb);
    printf("comp=%.17g, threshold = %f, varepsilon_a = %.17g \\varepsilon = %f, data_complexity %f, time_complexity %f.\n", comp, threshold, varepsilon_a, e, data_complexity, time_complexity);
}
// 134 0, 7, 8, 20, 21, 31, 35, 39, 44, 45, 46, 47, 51, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 71, 72, 73, 74, 77, 80, 83, 84, 85, 86, 89, 90, 91, 95, 99, 100, 104, 108, 109, 110, 115, 123, 124, 125, 126, 127, 128, 129, 130, 135, 140, 141, 142, 147, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 198, 199, 200, 201, 204, 205, 206, 207, 210, 211, 212, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 231, 244, 245, 246, 247, 248, 252, 255
// comp=3.814697265625e-06, threshold = 0.280000, varepsilon_a = 6.9297329901019111e-06 \varepsilon = 0.000000, data_complexity 86.237379, time_complexity 208.237379, alpha=48
void five_key_recovery_attack()
{
    printf("Five key-recovery attack\n");
    uint32_t ODmask[16] = { 0 };
    uint32_t ID[16] = { 0 };
    uint64_t N = (1 << 20);
    ALGORITHM alg;
    double neutrality_list[256];
    uint32_t pnb[256], number_of_pnb;
    double threshold = 0.28, varepsilon_a, varepsilon_d;

    define_alg(&alg, ALG_TYPE_CHACHA);

    //Go 4, back 2
    memset(ODmask, 0x00, sizeof(uint32_t) * 16);
    N = 1;
    N <<= 16;
    ID[15] = (1 << 29);
    ID[15] = (1 << 9);
    uint32_t finalListOfWords[5] = { 2, 6, 6, 10, 14 };
    uint32_t finalListOfBits[5] = { 0, 7, 19, 12, 0 };
    //uint32_t finalListOfWords[1] = { 2};
    //uint32_t finalListOfBits[1] = { 0};
    set_list_of_bits(ODmask, finalListOfWords, finalListOfBits, 5);

    compute_neutrality_for_every_key_bit(neutrality_list, N, ID, ODmask, 8, 3, alg);
    get_pnb_list(pnb, &number_of_pnb, threshold, neutrality_list);
    for (int i = 0; i < number_of_pnb; i++)
        printf("%d, ", pnb[i]);

    printf("\nNumber of pnb: %d.\n", number_of_pnb);

    N = 1;
    N<<=40;
    varepsilon_a = compute_mean_bias_of_g_cuda(N, ID, ODmask,
                                               8, 3, pnb, number_of_pnb, alg);

    varepsilon_d = pow(2,-34.15);
    double comp = (2 * (0.5 + 4 * sqrt(0.25 / N)) - 1);
    double time_complexity, data_complexity;
    double e = varepsilon_a * varepsilon_d;

    compute_complexity_of_the_attack(&data_complexity, &time_complexity, e, number_of_pnb);
    printf("comp=%.17g, threshold = %f, varepsilon_a = %.17g \\varepsilon = %f, data_complexity %f, time_complexity %f.\n", comp, threshold, varepsilon_a, e, data_complexity, time_complexity);
}


int main()
{
    srand(time(NULL));
    if (five_key_recovery == 0) {
        if (testChachaOnGPU())
        {
            printf("Error, aborting\n");
            return(1);
        }
    }

    // key-recovery attack using Differential-Linear Distinguisher 1
    first_key_recovery_attack();
    // key-recovery attack using Differential-Linear Distinguisher 4.
    second_key_recovery_attack();
    // key-recovery attack using Differential-Linear Distinguisher 2 (7 rounds)
    third_key_recovery_attack();
    // key-recovery attack using Differential-Linear Distinguisher 6 (7 rounds)
    fourth_key_recovery_attack();
    // key-recovery attack using Differential-Linear Distinguisher 2 (7.25 rounds)
    five_key_recovery_attack(); // ATTENTION:: to run this key recovery attack make sure five_key_recovery=1. Set five_key_recovery=0 for all another experiments


    return 0;
}