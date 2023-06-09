#ifndef CONSTANTS_CUH
#define CONSTANTS_CUH

#define RV_SUCESS 0
#define RV_ERROR 1

#define TRUE 1
#define FALSE 0

//CUDA
#define NUMBER_OF_DEVICES_PER_MACHINE 8
#define NUMBER_OF_TESTS_PER_THREAD (1<<15)//(1<<15)
#define NUMBER_OF_CUDA_THREADS (1<<8)//(1<<10) //(1<<8)
#define NUMBER_OF_CUDA_BLOCKS (1<<6)//68//(1<<6) //Must be >= (1<<7)!!!
#define TOTAL_EXECUTIONS_PER_KERNEL (NUMBER_OF_TESTS_PER_THREAD * NUMBER_OF_CUDA_THREADS * NUMBER_OF_CUDA_BLOCKS)

#define MAX_LEVEL 50 //We call LEVEL the log2(N), where N is the number of trials when computing a correlation.
#define MIN_LEVEL 34 //We call LEVEL the log2(N), where N is the number of trials when computing a correlation.
#define MAX_BITS_IN_LIST_OF_BITS 134

#define STATE_SIZE 16
#define STATE_SIZE_IN_BYTES (STATE_SIZE * sizeof(uint32_t))
#define STATE_SIZE_IN_BITS (STATE_SIZE_IN_BYTES * 8 )
#define KEY_SIZE 8
#define KEY_SIZE_IN_BITS (256)
#define NONCE_SIZE 2
#define CTR_SIZE 2
#define NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS 65536

#endif