/*
* This is an example of functionalities that might be useful
* to the researcher trying to find differentials or 
* linear approximations to ChaCha, Salsa or Forro.
*/
#include "arx_cryptanalysis.cuh"

int my_rank, num_procs;


void differential_linear_distinguisher_1()
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {2, 2, 6, 6, 10, 10, 14, 14},{17, 5, 24, 8, 5, 1, 25, 1},2, 8}, //id
            {{0}, {11},{0},7,1}, //od
            {"differential_linear_distinguisher_1",0.000002480046533, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 42;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_1() //0.008051972105325
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    { 8, 8 },
                    {20, 19},
                    6, 2

            }, //od
            {"differential_linear_distinguisher_2_1",2^-33, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 41;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_1_35() //0.008051972105325
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {2,7},
                    {0,0},
                    6, 2

            }, //od
            {"differential_linear_distinguisher_2_1",2^-33, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 47;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_1_3() //-0.000793696319306
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {8, 8, 2, 2, 2},
                    {20, 19, 3, 0, 4},
                    6, 5

            }, //od
            {"differential_linear_distinguisher_2_1_3",2^-13, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<=39;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_1_2() //0.001327190414789
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    { 8, 8, 13},
                    {20, 19, 4},
                    6, 3

            }, //od
            {"differential_linear_distinguisher_2_1_2",2^-11, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 38;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}
void differential_linear_distinguisher_2_2() //  0.109562245215178
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {13},
                    {4},
                    6, 1

            }, //od
            {"differential_linear_distinguisher_2_2",2^-33, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 41;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_3() //-0.026363534583955
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {2, 2, 2},
                    {3, 0, 4},
                    6, 3

            }, //od
            {"differential_linear_distinguisher_2_3",0.026363534583955, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 41;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_1_2_3() // 0.000144687348438
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {8, 8, 13, 2, 2, 2}, //7, 7, 7},
                    {20, 19, 4, 3, 0, 4}, //0, 4, 20}, 6, 9
                    //{13, 8},
                    //{4, 20},
                    6, 6

            }, //od
            {"differential_linear_distinguisher_2_2_3",2^-16, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 41;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_2_3() //0.005073826096426
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {13, 2, 2, 2}, //7, 7, 7},
                    {4, 3, 0, 4}, //0, 4, 20}, 6, 9
                    //{13, 8},
                    //{4, 20},
                    6, 4

            }, //od
            {"differential_linear_distinguisher_2_2_3",0.005073826096426, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));
    //printf("Aqui=%d\n", diff.output.words[0]);
    //printf("Aqui=%d, size=%d\n", diff.output.words[1], diff.output.number_of_bits);

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 41;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_3_4() //
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {2, 2, 2, 7, 7, 7},
                    {3, 0, 4, 0, 4, 20}, //6, 9
                    //{13, 8},
                    //{4, 20},
                    6, 6

            }, //od
            {"differential_linear_distinguisher_2_3_4",2^-23, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 47;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_2_3_4() //
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {13, 2, 2, 2, 7, 7, 7},
                    {4, 3, 0, 4, 0, 4, 20}, //6, 9
                    //{13, 8},
                    //{4, 20},
                    6, 7

            }, //od
            {"differential_linear_distinguisher_2_3_4",2^-23, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 51;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_4() // 0.000004901891164
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {7, 7, 7},
                    {0, 4, 20}, //6, 9
                    //{13, 8},
                    //{4, 20},
                    6, 3

            }, //od
            {"differential_linear_distinguisher_2_4",0.000004767074643, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 42;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void differential_linear_distinguisher_2_2_4() // 0.000000683569889
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    {13, 7, 7, 7},
                    {4, 0, 4, 20}, //6, 9
                    //{13, 8},
                    //{4, 20},
                    6, 4

            }, //od
            {"differential_linear_distinguisher_2_2_4",2^-20, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 44;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}


void differential_linear_distinguisher_2_1_4() //0.000000158730987
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    { 8, 8,7, 7, 7}, //13, 2, 2, 2, 7, 7, 7},
                    {20, 19,0, 4, 20}, //4, 3, 0, 4, 0, 4, 20}, 6, 9
                    //{13, 8},
                    //{4, 20},
                    6, 5

            }, //od
            {"differential_linear_distinguisher_2_1_4",2^-24, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 49;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}


void differential_linear_distinguisher_2_1_2_4() //
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {3, 3, 7, 7, 11, 11, 15, 15},{25, 5, 28, 12, 25, 21, 21, 13}, 2, 8}, //id
            {
                    {0},
                    { 8, 8, 13, 7, 7, 7}, //13, 2, 2, 2, 7, 7, 7},
                    {20, 19, 4, 0, 4, 20}, //4, 3, 0, 4, 0, 4, 20}, 6, 9
                    //{13, 8},
                    //{4, 20},
                    6, 6

            }, //od
            {"differential_linear_dist_2_1_2_4",2^-27, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));

    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 55;
    diff.correlation.number_of_trials = 3*diff.correlation.number_of_trials;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void example_differential_correlation()
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {14},{6},0, 1}, //id
            {{0}, {3,4},{0,0},6,2}, //od
            {"Sec 3.1 of [Coutinho 20]",0.00048, 0, 0, 0}
        };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));
    
    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 38;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

void lemma_1() // -6.5850
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {12,11},{0,0},0, 2}, //id
            {
                    {0},
                    {0, 0, 1, 1, 1, 2, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 14, 14, 15, 15, 15},
                    {8, 24, 0, 8, 24, 0, 0, 7, 19, 26, 2, 3, 14, 15, 19, 22, 23, 30, 31, 0, 6, 7, 12, 19, 0, 7, 8, 12, 27, 28, 0, 23, 24, 0, 8, 16, 0, 8, 0, 7, 14, 16, 24},
                    4,43
            }, //od
            {"Lemma 1",0.010416179111417, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 40;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void lemma_2_qr_1() //-7.58
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {0,0,4,4,4,8,12,12,12},{8,24, 7, 19, 26, 0, 0, 8, 16},12, 9}, //id
            {
                    {0},
                    {0,4,4,4,4,4,4,4,4,4,4,8,8,8,8,8,8,8,8,8,8,12,12,12,12,12,12,12},
                    {0,6,7,10,11,13,22,23,27,30,31,3,4,6,8,15,16,19,20,26,31,7,8,18,19,23,25,26 },
                    14,28
            }, //od
            {"Lemma 2 Quarter Round 1", 0.005210383681212, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 42;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void lemma_2_qr_2() // -11.76
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {1,1,1,5,5,5,5,5,5,9,9,9,13,13},{0,8,24,2,3,14,15,22,23,7,12,19,0,8},12, 14}, //id
            {
                    {0},
                    {1,1,1,1,1,1,1,5,5,5,5,5,5,5,5,5,5,5,9,9,9,9,9,9,9,9,9,13,13,13,13,13,13,13,13,13,13,13,13,13},
                    {6,7,11,12,16,18,19,1,2,9,11,19,21,22,26,27,30,31,0,3,4,8,20,22,23,26,27,2,3,7,8,11,12,18,20,21,22,23,26,27},
                    14,40
            }, //od
            {"Lemma 2 Quarter Round 2", 0.000288356850736, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 42;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void lemma_2_qr_3() // -6.74
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {2,6,6,10,10,10,10,10,10,14,14,14},{0,30,31,0,7,8,12,27,28,0,6,7},12, 12}, //id
            {
                    {0},
                    {2,2,2,2,2,2,2,2,2,2,6,6,6,6,6,10,10,10,10,10,10,10,10,14,14,14,14,14,14,14,14,14,14,14,14},
                    {0,6,8,11,12,16,22,23,27,28,13,14,17,18,19,6,8,10,11,27,28,30,31,3,4,7,11,12,15,16,19,20,24,27,28},
                    14,35
            }, //od
            {"Lemma 2 Quarter Round 3", 0.009331808910247, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 42;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void lemma_2_qr_4() // -4.58
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {3,11,11,11,15,15,15},{0,0,23,24,14,16,24},12, 7}, //id
            {
                    {0},
                    {3,3,3,3,3,3,7,7,7,7,7,7,7,7,11,11,11,11,11,11,15,15,15,15,15,15},
                    {0,8,14,16,23,30,7,19,20,21,22,23,30,31,0,12,13,14,15,16,0,6,16,23,24,31},
                    14,26
            }, //od
            {"Lemma 2 Quarter Round 4", 0.041666479221021, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 42;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}


void lemma_4_qr_1() //-
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {0,0,4,8,12,12},{11,12,7,12,7,19},12, 6}, //id
            {
                    {0},
                    {0,0,0,0,4,4,4,4,4,4,4,8,8,8,8,8,12,12,12,12,12,12},
                    {3,7,19,23,13,14,18,19,25,30,31,6,12,18,23,24,6,7,10,19,20,31},
                    14,22
            }, //od
            {"Lemma 4 Quarter Round 1", 0.010406723426058, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 41;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void lemma_4_qr_2() // -
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {9,13},{0,0},12, 2}, //id
            {
                    {0},
                    {1,5,13,13,13},
                    {16,7,0,8,24},
                    14,5
            }, //od
            {"Lemma 4 Quarter Round 2", 1, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 40;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void lemma_4_qr_3() // -
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {2,6,6,10,14},{0,6,26,12,24},12, 5}, //id
            {
                    {0},
                    {2,2,2,2,2,6,6,6,6,6,6,10,10,10,10,14,14,14,14,14,14,14,14,14,14},
                    {0,8,11,12,24,7,13,19,25,30,31,18,23,24,26,0,5,6,11,12,16,19,20,25,26},
                    14,25
            }, //od
            {"Lemma 4 Quarter Round 3", 0.041657735142508, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 42;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void lemma_4_qr_4() // -10.17
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {3,3,7,7,11,11,11,11,15,15,15,15},{0,16,7,19,6,7,18,31,11,12,19,20},12, 12}, //id
            {
                    {0},
                    {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,7,7,7,7,7,7,7,7,11,11,11,11,11,11,11,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15},
                    {0,3,4,6,7,11,12,16,17,18,19,20,27,28,30,31,2,3,6,7,18,22,23,27,6,11,18,19,20,27,28,0,3,4,5,7,11,12,14,16,18,19,25,26,30,31},
                    14,46
            }, //od
            {"Lemma 4 Quarter Round 4", 0.000869511311976, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 42;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void linear_correlation_4() // -23.199105672195753 3*(1<<49)
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {2,2,2,7,7,7,8,8,13},{4,3,0,20,4,0,20,19,4},2, 9}, //id
            {
                    {0},
                    {0,0,0,0,0,0,1,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,4,4,4,5,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,10,10,10,10,11,11,11,12,12,12,12,12,12,12,12,13,13,13,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15},
                    {23,22,19,7,3,2,16,24,23,12,8,0,31,28,20,16,12,11,7,6,4,3,0,31,19,14,7,31,25,19,13,7,27,26,23,22,7,6,3,2,24,12,26,25,24,18,28,27,20,31,30,20,19,11,10,7,6,24,8,0,26,20,16,12,11,6,5,0,31,19,18,16,14,12,11,7,4,0},
                    10,78
            }, //od
            {"linear_correlation_4 A'->B'",1.4901161193847656e-08, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 49;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void linear_correlation_7_second_subcipher() //
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {
                {0},
                {0,0,0,0,0,0,1,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,4,4,4,5,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,10,10,10,10,11,11,11,12,12,12,12,12,12,12,12,13,13,13,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15},
                {23,22,19,7,3,2,16,24,23,12,8,0,31,28,20,16,12,11,7,6,4,3,0,31,19,14,7,31,25,19,13,7,27,26,23,22,7,6,3,2,24,12,26,25,24,18,28,27,20,31,30,20,19,11,10,7,6,24,8,0,26,20,16,12,11,6,5,0,31,19,18,16,14,12,11,7,4,0},
                14, 78
             }, //id
            {
                    {0},
                    {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,5,5,5,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15},
                    {31,23,16,14,12,11,6,4,3,0,31,30,20,19,16,11,10,7,6,12,31,28,26,11,7,6,5,4,28,27,26,24,19,16,12,10,8,7,0,31,15,3,31,28,27,25,19,11,5,24,23,20,18,15,14,12,7,6,4,3,2,27,26,23,22,8,6,3,2,0,30,28,27,20,16,15,14,12,7,4,0,26,24,23,19,18,3,31,28,25,20,19,16,15,13,7,28,27,26,23,22,20,19,15,14,4,3,23,16,12,8,28,27,22,21,16,10,4,0,30,28,27,26,24,20,18,17,16,15,3,2,0},
                    15,132
            }, //od
            {"LT7_second_subcipher 6.5 -> 7  (B'->C')", 0.000000226505605, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 49;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}


void lemma_3() // -5.5850
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {2,2,2,7,7,7,8,8,13},{4,3,0,20,4,0,20,19,4},6, 9}, //id
            {
                    {0},
                    {0, 0, 2, 3, 3, 4, 6, 6, 7, 7, 8, 9, 10, 11, 11, 11, 11, 12, 12, 13, 14, 15, 15, 15, 15},
                    {11, 12, 0, 0, 16, 7, 6, 26, 7, 19, 12, 0, 12, 6, 7, 18, 31, 7, 19, 0, 24, 11, 12, 19, 20},
                    12,25
            }, //od
            {"Lemma 3", 0.020832863796992, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 41;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void lemma_2_partition_1_step_by_step()
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {8},{0},0, 1}, //id
            {
                    {0},
                    {0,8,12,12},
                    {0,0,0,8},
                    2,4
            }, //od
            {"Lemma 2 Partition 1",0.0078125, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 38;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void linear_trail_3_subcipher_1() //15.34
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {11},{0},7, 1}, //id
            {
                    {0},
                    {0,0,0,1,2,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,5,6,6,6,7,8,8,8,8,8,8,8,9,9,9,9,9,9,10,10,10,10,10,10,10,11,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,15,15,15},
                    {24,16,0,24,7,24,16,15,0,31,20,19,6,4,3,31,27,20,19,15,12,4,12,11,10,12,26,24,23,19,8,7,0,24,15,12,8,3,0,31,30,28,27,12,11,8,24,24,16,24,19,18,16,12,11,7,28,23,16,12,8,7,0,31,24,8},
                    13,66 //100
            }, //od
            {"linear trail 3 rounds 3.5 to 6.5", 0.000024167308766, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 41;
    //lin.correlation.number_of_trials = lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void linear_trail_3_subcipher_2() //20.17
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0},
             {0,0,0,1,2,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,5,6,6,6,7,8,8,8,8,8,8,8,9,9,9,9,9,9,10,10,10,10,10,10,10,11,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,15,15,15},
             {24,16,0,24,7,24,16,15,0,31,20,19,6,4,3,31,27,20,19,15,12,4,12,11,10,12,26,24,23,19,8,7,0,24,15,12,8,3,0,31,30,28,27,12,11,8,24,24,16,24,19,18,16,12,11,7,28,23,16,12,8,7,0,31,24,8},
             13, 66}, //id
            {
                    {0},
                    {0,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,11,11,11,12,12,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,15,15,15,15},
                    {0,19,18,16,12,11,7,28,23,16,12,8,6,0,31,16,8,0,31,30,27,26,23,22,13,11,10,7,6,31,30,27,26,22,19,11,6,2,19,18,17,14,31,30,23,19,7,31,26,23,20,18,16,15,8,6,4,3,31,27,20,19,14,11,8,7,4,3,2,0,31,28,11,10,8,16,12,0,26,25,19,8,27,26,20,19,12,8,3,28,24,20,16,15,12,4,24,16,7,0},

                    14,100
            }, //od
            {"linear trail 3 rounds 6.5 to 7", 9.5367431640625e-07, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 43;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void linear_trail_4()
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {2,2,2,7,7,7,8,8,13},{4,3,0,20,4,0,20,19,4},6, 9}, //id
            {
                    {0},
                    {0,0,0,0,0,0,1,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,4,4,4,5,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,10,10,10,10,11,11,11,12,12,12,12,12,12,12,12,13,13,13,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15},
                    {23,22,19,7,3,2,16,24,23,12,8,0,31,28,20,16,12,11,7,6,4,3,0,31,19,14,7,31,25,19,13,7,27,26,23,22,7,6,3,2,24,12,26,25,24,18,28,27,20,31,30,20,19,11,10,7,6,24,8,0,26,20,16,12,11,6,5,0,31,19,18,16,14,12,11,7,4,0},
                    14,78
            }, //od
            {"Linear Trail 4",1.1920928955078125e-07, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 53;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void linear_trail_5() //-11.19
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {0,0,1,1,1,2,2,3,3,4,4,5,5,5,5,5,5,6,7,9,9,9,10,13,13,13,13,13,13,14,15},
             {0, 4, 5, 17, 21, 11, 12, 27, 28, 0, 4, 5, 13, 17, 21, 25, 29, 12, 28, 13, 25, 29, 20, 4, 5, 16, 17, 20, 21, 12, 28}, 5, 31}, //id
            {
                    {0},
                    {0},
                    {0},
                    8,1
            }, //od
            {"Linear Trail 5",0.000426475446754, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 40;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

void linear_trail_6()
{
    linear_approximation_t lin = {
            ALG_TYPE_CHACHA,
            {{0}, {0,0,0,1,2,3,3,4,4,4,4,4,4,5,5,6,6,7,7,8,8,8,9,12,12,12,12,12,12,13,14},
             {5, 17, 21, 12, 28, 0, 4, 5, 13, 17, 21, 25, 29, 11, 12, 27, 28, 0, 4, 13, 25, 29, 20, 4, 5, 16, 17, 20, 21, 12, 28}, 5, 31}, //id
            {
                    {0},
                    {3},
                    {0},
                    8,1
            }, //od
            {"Linear Trail 6",0.000426475446754, 0, 0, 0}
    };

    lob_compute_mask_from_list_of_bits(&(lin.input));
    lob_compute_mask_from_list_of_bits(&(lin.output));

    lin.correlation.number_of_trials = 1;
    lin.correlation.number_of_trials <<= 40;
    lin.correlation.number_of_trials = 3*lin.correlation.number_of_trials;
    compute_differential_or_linear_correlation_original(&lin, TYPE_LINEAR);
    if(my_rank == 0)
        la_print(NULL, lin);
}

int main()
{
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    //differential_linear_distinguisher_1();
    //differential_linear_distinguisher_2_1();
    //differential_linear_distinguisher_2_2();
    differential_linear_distinguisher_2_3();
    differential_linear_distinguisher_2_4();
    //differential_linear_distinguisher_2_1_2();
    //differential_linear_distinguisher_2_1_3();
    //differential_linear_distinguisher_2_1_4();
    //differential_linear_distinguisher_2_2_4();
    differential_linear_distinguisher_2_2_3();
    //differential_linear_distinguisher_2_3_4();
    //differential_linear_distinguisher_2_1_2_3();
    //differential_linear_distinguisher_2_2_3_4();
    //differential_linear_distinguisher_2_1_2_4();

    lemma_1();
    lemma_2_qr_1();
    //lemma_2_qr_2();
    //lemma_2_qr_3();
    //lemma_2_qr_4();
    //lemma_3();
    //lemma_4_qr_1();
    //lemma_4_qr_2();
    lemma_4_qr_3();
    //lemma_4_qr_4();
    //linear trail 1 comes from lemma_1 and lemma 2
    //linear trail 2 comes from lemma 3 and lemma 4
    //linear_trail_3_subcipher_1();
    //linear_trail_3_subcipher_2();
    //linear_trail_4();
    linear_trail_5();
    //linear_trail_6();
    //linear_correlation_7_first_subcipher is linear_trail_4
   // linear_correlation_7_second_subcipher();
    MPI_Finalize();
    return 0;
}