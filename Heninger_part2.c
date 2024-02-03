#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>


// Define the ValidSolution struct before using it
typedef struct {
    char slice[5];
} ValidSolution;

// Read files that store known bits
void readKnownBits(const char* filename, int known_bits[32]) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return;
    }

    memset(known_bits, -1, sizeof(int) * 32); // Initialize all bits to -1

    int index, value;
    while (fscanf(file, "%d %d", &index, &value) == 2) {
        if (index < 1024) {
            known_bits[index] = value;
        }
    }

    fclose(file);
}

void read_and_convert_component(mpz_t component, FILE *file) {
    char buffer[2048]; 
    if (fgets(buffer, sizeof(buffer), file) != NULL) {

        mpz_set_str(component, strchr(buffer, '=') + 1, 10);
    }
}




// Utility function to get the ith bit of a GMP number
char getGmpBit(mpz_t number, int i) {
    return mpz_tstbit(number, i) ? 1 : 0;
}

// Equation 8
bool equation8(mpz_t N, mpz_t p0, mpz_t q0, char p_i, char q_i , int i) {
    mpz_t term;
    mpz_init(term);
    mpz_mul(term, p0 , q0);          // term = p0 * q0
    mpz_sub(term, N, term);          // term = N - p0 * q0
    bool result = getGmpBit(term, i) == (p_i ^ q_i);
    mpz_clear(term);
    return result;
}

// Equation 9
bool equation9(mpz_t N, mpz_t e, mpz_t k, int tau_k, mpz_t p0, mpz_t q0, mpz_t d0,
            char p_i, char q_i, char d_tau_k_i , int i) {
    mpz_t term, temp;
    mpz_init(term);
    mpz_init(temp);

    // Calculate k * (N + 1)
    mpz_add_ui(temp, N, 1);          // temp = N + 1
    mpz_mul(term, k, temp);          // term = k * (N + 1)

    // Calculate k * (p0 + q0)
    mpz_add(temp, p0 , q0);          // temp = p0 + q0
    mpz_mul(temp, k, temp);          // temp = k * (p0 + q0)
    mpz_sub(term, term, temp);       // term = k * (N + 1) - (k * (p0 + q0))

    // Calculate e * d0
    mpz_mul(temp, e, d0);            // temp = e * d0
    mpz_sub(term, term, temp);       // term = k * (N + 1) - (k * (p0 + q0)) - e * d0
    mpz_add_ui(term, term, 1);       // term += 1
              
    bool result = getGmpBit(term, i + tau_k) == (d_tau_k_i ^ p_i ^ q_i);
    mpz_clear(term);
    mpz_clear(temp);

    return result;
}


// Equation 10
bool equation10(mpz_t e, mpz_t kp, int tau_kp, mpz_t p0, mpz_t dp0, char p_i, char dp_i,  int i) {
    mpz_t term, e_term, temp;
    mpz_init(term);
    mpz_init(temp);
    mpz_init(e_term);

    // Calculate term focused on (i + tau_kp)-th bit
    mpz_mul(e_term, e, dp0);
    mpz_sub_ui(temp, p0, 1);
    mpz_mul(term, kp, temp);
    mpz_add_ui(term, term, 1);
    mpz_sub(term, term, e_term);

    bool result = getGmpBit(term, i + tau_kp) == ((dp_i ^ p_i) );

    mpz_clear(term);
    mpz_clear(e_term);
    mpz_clear(temp);
    return result;
}


// Equation 11
bool equation11(mpz_t e, mpz_t kq, int tau_kq, mpz_t q0, mpz_t dq0, char q_i, char dq_i, int i) {

    mpz_t term, e_term, temp;
    mpz_init(term);
    mpz_init(e_term);
    mpz_init(temp);

    // Calculate term focused on (i + tau_kq)-th bit
    mpz_mul(e_term, e, dq0);
    mpz_sub_ui(temp, q0, 1);
    mpz_mul(term, kq, temp);
    mpz_add_ui(term, term, 1);
    mpz_sub(term, term, e_term);

    bool result = getGmpBit(term, i + tau_kq) == ((dq_i ^ q_i));

    mpz_clear(temp);
    mpz_clear(term);
    mpz_clear(e_term);

    return result;
}


// Compute τ(x), the exponent of the largest power of 2 dividing x
int tau(const mpz_t x) {
    mpz_t temp;
    mpz_init_set(temp, x);
    int count = 0;

    while (mpz_even_p(temp)) {
        mpz_fdiv_q_2exp(temp, temp, 1);  // temp = temp / 2
        count++;
    }

    mpz_clear(temp);
    return count;
}


// Function to correct the least significant bits for dp, dq, d
void correct_lsb(mpz_t result, const mpz_t e, const mpz_t k, int tau_k, bool is_d) {
    mpz_t modulus;

    mpz_init(modulus);
    int power = is_d ? 2 + tau_k : 1 + tau_k;
    mpz_ui_pow_ui(modulus, 2, power);  // modulus = 2^(power)
    mpz_invert(result, e, modulus);
    mpz_clear(modulus);
}

// Core function for the first phase of the RSA key reconstruction algorithm
void reconstruct_rsa_key_first_phase(mpz_t p, mpz_t q, mpz_t d, mpz_t dp, mpz_t dq, const mpz_t e, const mpz_t k, const mpz_t kp, const mpz_t kq) {
    mpz_setbit(p, 0);  // p is odd
    mpz_setbit(q, 0);  // q is odd

    int tau_kp = tau(kp);
    int tau_kq = tau(kq);
    int tau_k = tau(k);  
    correct_lsb(dp, e, kp, tau_kp, false);
    correct_lsb(dq, e, kq, tau_kq, false);
    correct_lsb(d, e, k, tau_k, true);
}

// Calculating p and q from the my_dp and my_dq
void compute_qp_from_dpq( mpz_t result, mpz_t e,  mpz_t k, mpz_t temp_dp){

    mpz_mul(result, e , temp_dp);
    mpz_sub_ui(result, result, 1);
    mpz_divexact(result, result, k);
    mpz_add_ui(result,result, 1);

}

// The core of the Heninger-Shacham algorithm
void BranchAndPrune(
    mpz_t result_p, mpz_t result_q, mpz_t my_p, mpz_t my_q, mpz_t my_d, mpz_t my_dp, mpz_t my_dq,
    mpz_t e, mpz_t k, mpz_t kp, mpz_t kq, mpz_t N,
    int tau_k, int tau_kp, int tau_kq,
    char possibilities[][5], int num_possibilities,
    int known_bits_p[], int known_bits_q[], int known_bits_d[], int known_bits_dp[], int known_bits_dq[], bool verbose, 
    int counter) {
    ValidSolution validSolutions[2]; 
    int validSolutionsCount = 0;
    mpz_t result_n;
    mpz_init(result_n);


    for (int i = 0; i <= num_possibilities; i++) {
        char* possibility = possibilities[i];

        if (equation8(N, my_p, my_q, possibility[0], possibility[1], counter) &&
            equation9(N, e, k, tau_k, my_p, my_q, my_d, possibility[0], possibility[1], possibility[2], counter) &&
            equation10(e, kp, tau_kp, my_p, my_dp, possibility[0], possibility[3], counter) &&
            equation11(e, kq, tau_kq, my_q, my_dq, possibility[1], possibility[4], counter)) {
            
            bool matchesKnownBits = true;
            if ((known_bits_p[counter] != -1 && known_bits_p[counter] != possibility[0]) ||
                (known_bits_q[counter] != -1 && known_bits_q[counter] != possibility[1]) ||
                (known_bits_d[counter + tau_k] != -1 && known_bits_d[counter + tau_k] != possibility[2]) ||
                (known_bits_dp[counter + tau_kp] != -1 && known_bits_dp[counter + tau_kp] != possibility[3]) ||
                (known_bits_dq[counter + tau_kq] != -1 && known_bits_dq[counter + tau_kq] != possibility[4])) {
                matchesKnownBits = false;
            }

            if (matchesKnownBits) {


                for (int j = 0; j < 5; ++j) {
                    validSolutions[validSolutionsCount].slice[j] = possibility[j] + '0';
                }
                
                validSolutionsCount++;
                
                if (verbose){
                    printf("Valid combination for Slice(%d): p[%d]=%d, q[%d]=%d, d[%d]=%d, dp[%d]=%d, dq[%d]=%d\n",
                        counter,counter, possibility[0], counter, possibility[1], tau_k + counter, possibility[2], tau_kp+counter, possibility[3], tau_kq+counter, possibility[4]);
                    gmp_printf("------------------------------------------------------------\nThe value of my p is : %Zd\n", my_p);
                    gmp_printf("The value of my q is : %Zd\n", my_q);
                    gmp_printf("The value of my d is : %Zd\n", my_d);
                    gmp_printf("The value of my dp is : %Zd\n", my_dp);
                    gmp_printf("The value of my dq is : %Zd\n------------------------------------------------------------\n", my_dq);
                }
                // Correctly access the last valid solution for printing
                ValidSolution lastValidSolution = validSolutions[validSolutionsCount - 1];

                compute_qp_from_dpq(result_p, e, kp, my_dp );
                compute_qp_from_dpq(result_q, e, kq, my_dq );
                mpz_mul(result_n, result_p, result_q);
                if (mpz_cmp(result_n, N) == 0 ) {

                    gmp_printf("The resulted P is : %Zu\n", result_p);
                    gmp_printf("The resulted Q is : %Zu\n", result_q);
                    goto finish;
                }

            }
        }

    }

    for (int i = 0; i < validSolutionsCount; i++) {
        // Clone the mpz_t variables to avoid modifying the originals prematurely
        mpz_t cloned_my_p, cloned_my_q, cloned_my_d, cloned_my_dp, cloned_my_dq;
        mpz_init_set(cloned_my_p, my_p);
        mpz_init_set(cloned_my_q, my_q);
        mpz_init_set(cloned_my_d, my_d);
        mpz_init_set(cloned_my_dp, my_dp);
        mpz_init_set(cloned_my_dq, my_dq);

        // Process each valid solution
        char* slice = validSolutions[i].slice;
        if (slice[0] == '1') mpz_setbit(cloned_my_p, counter);
        if (slice[1] == '1') mpz_setbit(cloned_my_q, counter);
        if (slice[2] == '1') mpz_setbit(cloned_my_d, counter + tau_k);
        if (slice[3] == '1') mpz_setbit(cloned_my_dp, counter + tau_kp);
        if (slice[4] == '1') mpz_setbit(cloned_my_dq, counter + tau_kq);

        // Recursive call to BranchAndPrune with incremented counter
        BranchAndPrune(result_p, result_q, cloned_my_p, cloned_my_q, cloned_my_d, cloned_my_dp, cloned_my_dq,
                     e, k, kp, kq, N, tau_k, tau_kp, tau_kq,
                     possibilities, num_possibilities,
                     known_bits_p, known_bits_q, known_bits_d, known_bits_dp, known_bits_dq, verbose,
                     counter + 1); 

        // Clear the cloned mpz_t variables to avoid memory leaks
        mpz_clear(cloned_my_p);
        mpz_clear(cloned_my_q);
        mpz_clear(cloned_my_d);
        mpz_clear(cloned_my_dp);
        mpz_clear(cloned_my_dq);
        
    }

    finish:
        return;
    
}


int main(int argc, char *argv[]) {

    char valid_slice1[32][5];  // Array to store valid combinations for slice 1
    int valid_count = 0;       // Number of valid combinations for slice 1
    bool verbose, quit = false;

    char possibilities[32][5] = {
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 1}, {0, 0, 0, 1, 0}, {0, 0, 0, 1, 1},
    {0, 0, 1, 0, 0}, {0, 0, 1, 0, 1}, {0, 0, 1, 1, 0}, {0, 0, 1, 1, 1},
    {0, 1, 0, 0, 0}, {0, 1, 0, 0, 1}, {0, 1, 0, 1, 0}, {0, 1, 0, 1, 1},
    {0, 1, 1, 0, 0}, {0, 1, 1, 0, 1}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 1},
    {1, 0, 0, 0, 0}, {1, 0, 0, 0, 1}, {1, 0, 0, 1, 0}, {1, 0, 0, 1, 1},
    {1, 0, 1, 0, 0}, {1, 0, 1, 0, 1}, {1, 0, 1, 1, 0}, {1, 0, 1, 1, 1},
    {1, 1, 0, 0, 0}, {1, 1, 0, 0, 1}, {1, 1, 0, 1, 0}, {1, 1, 0, 1, 1},
    {1, 1, 1, 0, 0}, {1, 1, 1, 0, 1}, {1, 1, 1, 1, 0}, {1, 1, 1, 1, 1}
};


    // Check command-line arguments for "-v"
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0) {
            verbose = true;
        }
    }
    
    // Initialize GMP variables 
    mpz_t e, k, kp, kq, N, p, q, d, dp, dq, my_p, my_q, my_d, my_dp, my_dq,result_p, result_q;
    mpz_inits(N, k, kp, kq, p, q, d, dp, dq, result_p, result_q, NULL);

    mpz_init2(my_p, 512);   mpz_set_ui(my_p, 0);
    mpz_init2(my_q, 512);   mpz_set_ui(my_q, 0);
    mpz_init2(my_d, 1024);   mpz_set_ui(my_d, 0);
    mpz_init2(my_dp, 512);   mpz_set_ui(my_dp, 0);
    mpz_init2(my_dq, 512);   mpz_set_ui(my_dq, 0);

    // Set values for k, kp, kq
    mpz_set_ui(k, 40023); 
    mpz_set_ui(kp, 57300);  
    mpz_set_ui(kq, 48356); 

    // Set value for N (example value, replace with actual value)
    FILE *file = fopen("keys/rsakey.txt", "r");
    read_and_convert_component(N, file);
    //read_and_convert_component(e, file);

    // Set the value for e (common RSA public exponent)
    mpz_init_set_ui(e, 65537);


    int known_bits_p[1024], known_bits_q[1024], known_bits_d[2048], known_bits_dp[1024], known_bits_dq[1024];
    readKnownBits("keys/known_bits_p.txt", known_bits_p);
    readKnownBits("keys/known_bits_q.txt", known_bits_q);
    readKnownBits("keys/known_bits_d.txt", known_bits_d);
    readKnownBits("keys/known_bits_dp.txt", known_bits_dp);
    readKnownBits("keys/known_bits_dq.txt", known_bits_dq);

    // First phase of the RSA key reconstruction
    reconstruct_rsa_key_first_phase(my_p, my_q, my_d, my_dp, my_dq, e, k, kp, kq);
    // Compute τ values for k, kp, and kq
    int tau_k = tau(k);
    int tau_kp = tau(kp);
    int tau_kq = tau(kq);
    

    // Prepare slice0 array to hold initial bits
    char slice0[5];

    // getGmpBit(mpz_t number, int i)
    slice0[0] = getGmpBit(my_p, 0); 
    slice0[1] = getGmpBit(my_q, 0); 
    slice0[2] = getGmpBit(my_d, tau_k); 
    slice0[3] = getGmpBit(my_dp, tau_kp); 
    slice0[4] = getGmpBit(my_dq, tau_kq); 

    // Print initial bits
    gmp_printf("Correction: p[0]=%Zu, q[0]=%Zu, d[%d]=%Zu, dp[%d]=%Zu, dq[%d]=%Zu\n", 
           my_p, my_q,tau_k + 2 , my_d,tau_kp + 1, my_dp,tau_kq+1, my_dq);
    
    gmp_printf("Slice(0): p[0]=%d, q[0]=%d, d[%d]=%d, dp[%d]=%d, dq[%d]=%d\n", 
           slice0[0],slice0[1], tau_k , slice0[2],tau_kp, slice0[3],tau_kq, slice0[4]);
    
    // Start the core algorithm of Heninger and Shacham
    BranchAndPrune(result_p, result_q, my_p, my_q, my_d, my_dp, my_dq, e, k, kp, kq, N, tau_k, tau_kp, tau_kq, possibilities, 32, known_bits_p,known_bits_q,known_bits_d,known_bits_dp,known_bits_dq, verbose, 1);

    // Switching Kp and Kq for the second execution    
    BranchAndPrune(result_p, result_q, my_p, my_q, my_d, my_dp, my_dq, e, k, kq, kp, N, tau_k, tau_kq, tau_kp, possibilities, 32, known_bits_p,known_bits_q,known_bits_d,known_bits_dp,known_bits_dq, verbose, 1);

    // Clean up
    mpz_clears(e, k, kp, kq, N, p, q, d, dp, dq, result_p, result_q, NULL);

    return 0;
}
