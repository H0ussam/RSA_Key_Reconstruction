#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <string.h>

#define KNOWN_BIT_PERCENTAGE 80

typedef struct {
    int position;
    int value;
} KnownBit;



void read_and_convert_component(mpz_t component, FILE *file) {
    char buffer[2048]; // Adjust buffer size if necessary
    if (fgets(buffer, sizeof(buffer), file) != NULL) {
        // Skip the identifier (e.g., "N=") and read the number in hex
        mpz_set_str(component, strchr(buffer, '=') + 1, 16);
    }
}

// The is fucntion write a given component written in the format of known and unknow bits to the actual integer value
void write_component() {

}


void degrade_component(mpz_t component, const char *filename) {
    unsigned int total_bits = mpz_sizeinbase(component, 2);
    unsigned int known_bit_count = (total_bits * KNOWN_BIT_PERCENTAGE) / 100;
    KnownBit all_bits[total_bits];

    // Initialize the array with all bits set to unknown (-1)
    for (unsigned int i = 0; i < total_bits; ++i) {
        all_bits[i].position = i;
        all_bits[i].value = -1; // Initially mark all as unknown
    }

    // Generate an array of indices and shuffle it
    unsigned int indices[total_bits];
    for (unsigned int i = 0; i < total_bits; ++i) {
        indices[i] = i; // Fill the array with all possible indices
    }

    // Randomly shuffle the indices array
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    for (unsigned int i = 0; i < total_bits; ++i) {
        unsigned int j = gmp_urandomm_ui(state, total_bits);
        // Swap elements i and j
        unsigned int temp = indices[i];
        indices[i] = indices[j];
        indices[j] = temp;
    }

    // Now select the first known_bit_count bits as known
    for (unsigned int i = 0; i < known_bit_count; ++i) {
        all_bits[indices[i]].value = mpz_tstbit(component, indices[i]);
    }

    gmp_randclear(state);

    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        exit(1);
    }

    // Write the bit status (known or unknown) to the file
    for (unsigned int i = 0; i < total_bits; ++i) {
        fprintf(file, "%d %d\n", all_bits[i].position, all_bits[i].value);
    }

    fclose(file);
}

int main() {
    mpz_t n, e, d, p, q, dp, dq, qp;
    mpz_inits(n, e, d, p, q, dp, dq, qp, NULL);

    FILE *rsakey = fopen("hexaRSAkey.txt", "r");
    if (!rsakey) {
        perror("Error opening hexaRSAkey.txt");
        return 1;
    }

    // Read and convert components from hexaRSAkey.txt
    read_and_convert_component(n, rsakey);
    read_and_convert_component(e, rsakey);
    read_and_convert_component(d, rsakey);
    read_and_convert_component(p, rsakey);
    read_and_convert_component(q, rsakey);
    read_and_convert_component(dp, rsakey);
    read_and_convert_component(dq, rsakey);
    read_and_convert_component(qp, rsakey);
    fclose(rsakey);

    // Degrade each component separately
    degrade_component(n, "known_bits_n.txt");
    degrade_component(e, "known_bits_e.txt");
    degrade_component(d, "known_bits_d.txt");
    degrade_component(p, "known_bits_p.txt");
    degrade_component(q, "known_bits_q.txt");
    degrade_component(dp, "known_bits_dp.txt");
    degrade_component(dq, "known_bits_dq.txt");
    degrade_component(qp, "known_bits_qp.txt");

    mpz_clears(n, e, d, p, q, dp, dq, qp, NULL);
    return 0;
}
