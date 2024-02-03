#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <time.h>

void generate_distinct_primes(mpz_t p, mpz_t q, unsigned int bit_size) {
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL)); // Seed with current time

    do {
        mpz_urandomb(p, state, bit_size); // Generate a random number and find next prime
        mpz_nextprime(p, p);
    } while (mpz_probab_prime_p(p, 25) == 0);

    do {
        mpz_urandomb(q, state, bit_size); // Generate a random number and find next prime
        mpz_nextprime(q, q);
    } while ((mpz_probab_prime_p(q, 25) == 0) || (mpz_cmp(p, q) == 0)); // Ensure q is distinct from p

    gmp_randclear(state);
}

void save_to_bit_file(FILE *file, mpz_t num) {
    size_t size = mpz_sizeinbase(num, 2) + 2; // +2 for null-terminator and potential '-' sign
    char *str = malloc(size);
    if (str == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    mpz_get_str(str, 2, num); // Get number in binary format
    fprintf(file, "%s\n", str);
    free(str);
}

int main(int argc, char *argv[]) {
    mpz_t p, q, n, e, d, phi, dp, dq, qp;
    mpz_inits(p, q, n, e, d, phi, dp, dq, qp, NULL);

    // Check if at least one additional argument is provided
    if (argc < 2) {
        printf("Usage: %s <number of bits>\n", argv[0]);
        return 1;
    }

    // Convert the first argument from string to integer
    int bits = atoi(argv[1]);

    // Check if the conversion was successful
    if (bits <= 0) {
        printf("Please enter a valid positive integer for the number of bits.\n");
        return 1;
    }

    // Generate two distinct 512-bit primes p and q
    generate_distinct_primes(p, q, bits);
    
    // Compute n = p * q
    mpz_mul(n, p, q);

    // Choose e, typically 65537
    mpz_set_ui(e, 65537);

    // Compute phi(n) = (p-1)(q-1)
    mpz_t p_1, q_1;
    mpz_inits(p_1, q_1, NULL);
    mpz_sub_ui(p_1, p, 1);
    mpz_sub_ui(q_1, q, 1);
    mpz_mul(phi, p_1, q_1);

    // Compute d = e^-1 mod phi
    if (mpz_invert(d, e, phi) == 0) {
        fprintf(stderr, "Failed to compute d. e and phi are not coprime.\n");
        return 1;
    }

    // Compute dp = d mod (p-1) and dq = d mod (q-1)
    mpz_mod(dp, d, p_1);
    mpz_mod(dq, d, q_1);

    // Compute the inverse of q modulo p
    if (mpz_invert(qp, q, p) == 0) {
        fprintf(stderr, "Failed to compute inverse of q mod p.\n");
        return 1;
    }

    // Clean up
    mpz_clears(p_1, q_1, NULL);

    // Open files for saving in decimal, hexadecimal, and bit formats
    FILE *file = fopen("rsakey.txt", "w");
    FILE *hexFile = fopen("hexaRSAkey.txt", "w");
    FILE *bitsFile = fopen("rsabitskey.txt", "w");

    if (!file || !hexFile || !bitsFile) {
        fprintf(stderr, "Error opening files\n");
        return 1;
    }

    // Save RSA components in decimal format
    fprintf(file, "N=");
    mpz_out_str(file, 10, n);
    fprintf(file, "\ne=");
    mpz_out_str(file, 10, e);
    fprintf(file, "\nd=");
    mpz_out_str(file, 10, d);
    fprintf(file, "\np=");
    mpz_out_str(file, 10, p);
    fprintf(file, "\nq=");
    mpz_out_str(file, 10, q);
    fprintf(file, "\ndp=");
    mpz_out_str(file, 10, dp);
    fprintf(file, "\ndq=");
    mpz_out_str(file, 10, dq);
    fprintf(file, "\nqp=");
    mpz_out_str(file, 10, qp);

    // Save RSA components in hexadecimal format
    fprintf(hexFile, "N=");
    mpz_out_str(hexFile, 16, n);
    fprintf(hexFile, "\ne=");
    mpz_out_str(hexFile, 16, e);
    fprintf(hexFile, "\nd=");
    mpz_out_str(hexFile, 16, d);
    fprintf(hexFile, "\np=");
    mpz_out_str(hexFile, 16, p);
    fprintf(hexFile, "\nq=");
    mpz_out_str(hexFile, 16, q);
    fprintf(hexFile, "\ndp=");
    mpz_out_str(hexFile, 16, dp);
    fprintf(hexFile, "\ndq=");
    mpz_out_str(hexFile, 16, dq);
    fprintf(hexFile, "\nqp=");
    mpz_out_str(hexFile, 16, qp);

    // Save RSA components in bit format
    fprintf(bitsFile, "N=");
    save_to_bit_file(bitsFile, n);
    fprintf(bitsFile, "e=");
    save_to_bit_file(bitsFile, e);
    fprintf(bitsFile, "d=");
    save_to_bit_file(bitsFile, d);
    fprintf(bitsFile, "p=");
    save_to_bit_file(bitsFile, p);
    fprintf(bitsFile, "q=");
    save_to_bit_file(bitsFile, q);
    fprintf(bitsFile, "dp=");
    save_to_bit_file(bitsFile, dp);
    fprintf(bitsFile, "dq=");
    save_to_bit_file(bitsFile, dq);
    fprintf(bitsFile, "qp=");
    save_to_bit_file(bitsFile, qp);

    // Close all files
    fclose(file);
    fclose(hexFile);
    fclose(bitsFile);

    // Clean up GMP variables
    mpz_clears(p, q, n, e, d, phi, dp, dq, qp, NULL);

    return 0;
}
