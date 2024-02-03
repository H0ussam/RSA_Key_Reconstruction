#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

void flip_bits_and_compute_d(const char *orig_filename, const char *mod_filename) {
    FILE *orig_file = fopen(orig_filename, "r");
    if (!orig_file) {
        perror("Error opening original file");
        exit(1);
    }

    char line[1024];
    mpz_t d;
    mpz_init(d);
    int found_d = 0;

    // Search for 'd=' in the original file and extract its value
    while (fgets(line, sizeof(line), orig_file)) {
        if (strncmp(line, "d=", 2) == 0) {
            char *d_str = line + 2; // Skip the "d=" part to get the number
            mpz_set_str(d, d_str, 10); // Assume 'd' is in base 10
            found_d = 1;
            break;
        }
    }
    fclose(orig_file);

    if (!found_d) {
        fprintf(stderr, "Failed to find 'd' in the original file.\n");
        exit(1);
    }

    // Open the modification file to read the bits to be flipped
    FILE *mod_file = fopen(mod_filename, "r");
    if (!mod_file) {
        perror("Error opening modification file");
        exit(1);
    }

    unsigned int position;
    int value;
    while (fscanf(mod_file, "%u %d\n", &position, &value) == 2) {
        if (value == -1) { // If the bit needs to be flipped
            if (mpz_tstbit(d, position) == 1) {
                mpz_clrbit(d, position); // Flip from 1 to 0
            } else {
                mpz_setbit(d, position); // Flip from 0 to 1
            }
        }
    }
    fclose(mod_file);

    // Print the final value of 'd'
    gmp_printf("Final value of d: %Zd\n", d);

    mpz_clear(d);
}


int main() {

    flip_bits_and_compute_d("rsakey.txt", "known_bits_d.txt");

    return 0;
}