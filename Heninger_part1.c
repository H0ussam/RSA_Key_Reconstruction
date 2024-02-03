#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>


void calculate_k(mpz_t N, mpz_t e, mpz_t degraded_d, mpz_t final_k, mpz_t final_d) {
    int max_matching_bits = 0;
    int matching_bits;
    mpz_t k, d_tilde, temp, N_plus_1, new_d_tilde, new_k;
    mpz_inits(k, d_tilde, temp, N_plus_1, new_d_tilde, new_k, NULL);


    size_t total_bits, bits_to_compare = 0;

    mpz_add_ui(N_plus_1, N, 1); // N_plus_1 = N + 1

    // Determine the number of bits in degraded_d
    total_bits = mpz_sizeinbase(degraded_d, 2);
    bits_to_compare = (total_bits / 2) - 2; // Compare down to (n/2) - 2

    // Loop to find the possible values of k between "1 and e-1" and calculate d_tilde
    for (mpz_set_ui(k, 1); mpz_cmp(k, e) < 0; mpz_add_ui(k, k, 1)) {
        // d_tilde(k) = (k(N + 1) + 1) / e
        mpz_mul(temp, k, N_plus_1); // temp = k * (N + 1)
        mpz_add_ui(temp, temp, 1);  // temp = k * (N + 1) + 1
        mpz_fdiv_q(d_tilde, temp, e); // d_tilde = temp / e

        matching_bits = 0;
        for (size_t i = 0; i < bits_to_compare; ++i) {
            size_t bit_index = total_bits - i - 1; // Calculate the index of the bit to compare
            if (mpz_tstbit(d_tilde, bit_index) == mpz_tstbit(degraded_d, bit_index)) {
                matching_bits++; // Increment match count if the bits are the same
            }
        }

        if (matching_bits > max_matching_bits) {
            max_matching_bits = matching_bits;
            matching_bits = 0;
            mpz_set(new_k, k);
            mpz_set(new_d_tilde, d_tilde);
        }       
    }
    
    mpz_set(final_k, new_k);
    mpz_set(final_d, new_d_tilde);

    gmp_printf("Matching d_tilde for k=%Zd: %Zd\n", new_k, new_d_tilde);

    mpz_clears(k, d_tilde, temp, N_plus_1, new_d_tilde, new_k, NULL);
}

//Once we find the correct value of k we can correct
// The most signifcant half of d_degraded
void correct_msb_half_d(mpz_t A, mpz_t B, mpz_t C) {
    // Initialize temporary variables
    mpz_t mask;
    mpz_init(mask);
    mpz_init(C);

    // Calculate the total number of bits in A
    size_t total_bits = mpz_sizeinbase(A, 2);

    // Calculate the number of bits to take from A (MSBs)
    size_t msb_bits = (total_bits / 2) - 2;

    // Create a bitmask for the MSB part of A
    mpz_set_ui(mask, 1);              // Set mask to 1
    mpz_mul_2exp(mask, mask, msb_bits); // Shift mask to create the MSB part (leave space for LSBs)
    mpz_sub_ui(mask, mask, 1);        // Create mask with msb_bits 1's
    mpz_mul_2exp(mask, mask, total_bits - msb_bits); // Align the mask with the MSBs of A

    // Isolate MSBs from A and store in C
    mpz_and(C, A, mask);

    // Create a bitmask for the LSB part of B
    mpz_set_ui(mask, 1);                 // Reset mask to 1
    mpz_mul_2exp(mask, mask, total_bits - msb_bits); // Shift mask for LSBs
    mpz_sub_ui(mask, mask, 1);           // Create mask with LSBs

    // Isolate LSBs from B and store in temp
    mpz_and(mask, B, mask);

    // Combine LSBs from B into C
    mpz_ior(C, C, mask);

    // Clean up
    mpz_clear(mask);
}

void Tonelli_shanks(const mpz_t n, const mpz_t p, mpz_t firstSolution, mpz_t secondSolution) {
    mpz_t Q, S, z, M, c, t, R, temp, b, two, exponent;
    mpz_inits(Q, S, z, M, c, t, R, temp, b, two, exponent, NULL);

    // p - 1 = Q * 2^S
    mpz_sub_ui(Q, p, 1); // Q = p - 1
    mpz_set_ui(S, 0);
    mpz_set_ui(two, 2);

    // Factor out powers of 2 from Q
    while (mpz_even_p(Q)) {
        mpz_fdiv_q_2exp(Q, Q, 1); // Q /= 2
        mpz_add_ui(S, S, 1); // S += 1
    }

    // Find z, a quadratic non-residue modulo p
    mpz_set_ui(z, 2);
    while (mpz_legendre(z, p) != -1) {
        mpz_add_ui(z, z, 1);
    }

    mpz_powm(c, z, Q, p); // c = z^Q % p
    mpz_powm(t, n, Q, p); // t = n^Q % p
    mpz_add_ui(temp, Q, 1);
    mpz_fdiv_q_2exp(temp, temp, 1); // temp = (Q + 1) / 2
    mpz_powm(R, n, temp, p); // R = n^((Q + 1)/2) % p

    // Main loop
    while (mpz_cmp_ui(t, 0) != 0 && mpz_cmp_ui(t, 1) != 0) {
        mpz_set_ui(M, 0);
        mpz_set(temp, t);

        while (mpz_cmp_ui(temp, 1) != 0) {
            mpz_powm_ui(temp, temp, 2, p); // temp = temp^2 % p
            mpz_add_ui(M, M, 1);
        }

        if (mpz_cmp(M, S) == 0) {
            mpz_clears(Q, S, z, M, c, t, R, temp, b, two, exponent, NULL);
            return; // No solution
        }

        mpz_sub(exponent, S, M);
        mpz_sub_ui(exponent, exponent, 1);
        mpz_powm_ui(b, two, mpz_get_ui(exponent), p); // b = 2^(S-M-1)
        mpz_powm(b, c, b, p); // b = c^b % p

        mpz_powm_ui(c, b, 2, p); // c = b^2 % p
        mpz_mul(t, t, c);
        mpz_mod(t, t, p); // t = t*c % p
        mpz_mul(R, R, b);
        mpz_mod(R, R, p); // R = R*b % p
        mpz_set(S, M); // S = M
    }

    // Set the solutions
    if (mpz_cmp_ui(t, 1) == 0) {
        mpz_set(firstSolution, R);
        mpz_neg(secondSolution, R); // secondSolution = -R
        mpz_mod(secondSolution, secondSolution, p);
    } else {
        printf("No solutions found \n");
    }

    mpz_clears(Q, S, z, M, c, t, R, temp, b, two, exponent, NULL);
}

// (kp^2) − [k(N − 1) + 1] * kp − k ≡ 0 (mod e)
// [kp - ([[k(N − 1) + 1] + e]/2)] ^ 2 ≡ k + ([[k(N − 1) + 1] + e]/2) ^2

//([[k(N − 1) + 1] + e]/2) ^2
void solve_quad(const mpz_t N, const mpz_t e, const mpz_t k, mpz_t potential_kp, mpz_t potential_kq) {

    mpz_t A, B, tmpi, tmpo;
    mpz_inits(A, B, tmpi, tmpo, NULL);

    mpz_sub_ui(A, N, 1);
    mpz_mul(A, A, k);
    mpz_add_ui(A, A, 1);
    mpz_add(A, A, e);
    mpz_divexact_ui(A, A, 2); // This is the first term
    mpz_mod(A, A, e);

    mpz_mul(B, A, A);
    mpz_add(B, B, k);
    mpz_mod(B, B, e);

    Tonelli_shanks(B, e, tmpi, tmpo);

    mpz_add(potential_kp, tmpi, A);
    mpz_mod(potential_kp, potential_kp, e);
    mpz_add(potential_kq, tmpo, A);
    mpz_mod(potential_kq, potential_kq, e);

    mpz_clears(A, B, tmpi, tmpo, NULL);
}


int main() {
    mpz_t N, e, d_degraded, d_of_k, k, d_corrected_msb, kp, kq;
    
    // Initialize GMP variables
    mpz_init(N);
    mpz_init(e);
    mpz_init(d_degraded);
    mpz_init(d_of_k);
    mpz_init(k);
    mpz_init(d_corrected_msb);
    mpz_init(kp);
    mpz_init(kq);

    // Set the values of N and e (you need to set these values)
    mpz_set_str(N, "1779131295899651367704488804439300022159511204362688397204824497686024001544484872399088904033540065264683756731185268415775007794293495299602662239982195644719918039527453467728888708456223781573027854296707137645916161497343260997632507957391442835006592938087091954737741462030683760980944679650127689113", 10); 
    mpz_set_str(e, "65537", 10);
    mpz_set_str(d_degraded, "1093102732592072442165908446218389646200696268335976957206521089450096288445014241662663040189285875185365234182151140230717904288971619114224494660404708171128903115391894886409890093795141138372717908472774390489979047095066776569674359909342168983428691146160076553255412183602252810955637382610314935718", 10);

    
    calculate_k(N, e, d_degraded, k, d_of_k);

    gmp_printf("k=%Zd\n",k);

    correct_msb_half_d(d_of_k, d_degraded, d_corrected_msb);

    gmp_printf("d_corrected_msb=%Zd\n",d_corrected_msb);

    printf("Now calculating kp and kq, hold still...\n");

    solve_quad(N, e, k, kp, kq);

    gmp_printf("kp=%Zd kq=%Zd\n", kp, kq);

    mpz_clears(N, e, d_degraded, d_of_k, k, d_corrected_msb, kp, kq, NULL);

    return 0;
}