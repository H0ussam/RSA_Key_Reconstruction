We aimed to implement the Heninger-shacham algorithm [Reconstructing RSA Private Keys from Random Key Bits](https://link.springer.com/chapter/10.1007/978-3-642-03356-8_1), it is partially based on the premise that an attacker can exploit a phyicaly property of memories "type DRAM" by conducting a [Cold Boot Attack](https://citp.princeton.edu/our-work/memory/) and preforming a dump of the memory.

It turned out that DRAMs don't earease bits stored within it immediatly after shutdown, instead the capacitors leak charges over time until reaching thier ground states, this open a wide window for exploitation.

you find mainly 2 parts of the whole algorithm:

the file named Heninger_part1.c, deals with finding the value of K, correcting the most significant half of d, finally computing kp & kq, we used the [Tonelli-Shanks algorithm](https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm), to compute the roots of the equation.
For exhaustive details please read the mentioned paper above.

The file named  Heninger_part2.c deals with the branch and prune part, it verifies all possible combination of bits for the coponents (d, p, q, d_p, d_q) at a given position "bit i", if they verifiy certain conditions mentioned in the paper, they are added to in the position i otherwise pruned; again for exhaustive details check the paper x)

The keys degraded should be in a certain format for Heninger_part2.c to work properly.

Compile the file generate_key.c:

`$gcc -o outputName generate_key.c -lgmp`

to run it:

`./outputName 512` or `./outputName 1024`, replace the inetegrs with whatever size you wish to generate.

you ll get multiple file, we will let you exploare.

compile and run the file degradation.c

`gcc -o outputName degradation.c -lgmp` && `./outputName`

this will preform a degradation on the keys with a given probability that you should set manually within the file itself, it is a public variable named `KNOWN_BIT_PERCENTAGE`.

As we said before, it will yeild the degraded values of the key in a specific form, but you shall need the degraded version of d in it's decimal form to run the first phase of the algorithm Heninger_part2.c.

To tackle the above and get the value of d, compile and run the file named print_value_degraded.c.

`gcc -o outputName print_value_degraded.c -lgmp` && `./outputName`

At this point you are all set with your values, save your degraded version of d and off to Heninger_part1.c, where you should set manually the values of your N, and degraded version of d.
In the same way compile & run your prrogram then analyse the results.

note:`We implemented the first phase of the algorithm on a more real life scenario, where the attacker will have no idea which bits are correct and which ones have been flipped, thus we have found by experiment, any key with error probability superior than 45% the program Heninger_part1.c has a very low chance of yeilding the correct values, we advise you to degrade the keys with a fixed probability between 15% and 45%, until later when we add the exact implementation of Heninsger-Shacham where you know the position of the correct bits exactly`.

Finally from the output of Heninger_part1.c, get the values of k, kp, kq, and set them manually; Compile then run in the same way.

note`If you get inccorect results, you might want to interchange the values of kp and kq, because we are not sure which one is kp and which one is kq, it is only a matter of testing, and for further details please check the paper.`

After commpiling, you can run your program in 2 ways:

the normal way -----> `./outputName`; Will yeild directly the values of p and q thus you have factorized N and solved the mathematical problem linked to the security of RSA.
Verbosity ----------> `./outputName -v`, which will show you the process of branching and pruning up until the very last value.

Thank you for reading.













