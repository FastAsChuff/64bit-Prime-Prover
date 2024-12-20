# 64bit-Prime-Prover
Trial division based deterministic primality test of 64 bit integers.

Test if argument is definitely prime using naive trial division. For large arguments, the program may use trial division by primes if a precomputed file of constants is available. This file is generated from primes.txt (all 203280221 32bit primes) if present. See https://github.com/FastAsChuff/Primes-List to create primes.txt = all 32bit primes = output of primeslist.bin 203280221

The main motive for writing these programs was out of curiosity about just how fast, or slow, trial division can be on 64 bit inputs. It was not intended to be practical, but may be useful to anyone who doesn't want to rely on deep number theory, opaque library functions, or claims of exhaustive testing which are extrememly difficult to verify.

isprime2.bin
------------

Memory/Storage Requirement: 250MB (if files used), <1MB otherwise.

Usage:- ./isprime2.bin number

0 <= number < 2^64

e.g. ./isprime2.bin 18361375334787048247

18361375334787048247 is proven prime.

Next proven prime is 18361375334787048287

isprime4.bin
------------

Memory/Storage Requirement: 4GB (if files used), <1MB otherwise.

Usage:- ./isprime4.bin number

0 <= number < 2^64

e.g. ./isprime4.bin 18446744073680551637

18446744073680551637 is proven prime.

Time of test = 625 ms

Next proven prime is 18446744073680551639

===============================================================================================================

The program isprime5.c is a deterministic primality test on 64bit unsigned integers too, but uses some trial division up to a point and then uses a look-up table of exceptions. The table is generated by primetestmakeWLfile*.c programs. The *Feitsma.c one uses a file of all base 2 Fermat pseudo-primes generated by Jan Feitsma which was downloaded from http://www.cecm.sfu.ca/Pseudoprimes/. The other can be generated without relying on external sources, except a list of primes, but takes some time. The file for w = 2, L = 26 took 85 minutes on a 6th gen Core i7 desktop with 8 threads running. The runtime approximately doubles on each decrement of L. To clarify, an odd composite 64 bit integer N > 2 is in a wLfile if it has no prime factors below 2^L and w^(N-1)/2 = +- 1 mod N. For testing, the exceptions files with w = 2 can be compared with the Feitsma pseudoprimes using pseudoprimecompare.c. The isprime5 timings presented below were done on an i7-6700 @3.4GHz.

isprime5.bin
------------

./isprime5.bin 18176543219876544409 primes.txt primetestwLfile3_2_20.txt 

18176543219876544409 is proven prime.

Time of test = 65 us

Next proven prime is 18176543219876544421

Time to find nextprime = 115 us


./isprime5.bin 18176543219876544409 primes.txt primetestwLfileF_2_8.txt 

18176543219876544409 is proven prime.

Time of test = 3 us

Next proven prime is 18176543219876544421

Time to find nextprime = 5 us

Exception Files
---------------

	primetestwLfile3_2_24.txt			https://drive.google.com/file/d/11JiVQsPSFxvAriK07apgrfXT2oR1UbdA/view?usp=sharing
	primetestwLfile3_2_22.txt			https://drive.google.com/file/d/1zGFJtTVEbSM4ZqGhsqz1xfe6k7h24GfX/view?usp=sharing
	primetestwLfile3_2_20.txt			https://drive.google.com/file/d/13LR5ji0WvaoEcOz3SrmkXdryyxqQV6wR/view?usp=sharing
	primetestwLfileF_2_8.txt			https://drive.google.com/file/d/1pQK0o8oZQfP_MrgY5p80DxG4P8rCy8Rl/view?usp=sharing


