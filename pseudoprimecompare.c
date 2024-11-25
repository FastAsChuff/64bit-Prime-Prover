#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include </home/simon/fastmodinvpow2fns.c> // See https://github.com/FastAsChuff/Fast-Modular-Inverse-Modulo-Powers-Of-2
#include </home/simon/modpowu64.c> // See https://github.com/FastAsChuff/Fast-Modular-Exponentiation/tree/main
#include </home/simon/isprime5fns.c> 

// gcc pseudoprimecompare.c -o pseudoprimecompare.bin -lm -O3 -march=native -Wall

void reportcase1(uint64_t n, char*filename) {
  printf("%lu not in %s\n", n, filename);
  return;
}

void reportcase2(uint64_t n, char*filename, primedivcheck_t primedivcheck) {
  uint64_t res = modpowu64(2, (n-1)/2, n);
  if ((res != 1) && (res != n-1)) return;
  _Bool divisorfound = false;   
  uint64_t pmax = isqrt(n);     
  for (uint32_t i=0; (i < primedivcheck.primedivchecksize); i++) {
    if (primedivcheck.p[i] > pmax) break;
    if (primedivcheck.a[i]*n <= primedivcheck.m[i]) {
      divisorfound = true;
      break;
    }
  }
  if (!divisorfound) {
    printf("%lu not in %s\n", n, filename);
  }
  return;
}

int main(int argc, char** argv) {
  if (argc < 4) {
    printf("This program compares pseudoprimes in an isprime5.c exceptions list for w = 2 in exceptionsfilename with integers in the list of all 64 bit fermat-pseudoprimes to base 2 downloaded from http://www.cecm.sfu.ca/Pseudoprimes/ computed by Jan Feitsma.\nNumbers found in one and not in the other, but should be, are reported.\nUsage %s primesfilename exceptionsfilename Feitsmafilename\n", argv[0]);
    exit(1);
  }
  isprime5exceptions_t ex;
  ex.array = NULL;
  char *primesfilename = NULL;
  char *exceptionsfilename = NULL;
  char *Feitsmafilename = NULL;
  if (argc >= 4) {
    primesfilename = argv[1];
    exceptionsfilename = argv[2];
    Feitsmafilename = argv[3];
  }
  primedivcheck_t primedivcheck;
  if (!loadexceptions(&ex, exceptionsfilename)) {
    if (argc >= 4) {
      printf("Failed to load %s\n", exceptionsfilename);
      exit(1);
    }
  } else {
    if (!loadprimes(ex, &primedivcheck, primesfilename)) {
      printf("Load of %s failed.\n", primesfilename);
      free(ex.array);
      exit(1);
    }
  }
  // Case 1: N in exceptionsfilename, and not in Feitsmafilename
  //   Report N not in Feitsmafilename
  // Case 2: N not in exceptionsfilename, and is in Feitsmafilename
  //   If N has no prime factor < 2^L, 
  //     and 2^(N-1)/2 != +- 1 mod N,
  //     Report N not in exceptionsfilename
  // Case 3: N in exceptionsfilename, and is in Feitsmafilename
  //   Do nothing
  uint32_t exix = 0;
  uint64_t prevnum, num, twoto2L;
  FILE *fp = fopen(Feitsmafilename, "r");
  if (fp == NULL) {
    printf("Failed to open %s\n", Feitsmafilename);
    primedivcheckfree(&primedivcheck);
    free(ex.array);
    exit(1);
  }
  twoto2L = 1ULL << (2*ex.L);
  prevnum = 0;
  while (true) {
    if (fscanf(fp, "%lu", &num) != 1) break;
    if (num <= prevnum) {
      printf("%s not in strictly ascending order at %lu\n", Feitsmafilename, num);
      primedivcheckfree(&primedivcheck);
      free(ex.array);
      exit(1);
    }
    if (exix < ex.count) {
      while (ex.array[exix] < num) reportcase1(ex.array[exix++], Feitsmafilename);
      if (ex.array[exix] == num) {
        exix++;
      } else {
        if (num > twoto2L) reportcase2(num, exceptionsfilename, primedivcheck);
      }
    } else {
      if (num > twoto2L) reportcase2(num, exceptionsfilename, primedivcheck);
    }
    prevnum = num;
  }
  while (exix < ex.count) reportcase1(ex.array[exix++], Feitsmafilename);
}
