#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include </home/simon/fastmodinvpow2fns.c> // See https://github.com/FastAsChuff/Fast-Modular-Inverse-Modulo-Powers-Of-2
#include </home/simon/isprime5fns.c> 

// gcc isprime5.c -o isprime5.bin -lm -O3 -march=native -Wall


uint64_t isprime5_gettimeus(void) {
  struct timeval thistime;
  gettimeofday(&thistime, NULL);
  return thistime.tv_sec*1000000ULL + thistime.tv_usec;
}

_Bool isprime64(uint64_t n) {
  // Deterministic primality test using trial division.
  if (n<=1) return false;
  if (n <= 3) return true;
  if (n == 5) return true;
  if ((n & 1) == 0) return false;
  if ((n % 3) == 0) return false;
  if ((n % 5) == 0) return false;
  uint32_t imax = isqrt(n);
  // 1,7,11,13,17,19,23,29 mod 30 only
  for (uint64_t i = 7; ; i+=6) {
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 4;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 2;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 4;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 2;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 4;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 6;
    if (n%i == 0) return false;
    if (i >= imax) break;
    i += 2;
    if (n%i == 0) return false;
    if (i >= imax) break;
  }
  return true;
}

_Bool nextprime64(uint64_t n, uint64_t *nextprime) {
  while (n < 0xffffffffffffffffULL) {
    n++;
    if (isprime64(n)) {
      *nextprime = n;
      return true;
    }
  }
  return false;
}

uint64_t atou64(char *in) {
  uint64_t res = 0;
  while (*in) {
    res *= 10;
    res += *in - '0';
    in++;
  }
  return res;
}

int main(int argc, char** argv) {
  if ((argc < 2) || (argc == 3)) {
    printf("Test if argument N is definitely prime using a deterministic primality test utilising naive trial division to search for factors below 2^L, a test using 'wL' exceptions file to see if w^(N-1)/N = +-1 mod N or, if N is in list of composite exceptions.\nMemory/Storage Requirement: <5GB (if files used), <1MB otherwise.\nUsage:- %s number [primesfilename wLfilename]\n0 <= number < 2^64\n", argv[0]);
    exit(1);
  }
  isprime5exceptions_t ex;
  ex.array = NULL;
  char *primesfilename = NULL;
  char *exceptionsfilename = NULL;
  if (argc >= 4) {
    primesfilename = argv[2];
    exceptionsfilename = argv[3];
  }
  primedivcheck_t primedivcheck;
  uint64_t n = atou64(argv[1]);
  if (!loadexceptions(&ex, exceptionsfilename)) {
    if (argc >= 4) printf("Failed to load %s\n", exceptionsfilename);
  } else {
    if (!loadprimes(ex, &primedivcheck, primesfilename)) {
      printf("Load of %s failed.\n", primesfilename);
      free(ex.array);
      ex.array = NULL;
    }
  }
  uint64_t next;
  _Bool nisprime;
  uint64_t starttime, endtime;
  starttime = isprime5_gettimeus();
  if (primedivcheck.p) {
    nisprime = isprime64dce(n, primedivcheck, ex);
  } else {
    nisprime = isprime64(n);
  }
  endtime = isprime5_gettimeus();
  if (nisprime) {
    printf("%lu is proven prime.\n", n);
  } else {
    printf("%lu is not prime.\n", n);
  }
  printf("Time of test = %lu us\n", endtime - starttime);
  starttime = isprime5_gettimeus();
  if (primedivcheck.p) {
    if (nextprime64dce(n, &next, primedivcheck, ex)) {
      endtime = isprime5_gettimeus();
      printf("Next proven prime is %lu\n", next);
      printf("Time to find nextprime = %lu us\n", endtime - starttime);
    }
  } else {
    if (nextprime64(n, &next)) {
      endtime = isprime5_gettimeus();
      printf("Next proven prime is %lu\n", next);
      printf("Time to find nextprime = %lu us\n", endtime - starttime);
    }
  }
  if (primedivcheck.p) {
    primedivcheckfree(&primedivcheck);
    free(ex.array);
  }
}
