#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include </home/simon/fastmodinvpow2fns.c>

// gcc isprime4.c -o isprime4.bin -lm -O3 -march=native -Wall
#define NUMOF32BITPRIMES 203280221
#define SIZEOFDIVCHECKARRAYP (NUMOF32BITPRIMES-1)*sizeof(uint32_t)
#define SIZEOFDIVCHECKARRAY (NUMOF32BITPRIMES-1)*sizeof(uint64_t)

typedef struct {
  uint32_t *p;
  uint64_t *a, *m;
} primedivcheck_t;

char* primesfilename = "primes.txt";
char* primedivcheckfilename = "primedivcheck4.dat";
// sha256(primedivcheckfilename) = a419a7db179f4190b4c2380a900009b119ddef23ccfb0c90cdf9811d0c48d910

uint64_t isprime3_gettimems(void) {
  struct timeval thistime;
  gettimeofday(&thistime, NULL);
  return thistime.tv_sec*1000ULL + thistime.tv_usec/1000ULL;
}

void primedivcheckfree(primedivcheck_t *primedivcheck) {
  free(primedivcheck->p);
  free(primedivcheck->a);
  free(primedivcheck->m);
} 

_Bool loadprimedivcheck(primedivcheck_t *primedivcheck, char *filename) {
  FILE *fp = fopen(filename, "rb");
  if (fp == NULL) return false;
  primedivcheck->p = malloc(SIZEOFDIVCHECKARRAYP);
  if (primedivcheck->p == NULL) {
    fclose(fp);
    return false;
  }
  size_t bytes = fread(primedivcheck->p, 1, SIZEOFDIVCHECKARRAYP, fp);
  if (bytes != SIZEOFDIVCHECKARRAYP) {
    free(primedivcheck->p);
    fclose(fp);
    return false;
  }
  primedivcheck->a = malloc(SIZEOFDIVCHECKARRAY);
  if (primedivcheck->a == NULL) {
    free(primedivcheck->p);
    fclose(fp);
    return false;
  }
  bytes = fread(primedivcheck->a, 1, SIZEOFDIVCHECKARRAY, fp);
  if (bytes != SIZEOFDIVCHECKARRAY) {
    free(primedivcheck->p);
    free(primedivcheck->a);
    fclose(fp);
    return false;
  }
  primedivcheck->m = malloc(SIZEOFDIVCHECKARRAY);
  if (primedivcheck->m == NULL) {
    free(primedivcheck->p);
    free(primedivcheck->a);
    fclose(fp);
    return false;
  }
  bytes = fread(primedivcheck->m, 1, SIZEOFDIVCHECKARRAY, fp);
  if (bytes != SIZEOFDIVCHECKARRAY) {
    free(primedivcheck->p);
    free(primedivcheck->a);
    free(primedivcheck->m);
    fclose(fp);
    return false;
  }
  fclose(fp);
  return true;
}

_Bool saveprimedivcheck(primedivcheck_t primedivcheck, char *filename) {
  FILE *fp = fopen(filename, "wb");
  if (fp == NULL) return false;
  size_t bytes = fwrite(primedivcheck.p, 1, SIZEOFDIVCHECKARRAYP, fp);
  if (bytes != SIZEOFDIVCHECKARRAYP) {
    fclose(fp);
    return false;
  }
  bytes = fwrite(primedivcheck.a, 1, SIZEOFDIVCHECKARRAY, fp);
  if (bytes != SIZEOFDIVCHECKARRAY) {
    fclose(fp);
    return false;
  }
  bytes = fwrite(primedivcheck.m, 1, SIZEOFDIVCHECKARRAY, fp);
  if (bytes != SIZEOFDIVCHECKARRAY) {
    fclose(fp);
    return false;
  }
  fclose(fp);
  return true;
}

_Bool loadprimedivcheck0(primedivcheck_t *primedivcheck, char *filename) {
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) return false;
  uint64_t primesum = 2;
  uint32_t prime;
  if (fscanf(fp, "%u", &prime) != 1) {
    fclose(fp);
    return false;
  }
  if (prime != 2) {
    fclose(fp);
    return false;
  }
  primedivcheck->p = malloc(SIZEOFDIVCHECKARRAYP);
  if (primedivcheck->p == NULL) {
    fclose(fp);
    return false;
  }
  primedivcheck->a = malloc(SIZEOFDIVCHECKARRAY);
  if (primedivcheck->a == NULL) {
    free(primedivcheck->p);
    fclose(fp);
    return false;
  }
  primedivcheck->m = malloc(SIZEOFDIVCHECKARRAY);
  if (primedivcheck->m == NULL) {
    free(primedivcheck->p);
    free(primedivcheck->a);
    fclose(fp);
    return false;
  }
  uint32_t i = 0;
  while (fscanf(fp, "%u", &prime) == 1) {
    primedivcheck->p[i] = prime;
    primedivcheck->a[i] = modinv64x(prime);
    primedivcheck->m[i] = 0xffffffffffffffffULL / prime;
    primesum += prime;
    i++;
  }
  fclose(fp);
  //printf("prime sum is %lu\n", primesum);
  if (primesum != 425649736193687430ULL) {
    primedivcheckfree(primedivcheck);
    return false;
  }
  return true;
}

uint32_t isqrt(uint64_t n) {
  if (n < 2) return n;
  uint64_t ai = sqrt(n);
  while (!((ai <= n/ai) && ((ai+1) > n/(ai+1)))) {    
    ai = (ai + n/ai)/2;
  }
  return ai;
}

_Bool isprime64dc(uint64_t n, primedivcheck_t primedivcheck) {
  // Deterministic primality test using trial division by primes.
  if (n<=1) return false;
  if (n <= 3) return true;
  if ((n & 1) == 0) return false;
  uint64_t pmax = isqrt(n);
  if (n == pmax*pmax) return false; 
  for (uint32_t i = 0; i < NUMOF32BITPRIMES-1; i++) {
    if (primedivcheck.p[i] > pmax) break;
    if (primedivcheck.a[i]*n <= primedivcheck.m[i]) return false;
  }
  return true;
}

_Bool nextprime64dc(uint64_t n, uint64_t *nextprime, primedivcheck_t primedivcheck) {
  while (n < 0xffffffffffffffffULL) {
    n++;
    if (isprime64dc(n, primedivcheck)) {
      *nextprime = n;
      return true;
    }
  }
  return false;
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
  if (argc < 2) {
    printf("Test if argument is definitely prime using naive trial division. For large arguments, the program may use trial division by primes if %s is available. This file is generated from %s (all 203280221 32bit primes) if present.\nMemory/Storage Requirement: 4GB (if files used), <1MB otherwise.\nUsage:- %s number\n0 <= number < 2^64\n", primedivcheckfilename, primesfilename, argv[0]);
    exit(1);
  }
  primedivcheck_t primedivcheck;
  primedivcheck.p = NULL;
  uint64_t n = atou64(argv[1]);
  if (n > 1000000000000000ULL) {
    if (!loadprimedivcheck(&primedivcheck, primedivcheckfilename)) {
      printf("Load of %s failed.\n", primedivcheckfilename);
      // See https://github.com/FastAsChuff/Primes-List
      // Create primes.txt = all 32bit primes = output of primeslist.bin 203280221
      if (loadprimedivcheck0(&primedivcheck, primesfilename)) {
        if (!saveprimedivcheck(primedivcheck, primedivcheckfilename)) {
          printf("Save of %s failed.\n", primedivcheckfilename);          
        }
      } else {
        printf("Load of %s failed.\n", primesfilename);  
      }
    }
  }
  uint64_t next;
  _Bool nisprime;
  uint64_t starttime, endtime;
  starttime = isprime3_gettimems();
  if (primedivcheck.p) {
    nisprime = isprime64dc(n, primedivcheck);
  } else {
    nisprime = isprime64(n);
  }
  endtime = isprime3_gettimems();
  if (nisprime) {
    printf("%lu is proven prime.\n", n);
  } else {
    printf("%lu is not prime.\n", n);
  }
  printf("Time of test = %lu ms\n", endtime - starttime);
  if (primedivcheck.p) {
    if (nextprime64dc(n, &next, primedivcheck)) printf("Next proven prime is %lu\n", next);
  } else {
    if (nextprime64(n, &next)) printf("Next proven prime is %lu\n", next);
  }
  if (primedivcheck.p) primedivcheckfree(&primedivcheck);
}
// largest 64 bit semi-prime = 4294967279*4294967291 = 18446743979220271189
