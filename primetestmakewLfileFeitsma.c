#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include </home/simon/fastmodinvpow2fns.c> // See https://github.com/FastAsChuff/Fast-Modular-Inverse-Modulo-Powers-Of-2
#include </home/simon/isprime5fns.c> 

// gcc primetestmakewLfileFeitsma.c -o primetestmakewLfileFeitsma.bin -lm -O3 -march=native -Wall

#define NUMOFFEITSMANUMS 118968378

int main(int argc, char** argv) {
  if (argc < 5) {
    printf("This program generates an isprime5.c exception file using Jan Feitsma's file of all base 2 Fermat pseudo-primes.\nUsage %s L primesfilename exceptionsfilename Feitsmafilename\n", argv[0]);
    exit(1);
  }
  isprime5exceptions_t ex;
  ex.array = NULL;
  char *primesfilename = NULL;
  char *exceptionsfilename = NULL;
  char *Feitsmafilename = NULL;
  uint64_t L;
  if (argc >= 4) {
    L = atoi(argv[1]);
    primesfilename = argv[2];
    exceptionsfilename = argv[3];
    Feitsmafilename = argv[4];
  }
  if ((L < 2) || (L > 31)) {
    printf("2 < L < 32\n");
    exit(1);
  }
  primedivcheck_t primedivcheck;
  ex.array = malloc(NUMOFFEITSMANUMS*sizeof(uint64_t));
  if (ex.array == NULL) {
    printf("Failed to allocate enough memory.\n");
    exit(1);
  }
  ex.count = 0;
  ex.sum = 0;
  ex.w = 2;
  ex.L = L;
  if (!loadprimes(ex, &primedivcheck, primesfilename)) {
    printf("Load of %s failed.\n", primesfilename);
    exit(1);
  }
  FILE *fp = fopen(Feitsmafilename, "r");
  if (fp == NULL) {
    printf("Failed to open %s\n", Feitsmafilename);
    primedivcheckfree(&primedivcheck);
    free(ex.array);
    exit(1);
  }
  FILE *fpex = fopen(exceptionsfilename, "w");
  if (fpex == NULL) {
    printf("Failed to open %s\n", exceptionsfilename);
    primedivcheckfree(&primedivcheck);
    free(ex.array);
    fclose(fp);
    exit(1);
  }
  uint64_t n, twoto2L;
  twoto2L = 1ULL << (2*ex.L);
  while (true) {
    if (fscanf(fp, "%lu", &n) != 1) break;
    if ((n & 1ULL) && (n > twoto2L)) {
      uint64_t res = modpowu64(2, (n-1)/2, n);
      if ((res == 1) || (res == n-1)) {
        _Bool divisorfound = false;   
        for (uint32_t i=0; (i < primedivcheck.primedivchecksize); i++) {
          if (primedivcheck.a[i]*n <= primedivcheck.m[i]) {
            divisorfound = true;
            break;
          }
        }
        if (!divisorfound) {
          ex.array[ex.count] = n;
          ex.count++;
          ex.sum += n;
        }
      }
    }
  }
  fprintf(fpex, "2 %lu %lu %lu\n", ex.L, ex.count, ex.sum);
  for (uint32_t i=0; i<ex.count; i++) {
    fprintf(fpex, "%lu\n", ex.array[i]);
  }
  fclose(fp);
  fclose(fpex);
  primedivcheckfree(&primedivcheck);
  free(ex.array);
  printf("OK\n");
}
