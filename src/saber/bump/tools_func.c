#include <unistd.h>
#include <stdio.h>
#include <stdint.h>

uint32_t fletcher32( uint32_t *n, uint16_t const *var )
{
   uint32_t sum1 = 0xffff, sum2 = 0xffff;
   size_t tlen; uint32_t words = *n ;
   while (words) {
      tlen = ((words >= 359) ? 359 : words);
      words -= tlen;
      do {
         sum2 += sum1 += *var++;
         tlen--;
      } while (tlen);
      sum1 = (sum1 & 0xffff) + (sum1 >> 16);
      sum2 = (sum2 & 0xffff) + (sum2 >> 16);
   }
   sum1 = (sum1 & 0xffff) + (sum1 >> 16);
   sum2 = (sum2 & 0xffff) + (sum2 >> 16);
   return (sum2 << 16) | sum1;
}
