#include "interval.hpp"
#include <random>
int main(){std::mt19937_64 rng(150022);long checks=0;
 for(int precision:{64,256,384}){PREC=precision;for(int k=0;k<10000;k++){
  I a=I(long(rng()%20001)-10000)*I::pow2(long(rng()%401)-200);I b=I(long(rng()%20001)-10000)*I::pow2(long(rng()%401)-200);
  a=hull(a,I(long(rng()%20001)-10000)*I::pow2(long(rng()%401)-200));b=hull(b,I(long(rng()%20001)-10000)*I::pow2(long(rng()%401)-200));
  if(k%11==0)a=I(0);if(k%13==0)b=I(0);
  I old=mul_four_corners(a,b),fast=mul_by_sign(a,b);
  if(!mpfr_equal_p(old.l,fast.l)||!mpfr_equal_p(old.u,fast.u)){std::cerr<<"mismatch "<<precision<<" "<<k<<"\n";return 1;}checks++;
 }}std::cout<<"{\"status\":\"PASS\",\"seed\":150022,\"checks\":"<<checks<<",\"precisions\":[64,256,384],\"comparison\":\"exact MPFR endpoints; four-corner versus sign-selected multiplication\"}\n";
}
