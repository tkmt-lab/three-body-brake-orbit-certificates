#pragma once
#include <mpfr.h>
#include <gmp.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
inline mpfr_prec_t PREC=256;
struct I {
 mpfr_t l,u;
 I(){mpfr_init2(l,PREC);mpfr_init2(u,PREC);mpfr_set_zero(l,0);mpfr_set_zero(u,0);}
 I(long n):I(){mpfr_set_si(l,n,MPFR_RNDD);mpfr_set_si(u,n,MPFR_RNDU);}
 I(const I& a):I(){mpfr_set(l,a.l,MPFR_RNDD);mpfr_set(u,a.u,MPFR_RNDU);}
 I(I&& a) noexcept:I(){mpfr_swap(l,a.l);mpfr_swap(u,a.u);}
 ~I(){mpfr_clear(l);mpfr_clear(u);}
 I& operator=(const I&a){if(this!=&a){mpfr_set(l,a.l,MPFR_RNDD);mpfr_set(u,a.u,MPFR_RNDU);}return *this;}
 I& operator=(I&&a) noexcept {mpfr_swap(l,a.l);mpfr_swap(u,a.u);return *this;}
 static I pow2(long e){I a;mpfr_set_ui_2exp(a.l,1,e,MPFR_RNDD);mpfr_set(a.u,a.l,MPFR_RNDU);return a;}
 static I decimal(std::string s){ // Finite decimal -> exact GMP rational -> directed endpoints.
  long ex=0;auto p=s.find_first_of("eE");if(p!=std::string::npos){ex=std::stol(s.substr(p+1));s.resize(p);}
  p=s.find('.');if(p!=std::string::npos){ex-=s.size()-p-1;s.erase(p,1);}
  mpq_t q;mpq_init(q);if(mpz_set_str(mpq_numref(q),s.c_str(),10))throw std::runtime_error("decimal parse");
  mpz_t ten;mpz_init(ten);mpz_ui_pow_ui(ten,10,std::abs(ex));
  if(ex<0)mpz_set(mpq_denref(q),ten);else mpz_mul(mpq_numref(q),mpq_numref(q),ten);
  mpq_canonicalize(q);I a;mpfr_set_q(a.l,q,MPFR_RNDD);mpfr_set_q(a.u,q,MPFR_RNDU);mpz_clear(ten);mpq_clear(q);return a;
 }
 bool zero()const{return mpfr_sgn(l)<=0&&mpfr_sgn(u)>=0;}
 bool point()const{return mpfr_equal_p(l,u);}
 bool finite()const{return mpfr_number_p(l)&&mpfr_number_p(u);}
 double upper()const{return mpfr_get_d(u,MPFR_RNDU);}
 double lower()const{return mpfr_get_d(l,MPFR_RNDD);}
};
inline I operator+(const I&a,const I&b){I c;mpfr_add(c.l,a.l,b.l,MPFR_RNDD);mpfr_add(c.u,a.u,b.u,MPFR_RNDU);return c;}
inline I operator-(const I&a){I c;mpfr_neg(c.l,a.u,MPFR_RNDD);mpfr_neg(c.u,a.l,MPFR_RNDU);return c;}
inline I operator-(const I&a,const I&b){I c;mpfr_sub(c.l,a.l,b.u,MPFR_RNDD);mpfr_sub(c.u,a.u,b.l,MPFR_RNDU);return c;}
inline I mul_four_corners(const I&a,const I&b){I c;mpfr_t t;mpfr_init2(t,PREC);mpfr_set_inf(c.l,1);mpfr_set_inf(c.u,-1);
 for(auto x:{a.l,a.u})for(auto y:{b.l,b.u}){mpfr_mul(t,x,y,MPFR_RNDD);if(mpfr_less_p(t,c.l))mpfr_set(c.l,t,MPFR_RNDD);mpfr_mul(t,x,y,MPFR_RNDU);if(mpfr_greater_p(t,c.u))mpfr_set(c.u,t,MPFR_RNDU);}mpfr_clear(t);return c;}
inline I mul_by_sign(const I&a,const I&b){
 I c;mpfr_srcptr xl=nullptr,yl=nullptr,xu=nullptr,yu=nullptr;
 if(mpfr_sgn(a.l)>=0){
  if(mpfr_sgn(b.l)>=0){xl=a.l;yl=b.l;xu=a.u;yu=b.u;}
  else if(mpfr_sgn(b.u)<=0){xl=a.u;yl=b.l;xu=a.l;yu=b.u;}
  else {xl=a.u;yl=b.l;xu=a.u;yu=b.u;}
 }else if(mpfr_sgn(a.u)<=0){
  if(mpfr_sgn(b.l)>=0){xl=a.l;yl=b.u;xu=a.u;yu=b.l;}
  else if(mpfr_sgn(b.u)<=0){xl=a.u;yl=b.u;xu=a.l;yu=b.l;}
  else {xl=a.l;yl=b.u;xu=a.l;yu=b.l;}
 }else{
  if(mpfr_sgn(b.l)>=0){xl=a.l;yl=b.u;xu=a.u;yu=b.u;}
  else if(mpfr_sgn(b.u)<=0){xl=a.u;yl=b.l;xu=a.l;yu=b.l;}
  else {mpfr_t t;mpfr_init2(t,PREC);mpfr_mul(c.l,a.l,b.u,MPFR_RNDD);mpfr_mul(t,a.u,b.l,MPFR_RNDD);mpfr_min(c.l,c.l,t,MPFR_RNDD);mpfr_mul(c.u,a.l,b.l,MPFR_RNDU);mpfr_mul(t,a.u,b.u,MPFR_RNDU);mpfr_max(c.u,c.u,t,MPFR_RNDU);mpfr_clear(t);return c;}
 }
 mpfr_mul(c.l,xl,yl,MPFR_RNDD);mpfr_mul(c.u,xu,yu,MPFR_RNDU);return c;
}
inline I operator*(const I&a,const I&b){
#ifdef FAST_INTERVAL
 return mul_by_sign(a,b);
#else
 return mul_four_corners(a,b);
#endif
}
inline I inv(const I&a){if(a.zero())throw std::runtime_error("division contains zero");I c;mpfr_ui_div(c.l,1,a.u,MPFR_RNDD);mpfr_ui_div(c.u,1,a.l,MPFR_RNDU);return c;}
inline I operator/(const I&a,const I&b){return a*inv(b);}
inline I sqr(const I&a){I c;if(a.zero())mpfr_set_zero(c.l,0);else {mpfr_sqr(c.l,a.l,MPFR_RNDD);mpfr_sqr(c.u,a.u,MPFR_RNDD);if(mpfr_less_p(c.u,c.l))mpfr_set(c.l,c.u,MPFR_RNDD);}mpfr_t t;mpfr_init2(t,PREC);mpfr_sqr(c.u,a.l,MPFR_RNDU);mpfr_sqr(t,a.u,MPFR_RNDU);if(mpfr_greater_p(t,c.u))mpfr_set(c.u,t,MPFR_RNDU);mpfr_clear(t);return c;}
inline I sqrt(const I&a){if(mpfr_sgn(a.l)<0)throw std::runtime_error("sqrt negative");I c;mpfr_sqrt(c.l,a.l,MPFR_RNDD);mpfr_sqrt(c.u,a.u,MPFR_RNDU);return c;}
inline I mid(const I&a){I c;mpfr_add(c.l,a.l,a.u,MPFR_RNDN);mpfr_div_2ui(c.l,c.l,1,MPFR_RNDN);mpfr_set(c.u,c.l,MPFR_RNDN);if(mpfr_less_p(c.l,a.l)||mpfr_greater_p(c.u,a.u))throw std::runtime_error("mid outside");return c;}
inline I absup(const I&a){I c;mpfr_abs(c.u,a.l,MPFR_RNDU);mpfr_t t;mpfr_init2(t,PREC);mpfr_abs(t,a.u,MPFR_RNDU);if(mpfr_greater_p(t,c.u))mpfr_set(c.u,t,MPFR_RNDU);mpfr_set(c.l,c.u,MPFR_RNDD);mpfr_clear(t);return c;}
inline I upperpoint(const I&a){I c;mpfr_set(c.l,a.u,MPFR_RNDD);mpfr_set(c.u,a.u,MPFR_RNDU);return c;}
inline I lowerpoint(const I&a){I c;mpfr_set(c.l,a.l,MPFR_RNDD);mpfr_set(c.u,a.l,MPFR_RNDU);return c;}
inline I sym(const I&a){I c=absup(a);mpfr_neg(c.l,c.u,MPFR_RNDD);return c;}
inline I hull(const I&a,const I&b){I c;mpfr_min(c.l,a.l,b.l,MPFR_RNDD);mpfr_max(c.u,a.u,b.u,MPFR_RNDU);return c;}
inline bool subset(const I&a,const I&b){return mpfr_lessequal_p(b.l,a.l)&&mpfr_lessequal_p(a.u,b.u);}
inline bool lt(const I&a,const I&b){return mpfr_less_p(a.u,b.l);}
inline I maxupper(const I&a,const I&b){return mpfr_greater_p(a.u,b.u)?upperpoint(a):upperpoint(b);}
inline I minlower(const I&a,const I&b){return mpfr_less_p(a.l,b.l)?lowerpoint(a):lowerpoint(b);}
inline std::string dyadic(mpfr_srcptr x){if(mpfr_zero_p(x))return "[\"0\",0]";mpz_t z;mpz_init(z);long e=mpfr_get_z_2exp(z,x);char*s=mpz_get_str(nullptr,10,z);std::string out="[\""+std::string(s)+"\","+std::to_string(e)+"]";free(s);mpz_clear(z);return out;}
inline std::string json(const I&a){I m=mid(a);I r=absup(a-m);return "{\"mid\":"+dyadic(m.l)+",\"rad\":"+dyadic(r.u)+"}";}
using V=std::vector<I>;using Mat=std::vector<V>;
inline Mat zeros(int n,int m){return Mat(n,V(m));}
inline Mat eye(int n){Mat a=zeros(n,n);for(int i=0;i<n;i++)a[i][i]=I(1);return a;}
inline V add(const V&a,const V&b){V c=a;for(size_t i=0;i<c.size();i++)c[i]=a[i]+b[i];return c;}
inline V sub(const V&a,const V&b){V c=a;for(size_t i=0;i<c.size();i++)c[i]=a[i]-b[i];return c;}
inline V scale(const V&a,const I&x){V c=a;for(I&z:c)z=z*x;return c;}
inline V mv(const Mat&a,const V&b){V c(a.size());for(size_t i=0;i<a.size();i++)for(size_t j=0;j<b.size();j++)c[i]=c[i]+a[i][j]*b[j];return c;}
inline Mat mm(const Mat&a,const Mat&b){Mat c=zeros(a.size(),b[0].size());for(size_t i=0;i<a.size();i++)for(size_t k=0;k<b.size();k++)for(size_t j=0;j<b[0].size();j++)c[i][j]=c[i][j]+a[i][k]*b[k][j];return c;}
inline Mat ma(const Mat&a,const Mat&b){Mat c=a;for(size_t i=0;i<a.size();i++)c[i]=add(a[i],b[i]);return c;}
inline Mat ms(const Mat&a,const Mat&b){Mat c=a;for(size_t i=0;i<a.size();i++)c[i]=sub(a[i],b[i]);return c;}
inline Mat scale(const Mat&a,const I&x){Mat c=a;for(V&v:c)v=scale(v,x);return c;}
inline I norm(const V&a){I x;for(const I&z:a)x=maxupper(x,absup(z));return x;}
inline I norm(const Mat&a){I x;for(const V&v:a){I s;for(const I&z:v)s=s+absup(z);x=maxupper(x,s);}return x;}
inline Mat midpoint(const Mat&a){Mat c=a;for(V&v:c)for(I&z:v)z=mid(z);return c;}
inline V midpoint(const V&a){V c=a;for(I&z:c)z=mid(z);return c;}
inline V symmetric(const V&a){V c=a;for(I&z:c)z=sym(z);return c;}
inline std::string json(const V&a){std::string s="[";for(size_t i=0;i<a.size();i++)s+=(i?",":"")+json(a[i]);return s+"]";}
inline std::string json(const Mat&a){std::string s="[";for(size_t i=0;i<a.size();i++)s+=(i?",":"")+json(a[i]);return s+"]";}
