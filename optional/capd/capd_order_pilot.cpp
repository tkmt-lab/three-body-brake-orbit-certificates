// R015 fresh local benchmark; existing research/verification code was not read.
// Compile and run on Modal only. This is not a whole-orbit certificate.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>
#include "capd/capdlib.h"
#include "capd/mpcapdlib.h"
#include "capd/dynsys/MpStepControl.h"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/dynsys/BasicOdeSolver.hpp"
#include "capd/poincare/TimeMap.hpp"
#include "capd/dynset/C1DoubletonSet.hpp"

using R = capd::MpFloat;
using I = capd::MpInterval;
using V = capd::MpIVector;
using M = capd::MpIMatrix;
using Solver = capd::dynsys::OdeSolver<capd::MpIMap, capd::dynsys::MpLastTermsStepControl>;
using TimeMap = capd::poincare::TimeMap<Solver>;
using Clock = std::chrono::steady_clock;

// Access the representation of a COPY through a public derived helper; no CAPD patch.
struct ExactReal : R {
  explicit ExactReal(const R& r):R(r) {}
  std::string json() const {
    if(!mpfr_number_p(mpfr_rep)) throw std::runtime_error("nonfinite endpoint");
    mpz_t n; mpz_init(n);
    auto e=mpfr_get_z_2exp(n,mpfr_rep);
    char* raw=mpz_get_str(nullptr,10,n);
    std::string value(raw);
    void (*release)(void*,size_t);
    mp_get_memory_functions(nullptr,nullptr,&release);
    release(raw,std::strlen(raw)+1); mpz_clear(n);
    return "{\"n\":\""+value+"\",\"e\":"+std::to_string(e)+"}";
  }
};
std::string exact(const R& r) { return ExactReal(r).json(); }
std::string exact(const I& r) { return "["+exact(r.leftBound())+","+exact(r.rightBound())+"]"; }
std::string vector_json(const V& v) {
  std::string s="[";
  for(int i=0;i<8;++i) { if(i) s+=","; s+=exact(v[i]); }
  return s+"]";
}
std::string matrix_json(const M& v) {
  std::string s="[";
  for(int i=0;i<8;++i) {
    if(i) s+=","; s+="[";
    for(int j=0;j<8;++j) { if(j) s+=","; s+=exact(v[i][j]); }
    s+="]";
  }
  return s+"]";
}
std::string quote(const std::string& s) {
  std::string out="\"";
  for(unsigned char c:s) {
    if(c=='\"' || c=='\\') { out+='\\'; out+=c; }
    else if(c=='\n') out+="\\n";
    else if(c=='\r') out+="\\r";
    else if(c=='\t') out+="\\t";
    else if(c<32) out+='?';
    else out+=c;
  }
  return out+"\"";
}
I decimal(const std::string& s) { return I(s,s); }
I dyadic(int e) {
  I r(1);
  while(e<0) { r/=I(2); ++e; }
  return r;
}
I rational(const std::string& n,const std::string& d) { return decimal(n)/decimal(d); }
std::vector<I> distances2(const V& y) {
  I x12=y[2]-y[0], y12=y[3]-y[1];
  I x13=-I(2)*y[0]-y[2], y13=-I(2)*y[1]-y[3];
  I x23=-y[0]-I(2)*y[2], y23=-y[1]-I(2)*y[3];
  return {capd::intervals::sqr(x12)+capd::intervals::sqr(y12),capd::intervals::sqr(x13)+capd::intervals::sqr(y13),
          capd::intervals::sqr(x23)+capd::intervals::sqr(y23)};
}
R width_upper(const I& x) { return (I(x.rightBound())-I(x.leftBound())).rightBound(); }
R abs_upper(const I& x) { return std::max(abs(x.leftBound()),abs(x.rightBound())); }
std::string field_string() {
  std::string x12="(x2-x1)",y12="(y2-y1)",x13="(-2*x1-x2)",y13="(-2*y1-y2)";
  std::string x23="(-x1-2*x2)",y23="(-y1-2*y2)";
  auto g=[](std::string x,std::string y,std::string numerator) {
    std::string s="("+x+"^2+"+y+"^2)";
    return "("+numerator+"/("+s+"*sqrt("+s+")))";
  };
  return "var:x1,y1,x2,y2,u1,w1,u2,w2;fun:u1,w1,u2,w2,"
    +g(x12,y12,x12)+"+"+g(x13,y13,x13)+","
    +g(x12,y12,y12)+"+"+g(x13,y13,y13)+",-"
    +g(x12,y12,x12)+"+"+g(x23,y23,x23)+",-"
    +g(x12,y12,y12)+"+"+g(x23,y23,y23)+";";
}

int scout(const char* path) {
  auto started=Clock::now(); auto cpu=std::clock();
  R::setDefaultPrecision(384); R::setDefaultRndMode(R::RoundNearest);
  std::ifstream input(path); std::string aa,bb,tt,unused;
  if(!(input>>aa>>bb>>tt>>unused)) throw std::runtime_error("missing row for scout");
  capd::MpVector y(8); y.clear();
  y[0]=R(-1); y[2]=R(1); y[4]=y[6]=R(aa); y[5]=y[7]=R(bb);
  capd::MpMap field(field_string());
  capd::MpOdeSolver solver(field,28,capd::dynsys::DLastTermsStepControl(2,1e-100));
  solver.setAbsoluteTolerance(1e-50); solver.setRelativeTolerance(1e-50);
  capd::poincare::TimeMap<capd::MpOdeSolver> map(solver); map.stopAfterStep(true);
  R target=R(tt)/R(4), length=std::min(dyadic(-8).leftBound(),R(tt)/R(256)), t(0);
  std::vector<R> times, ds;
  std::vector<capd::MpVector> states;
  auto save=[&] {
    V point(8); for(int i=0;i<8;++i) point[i]=I(y[i]);
    auto distances=distances2(point);
    R d=std::min(distances[0].leftBound(),std::min(distances[1].leftBound(),distances[2].leftBound()));
    times.push_back(t); states.push_back(y); ds.push_back(d);
  };
  save();
  do {
    if(std::chrono::duration<double>(Clock::now()-started).count()>55)
      throw std::runtime_error("384-bit point scout exceeded 55 seconds");
    R before=t;
    map(target,y,t);
    if(!(t>before)) throw std::runtime_error("384-bit point scout time stalled");
    save();
  } while(!map.completed());
  size_t closest=std::min_element(ds.begin(),ds.end())-ds.begin();
  R desired=std::min(std::max(R(0),times[closest]-length/R(2)),target-length);
  size_t chosen=0;
  R best=abs(times[0]-desired);
  for(size_t i=1;i<times.size();++i) if(times[i]<=target-length && abs(times[i]-desired)<best) {
    chosen=i; best=abs(times[i]-desired);
  }
  std::cout<<"{\"kind\":\"scout\",\"validated\":false,\"method\":\"CAPD BasicOdeSolver point Taylor\","
    <<"\"precision\":384,\"order\":28,\"tolerance\":\"1e-50\",\"complete\":true,\"steps\":"<<times.size()-1
    <<",\"closest_time\":"<<exact(times[closest])<<",\"min_sampled_distance2\":"<<exact(ds[closest])
    <<",\"fixture_start_time_exact\":"<<exact(times[chosen])<<",\"fixture_length_exact\":"<<exact(length)
    <<",\"fixture_center_exact\":[";
  for(int i=0;i<8;++i) { if(i) std::cout<<","; std::cout<<exact(states[chosen][i]); }
  std::cout<<"],\"elapsed_seconds\":"<<std::chrono::duration<double>(Clock::now()-started).count()
    <<",\"cpu_seconds\":"<<double(std::clock()-cpu)/CLOCKS_PER_SEC<<"}"<<std::endl;
  return 0;
}

int main(int argc,char** argv) {
  auto started=Clock::now(); auto cpu=std::clock();
  auto elapsed=[&] { return std::chrono::duration<double>(Clock::now()-started).count(); };
  std::cout.precision(17);
  try {
    if(argc!=4 && argc!=5) throw std::runtime_error("usage: capd_pilot fixture.txt precision order [c1_only]");
    if(std::string(argv[2])=="scout") return scout(argv[1]);
    bool only_c1=argc==5 && std::string(argv[4])=="c1_only";
    int precision=std::stoi(argv[2]), order=std::stoi(argv[3]);
    R::setDefaultPrecision(precision); R::setDefaultRndMode(R::RoundNearest);
    std::ifstream input(argv[1]); std::string n,d;
    if(!(input>>n>>d)) throw std::runtime_error("missing fixture time");
    I target=rational(n,d);
    V center(8);
    for(int i=0;i<8;++i) { if(!(input>>n>>d)) throw std::runtime_error("missing fixture state"); center[i]=rational(n,d); }
    capd::MpIMap field(field_string());
    // Elementary rounding and original Newton-equation checks, performed remotely.
    I tenth=decimal("0.1")*I(10);
    if(!(tenth.leftBound()<=R(1) && tenth.rightBound()>=R(1))) throw std::runtime_error("rounding smoke failure");
    V initial(8); initial.clear(); initial[0]=I(-1); initial[2]=I(1);
    V f=field(initial);
    if(!(f[4].leftBound()<=R(1.25) && f[4].rightBound()>=R(1.25)
      && f[6].leftBound()<=R(-1.25) && f[6].rightBound()>=R(-1.25)))
      throw std::runtime_error("Newton field smoke failure");
    M jac=field.derivative(initial);
    double expected[8][8] = {};
    for(int i=0;i<4;++i) expected[i][i+4]=1;
    expected[4][0]=expected[6][2]=4.25;
    expected[4][2]=expected[6][0]=1.75;
    expected[5][1]=expected[7][3]=-2.125;
    expected[5][3]=expected[7][1]=-.875;
    for(int i=0;i<8;++i) for(int j=0;j<8;++j)
      if(!(jac[i][j].leftBound()<=R(expected[i][j]) && jac[i][j].rightBound()>=R(expected[i][j])))
        throw std::runtime_error("Newton Jacobian smoke failure");
    auto integrate=[&](auto& set,const char* phase) {
      Solver solver(field,order,capd::dynsys::MpLastTermsStepControl(1,R("1e-100")));
      R tol=dyadic(-120).rightBound();
      solver.setAbsoluteTolerance(tol); solver.setRelativeTolerance(tol);
      TimeMap map(solver); map.stopAfterStep(true);
      unsigned steps=0; bool positive=true, complete=false;
      R min_d2(1), min_h(1);
      do {
        if(elapsed()>116) break;
        map(target,set); ++steps;
        auto pair=distances2(set.getLastEnclosure());
        for(const I& v:pair) { min_d2=std::min(min_d2,v.leftBound()); positive=positive && v.leftBound()>R(0); }
        min_h=std::min(min_h,solver.getStep().leftBound());
        std::cout<<"{\"kind\":\"step\",\"phase\":"<<quote(phase)<<",\"steps\":"<<steps
          <<",\"time\":"<<exact(map.getCurrentTime())<<",\"h\":"<<exact(solver.getStep())
          <<",\"min_distance2_lower\":"<<exact(min_d2)<<",\"elapsed_seconds\":"<<elapsed();
        if(std::string(phase)=="c1_box") std::cout<<",\"all_time_tube\":"<<vector_json(set.getLastEnclosure());
        std::cout<<"}"<<std::endl;
        complete=map.completed();
        if(!positive) break;
      } while(!complete);
      std::cout<<"{\"kind\":\"phase\",\"phase\":"<<quote(phase)<<",\"complete\":"<<(complete?"true":"false")
        <<",\"all_time_distances_positive\":"<<(positive?"true":"false")<<",\"steps\":"<<steps
        <<",\"reached_time\":"<<exact(map.getCurrentTime())<<",\"min_h\":"<<exact(min_h)
        <<",\"min_distance2_lower\":"<<exact(min_d2)<<",\"endpoint\":"<<vector_json(V(set))
        <<",\"elapsed_seconds\":"<<elapsed()<<"}"<<std::endl;
      return complete && positive;
    };
    capd::MpC0Rect2Set point(center);
    bool c0=only_c1 ? false : integrate(point,"c0_center");
    R max_point_width(0);
    V pv(point);
    for(int i=0;i<8;++i) max_point_width=std::max(max_point_width,width_upper(pv[i]));
    bool c0_quality=c0 && max_point_width<=dyadic(-100).leftBound();
    V box(center); for(int i=0;i<8;++i) box[i]+=I(-1,1)*dyadic(-80);
    capd::MpC1Rect2Set set(box);
    bool c1=false, c1_quality=false;
    if((c0 || only_c1) && elapsed()<=116) {
      c1=integrate(set,"c1_box"); c1_quality=c1;
      M matrix(set);
      for(int i=0;i<8;++i) for(int j=0;j<8;++j) {
        R scale=std::max(R(1),abs_upper(matrix[i][j]));
        R allowed=(dyadic(-60)*I(scale)).leftBound();
        c1_quality=c1_quality && width_upper(matrix[i][j])<=allowed;
      }
      std::cout<<"{\"kind\":\"c1_endpoint\",\"matrix\":"<<matrix_json(matrix)<<"}"<<std::endl;
    }
    struct rusage usage; getrusage(RUSAGE_SELF,&usage);
    std::cout<<"{\"kind\":\"final\",\"precision\":"<<precision<<",\"order\":"<<order
      <<",\"mode\":"<<quote(only_c1 ? "c1_only" : "combined")
      <<",\"c0_complete\":"<<(c0?"true":"false")<<",\"c1_complete\":"<<(c1?"true":"false")
      <<",\"c0_quality\":"<<(c0_quality?"true":"false")<<",\"c1_quality\":"<<(c1_quality?"true":"false")
      <<",\"point_width_upper\":"<<exact(max_point_width)<<",\"elapsed_seconds\":"<<elapsed()
      <<",\"cpu_seconds\":"<<double(std::clock()-cpu)/CLOCKS_PER_SEC
      <<",\"max_rss_kib_linux\":"<<usage.ru_maxrss<<"}"<<std::endl;
    return c0_quality && c1_quality ? 0 : 2;
  } catch(const std::exception& e) {
    std::cout<<"{\"kind\":\"error\",\"message\":"<<quote(e.what())<<",\"elapsed_seconds\":"<<elapsed()<<"}"<<std::endl;
    return 1;
  }
}
