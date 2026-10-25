// Batch2: unchanged field/integrator; retain full C1 doubleton at final record. Uses only this team's equations/CAPD.
#define main unused_pilot_main
#include "capd_order_pilot.cpp"
#undef main
#include <type_traits>

int radius_exp=-90;
std::string vjson(const V& v) {
  std::string s="[";for(int i=0;i<v.dimension();++i){if(i)s+=",";s+=exact(v[i]);}return s+"]";
}
std::string mjson(const M& m) {
  std::string s="[";for(int i=0;i<m.numberOfRows();++i){if(i)s+=",";s+="[";
    for(int j=0;j<m.numberOfColumns();++j){if(j)s+=",";s+=exact(m[i][j]);}s+="]";}return s+"]";
}
template<class Set>std::string c0_representation(const Set& set) {
  return "{\"x\":"+vjson(set.get_x())+",\"C\":"+mjson(set.get_C())+",\"r0\":"+vjson(set.get_r0())
    +",\"B\":"+mjson(set.get_B())+",\"r\":"+vjson(set.get_r())+"}";
}
std::string normalized_field() {
  std::string old=field_string(),functions=old.substr(old.find("fun:")+4);functions.pop_back();
  std::stringstream stream(functions);std::string term;
  std::string out="var:x1,y1,x2,y2,u1,w1,u2,w2,theta;fun:";
  for(int i=0;i<8;++i){if(!std::getline(stream,term,','))throw std::runtime_error("field split");
    if(i)out+=",";out+="(theta/512)*("+term+")";}
  return out+",0;";
}
M parameter_derivative(const M& matrix) {
  M d(9,3);for(int i=0;i<9;++i){d[i][0]=matrix[i][4]+matrix[i][6];d[i][1]=matrix[i][5]+matrix[i][7];d[i][2]=matrix[i][8];}return d;
}
void smoke(capd::MpIMap& field) {
  V y(9);y.clear();y[0]=I(-1);y[2]=I(1);y[8]=I(512);
  V f=field(y);M j=field.derivative(y);
  double a[9][9]={};for(int i=0;i<4;++i)a[i][i+4]=1;
  a[4][0]=a[6][2]=4.25;a[4][2]=a[6][0]=1.75;
  a[5][1]=a[7][3]=-2.125;a[5][3]=a[7][1]=-.875;
  a[4][8]=1.25/512;a[6][8]=-1.25/512;
  for(int i=0;i<9;++i)for(int k=0;k<9;++k)
    if(!(j[i][k].leftBound()<=R(a[i][k])&&j[i][k].rightBound()>=R(a[i][k])))throw std::runtime_error("9D Jacobian smoke");
  if(!(f[4].leftBound()<=R(1.25)&&f[4].rightBound()>=R(1.25)&&f[8]==I(0)))throw std::runtime_error("9D field smoke");
  std::cout<<"{\"kind\":\"smoke\",\"status\":\"PASS\",\"clock_target\":128,\"physical_clock\":\"theta*s/512\"}"<<std::endl;
}

template<class Set>
int flow(Set& set,capd::MpIMap& field,const std::string& phase,int order,int tol_exp,int max_steps,int seconds) {
  auto start=Clock::now();auto cpu=std::clock();
  Solver solver(field,order,capd::dynsys::MpLastTermsStepControl(1,R("1e-100")));
  R tol=dyadic(tol_exp).rightBound();solver.setAbsoluteTolerance(tol);solver.setRelativeTolerance(tol);
  TimeMap map(solver);map.stopAfterStep(true);unsigned count=0;bool complete=false,positive=true;R minimum(1);
  std::cout<<"{\"kind\":\"initial\",\"phase\":"<<quote(phase)<<",\"box\":"<<vjson(V(set))<<"}"<<std::endl;
  do {
    if(std::chrono::duration<double>(Clock::now()-start).count()>seconds)break;
    if(max_steps&&count>=unsigned(max_steps))break;
    V before(set);map(I(128),set);++count;
    V tube=set.getLastEnclosure();auto ds=distances2(tube);
    for(const I& d:ds){minimum=std::min(minimum,d.leftBound());positive=positive&&d.leftBound()>R(0);}
    std::cout<<"{\"kind\":\"step\",\"phase\":"<<quote(phase)<<",\"index\":"<<count
      <<",\"time\":"<<exact(map.getCurrentTime())<<",\"h\":"<<exact(solver.getStep())
      <<",\"input_box\":"<<vjson(before)<<",\"all_time_tube\":"<<vjson(tube)
      <<",\"endpoint\":"<<vjson(V(set))<<",\"doubleton\":"<<c0_representation(set)
      <<",\"distance2_running_lower\":"<<exact(minimum);
    if constexpr(std::is_same<Set,capd::MpC1Rect2Set>::value) {
      // CAPD C1DoubletonSet publicly exposes jacPhi AFTER jacRem was added.
      // encloseC1Map evaluated this derivative over its full input axis box.
      std::cout<<",\"local_jacobian\":"<<mjson(set.jacPhi)
        <<",\"global_parameter_derivative\":"<<mjson(parameter_derivative(M(set)));
    }
    std::cout<<"}"<<std::endl;
    if(count<=4||count%16==0)std::cout<<"{\"kind\":\"progress\",\"phase\":"<<quote(phase)<<",\"steps\":"<<count
      <<",\"time\":"<<exact(map.getCurrentTime())<<",\"elapsed_seconds\":"<<std::chrono::duration<double>(Clock::now()-start).count()<<"}"<<std::endl;
    complete=map.completed();if(!positive)break;
  }while(!complete);
  struct rusage usage;getrusage(RUSAGE_SELF,&usage);
  std::cout<<"{\"kind\":\"final\",\"phase\":"<<quote(phase)<<",\"complete\":"<<(complete?"true":"false")
    <<",\"positive\":"<<(positive?"true":"false")<<",\"steps\":"<<count
    <<",\"time\":"<<exact(map.getCurrentTime())<<",\"endpoint\":"<<vjson(V(set))
    <<",\"nondegeneracy\":"<<exact(V(set)[0]+I(2)*V(set)[2])
    <<",\"distance2_lower\":"<<exact(minimum)<<",\"elapsed_seconds\":"<<std::chrono::duration<double>(Clock::now()-start).count()
    <<",\"cpu_seconds\":"<<double(std::clock()-cpu)/CLOCKS_PER_SEC<<",\"max_rss_kib_linux\":"<<usage.ru_maxrss;
  std::cout<<",\"checkpoint_c0\":"<<c0_representation(set)<<",\"checkpoint_invB\":"<<mjson(set.get_invB());
  if constexpr(std::is_same<Set,capd::MpC1Rect2Set>::value)std::cout<<",\"global_parameter_derivative\":"<<mjson(parameter_derivative(M(set)))
    <<",\"checkpoint_full_M\":"<<mjson(M(set))
    <<",\"checkpoint_c1\":{\"D\":"<<mjson(set.get_D())<<",\"Cjac\":"<<mjson(set.get_Cjac())<<",\"R0\":"<<mjson(set.get_R0())
    <<",\"Bjac\":"<<mjson(set.get_Bjac())<<",\"R\":"<<mjson(set.get_R())<<",\"invBjac\":"<<mjson(set.get_invBjac())<<"}";
  std::cout<<"}"<<std::endl;
  return positive&&(complete||max_steps)?0:2;
}

struct BinaryReal:R {
  BinaryReal(const std::string& n,int e):R(n){mpfr_mul_2si(mpfr_rep,mpfr_rep,e,MPFR_RNDN);}
};
I read_interval(std::istream& input) {
  std::string a,b;int e,f;if(!(input>>a>>e>>b>>f))throw std::runtime_error("missing interval");
  return I(BinaryReal(a,e),BinaryReal(b,f));
}
int terminal(const char* path) {
  R::setDefaultPrecision(1024);R::setDefaultRndMode(R::RoundNearest);
  std::ifstream input(path);V f(3);M j(3,3);
  for(int i=0;i<3;++i)f[i]=read_interval(input);
  for(int i=0;i<3;++i)for(int k=0;k<3;++k)j[i][k]=read_interval(input);
  I nd=read_interval(input),radius=dyadic(radius_exp);R aug[3][6];
  for(int i=0;i<3;++i)for(int k=0;k<6;++k)aug[i][k]=k<3?(j[i][k].leftBound()+j[i][k].rightBound())/R(2):R(k-3==i?1:0);
  for(int k=0;k<3;++k){int pivot=k;for(int i=k+1;i<3;++i)if(abs(aug[i][k])>abs(aug[pivot][k]))pivot=i;
    if(aug[pivot][k]==R(0))throw std::runtime_error("singular midpoint Jacobian");
    for(int q=0;q<6;++q)std::swap(aug[pivot][q],aug[k][q]);R d=aug[k][k];
    for(int q=0;q<6;++q)aug[k][q]/=d;
    for(int i=0;i<3;++i)if(i!=k){R v=aug[i][k];for(int q=0;q<6;++q)aug[i][q]-=v*aug[k][q];}
  }
  M c(3,3);for(int i=0;i<3;++i)for(int k=0;k<3;++k)c[i][k]=I(aug[i][k+3]);
  M w=M::Identity(3)-c*j;V delta(3),box(3);for(int i=0;i<3;++i)box[i]=I(-1,1)*radius;
  delta=-c*f+w*box;R q(0);bool pass=true;
  for(int i=0;i<3;++i){I row(0);for(int k=0;k<3;++k)row+=I(abs_upper(w[i][k]));q=std::max(q,row.rightBound());pass=pass&&abs_upper(delta[i])<radius.leftBound();}
  I det=c[0][0]*(c[1][1]*c[2][2]-c[1][2]*c[2][1])-c[0][1]*(c[1][0]*c[2][2]-c[1][2]*c[2][0])+c[0][2]*(c[1][0]*c[2][1]-c[1][1]*c[2][0]);
  pass=pass&&q<R(1)&&(det.leftBound()>R(0)||det.rightBound()<R(0))&&(nd.leftBound()>R(0)||nd.rightBound()<R(0));
  std::cout<<"{\"kind\":\"terminal\",\"precision\":1024,\"F0\":"<<vjson(f)<<",\"J\":"<<mjson(j)
    <<",\"C\":"<<mjson(c)<<",\"W\":"<<mjson(w)<<",\"K_delta\":"<<vjson(delta)
    <<",\"q_upper\":"<<exact(q)<<",\"det_C\":"<<exact(det)<<",\"nondegeneracy\":"<<exact(nd)
    <<",\"radius\":"<<exact(radius)<<",\"pass\":"<<(pass?"true":"false")<<"}"<<std::endl;
  return pass?0:2;
}
int main(int argc,char** argv) {
  std::cout.precision(17);
  try {
    if((argc==3||argc==4)&&std::string(argv[1])=="terminal"){if(argc==4)radius_exp=std::stoi(argv[3]);return terminal(argv[2]);}
    if(argc!=9)throw std::runtime_error("usage: batch1_solver row phase precision order tol_exp max_steps seconds");
    radius_exp=std::stoi(argv[8]);if(radius_exp!=-130&&radius_exp!=-140)throw std::runtime_error("unauthorized radius");
    std::string phase=argv[2];int p=std::stoi(argv[3]),order=std::stoi(argv[4]),tol=std::stoi(argv[5]),limit=std::stoi(argv[6]),seconds=std::stoi(argv[7]);
    R::setDefaultPrecision(p);R::setDefaultRndMode(R::RoundNearest);
    std::ifstream input(argv[1]);std::string a,b,t,unused;if(!(input>>a>>b>>t>>unused))throw std::runtime_error("missing source row");
    V center(9);center.clear();center[0]=I(-1);center[2]=I(1);center[4]=center[6]=decimal(a);center[5]=center[7]=decimal(b);center[8]=decimal(t);
    capd::MpIMap field(normalized_field());smoke(field);
    std::cout<<"{\"kind\":\"setup\",\"phase\":"<<quote(phase)<<",\"precision\":"<<p<<",\"order\":"<<order
      <<",\"tol_exp\":"<<tol<<",\"center\":"<<vjson(center)<<",\"parameter_radius\":"<<exact(dyadic(radius_exp))<<"}"<<std::endl;
    if(phase=="c0"){capd::MpC0Rect2Set set(center);return flow(set,field,phase,order,tol,limit,seconds);}
    if(phase!="c1")throw std::runtime_error("invalid phase");
    M basis(9,9);basis.clear();basis[4][0]=basis[6][0]=I(1);basis[5][1]=basis[7][1]=I(1);basis[8][2]=I(1);
    V radii(9);radii.clear();for(int i=0;i<3;++i)radii[i]=I(-1,1)*dyadic(radius_exp);
    capd::MpC1Rect2Set set(center,basis,radii);return flow(set,field,phase,order,tol,limit,seconds);
  }catch(const std::exception& e){std::cout<<"{\"kind\":\"error\",\"message\":"<<quote(e.what())<<"}"<<std::endl;return 1;}
}
