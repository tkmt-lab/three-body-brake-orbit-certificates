#include "dynamics.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cstdio>
using Clock=std::chrono::steady_clock;
inline void save(const std::string&path,const std::string&data){std::ofstream f(path+".tmp");f<<data<<"\n";f.close();if(std::rename((path+".tmp").c_str(),path.c_str()))throw std::runtime_error("rename output");}
inline std::string fraction(const I&a){if(!a.point())throw std::runtime_error("fraction requires point");mpq_t q;mpq_init(q);mpfr_get_q(q,a.l);char*n=mpz_get_str(nullptr,10,mpq_numref(q));char*d=mpz_get_str(nullptr,10,mpq_denref(q));std::string s="[\""+std::string(n)+"\",\""+d+"\"]";free(n);free(d);mpq_clear(q);return s;}
int main(int argc,char**argv){
 if(argc!=10){std::cerr<<"certify data row precision order radius_exp tolerance_exp max_steps seconds output_prefix\n";return 2;}
 int row=std::stoi(argv[2]);PREC=std::stoi(argv[3]);int N=std::stoi(argv[4]),re=std::stoi(argv[5]),te=std::stoi(argv[6]),maxsteps=std::stoi(argv[7]);double cap=std::stod(argv[8]);std::string prefix=argv[9];
 auto start=Clock::now();auto elapsed=[&](){return std::chrono::duration<double>(Clock::now()-start).count();};
 std::ifstream input(argv[1]);std::string line;for(int i=0;i<row;i++)if(!std::getline(input,line))throw std::runtime_error("input row missing");std::istringstream in(line);std::vector<std::string> xs(4);for(auto&s:xs)in>>s;
 I r=I::pow2(re),tol=I::pow2(te),a0=I::decimal(xs[0]),b0=I::decimal(xs[1]),T0=I::decimal(xs[2]);
 std::string metadata="\"schema\":\"R015.validated-quarter.v1\",\"implementation\":\"MPFR-directed-endpoints-second\",\"source_sha256\":\"e93d789a7c1dff2931ff2772754d4072d464704ec35cb9f466274bfc5d2977d3\",\"orbit\":"+std::to_string(row)+",\"precision_bits\":"+std::to_string(PREC)+",\"order\":"+std::to_string(N)+",\"box_radius_power_two\":"+std::to_string(re)+",\"local_remainder_tolerance_power_two\":"+std::to_string(te)+",\"center\":[\""+xs[0]+"\",\""+xs[1]+"\",\""+xs[2]+"\"],\"Tstar\":\""+xs[3]+"\",\"mpfr_version\":\""+mpfr_get_version()+"\",\"gmp_version\":\""+gmp_version+"\"";
 V initial{I(-1),I(),I(1),I(),a0,b0,a0,b0};V c=midpoint(initial),p=symmetric(sub(initial,c)),b=p;for(int i=4;i<8;i++)b[i]=b[i]+sym(r);
 Mat Q=eye(8),R=zeros(8,2);R[4][0]=R[6][0]=I(1);R[5][1]=R[7][1]=I(1);
 I t,tgoal=mid(T0/I(4));mpfr_prec_round(tgoal.l,PREC/2,MPFR_RNDN);mpfr_set(tgoal.u,tgoal.l,MPFR_RNDN);mpfr_prec_round(tgoal.l,PREC,MPFR_RNDN);
 int he=-4,count=0,rejects=0;I d2=I::pow2(10000),maxstate,maxvar,maxha;std::vector<I>steps;std::string reason,status="PARTIAL_INTEGRATION";
 auto progress=[&](){std::ostringstream s;s<<"{"<<metadata<<",\"status\":\"RUNNING\",\"steps_completed\":"<<count<<",\"rejections\":"<<rejects<<",\"elapsed_seconds\":"<<elapsed()<<",\"time\":"<<fraction(t)<<",\"target\":"<<fraction(tgoal)<<",\"next_step_exponent\":"<<he<<",\"distance_squared_lower\":"<<json(d2)<<",\"parameter_radius_max\":"<<json(norm(b))<<",\"center_radius_max\":"<<json(norm(p))<<",\"R_norm\":"<<json(norm(R))<<"}";save(prefix+".progress.json",s.str());};
 try {
  if(!lt(I(),T0-r)||(!lt(r,absup(a0))&&!lt(r,absup(b0))))throw std::runtime_error("initial positive period or nonzero speed unproved");
  while(lt(t,tgoal)){
   if(elapsed()>cap){reason="internal wall-time cap";break;}if(count>=maxsteps){reason="requested step cap";break;}
   I h=I::pow2(he),remaining=tgoal-t;if(!lt(h,remaining))h=remaining;if(!h.point()||mpfr_sgn(h.l)<=0)throw std::runtime_error("step not positive exact dyadic");
   V Y=add(c,mv(Q,b));for(int i=0;i<8;i++)Y[i]=hull(Y[i],add(c,mv(Q,p))[i]);
   try {
    Slab s=slab(Y,h);I E=s.contraction/(I(1)-s.contraction);Mat Mbox=eye(8);for(V&v:Mbox)for(I&z:v)z=z+sym(E);
    Jet jb=jet(s.box,h,N,true);Mat remD=mm(jb.M[N],Mbox);I ns=norm(jb.y[N]),nv=norm(remD);
    if(!lt(ns,tol)||!lt(nv,tol))throw std::runtime_error("Taylor remainder exceeds local tolerance");
    Jet jy=jet(Y,h,N-1,true),jc=jet(c,h,N-1,false);Mat D=zeros(8,8);V C(8);for(int k=0;k<N;k++){D=ma(D,jy.M[k]);C=add(C,jc.y[k]);}D=ma(D,remD);C=add(C,jb.y[N]);
    Mat DQ=mm(D,Q);auto [Qnew,Qinv]=basis(DQ);Mat L=mm(Qinv,DQ);V cnew=midpoint(C),err=mv(Qinv,sub(C,cnew));V bnew=symmetric(add(err,mv(L,b))),pnew=symmetric(add(err,mv(L,p)));for(int i=0;i<8;i++)bnew[i]=hull(bnew[i],pnew[i]);Mat Rnew=mm(L,R);
    I tnew=t+h;if(!tnew.point())throw std::runtime_error("time sum is not exact");
    c=std::move(cnew);Q=std::move(Qnew);b=std::move(bnew);p=std::move(pnew);R=std::move(Rnew);t=tnew;steps.push_back(h);count++;d2=minlower(d2,s.ev.d2);maxstate=maxupper(maxstate,ns);maxvar=maxupper(maxvar,nv);maxha=maxupper(maxha,s.contraction);
    if(count%5==0){progress();std::cout<<"row "<<row<<" step "<<count<<" time "<<std::setprecision(10)<<t.upper()<<"/"<<tgoal.upper()<<" h "<<h.upper()<<" rb "<<norm(b).upper()<<" elapsed "<<elapsed()<<std::endl;}
    if(count%8==0&&he< -2)he++;
   }catch(const std::runtime_error&e){rejects++;he--;reason=e.what();if(he< -100)throw std::runtime_error("step floor: "+reason);}
  }
  progress();
  std::string extra;
  if(!lt(t,tgoal)){
   V Y=add(c,mv(Q,b));I delta0=T0/I(4)-t,Delta=upperpoint(absup(delta0)+r/I(4));Slab se=slab(Y,Delta,true);d2=minlower(d2,se.ev.d2);
   V Fc=add(add(c,mv(Q,p)),scale(se.ev.f,delta0));V F{Fc[4],Fc[5],Fc[6]};Mat W=mm(Q,R);I expansion=se.contraction/(I(1)-se.contraction)*norm(W);Mat We=W;for(V&v:We)for(I&z:v)z=z+sym(expansion);Mat Wend=ma(W,scale(mm(se.ev.A,We),sym(Delta)));Mat J=zeros(3,3);for(int i=0;i<3;i++){J[i][0]=Wend[i+4][0];J[i][1]=Wend[i+4][1];J[i][2]=se.ev.f[i+4]/I(4);}
   Mat C0=inverse_proposal(J),Z=ms(eye(3),mm(C0,J));V K=add(scale(mv(C0,F),I(-1)),mv(Z,V(3,sym(r))));I contraction=norm(Z);V ratios;bool inside=true;for(const I&v:K){ratios.push_back(upperpoint(absup(v)/r));inside&=lt(absup(v),r);}I nondeg=se.box[0]+I(2)*se.box[2];bool pass=inside&&lt(contraction,I(1))&&!det3(C0).zero()&&!nondeg.zero()&&mpfr_sgn(d2.l)>0;
   status=pass?"CERTIFIED":"KRAWCZYK_NOT_PROVED";reason=pass?"all specified sufficient conditions passed":"endpoint internal inclusion, contraction, invertibility or nondegeneracy not proved";
   extra=",\"krawczyk\":{\"F\":"+json(F)+",\"J\":"+json(J)+",\"C\":"+json(C0)+",\"K\":"+json(K)+",\"Z\":"+json(Z)+",\"endpoint\":"+json(se.box)+"},\"inclusion_ratios\":"+json(ratios)+",\"contraction_upper\":"+json(upperpoint(contraction))+",\"nondegeneracy_interval\":"+json(nondeg)+",\"endpoint_delta\":"+json(Delta)+",\"center_time_error\":"+json(delta0);
  }
  std::string ss="[";for(size_t i=0;i<steps.size();i++)ss+=(i?",":"")+fraction(steps[i]);ss+="]";
  std::ostringstream out;out<<"{"<<metadata<<",\"status\":\""<<status<<"\",\"reason\":\""<<reason<<"\",\"elapsed_seconds\":"<<elapsed()<<",\"steps\":"<<ss<<",\"integrated_time\":"<<fraction(t)<<",\"target_time\":"<<fraction(tgoal)<<",\"rejections\":"<<rejects<<",\"distance_squared_lower\":"<<json(d2)<<",\"distance_lower\":"<<json(lowerpoint(sqrt(d2)))<<",\"max_state_remainder\":"<<json(maxstate)<<",\"max_variational_remainder\":"<<json(maxvar)<<",\"max_slab_contraction\":"<<json(maxha)<<",\"endpoint\":{\"c\":"<<json(c)<<",\"Q\":"<<json(Q)<<",\"rp\":"<<json(p)<<",\"rb\":"<<json(b)<<",\"R\":"<<json(R)<<"}"<<extra<<"}";save(prefix+".result.json",out.str());std::cout<<status<<" "<<reason<<" elapsed "<<elapsed()<<std::endl;return status=="CERTIFIED"?0:3;
 }catch(const std::exception&e){std::cerr<<"ERROR "<<e.what()<<std::endl;save(prefix+".error.json","{"+metadata+",\"status\":\"TECHNICAL_ERROR\",\"reason\":\""+e.what()+"\",\"steps_completed\":"+std::to_string(count)+",\"integrated_time\":"+fraction(t)+"}");return 2;}
}
