#pragma once
#include "interval.hpp"
struct Eval{V f;Mat A;I d2;};
inline Eval eval(const V&y){
 V dx{y[2]-y[0],-I(2)*y[0]-y[2],-y[0]-I(2)*y[2]};
 V dy{y[3]-y[1],-I(2)*y[1]-y[3],-y[1]-I(2)*y[3]};
 V gx(3),gy(3);std::vector<Mat>B(3,zeros(2,2));I d2=I::pow2(10000);
 for(int p=0;p<3;p++){I rho=sqr(dx[p])+sqr(dy[p]);if(mpfr_sgn(rho.l)<=0)throw std::runtime_error("distance box contains collision");d2=minlower(d2,rho);I z=inv(rho*sqrt(rho)),w=z/rho;gx[p]=dx[p]*z;gy[p]=dy[p]*z;V d{dx[p],dy[p]};for(int i=0;i<2;i++)for(int j=0;j<2;j++)B[p][i][j]=(i==j?z:I())-I(3)*w*(i==j?sqr(d[i]):d[i]*d[j]);}
 V f{y[4],y[5],y[6],y[7],gx[0]+gx[1],gy[0]+gy[1],-gx[0]+gx[2],-gy[0]+gy[2]};Mat A=zeros(8,8);for(int i=0;i<4;i++)A[i][i+4]=I(1);
 for(int i=0;i<2;i++)for(int j=0;j<2;j++){A[i+4][j]=-B[0][i][j]-I(2)*B[1][i][j];A[i+4][j+2]=B[0][i][j]-B[1][i][j];A[i+6][j]=B[0][i][j]-B[2][i][j];A[i+6][j+2]=-B[0][i][j]-I(2)*B[2][i][j];}
 return {f,A,d2};
}
// Coefficients in u=t/h, so every kth coefficient already contains h^k.
struct Jet{std::vector<V> y;std::vector<Mat>M;};
inline Jet jet(const V&init,const I&h,int N,bool variational){
 std::vector<V> y(N+1,V(8));y[0]=init;
 std::vector<Mat>H(N,zeros(4,4));
 std::vector<V> dx(3,V(N)),dy=dx,rho=dx,z=dx,w=dx;
 for(int k=0;k<N;k++){
  V gx(3),gy(3);std::vector<Mat>bk(3,zeros(2,2));
  for(int p=0;p<3;p++){
   dx[p][k]=p==0?y[k][2]-y[k][0]:(p==1?-I(2)*y[k][0]-y[k][2]:-y[k][0]-I(2)*y[k][2]);
   dy[p][k]=p==0?y[k][3]-y[k][1]:(p==1?-I(2)*y[k][1]-y[k][3]:-y[k][1]-I(2)*y[k][3]);
   for(int j=0;j<=k;j++)rho[p][k]=rho[p][k]+(j==k-j?sqr(dx[p][j])+sqr(dy[p][j]):dx[p][j]*dx[p][k-j]+dy[p][j]*dy[p][k-j]);
   if(k==0){if(mpfr_sgn(rho[p][0].l)<=0)throw std::runtime_error("jet collision");z[p][0]=inv(rho[p][0]*sqrt(rho[p][0]));if(variational)w[p][0]=z[p][0]/rho[p][0];}
   else{I sz,sw;for(int j=1;j<=k;j++){sz=sz+I(-j-2*k)*rho[p][j]*z[p][k-j]/I(2);if(variational)sw=sw+I(-3*j-2*k)*rho[p][j]*w[p][k-j]/I(2);}z[p][k]=sz/(I(k)*rho[p][0]);if(variational)w[p][k]=sw/(I(k)*rho[p][0]);}
   for(int j=0;j<=k;j++){gx[p]=gx[p]+dx[p][j]*z[p][k-j];gy[p]=gy[p]+dy[p][j]*z[p][k-j];}
   if(variational){for(int i=0;i<2;i++)for(int j=i;j<2;j++){I t;const V&di=i==0?dx[p]:dy[p];const V&dj=j==0?dx[p]:dy[p];for(int l=0;l<=k;l++){I dd;for(int m=0;m<=l;m++)dd=dd+(i==j&&m==l-m?sqr(di[m]):di[m]*dj[l-m]);t=t+dd*w[p][k-l];}bk[p][i][j]=bk[p][j][i]=(i==j?z[p][k]:I())-I(3)*t;}}
  }
  V fk{y[k][4],y[k][5],y[k][6],y[k][7],gx[0]+gx[1],gy[0]+gy[1],-gx[0]+gx[2],-gy[0]+gy[2]};y[k+1]=scale(fk,h/I(k+1));
  if(variational)for(int i=0;i<2;i++)for(int j=0;j<2;j++){H[k][i][j]=-bk[0][i][j]-I(2)*bk[1][i][j];H[k][i][j+2]=bk[0][i][j]-bk[1][i][j];H[k][i+2][j]=bk[0][i][j]-bk[2][i][j];H[k][i+2][j+2]=-bk[0][i][j]-I(2)*bk[2][i][j];}
 }
 std::vector<Mat>M;
 if(variational){M.assign(N+1,zeros(8,8));M[0]=eye(8);for(int k=0;k<N;k++){for(int i=0;i<4;i++)for(int j=0;j<8;j++){M[k+1][i][j]=M[k][i+4][j]*h/I(k+1);I s;for(int l=0;l<=k;l++)for(int m=0;m<4;m++)s=s+H[l][i][m]*M[k-l][m][j];M[k+1][i+4][j]=s*h/I(k+1);}}}
 return {std::move(y),std::move(M)};
}
struct Slab{V box;Eval ev;I contraction;};
inline bool contains(const V&big,const V&small){for(size_t i=0;i<big.size();i++)if(!subset(small[i],big[i]))return false;return true;}
inline Slab slab(const V&Y,const I&h,bool symmetric_time=false){
 I time=h;if(symmetric_time)time=sym(h);else mpfr_set_zero(time.l,0);
 V B=Y;Eval ev=eval(Y);
 for(int tries=0;tries<14;tries++){
  V next=add(Y,scale(ev.f,time));
  if(contains(B,next)){I a=norm(ev.A),ha=absup(h)*a;if(!lt(ha,I(1)/I(4)))throw std::runtime_error("Picard contraction >= 1/4");return {B,ev,ha};}
  for(int i=0;i<8;i++){I z=hull(Y[i],next[i]);I m=mid(z);B[i]=m+sym(absup(z-m)*I(9)/I(8)+I::pow2(-PREC+16));}
  ev=eval(B);
 }
 throw std::runtime_error("Picard inclusion not obtained");
}
// Modified Gram-Schmidt only proposes a point basis. Orthogonality is NOT assumed.
inline std::pair<Mat,Mat> basis(const Mat&input){
 Mat a=midpoint(input),q=zeros(8,8);
 for(int j=0;j<8;j++){V v(8);for(int i=0;i<8;i++)v[i]=a[i][j];for(int pass=0;pass<2;pass++)for(int k=0;k<j;k++){I dot;for(int i=0;i<8;i++)dot=mid(dot+q[i][k]*v[i]);for(int i=0;i<8;i++)v[i]=mid(v[i]-dot*q[i][k]);}I n;for(const I&z:v)n=mid(n+sqr(z));n=mid(sqrt(n));if(n.zero())throw std::runtime_error("QR proposal singular");for(int i=0;i<8;i++)q[i][j]=mid(v[i]/n);}
 Mat approx=zeros(8,8);for(int i=0;i<8;i++)for(int j=0;j<8;j++)approx[i][j]=q[j][i];
 I defect=norm(ms(eye(8),mm(approx,q)));if(!lt(defect,I(1)))throw std::runtime_error("basis inverse not proved");
 I error=upperpoint(defect/(I(1)-defect)*norm(approx));Mat qi=approx;for(V&v:qi)for(I&z:v)z=z+sym(error);
 return {q,qi};
}
inline Mat inverse_proposal(const Mat&input){int n=input.size();Mat a=midpoint(input),b=eye(n);for(int k=0;k<n;k++){int pivot=k;for(int i=k+1;i<n;i++)if(lt(absup(a[pivot][k]),absup(a[i][k])))pivot=i;std::swap(a[k],a[pivot]);std::swap(b[k],b[pivot]);I d=a[k][k];if(d.zero())throw std::runtime_error("inverse proposal singular pivot");for(int j=0;j<n;j++){a[k][j]=mid(a[k][j]/d);b[k][j]=mid(b[k][j]/d);}for(int i=0;i<n;i++)if(i!=k){I f=a[i][k];for(int j=0;j<n;j++){a[i][j]=mid(a[i][j]-f*a[k][j]);b[i][j]=mid(b[i][j]-f*b[k][j]);}}}return b;}
inline I det3(const Mat&a){return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);}
