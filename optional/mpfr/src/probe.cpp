#include "dynamics.hpp"
#include <fstream>
int main(){PREC=256;
 for(int k=0;k<120;k++){I a=I((k*37)%83-41)/I(8),b=I((k*19)%61-30)/I(16);a=hull(a,a+I((k%7)+1)/I(32));b=hull(b,b+I((k%5)+1)/I(32));std::cout<<"{\"kind\":\"arithmetic\",\"a\":"<<json(a)<<",\"b\":"<<json(b)<<",\"sum\":"<<json(a+b)<<",\"difference\":"<<json(a-b)<<",\"product\":"<<json(a*b)<<",\"square\":"<<json(sqr(a));if(!b.zero())std::cout<<",\"quotient\":"<<json(a/b);I pos=sqr(a)+I(1)/I(10);std::cout<<",\"positive\":"<<json(pos)<<",\"sqrt\":"<<json(sqrt(pos))<<"}\n";}
 int negative_tests=0;try{inv(sym(I(1)));}catch(const std::runtime_error&){negative_tests++;}try{sqrt(I(-1));}catch(const std::runtime_error&){negative_tests++;}
 V collision(8);try{eval(collision);}catch(const std::runtime_error&){negative_tests++;}
 std::ifstream f("studies/R004/r015-second/euler_100_4860.txt");std::string line;for(int i=0;i<4;i++)std::getline(f,line);std::istringstream ss(line);std::string a,b;ss>>a>>b;V y{I(-1),I(),I(1),I(),I::decimal(a),I::decimal(b),I::decimal(a),I::decimal(b)};I h=I::pow2(-7);auto j=jet(y,h,12,true);auto e=eval(y);
 std::cout<<"{\"kind\":\"jets\",\"negative_tests\":"<<negative_tests<<",\"initial\":"<<json(y)<<",\"h\":"<<json(h)<<",\"f\":"<<json(e.f)<<",\"A\":"<<json(e.A)<<",\"coefficients\":[";for(int k=0;k<=12;k++)std::cout<<(k?",":"")<<json(j.y[k]);std::cout<<"],\"variational_coefficients\":[";for(int k=0;k<=12;k++)std::cout<<(k?",":"")<<json(j.M[k]);std::cout<<"]}\n";
}
