#include "dynamics.hpp"
#include <fstream>
#include <iomanip>
int main(int argc,char**argv){PREC=256;int row=argc>1?std::stoi(argv[1]):4;std::ifstream f("studies/R004/r015-second/euler_100_4860.txt");std::string line;for(int i=0;i<row;i++)std::getline(f,line);std::istringstream ss(line);std::string a,b;ss>>a>>b;V y{I(-1),I(),I(1),I(),I::decimal(a),I::decimal(b),I::decimal(a),I::decimal(b)};
 for(int he=-4;he>=-17;he--){try{I h=I::pow2(he);Slab s=slab(y,h);Jet j=jet(s.box,h,24,true);I E=s.contraction/(I(1)-s.contraction);Mat M=eye(8);for(V&v:M)for(I&z:v)z=z+sym(E);std::cout<<he<<" ha="<<s.contraction.upper()<<" rem="<<norm(j.y[24]).upper()<<","<<norm(mm(j.M[24],M)).upper()<<"\n";if(he==-8){for(int k=0;k<=24;k++)std::cout<<"coeff "<<k<<" "<<norm(j.y[k]).upper()<<" "<<norm(j.M[k]).upper()<<"\n";}}catch(const std::exception&e){std::cout<<he<<" "<<e.what()<<"\n";}}
}
