
//This class will the class of functions from R3 or R2 to the same space that
// are polynomial by coordinates along with all the interactions
//that can occur between them like the dot product.
#include "Polynomial.hh"


class PolynomialByCoordinates { public:
       
    PolynomialByCoordinates(int dim,int deg){dim_ = dim;
                                             deg_ = deg;
                                             pol1_.Reshape(dim_,deg_,true);
                                             pol2_.Reshape(dim_,deg_,true);
                                             pol3_.Reshape(dim_,deg_,true);
                                             };
    PolynomialByCoordinates(int dim,int deg, 
                            Amanzi::WhetStone::Polynomial& pol1,
                            Amanzi::WhetStone::Polynomial& pol2,
                            Amanzi::WhetStone::Polynomial& pol3)
                            {dim_ = dim, deg_ = deg;
                             pol1_ = pol1;
                             pol2_ = pol2;
                             pol3_ = pol3;
                            };                                 
    ~PolynomialByCoordinates(){};
    void SetDimDeg(int dim,int deg){dim_=dim, deg_=deg;
                                    pol1_.Reshape(dim_,deg_,true);
                                    pol2_.Reshape(dim_,deg_,true);
                                    pol3_.Reshape(dim_,deg_,true);
                                    };
    void SetPol1(Amanzi::WhetStone::Polynomial pol1){pol1_=pol1;};
    void SetPol2(Amanzi::WhetStone::Polynomial pol2){pol2_=pol2;};
    void SetPol3(Amanzi::WhetStone::Polynomial pol3){pol3_=pol3;};
    //Accessors
    int GetDim(){return dim_;};
    int GetDeg(){return deg_;};
    Amanzi::WhetStone::Polynomial GetPol1(){return pol1_;};
    Amanzi::WhetStone::Polynomial GetPol2(){return pol2_;};
    Amanzi::WhetStone::Polynomial GetPol3(){return pol3_;};
    
    //Dot Product between Polynomials by coodinates
    Amanzi::WhetStone::Polynomial dot(PolynomialByCoordinates b){
        Amanzi::WhetStone::Polynomial apol1 = this->GetPol1();
        Amanzi::WhetStone::Polynomial apol2 = this->GetPol2();
        Amanzi::WhetStone::Polynomial apol3 = this->GetPol3();
        Amanzi::WhetStone::Polynomial bpol1 = b.GetPol1();
        Amanzi::WhetStone::Polynomial bpol2 = b.GetPol2();
        Amanzi::WhetStone::Polynomial bpol3 = b.GetPol3();

        return (apol1*bpol1)+(apol2*bpol2)+(apol3*bpol3);
    }
    //The degree of the polynomials
    int deg_;
    //The dimensionality of the mapping (2D or 3D) 
    int dim_;
    //The lists of coefficients of the polynomials.
    //The maximum degree of the polynomials supported is 2
    //Thus we have 
    Amanzi::WhetStone::Polynomial pol1_;
    Amanzi::WhetStone::Polynomial pol2_;
    Amanzi::WhetStone::Polynomial pol3_;
};