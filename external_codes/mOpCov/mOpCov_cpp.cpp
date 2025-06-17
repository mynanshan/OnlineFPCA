#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <stdlib.h>
#include <chrono>
using namespace std::chrono;
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


vec K_cos_cpp2(subview_col<double>  s, double t){
    vec k1 = abs(s-t)/2.0;
    vec k2 = (s+t)/2.0;
    vec s1 = pow(k1,4) - 2* pow(k1,3) + pow(k1,2) - 1.0/30.0;
    vec s2 = pow(k2,4) - 2* pow(k2,3) + pow(k2,2) - 1.0/30.0;
    return (-1.0/3.0)* (s1 + s2);
}
// [[Rcpp::export]]
vec K_cos_cpp3(double t, const vec & s){
    vec k1 = abs(s-t)/2.0;
    vec k2 = (s+t)/2.0;
    vec s1 = pow(k1,4) - 2* pow(k1,3) + pow(k1,2) - 1.0/30.0;
    vec s2 = pow(k2,4) - 2* pow(k2,3) + pow(k2,2) - 1.0/30.0;
    return (-1.0/3.0)* (s1 + s2);
}


// [[Rcpp::export]]
vec K_sob_cpp(const vec & s, const vec & t){
    vec k1s = s-0.5;
    vec k1t = t-0.5;
    vec k1abst = abs(s-t)-0.5;
    
    vec k2s = (k1s % k1s - 1.0/12)/2;
    vec k2t = (k1t % k1t - 1.0/12)/2;
    
    vec k4abst = (pow(k1abst,4) - k1abst%k1abst/2 + 7.0/240)/24;
    
    return 1.0 + k1s % k1t + k2s % k2t - k4abst;
}

vec K_sob_cpp2(subview_col<double>  s, double t){
    vec k1s = s-0.5;
    double k1t = t-0.5;
    vec k1abst = abs(s-t)-0.5;
    
    vec k2s = (k1s % k1s - 1.0/12)/2;
    double k2t = (k1t * k1t - 1.0/12)/2;
    
    vec k4abst = (pow(k1abst,4) - k1abst%k1abst/2 + 7.0/240)/24;
    
    return 1.0 + k1s * k1t + k2s * k2t - k4abst;
}

// [[Rcpp::export]]
vec K_sob_cpp3(double t,const vec & s){
    vec k1s = s-0.5;
    double k1t = t-0.5;
    vec k1abst = abs(s-t)-0.5;
    
    vec k2s = (k1s % k1s - 1.0/12)/2;
    double k2t = (k1t * k1t - 1.0/12)/2;
    
    vec k4abst = (pow(k1abst,4) - k1abst%k1abst/2 + 7.0/240)/24;
    
    return 1.0 + k1s * k1t + k2s * k2t - k4abst;
}


// [[Rcpp::export]]
mat getK_cos(const vec & t){
    int n = t.n_elem, i;
    mat out(n,n);
    for (i=0; i<n; i++){
        out.submat(i, i, n-1, i) = K_cos_cpp2(t.subvec(i, n-1), t[i]);
    }
    return symmatl(out);
}

// [[Rcpp::export]]
mat getK_sob(const vec & t){
    int n = t.n_elem, i;
    mat out(n,n);
    for (i=0; i<n; i++){
        out.submat(i, i, n-1, i) = K_sob_cpp2(t.subvec(i, n-1), t[i]);
    }
    return symmatl(out);
}




#define CLICKJ                                                   \
for (itmp = 0; itmp < d; itmp++)                        \
if (iip(itmp) == isr(itmp)-1) iip(itmp) = 0;    \
else {                                                         \
iip(itmp)++;                                           \
break;                                                    \
}                                                              \
for (lj = 0, itmp = 0; itmp < d; itmp++)           \
lj += iip(itmp) * stride(itmp);



//permutation functions for tensors
// [[Rcpp::export]]
vec unfold(const vec &a, const vec &perm, const vec &dims) {
    int d = dims.n_elem;
    int i,itmp;
    vec isr = zeros(d);
    vec pp = zeros(d);
    for (i = 0; i < d; i++) {
        pp(i) = perm(i) - 1;
        isr(i) = dims(pp(i));
    }
    
    int n = a.n_elem;
    vec r = zeros(n);
    
    vec iip = zeros(d);
    for (i = 0; i < d; i++)
    if (perm(i) >= 0 && pp(i) < d) iip(pp(i))++;
    
    vec stride = zeros(d);
    for (iip(0) = 1, i = 1; i<d; i++) iip(i) = iip(i-1) * dims(i-1);
    for (i = 0; i < d; i++) stride(i) = iip(pp(i));
    
    for (i = 0; i < d; iip(i++) = 0);
    
    int li,lj;
    for (lj = 0, li = 0; li < n; li++) {
        r(li) = a(lj);
        CLICKJ;
    }
    return r;
    
}


//khatri-rao
// [[Rcpp::export]]
mat rowkhatrirao(const mat A, const mat B){
    int i;
    int m=A.n_rows, q=A.n_cols*B.n_cols;
    mat C(m,q);
    for (i=0;i<m;i++){
        C.row(i)=kron(A.row(i),B.row(i));
    }
    return C;
}


void Qlossprep(mat &G, vec &h, double *c, List &Xs, List &Ms,vec q,int d, const ivec & include){
    // R, Qv and c are outputs
    double lconst;
    mat Z;
    mat katriM;
    vec cumpq = cumprod(q);
    vec cumsq = cumsum(q);
    G.zeros();
    h.zeros();
    int n = include.n_elem, i, j,l;
    
    for (j = 0; j < include.n_elem; ++j)
    {
        i = include[j] - 1; // p.s. R indexing to C indexing
        int* Mdim = INTEGER(Rf_getAttrib(Ms[i], R_DimSymbol));
        int nX = Rf_length(Xs[i]);
        mat M(REAL(Ms[i]), Mdim[0], Mdim[1], false);
        
        katriM.set_size(nX,prod(q));
        katriM.head_cols(q(0)) = M.head_cols(q(0));
        for(l=1;l<d;l++){
            katriM.head_cols(cumpq(l)) = rowkhatrirao(katriM.head_cols(cumpq(l-1)), M.cols(cumsq(l-1),cumsq(l)-1));
        }
        
        vec x(REAL(Xs[i]), nX, false);
        lconst = 1.0 / (n * nX * (nX - 1.0) );
        mat MM = kron(katriM, katriM);
        G += MM.t() * diagmat((vectorise(ones(nX,nX) - eye(nX, nX))) *lconst) * MM;
        Z = x * x.t();
        h += MM.t() * diagmat((vectorise(ones(nX,nX) - eye(nX, nX))) *lconst)*  vectorise(Z);
        *c += pow(norm(Z - diagmat(Z), "fro"), 2.0) * lconst;
    }
}


// [[Rcpp::export]]
List Qlossprep_cpp(List Xs, List Ms, vec q, const ivec & include){

    int d = q.n_elem;
    int dim = prod(q)*prod(q);
    mat G = zeros(dim, dim);
    vec h = zeros(dim);
    double c=0.0;
    
    Qlossprep(G, h, &c, Xs, Ms,q,d, include);
    
    List out;
    out["G"] = wrap(G);
    out["h"] = wrap(h);
    out["c"] = wrap(c);
    return out;
}


mat vec2m(const vec v, int row){
    int col=v.n_elem/row;
    mat O = reshape(v,row,col);
    return O;
}


void Z0UP(const vec &v,const vec &u0,vec &z0,vec q, const double lam1,const double eta){
    mat Y0=vec2m(v+u0,prod(q));
    vec sv;
    mat U;
    double temp=lam1/eta;
    Y0=0.5*Y0+0.5*Y0.t();
    eig_sym(sv,U,Y0);
    uvec ind = find(sv >= temp);
    vec  tv=temp*ones(ind.n_elem);
    z0=vectorise(U.cols(ind)*diagmat(sv.elem(ind)-tv)*U.cols(ind).t());
}


vec ZUP(const vec &v,int q, const double rho){
    mat Y=vec2m(v,q);
    vec sv;
    mat U;
    mat V;
    svd_econ(U,sv,V,Y);
    uvec ind = find(sv >= rho);
    vec  tv=rho*ones(ind.n_elem);
    return vectorise(U.cols(ind)*diagmat(sv.elem(ind)-tv)*V.cols(ind).t());
}

//objective function:
// [[Rcpp::export]]
double obj(vec &v,double lam1,double lam2,const vec q,const mat &G, const vec &h){
    vec v0;
    int d = q.n_elem;
    vec permind=linspace(d,1,d);
    vec vemp;
    vec dims=join_cols(reverse(q), reverse(q)),  perm=join_cols(permind, permind+d);
    eig_sym(v0,vec2m(v,prod(q)));
    int i;
    double s1=0.0;
    for (i=0;i<d;i++){
        if(i!=0){
            perm.swap_rows(0,i);
        }
        svd(vemp,vec2m(unfold(v,perm,dims),q(i)));
        s1 += sum(vemp);
    }
    double value;
    value=as_scalar(v.t()*G*v-2*v.t()*h)+lam1*sum(v0)+lam2*s1;
    return value;
}


void acadmm(const vec &q,vec &b, vec &z0,vec &z_0, mat &Z, mat &Z_,vec &u0,vec &u_0,mat &U, mat &U_, vec &Obj, const int k, const mat &GU, const vec &vG, mat &A, const mat &G, const vec &h ,double *eta,const double lam1,const double lam2, const double tol){
    int i,j ;
    double r;
    double s;
    int count=0;
    int d=q.n_elem;
    double a = 1.0;
    double a_;
    vec temp;
    vec permind=linspace(d,1,d);
    vec foldb;
    vec dims=join_cols(q, q),  perm=join_cols(permind, permind+d);
    List out;
    
    Obj(0) = obj(z0,lam1,lam2,q,G,h);

    for (i=0;i<k;i++){
        b = z0-u0;
        
        for(j=0;j<d;j++) {
            if(j!=0) { perm.swap_rows(d-j-1,d-j);
                dims.swap_rows(0,j);
            }
            temp = Z.col(j)-U.col(j);
            
            b = b + unfold(temp,perm,dims);
        }
        
        
        b = A * ((*eta) * b + 2*h);
        
        
        dims.head(d) = reverse(q);
        dims.tail(d) = reverse(q);
        perm.head(d) = permind;
        
        
        r=0.0;
        for (j=0; j<d ; j++){
            if(j!=0){perm.swap_rows(0,j);}
            foldb = unfold(b,perm,dims);
            
            Z.col(j) = ZUP(foldb+U.col(j), q(j),lam2/(*eta));
            r+= sum(square(Z.col(j)-foldb));
            U.col(j)=U.col(j) + foldb - Z.col(j);
        }
        
        
        dims.head(d) = q;
        dims.tail(d) = q;
        perm.head(d) = permind;
        
        
        Z0UP(b,u0,z0,q,lam1,(*eta));
       
        
        u0=u0+b-z0;
        
            r = r + sum(square(z0-b));
            r = pow(r,0.5);
            s = (*eta) * pow((accu(square(Z-Z_)) +  sum(square(z0-z_0))),0.5);

      
       
     
            Obj(i+1) = obj(z0,lam1,lam2,q,G,h);
            if(((Obj(i)-Obj(i+1))/abs(Obj(i))<=tol) && (Obj(i)>=Obj(i+1))) { break;}
    
        
        count += 1;
        
        if(i>500 && count>300 && i<k/2){
            if(r>s){*eta = 2.0*(*eta);
                A = GU * diagmat(1.0/(2.0*(vG + 0.5* (q.n_elem+1)*(*eta)*ones(b.n_elem)))) * GU.t();
                u0 = 0.5 * u0;
                u_0 = 0.5 * u_0;
                U = 0.5 * U;
                U_ = 0.5 * U_;
                count = 0;
            }
            if(s>100*r){*eta = (*eta)/2.0;
                A = GU *  diagmat(1.0/(2.0*(vG + 0.5* (q.n_elem+1)*(*eta)*ones(b.n_elem))))  * GU.t();
                u0 = 2.0 * u0;
                u_0 = 2.0 * u_0;
                U = 2.0 * U;
                U_ = 2.0 * U_;
                count = 0;
            
            }
           
        }
        
        a_=1.0 + sqrt(1.0+4.0*pow(a,2))/2.0 ;
        
        z0 = z0 + ((a-1.0)/a_) * (z0 - z_0);
        u0 = u0 + ((a-1.0)/a_) * (u0 - u_0);
        
        
        Z = Z + ((a-1.0)/a_) * (Z - Z_);
        U = U + ((a-1.0)/a_) * (U - U_);
        
        
        z_0 = z0;
        Z_ = Z;
        u_0 = u0;
        U_=U;
        a=a_;
        
    }
    
    b = z0-u0;
    for(j=0;j<d;j++) {
        if(j!=0) { perm.swap_rows(d-j-1,d-j);
            dims.swap_rows(0,j);
        }
        temp = Z.col(j)-U.col(j);
        b = b + unfold(temp,perm,dims);
    }
    b = A * ((*eta) * b + 2*h);
    
    Z0UP(b,u0,z0,q,lam1,(*eta));

}




// [[Rcpp::export]]
List ACADMM(const vec q, const int k, const mat &G, const mat &GU, const vec &vG, const vec &h, vec &b, double eta,const double lam1,const double lam2,const double tol){
    int dim=prod(q) * prod(q);
    List OUT;
    vec z0=1.0*b;
    vec z_0=1.0*z0;
    mat Z = zeros<mat>(dim,2);
    mat Z_ = 1.0*Z;
    mat U= zeros<mat>(dim,2);
    mat U_ = 1.0*U;
    
    vec Obj(k); Obj.zeros();
    vec u0=zeros(dim);
    vec u_0=u0;
    
    mat A=GU * diagmat(1.0/(2.0*(vG + 0.5* (q.n_elem+1)*(eta)*ones(dim))))  * GU.t();
    
    
    acadmm(q,b,z0,z_0,Z,Z_,u0,u_0,U,U_,Obj, k,GU, vG, A, G, h, &eta,lam1,lam2,tol);
    
    OUT["b"]=z0;
    OUT["obj"]=Obj;
    return OUT;
}


// [[Rcpp::export]]
vec mOpCov_cv_cpp(List &Xs, List &Ms, const vec q,ivec traingroup, ivec testgroup, List lambda,
                   vec &b0,  double eta=1e-07,int
                   maxiter=7000, double tol=1e-8){
    int dim = prod(q) * prod(q);
    mat Gtrain(dim, dim), Gtest(dim, dim);
    vec htrain(dim), htest(dim);
    double ctrain=0.0, ctest=0.0;
    List out;
    List lamb;
    double lam,alpha;
    vec z0=1.0*b0;
    vec z_0=1.0*z0;
    mat Z = zeros<mat>(dim,q.n_elem);
    mat Z_ = 1.0*Z;
    mat U= zeros<mat>(dim,q.n_elem);
    mat U_ = 1.0*U;
    vec u0=zeros(dim);
    vec u_0=1.0*u0;
    vec Obj(maxiter);
    Qlossprep(Gtrain, htrain, &ctrain, Xs,Ms,q,q.n_elem, traingroup);
    mat GU;
    vec vG;
    vec b;
    eig_sym(vG, GU, Gtrain);
    mat A=GU * diagmat(1.0/(2.0*(vG + 0.5* (q.n_elem+1)*(eta)*ones(dim))))  * GU.t();

    
    Qlossprep(Gtest, htest, &ctest, Xs,Ms,q,q.n_elem, testgroup);
    
    int nlam = lambda.size();
    int i;
    vec error = zeros(nlam);
    
    for(i=0;i<nlam;i++){
        lamb = lambda[i];
        lam = lamb["lam"];
        alpha = lamb["alpha"];
        if(eta>10*lam) {eta = 10*lam;
           A = GU * diagmat(1.0/(2.0*(vG + 0.5* (q.n_elem+1)*(eta)*ones(dim))))  * GU.t();
        }
        Obj.zeros();

        acadmm(q,b0,z0,z_0,Z,Z_,u0,u_0,U,U_,Obj, maxiter,GU, vG, A, Gtrain, htrain, &eta,lam*alpha,lam*(1.0-alpha)/2.0,tol);
        

        error(i) = as_scalar(z0.t() * Gtest * z0) - 2* sum(z0%htest) + ctest;
    }

    return error;
}

