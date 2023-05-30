#include<bits/stdc++.h>
#include"TDMA.cpp"



int main(){
    auto start = std::chrono::high_resolution_clock::now();
    double beta=2e-13,  e=1.6e-19, eps0=8.854e-12,  gamma=0.07, L=0.01, p=1,    mue=30.0/p,  mui=0.8/p,  Te=1.0,   Ti=0.025,   Ud=400.0, De=mue*Te,  Di=mui*Ti,  C=0.1; double err=3.0;
    int N=100;
    double h=L/(N-1);
    std::vector<double>x; std::vector<double>Phi;
    Phi.assign(N, 0.0);
    for(int i=0; i<N; i++){
        double temp=0+i*h;
        x.push_back(temp);
    }
    std::vector<double>ni;  std::vector<double>ne;
    ni.assign(N,0.0);   ne.assign(N,0.0);
    for(int i=0; i<N; i++){
        ni[i]=-1e20*(x[i]*(x[i]-L));
    }
    for(int i=0; i<N; i++){
        ne[i]=-1e20*(x[i]*(x[i]-L));
    }
    std::vector<double>fe;  std::vector<double>fi;  std::vector<double>E;   std::vector<double>EE;  
    fe.assign(N-1,0);   fi.assign(N-1,0);   E.assign(N-1,0);    EE.assign(N-1,0);   

    std::vector<double>Pe;  std::vector<double>Pi;  std::vector<double>Ape;   std::vector<double>Api;  
    Pe.assign(N,0);   Pi.assign(N,0);   Ape.assign(N,0);    Api.assign(N,0);

    std::vector<double>Aphi;    std::vector<double>Bphi;    std::vector<double>Cphi;    std::vector<double>Fphi;
    Aphi.assign(N, 0.0);    Bphi.assign(N, 0.0);    Cphi.assign(N, 0.0);    Fphi.assign(N, 0.0);

    Bphi[0]=1;    Fphi[0]=Ud;    Bphi[N-1]=1;    Fphi[N-1]=0;
        
    std::vector<double>Ae;    std::vector<double>Be;    std::vector<double>Ce;    std::vector<double>Fe;    
    Ae.assign(N, 0.0);    Be.assign(N, 0.0);    Ce.assign(N, 0.0);    Fe.assign(N, 0.0);

    std::vector<double>Ai;    std::vector<double>Bi;    std::vector<double>Ci;    std::vector<double>Fi;   std::vector<double>fluxe;    std::vector<double>fluxi;
    Ai.assign(N, 0.0);    Bi.assign(N, 0.0);    Ci.assign(N, 0.0);    Fi.assign(N, 0.0);    fluxe.assign(N,1.0);  fluxi.assign(N,1.0);
    Bi[0]=1.0;    Ai[N-1]=1.0;  Bi[N-1]=1.0;  double tau=2e-15;   Be[0]=1.0;    Ce[0]=1.0;    Be[N-1]=1.0;  double time=0.0;



    std::vector<double>tempn;

    for(int i=1; i<N-1; i++){
        Aphi[i]=1/(h*h);
        Cphi[i]=1/(h*h);
        Bphi[i]=2/(h*h);
        Fphi[i]=-e/eps0*(ni[i]-ne[i]);
        //Fphi[i]=1.0;
    }
    TDMA(Aphi, Bphi, Cphi, Fphi, Phi);
    Fe[N-1]=(ni[N-1]*gamma*mui)/mue;

    for(int i=1; i<N-1; i++){
        EE[i]=(-Phi[i+1]+Phi[i-1])/(2*h);
    }
    for(int i=0; i<N-1; i++){
        E[i]=(-Phi[i+1]+Phi[i])/(h);
    }
    for(int i=0; i<N-1; i++){
        fe[i]=-mue*E[i];
        fi[i]=mui*E[i];
    }
    for(int i=0; i<N-1; i++){
        Pe[i]=fe[i]*h/De;
        Pi[i]=fi[i]*h/Di;
    }
    for(int i=0; i<N-1; i++){
        if(std::abs(Pe[i])<1e-5){
            Ape[i]=1;
        }else{
            Ape[i]=std::abs(Pe[i])/(std::exp(std::abs(Pe[i]))-1);
        }
    }
    for(int i=0; i<N-1; i++){
        if(std::abs(Pi[i])<1e-5){
            Api[i]=1;
        }else{
            Api[i]=std::abs(Pi[i])/(std::exp(std::abs(Pi[i]))-1);
        }
    }



    
    int it=0;
    while(it<1e5){
    //tempn.assign(ni.begin(),ni.end());
    it+=1;
    double a=0.0;
    for(int i=1; i<N-1; i++){
        fluxe[i]=-De*(ne[i+1]-ne[i])/h-ne[i]*EE[i]*mue;
        a+=fluxe[i];
    }

    for(int i=1; i<N-1; i++){
            fluxe[i]=a/(N-2); 
    }

    double b=0.0;
    for(int i=0; i<N-1; i++){
            fluxi[i]=-Di*(ni[i+1]-ni[i])/h+ni[i]*EE[i]*mui;
            b+=fluxi[i];
    }

    for(int i=0; i<N-1; i++){
            fluxi[i]=b/(N-1); 
    }

    
    for(int i=1; i<N-1; i++){
            Ai[i]= Di*Api[i-1]/h+ std::max(fi[i-1],0.0);
            Ci[i]= Di*Api[i]/h+ std::max(-fi[i],0.0);
            Bi[i]= Ai[i]+ Ci[i]+ (fi[i]-fi[i-1])+ h/tau;
            Fi[i]= std::abs(fluxe[i])*1200*p*std::exp(-18e3*p/std::abs(EE[i]))*h- beta*ni[i]*ne[i]*h+ ni[i]*h/tau;
    }

    TDMA(Ai, Bi, Ci, Fi, ni);


    //double tempy=0.0;

    //for(int i=1; i<N-1; i++){
            //tempy+=std::sqrt(1.0/N*std::pow((tempn[i]-ni[i]),2)/std::abs(ni[i]));
    //}
    //std::cout<<tempy<<std::endl;


    for(int i=1; i<N-1; i++){
            Ae[i]= De*Ape[i-1]/h+ std::max(fe[i-1], 0.0);
            Ce[i]= De*Ape[i]/h+ std::max(-fe[i], 0.0);
            Be[i]= Ae[i]+Ce[i]+(fe[i]-fe[i+1])+ h/tau;
            Fe[i]= std::abs(fluxe[i])*1200*p*std::exp(-18e3*p/std::abs(EE[i]))*h- beta*ni[i]*ne[i]*h+ ne[i]*h/tau;
    }

    TDMA(Ae, Be, Ce, Fe, ne);  
    double tempEE=*std::max_element(EE.begin(),EE.end());
    tau=C*h/(abs(tempEE)*mue);
    time+=tau;
    }
    std::cout << std::setprecision(10);
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout<<duration.count()<<std::endl;

    //for(int i=0; i<ni.size(); i++){
    //   std::cout<<ni[i]<<std::endl;
    //}
    std::cout<<time<<std::endl;


    return 0;
}