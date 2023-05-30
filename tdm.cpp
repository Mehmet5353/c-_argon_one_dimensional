#include<bits/stdc++.h>

void tdm(double A[], double B[], double C[], double F[], double result[], int n){
    C[0]=C[0]/B[0];
    F[0]=F[0]/B[0];
    for(int i=1; i<n-1; i++){
        double dummy= B[i]-A[i]*C[i-1];
        C[i]=C[i]/dummy;
        F[i]=(F[i]+A[i]*F[i-1])/dummy; 
    }

    F[n-1]=(F[n-1]+A[n-1]*F[n-2])/(B[n-1]-A[n-1]*C[n-2]);
    result[n-1]=F[n-1];
    for(int i=n-2; i>-1; i--){
        result[i]=F[i]+C[i]*result[i+1];
    }
}

int main(){
    auto start = std::chrono::high_resolution_clock::now();
    double beta=2e-13,  e=1.6e-19, eps0=8.854e-12,  gamma=0.07, L=0.01, p=1,    mue=30.0/p,  mui=0.8/p,  Te=1.0,   Ti=0.025,   Ud=400.0, De=mue*Te,  Di=mui*Ti,  C=0.1; double err=3.0;
    int N=100;
    double h=L/(N-1);
    double x[N];
    double Phi[N],  Aphi[N], Bphi[N], Cphi[N], Fphi[N];
    for(int i=0; i<N; i++){
        Phi[i]=0;
        Aphi[i]=0;
        Bphi[i]=0;
        Cphi[i]=0;
        Fphi[i]=0;
    }
    Bphi[0]=1;    Fphi[0]=Ud;    Bphi[N-1]=1;    Fphi[N-1]=0;
    double fe[N-1], fi[N-1],  E[N-1],  EE[N-1];  
    for(int i=0; i<N-1; i++){
        fe[i]=0;
        fi[i]=0;
        E[i]=0;
        EE[i]=0;
    }

    for(int i=0; i<N; i++){
        x[i]=i*h;
    }
    double ni[N],   ne[N];
    for(int i=0; i<N; i++){
        ni[i]=-1e20*(x[i]*(x[i]-L));
    }
    for(int i=0; i<N; i++){
        ne[i]=-1e20*(x[i]*(x[i]-L));
    }
    for(int i=1; i<N-1; i++){
        Aphi[i]=1/(h*h);
        Cphi[i]=1/(h*h);
        Bphi[i]=2/(h*h);
        Fphi[i]=-e/eps0*(ni[i]-ne[i]);
        //Fphi[i]=1.0;
    }

    double Pe[N],  Pi[N],  Ape[N],  Api[N];  
    for(int i=0; i<N-1; i++){
        Pe[i]=0;
        Pi[i]=0;
        Ape[i]=0;
        Api[i]=0;
    }
    double Ae[N],    Be[N],    Ce[N],    Fe[N];
    for(int i=0; i<N-1; i++){
        Ae[i]=0;
        Be[i]=0;
        Ce[i]=0;
        Fe[i]=0;
    }   
    double Ai[N],    Bi[N],    Ci[N],    Fi[N],   fluxe[N],    fluxi[N];
    for(int i=0; i<N-1; i++){
        Ai[i]=0;
        Bi[i]=0;
        Ci[i]=0;
        Fi[i]=0;
        fluxe[i]=1.0;
        fluxi[i]=1.0;
    } 
    Bi[0]=1.0;    Ai[N-1]=1.0;  Bi[N-1]=1.0;  double tau=2e-15;   Be[0]=1.0;    Ce[0]=1.0;    Be[N-1]=1.0;  double time=0.0; 

    tdm(Aphi, Bphi, Cphi, Fphi, Phi, N);
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
    while(it<1e7){
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

    tdm(Ai, Bi, Ci, Fi, ni, N);


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

    tdm(Ae, Be, Ce, Fe, ne, N);  
    double tempEE=*std::max_element(EE,EE+N);
    tau=C*h/(abs(tempEE)*mue);
    time+=tau;
    }

    std::cout<<time<<std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout<<duration.count()<<std::endl;



    return 0;
}