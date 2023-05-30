#include<bits/stdc++.h>

void TDMA(std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, std::vector<double>& F, std::vector<double>& result){
    size_t n=F.size();

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