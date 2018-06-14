//#include <stdlib.h>
//#include <stdio.h>
//#include <iostream>
//#include <math.h>
//#include <time.h>
//#include <vector>
//#include <algorithm>
//
//using namespace std;
//FILE  *datastream, *data_2;
//
//#define L 300
//#define N (L*L)
//
//
//int main(){
//  int seed=123456,i=0,n=0;
//
//  double  du=1.0,p=0.0, final=1.0,previous=1.0,alal=0.0,fog=0.0;
//    double  xini=0.0, xf=0.0,yini=0.0,asus=0.0,dell=0.0,samsung=0.0;
//    seed=time(NULL);
//    srand(seed);
//
//    datastream=fopen("L_0_entropy_slope(s_300_en_100000)_con.txt","w");
//    data_2=fopen("L_0_entropy_slope(s_300_en_100000).txt","r");
//    vector<double > microcanonical(N+2,0.0);
//
//for(i=1;i<=N;i++)
//   {
//     fscanf(data_2,"%lf     %lf\n",&xini,&xf);
//
//    microcanonical[i]=xf;
//   }
//
//  double binom[N+1];
//    for (int j = 1; j <=N; j++)
//    {
//        int n, i;
//        double prob = j*1.0/N;
//        binom[j] = 1;
//	//double prob = 0.526846;
//
//        for (n = j+1; n <= N; ++n)
//            binom[n] = binom[n-1]*(N-n+1)*1.0/n*prob/(1-prob);
//
//        for (n = j - 1; n >= 0; --n)
//            binom[n] = binom[n+1]*(n+1)*1.0/(N-n)*(1-prob)/prob;
//
//        double sum = 0;
//        for (i = 0; i <= N; ++i) sum += binom[i];
//
//        for (i = 0; i <= N; ++i) binom[i] /= sum;
//
//        sum = 0;
//        for (n = 1; n <= N; ++n) sum += microcanonical[n]*binom[n];
//        //printf("%12d%10.5f%18.8e%18.8e\n", j, prob, microcanonical[j]*1.0, sum);
//       fprintf(datastream," %lf        %lf        %lf\n",prob,microcanonical[j]*1.0, sum);
//     }
//
//
//    fclose(datastream);
//    fclose(data_2);
//    return 0;
//}
