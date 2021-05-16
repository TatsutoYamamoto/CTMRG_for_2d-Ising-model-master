# include <stdlib.h>
# include <stdio.h>
# include <iostream>
# include <fstream>
# include <math.h>
# include <time.h>
# include <vector>

# include "StdVector"
# include "Eigenvalues"

// This code calulates the thermodynamics proparties of 2-d Ising model (vertex model).

// "CTMRG"
// T.Nishino and K.Okunishi, J. Phys. Soc. Jpn. 65, 4 (1996).
// "2-d vertex model"
// R.J.Baxter, J. Stat. Phys. 19, 5 (1978).

// In this code, we use "Eigen" for linear algebla.
// Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
//
// and you compile this as...
// g++ -O -I (path-to-Eigen) CTMRG.cpp -o (proglam name)

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;

#define J 1.0//結合定数
#define q 2//スピン自由度。イジングであれば+1,-1で q=2。
#define small_beta 0.0003//数値微分(差分)で用いる、微小の逆温度の幅

#define N_max 100//実際に求める系の角転送行列の大きさ
#define T_i 2.0//最初の温度
#define T_f 3.0//最後の温度
#define T_step 0.01//温度上昇のステップ
#define dim 50//imension of CTM
#define g 0;// g=0:free boundary condition , g=1:fixed boundary condition (z=+1)

double calculate(double &);//全体の操作を行う関数。逆温度を与えると、求めたい系のサイズの分配関数の対数を返す
void calculate_thermo(vector<double> &,double &);//calculateに追加して、系の熱力学量をファイルに出力する。

void Initialize(double &, vector<double> &, vector<double> &);//最小の角転送行列,半列転送行列を作る
void expand_C(vector<double> &, vector<double> &, vector<double> &, int &,vector<double> &);//角転送行列の拡大
void expand_P(vector<double> &, vector<double> &, int &);//半列転送行列の拡大
MatrixXd densitymatrix(vector<double> &, int &);//密度行列の構成
void Renor_constract(MatrixXd &, vector<double> &, vector<double> &,int &);//繰り込み操作の行列を構成
void renormalize_C(vector<double> &, vector<double> &,vector<double> &, vector<double> &,int &);//角転送行列の繰り込み
void renormalize_P(vector<double> &, vector<double> &,vector<double> &, vector<double> &,int &);//半列転送行列の繰り込み
double log_partition(vector<double> &, int &, int &);//分配関数の対数を計算
void normalize_C(vector<double> &, int &);//規格化定数を与える。
void normalize_P(vector<double> &, int &);//規格化定数を与える。

double specific_heat(double &,double &,double &,double &);//分配関数の対数から比熱を計算する
double energy_cal(double &,double &,double &);//分配関数からエネルギーの計算
double magnetization_cal(vector<double> &, int &);//中心付近のスピンの期待値を計算する。対称性が保たれる限り、必ず0になる。
double energy_spin_cal(vector<double> &, int &);//分配関数を用いずにエネルギーの計算を行う

double W_B[q*q*q*q];//ボルツマン重み

double D_max[N_max] = {};//角転送行列の成分の最大値を計算するための配列。全て0で初期化
double Q_max[N_max] = {};//半列転送行列の成分の最大値を計算するための配列。全て0で初期化

int main(){
  clock_t start,end;

  double temp = T_i;//温度
  double beta;//逆温度
  cout << "start"<<endl;
  cout << "temp ="<<temp<<endl;


  //ファイルに出力するために使う関数。
  int name = N_max;
  int name_1 = dim;
  int boundary = g;
  char filename[100];
  sprintf(filename, "N=%d_dim=%d_g=%d.txt", name,name_1,boundary);
  ofstream outputfile(filename);

  char initial_condition[100];
  sprintf(initial_condition, "# N=%d, dim=%d, small_beta=%f,g=%d",N_max,dim,small_beta,boundary);
  outputfile << initial_condition <<endl;

  outputfile << "g=0:open boundary condition, g=1:fixed boundary condition (z=+1)" <<endl;

  char value[100];
  sprintf(value, "#temperature\tenergy\t\t\tspecific_heat\t\tmagnetization\t\tenergy_spin");
  outputfile << value <<endl;

  while(temp < T_f){
    start = clock();

    beta = 1/temp;//逆温度
    //差分を用いて数値微分を得るため、三つの温度での値が必要
    double log_Z;//逆温度 beta の分配関数の対数
    double log_Z_up;//逆温度 beta + Δbetaの分配関数の対数
    double log_Z_low;//逆温度 beta - Δbetaの分配関数の対数

    double beta_up = beta + small_beta;
    double beta_low = beta - small_beta;

    double magnetization;//磁化
    double energy_spin;//中心付近のスピン配位から計算したエネルギー

    int number_of_thermo_quantity = 3;
    vector<double> thermo(number_of_thermo_quantity);//熱力学量が色々入ったやつ
    calculate_thermo(thermo,beta);

    log_Z = thermo[0];
    magnetization = thermo[1];
    energy_spin = thermo[2];

    log_Z_up = calculate(beta_up);
    log_Z_low = calculate(beta_low);

    cout << "temperature = " << temp <<endl;

    double energy = 0;//internal energy
    energy = energy_cal(log_Z_up,log_Z_low,beta);
    cout << "enegy = " << energy <<endl;

    double spe_heat = 0;//系のスピン当たりの比熱
    spe_heat = specific_heat(log_Z, log_Z_up, log_Z_low, beta);
    cout << "spe_heat = " << spe_heat <<endl;

    char thermodynamics[1000];
    sprintf(thermodynamics, "%1.6E\t%1.16E\t%1.16E\t%1.16E\t%1.16E"
    ,temp,energy,spe_heat,magnetization,energy_spin);
    outputfile << thermodynamics <<endl;

    temp += T_step;

    end = clock();
    printf("一つの温度領域で、%.2f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
  }


  double cite = 8*pow(N_max,2) + 4*N_max;//系のスピンの数
  cout << "系全体でのスピン数 = "<<cite <<endl;

  return 0;
}

void calculate_thermo(vector<double> & thermo, double & beta){//与えた温度に対して、熱力学量を返す関数。

  int N = 1;//角転送行列の一辺のスピンの数、初めは１
  int dim_C = pow(q,2);//角転送行列の行列の大きさ (行列を配列で保持した場合の、配列の長さ
  int dim_P = pow(q,3);//半列転送行列の行列の大きさ (行列を配列で保持した場合の、配列の長さ

  double magnetization;//磁化
  double energy_spin;//スピン配位から計算したエネルギー

  vector<double> C_old(dim_C);
  vector<double> P_old(dim_P);
  Initialize(beta,C_old,P_old);//ボルツマン重みと、最小の角転送行列,最小の半列転送行列を作る。

  int l = q;//角転送行列の次元。

  while(N < N_max){//角転送行列が望みの大きさになるまで
    //   角転送行列
    //
    //          +
    // C   +--- + j,(0 <= j <=  l-1)
    //     |    +
    //     |
    //    +++ i,(0 <= i <=  l-1)
    //
    //  拡大前の角転送行列の、一辺のブロックスピンを、l次元のベクトルとする。
    //  半列転送行列についても、同様にブロックスピンを、l次元のベクトルとする。
    //
    //  繰り込み操作がなければ「l=2^N」 であり、　繰り込み操作によって「l = dim」 に保たれる。

    //dim_C,dim_P は、拡大前の行列の大きさ (行列を配列で保持した場合の、配列の長さ)
    //繰り込み操作がなければ「dim_C=(2^N)*(2^N)」 であり、　繰り込み操作によって「dim_C = l*l」 に保たれる。


    //半列転送行列の拡大。
    vector<double> P_new(q*q*dim_P);//拡大すると、自由度qのスピンが二個足される
    expand_P(P_new,P_old,l);//Nは、係数の規格化のために必要

    //角転送行列の拡大。
    vector<double> C_new(q*q*dim_C);//拡大すると、自由度qのスピンが二個足される
    expand_C(C_new,C_old,P_old,l,P_new);//ここではまだ規格化されていない。

    //ここで、繰り込みの有無の分岐
    if( q * l > dim){//繰り込む場合
      normalize_C(C_new,N);
      normalize_P(P_new,N);

      //密度行列の構成
      MatrixXd DM(q*l,q*l);// Eigen を用いて、(q*l)×(q*l)の行列を宣言
      DM = densitymatrix(C_new,l);//密度行列の構成

      if( N+1 == N_max){//拡大が完了した場合,ここで熱力学量を計算する
        magnetization = magnetization_cal(C_new,l);
        energy_spin = energy_spin_cal(C_new,l);
      }

      vector<double> Renor(dim*q*l);//繰り込み操作の行列の宣言
      vector<double> Renor_T(dim*q*l);//Renor の転置行列の宣言
      Renor_constract(DM,Renor,Renor_T,l);//繰り込み操作の行列と、その逆行列を構成

      vector<double> C_renor(dim*dim);//繰り込まれた角転送行列の宣言
      renormalize_C(C_renor,C_new,Renor,Renor_T,l);
      vector<double> P_renor(q*dim*dim);//繰り込まれた半列転送行列の宣言
      renormalize_P(P_renor,P_new,Renor,Renor_T,l);

      l = dim;
      dim_C = l * l;
      dim_P = q * dim_C;

      C_old = C_renor;
      P_old = P_renor;
    }
    else{//繰り込まない場合
      normalize_C(C_new,N);
      normalize_P(P_new,N);

      l = q * l;
      dim_C = l * l;
      dim_P = q * dim_C;

      C_old = C_new;
      P_old = P_new;
    }

    //一辺の長さを追加。
    N += 1;
  }

  //分配関数を計算する関数。
  double log_Z = 0;//分配関数の対数

  if(l < dim){//繰り込みを一度も行っていない場合
    int size = pow(2,N_max);
    log_Z = log_partition(C_old,N,size);
  }
  else{//繰り込みを行った場合
    int size = dim;
    log_Z = log_partition(C_old,N,size);
  }

  thermo[0] = log_Z;
  thermo[1] = magnetization;
  thermo[2] = energy_spin;
}

double calculate(double & beta){//与えた温度に対して、分配関数の対数を返す関数。

  int N = 1;//角転送行列の一辺のスピンの数、初めは１
  int dim_C = pow(q,2);//角転送行列の行列の大きさ (行列を配列で保持した場合の、配列の長さ
  int dim_P = pow(q,3);//半列転送行列の行列の大きさ (行列を配列で保持した場合の、配列の長さ

  vector<double> C_old(dim_C);
  vector<double> P_old(dim_P);
  Initialize(beta,C_old,P_old);//ボルツマン重みと、最小の角転送行列,最小の半列転送行列を作る。

  int l = q;//角転送行列の次元。

  while(N < N_max){//角転送行列が望みの大きさになるまで
    //   角転送行列
    //
    //          +
    // C   +--- + j,(0 <= j <=  l-1)
    //     |    +
    //     |
    //    +++ i,(0 <= i <=  l-1)
    //
    //  拡大前の角転送行列の、一辺のブロックスピンを、l次元のベクトルとする。
    //  半列転送行列についても、同様にブロックスピンを、l次元のベクトルとする。
    //
    //  繰り込み操作がなければ「l=2^N」 であり、　繰り込み操作によって「l = dim」 に保たれる。

    //dim_C,dim_P は、拡大前の行列の大きさ (行列を配列で保持した場合の、配列の長さ)
    //繰り込み操作がなければ「dim_C=(2^N)*(2^N)」 であり、　繰り込み操作によって「dim_C = l*l」 に保たれる。


    //半列転送行列の拡大。
    vector<double> P_new(q*q*dim_P);//拡大すると、自由度qのスピンが二個足される
    expand_P(P_new,P_old,l);//Nは、係数の規格化のために必要

    //角転送行列の拡大。
    vector<double> C_new(q*q*dim_C);//拡大すると、自由度qのスピンが二個足される
    expand_C(C_new,C_old,P_old,l,P_new);

    //ここで、繰り込みの有無の分岐
    if( q * l > dim){//繰り込む場合
      normalize_C(C_new,N);
      normalize_P(P_new,N);

      //密度行列の構成
      MatrixXd DM(q*l,q*l);// Eigen を用いて、(q*l)×(q*l)の行列を宣言
      DM = densitymatrix(C_new,l);//密度行列の構成

      vector<double> Renor(dim*q*l);//繰り込み操作の行列の宣言
      vector<double> Renor_T(dim*q*l);//Renor の転置行列の宣言
      Renor_constract(DM,Renor,Renor_T,l);//繰り込み操作の行列と、その逆行列を構成

      vector<double> C_renor(dim*dim);//繰り込まれた角転送行列の宣言
      renormalize_C(C_renor,C_new,Renor,Renor_T,l);
      vector<double> P_renor(q*dim*dim);//繰り込まれた半列転送行列の宣言
      renormalize_P(P_renor,P_new,Renor,Renor_T,l);

      l = dim;
      dim_C = l * l;
      dim_P = q * dim_C;

      //新しい関数への書き換え。型さえ合っていれば、サイズは別々でも大丈夫。
      C_old = C_renor;
      P_old = P_renor;
    }
    else{//繰り込まない場合
      normalize_C(C_new,N);
      normalize_P(P_new,N);

      l = q * l;
      dim_C = l * l;
      dim_P = q * dim_C;

      //新しい関数への書き換え。型さえ合っていれば、サイズは別々でも大丈夫。
      C_old = C_new;
      P_old = P_new;
    }

    //一辺の長さを追加。
    N += 1;
  }


  //分配関数を計算する関数。
  double log_Z = 0;//分配関数の対数

  if(l < dim){//繰り込みを一度も行っていない場合
    int size = pow(2,N_max);
    log_Z = log_partition(C_old,N,size);
  }
  else{//繰り込みを行った場合
    int size = dim;
    log_Z = log_partition(C_old,N,size);
  }

  return log_Z;
}


void Initialize(double & beta, vector<double> & C_old, vector<double> & P_old){//固定端境界条件の実装
  /*
  4次元配列を1次元に直す
  W_B[z1][z2][z3][z4]　→　W_B[q*q*q*z1 + q*q*z2 + q*z3 + z4]
  W_B[0]        W_B[1]
      z1          -1          -1
      |            |           |
  z2- * - z4   -1 -*- -1   -1 -*- +1
      |            |           |
      z3          -1           -1
  */
  /*
  2次元配列を1次元に直す
  C_1[z3][z4]　→　C_1[q*z3 + z4]
  C_1[0]        C_1[1]
     #            #           #
     |            |           |
  #- * - z4   #  -*- -1    # -*- +1
     |            |           |
     z3          -1           -1
  */
  /*
  3次元配列を1次元に直す
  P_1[z3][z2][z4] → P_1[q*q*z3 + q*z2 + z4]
  P_1[0]        P_1[1]
      #            #           #
      |            |           |
  z2- * - z4   -1 -*- -1   -1 -*- +1
      |            |           |
      z3          -1           -1
  */

int boundary = g;

if(boundary == 0){//自由端境界条件の場合
  for (int s4=0; s4<q; s4++)//s_1はスピン変数。Isingの場合はq=2のため、s_1 = 0 or 1のみ。
  {
    double z4 = 2*s4 - 1;

    for (int s3=0; s3<q; s3++)
    {
      double z3 = 2*s3 - 1;
      double sum_C=0;//角転送行列の中身C_1[z3,z4]

      for (int s2=0; s2<q; s2++)
      {
        double z2 = 2*s2 - 1;
        double sum_P=0;//半列転送行列の中身

        for (int s1=0; s1<q; s1++)
        {
          double z1 = 2*s1 - 1;

          W_B[q*q*q*s1 + q*q*s2 + q*s3 + s4] = exp(J*beta*( z1*z2 + z2*z3 + z3*z4 + z4*z1 ));
          sum_P +=  W_B[q*q*q*s1 + q*q*s2 + q*s3 + s4] ;
          sum_C +=  W_B[q*q*q*s1 + q*q*s2 + q*s3 + s4] ;
        }
        //半列転送行列
        P_old[q*q*s3 + q*s2 + s4] = sum_P;
      }

      //Cは転送行列のため、forの中でsumをとってそれを代入する
      C_old[q*s3 + s4] = sum_C;
    }
  }
}//自由端境界条件の場合
else if(boundary == 1){//固定端境界条件(z=+1)の場合
  for (int s4=0; s4<q; s4++){//s_1はスピン変数。Isingの場合はq=2のため、s_1 = 0 or 1のみ。
    double z4 = 2*s4 - 1;

    for (int s3=0; s3<q; s3++){
      double z3 = 2*s3 - 1;
      double sum_C=0;//角転送行列の中身C_1[z3,z4]

      for (int s2=0; s2<q; s2++){
        double z2 = 2*s2 - 1;

        for (int s1=0; s1<q; s1++){
          double z1 = 2*s1 - 1;
          W_B[q*q*q*s1 + q*q*s2 + q*s3 + s4] = exp(J*beta*( z1*z2 + z2*z3 + z3*z4 + z4*z1 ));
        }

        //半列転送行列
        P_old[q*q*s3 + q*s2 + s4] = W_B[q*q*q*1 + q*q*s2 + q*s3 + s4];
      }
      //Cは転送行列のため、forの中でsumをとってそれを代入する
      C_old[q*s3 + s4] = W_B[q*q*q*1 + q*q*1 + q*s3 + s4];
    }
  }
}//固定端境界条件の場合
else{
  cout << "boundary condition error" << endl;
  exit(1);
}

    //最小の大きさの角転送行列の最大の成分を保管する。(D_max,Q_maxは0で初期化済み)
    double D_max_val = 0;
    for (int i = 0; i < C_old.size(); i++) {
      if(D_max_val < C_old[i]){
        D_max_val = C_old[i];
      }
    }
    D_max[0] = D_max_val;

    //最小の大きさの半列転送行列の最大の成分を保管する。
    double Q_max_val = 0;
    for (int i = 0; i < P_old.size(); i++) {
      if(Q_max_val < P_old[i]){
        Q_max_val = P_old[i];
      }
    }
    Q_max[0] = Q_max_val;

    //ここで、角転送行列と半列転送行列を規格化しておく。
    for (int i = 0; i < C_old.size(); i++) {
      C_old[i] = C_old[i]/D_max_val;
    }
    for (int i = 0; i < P_old.size(); i++) {
      P_old[i] = P_old[i]/Q_max_val;
    }

}//Initializeここまで


void expand_C(vector<double> & C_new, vector<double> & C_old, vector<double> & P_old, int & l, vector<double> & P_new){
  //拡大前の角転送行列の、一辺のブロックスピンを、l次元のベクトルとする。

  vector<double> Theta(q*l*l);//計算の補助。C_old と P_old を先に和をとっておく。
  // Theta(sigma_c|eta_2|i) = sum_{j} P_N(sigma_c|eta_2|j) * C_N(j|i)

  // Theta[sigma_c][eta_2][i]
  //
  //     sigma_c
  //     +        +++ i
  //     |         |
  // +   |   +  +  |
  // + - + - +  +--+
  // +       +  +
  // eta_2    j
  //

 for (int sigma_c = 0; sigma_c < q; sigma_c++) {
   for (int eta_2 = 0; eta_2 < l; eta_2++) {
     for (int i = 0; i < l; i++) {

       double sum_theta = 0;
       for (int j = 0; j < l; j++) {
         sum_theta += P_old[l*l*sigma_c + l*eta_2 + j] * C_old[l*j + i];
       }
       Theta[l*l*sigma_c + l*eta_2 + i] = sum_theta;

     }
   }
 }
 //  C_new[eta_1][sigma_1][eta_2][sigma_2]  (sigma_1,sigma_2=0,1 / 0<i,j<l)
 //  これを一次元に直すと　C[ q*q*k*eta_1 + q*k*sigma_1 + q*eta_2 + sigma_2 ]
 //
 //    P_new[sigma_2][eta_1][sigma_1][i][sigma_c]
 //
 //          sigma_1    eta_1
 //            +        +++
 //            |         |
 // sigma_2 +--*---------*
 //            |         |
 //            +        +++
 //           sigma_c   i
 //            +        +++
 //            |         |
 //         +  |         |
 //  eta_2  + -*-- ------+
 //         +
 //
 //  Theta[sigma_c][eta_2][i]
 //
 //

 for (int eta_1 = 0; eta_1 < l; eta_1++) {
   for (int sigma_1 = 0; sigma_1 < q; sigma_1++) {
     for (int eta_2 = 0; eta_2 < l; eta_2++) {
       for (int sigma_2 = 0; sigma_2 < q; sigma_2++) {
         double sum_C = 0;

         for (int sigma_c = 0; sigma_c < q; sigma_c++) {
           for (int i = 0; i < l; i++) {
             sum_C += Theta[l*l*sigma_c + l*eta_2 + i]
             * P_new[q*q*l*l*sigma_2 + q*q*l*eta_1+ q*l*sigma_1 + q*i + sigma_c];
           }
         }

         C_new[q*q*l*eta_1 + q*l*sigma_1 + q*eta_2 + sigma_2] = sum_C;
       }
     }
   }
 }
}//expand_C ここまで

void expand_P(vector<double> & P_new, vector<double> & P_old, int & l){
// <P_1の場合>
// 3次元配列を1次元に直す
// P_1[z3][z2][z4] → P_1[q*q*z3 + q*z2 + z4]
// P_1[0]        P_1[1]
//     #            #           #
//     |            |           |
// z2- * - z4   -1 -*- -1   -1 -*- +1
//     |            |           |
//     z3          -1           -1
//
// P_N+1
//
//  j
// +++         +　sigma_2
//  |   delta  |
//  *---*   *--*-- + psi
//  |          |
// +++         + sigma_1
//  i
//    psi = 0,1 / i,j = 0 ~ 2^N
// P_new[psi][i][sigma_1][j][sigma_2] = sum { P_old[delta][i][j] * W[psi][sigma_2][delta][sigma_1] }
//
// k = 2^Nとして
//   P(2*k*kの三次元テンソル) →　成分は P[s][i][j] (s=0,1 / 0<i,j<k)
//  これを一次元の成分に直すと P[ k*k*s + k*i + j ]
//


//       P_new[psi][i][sigma_1][j][sigma_2] → P_new[q*q*k*k*psi + q*q*k*i + q*k*sigma_1 + q*j + sigma_2]
//拡大前の角転送行列の、一辺のブロックスピンを、l次元のベクトルとする。

    for (int psi = 0; psi < q; psi++) {
      for (int i = 0; i < l; i++) {
        for (int sigma_1 = 0; sigma_1 < q; sigma_1++) {
          for (int j = 0; j < l; j++) {
            for (int sigma_2 = 0; sigma_2 < q; sigma_2++) {
              double sum_P = 0;

              for (int delta = 0; delta < q; delta++) {
                sum_P += P_old[l*l*delta + l*i +j] * W_B[q*q*q*psi + q*q*sigma_2 + q*delta + sigma_1];
              }

              P_new[q*q*l*l*psi + q*q*l*i + q*l*sigma_1 + q*j + sigma_2] = sum_P;
            }
          }
        }
      }
    }//for文の中身ここまで

 }//Expand_Pここまで

 MatrixXd densitymatrix(vector<double> & C_new, int & l){
   //角転送行列の方で規格化をしているため、密度行列の規格化は考えなくて良い。
   //拡大前の角転送行列の、一辺のブロックスピンを、l次元のベクトルとする。
   // C_new[theta][eta] → C_new[k*theta + eta]
   //          +
   //     +--- + eta,(0 <= eta <=  k-1)
   //     |    +
   //     |
   //    +++ theta,(0 <= theta <=  k-1)
   //
   //
   // DM[i][j] → DM[k*i + j]
   //          +     +
   //     +--- + eta + ---+
   //     |    +     +    |
   //     |               |
   //    +++             +++
   //    theta           sigma
   //    +++             +++
   //     |               |
   //     |    +    +     |
   //     +----+    +-----+
   //          +    +
   //          i    j
   //
   //Eigen を用いて、直接行列として取る。

   int k = q * l;//拡大後の角転送行列のブロックスピンを、k次元のベクトルとする。
   MatrixXd DM(k,k);// Eigen を用いて、l*l の行列を宣言

   vector<double> C_2(k*k);// Squared C, for auxiliary

   for (int i = 0; i < k; i++) {
     for (int j = 0; j < k; j++) {
       double sum_C2 = 0;
       for (int theta = 0; theta < k; theta++) {
         sum_C2 += C_new[k*i+theta] * C_new[k*theta + j];
       }
       C_2[k*i + j] = sum_C2;
     }
   }

   for (int i = 0; i < k; i++) {
     for (int j = 0; j < k; j++) {
       double sum_DM = 0;
         for (int theta = 0; theta < k; theta++) {
           sum_DM += C_2[k*i + theta] * C_2[k*theta + j];
         }
       DM(i,j) = sum_DM;
     }
   }
   return DM;
  }

 void Renor_constract(MatrixXd & DM, vector<double> & Renor, vector<double> & Renor_T,int & l){
   // R[i][j]:繰り込みで用いる行列
   //
   //  R = [v_1,v_2,...,v_dim] (v_iは、l次元の固有ベクトル)
   //    →　R は l行 dim列 の行列(縦長の行列)
   //  R[i][j] →　R[ i*dim + j ]   ( 0 < i < l-1, 0 < j < dim-1 )
   //
     // step
     // . DMの固有ベクトルを求める
     // . 固有値の大きい順に、DMの固有ベクトルをdim個だけ取って、それを配列Renor[]に代入していく

     //DMの固有値と固有ベクトルをえる。
     SelfAdjointEigenSolver<MatrixXd> es(DM);//対称行列と信じれば、全ての固有値は実数
     // es.eigenvectors() は、固有値の小さい順に固有ベクトル(縦ベクトル)を並べた直交行列を生成している
     // これを後ろの列から順に成分を拾って行って、dim個のベクトルだけで Renor を構成する

     int k = q * l;//繰り込み前の角転送行列のサイズ

     for (int i=0; i<k; i++){
       for (int j=0; j<dim; j++){//dim は#define で与えている
           Renor[i*dim + j] = es.eigenvectors()(i, (k-1)-j);
         }
      }

     //転置行列の構成
      for (int i=0; i<dim; i++){
         for (int j=0; j<k; j++){
           Renor_T[i*k + j] = Renor[i+j*dim];
         }
       }

 }


 void renormalize_C(vector<double> & C_renor, vector<double> & C_new, vector<double> & Renor, vector<double> & Renor_T,int & l){
   //  拡大前の角転送行列の、一辺のブロックスピンを、l次元のベクトルとする。
   //  繰り込み行列(q*l)×(dim) R と 角転送行列(q*l)×(q*l) C_new を用いて
   //  C_renor = Renor_T C_new Renor
   //  を計算する。
   //  C_renor[i][j] = Renor_T[i][eta]*C_new[eta][theta]*Renor[theta][j] (0<=i,j<dim, 0<= eta,theta < l)
   //
   //  Renor[theta][j] → Renor[theta*dim + j]
   //  Renor_T[i][eta] → Renor_T[i*l + eta]

     int k = q * l;//繰り込み前の角転送行列のサイズ
     vector<double> CR(dim*k);//計算量の都合上、C_new と R で先に積をとる

     for (int eta = 0; eta < k; eta++) {
       for (int j = 0; j < dim; j++) {
         double sum_CR = 0;
         for (int theta = 0; theta < k; theta++) {
           sum_CR += C_new[k*eta + theta] * Renor[theta*dim + j];
         }
         CR[dim*eta +j] = sum_CR;
       }
     }

     for (int i = 0; i < dim; i++) {
       for (int j = 0; j < dim; j++) {
         double sum_RCR = 0;
           for (int eta = 0; eta < k; eta++) {
             sum_RCR += Renor_T[k*i + eta]*CR[dim*eta +j];
           }
         C_renor[i * dim + j] = sum_RCR;
       }
     }

 }//renormalize_Cの末尾

 void renormalize_P(vector<double> & P_renor, vector<double> & P_new, vector<double> & Renor, vector<double> & Renor_T,int & l){
   //  拡大前の角転送行列の、一辺のブロックスピンを、l次元のベクトルとする。
   //  繰り込み行列(q*l)×(dim) R と 半列転送行列 q*(q*l)×(q*l) P_new を用いて
   //  P_renor = Renor_T P_new Renor
   //  を計算する。
   //  P_renor[sigma][i][j] = Renor_T[i][eta]*P_new[sigma][eta][theta]*Renor[theta][j] (0<=i,j<dim, 0<= eta,theta < l, sigma=0,1)
   //
   //  Renor[theta][j] → Renor[theta*dim + j]
   //  Renor_T[i][eta] → Renor_T[i*l + eta]
   //  P_renor[sigma][i][j] → P_renor[sigma*dim*dim + i*dim + j]

   int k = q * l;//繰り込み前の角転送行列のサイズ
   vector<double> PR(q*dim*k);//計算量の都合上、P_new と Renor で先に積をとる

   for (int sigma = 0; sigma < q; sigma++) {
     for (int xi = 0; xi < k; xi++) {
       for (int j = 0; j < dim; j++) {
         double sum_PR = 0;
         for (int eta = 0; eta < k; eta++) {
           sum_PR += P_new[k*k*sigma + k*xi + eta] * Renor[dim*eta + j];
         }
         PR[k*dim*sigma + dim*xi + j] = sum_PR;
       }
     }
   }

   for (int sigma = 0; sigma < q; sigma++) {
     for (int i = 0; i < dim; i++) {
       for (int j = 0; j < dim; j++) {
         double sum = 0;
         for (int xi = 0; xi < k; xi++) {
           sum += Renor_T[i*k + xi] * PR[k*dim*sigma + dim*xi + j];
         }
         P_renor[dim*dim*sigma + i*dim + j] = sum;
       }
     }
   }
 }//renormalize_Pの末尾

void normalize_C(vector<double> & C_new,int & N){
  //この関数内の C_new は、C_new と C_renor のどちらもありうる。
   double C_new_max = 0;
   int size = C_new.size();
   for (int i = 0; i < size; i++) {
     if( C_new_max < C_new[i]){
        C_new_max = C_new[i];
     }
   }
   D_max[N] = C_new_max;

   for (int i = 0; i < size; i++) {
     C_new[i] = C_new[i]/C_new_max;
   }
}

void normalize_P(vector<double> & P_new, int & N){
  //この関数内の P_new は、P_new と P_renor のどちらもありうる。

    double P_new_max = 0;
    int size = P_new.size();
     for (int i = 0; i < size; i++) {
       if( P_new_max < P_new[i]){
          P_new_max = P_new[i];
       }
     }
     Q_max[N] = P_new_max;

     for (int i = 0; i < size; i++) {
       P_new[i] = P_new[i]/P_new_max;
     }
  }

double log_partition(vector<double> & C_old, int & N, int & size){//分配関数の対数を返す。
/*<分配関数の規格化>
構成する時の角転送行列の規格化定数 C_max を用いると
Z_nol = sum C * C * C * C
であり
log Z = log Z_nol + 4*log(C_max)
と書ける。

ここで、log(C_max)は
log(C_max) = sum_i {log(Di_max)} + 2*sum_j {(N-j)*log(Qj_max)}
である。
*/

  double Z_nol = 0;//規格化された分配関数
  double log_Z = 0;

  //size は、角転送行列のブロックスピンの次元。
  //繰り込みを行っていれば必ず size = dim で、繰り込みがなければ size = 2^N 。

  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      for(int l=0; l<size; l++){
        for(int n=0; n<size; n++){
          Z_nol += C_old[size*i +j] * C_old[size*j + l] * C_old[size*l + n] * C_old[size*n +i];
        }
      }
    }
  }

  if(Z_nol < 0){
    cout<<"error, Z_nol<0"<<endl;
  }

  double sum_C = 0;
  double sum_P = 0;

  for (int i = 0; i < N_max ; i++) {//配列 D_max[]の大きさの分だけ足す。
    sum_C += log(D_max[i]);
  }
  for (int i = 0; i < N_max-1 ; i++) {//配列 Q_max[]の大きさより「1つ少なく(!)」足す。半列転送行列は、最後の一回は無駄に構成しているため。
    sum_P += (N_max-(i+1))*log(Q_max[i]);
  }

  log_Z = log(Z_nol) + 4*sum_C + 8*sum_P;

  return log_Z;
}

double specific_heat(double & log_Z,double & log_Z_up,double & log_Z_low,double & beta){
  //中点差分を用いたlog_Zの二階微分
  double spe_heat = 0;//ボルツマン定数で規格化されている。
  double cite = 8*pow(N_max,2) + 4*N_max;//系のスピンの数
  spe_heat =  ( (log_Z_up - 2*log_Z + log_Z_low)*beta*beta / (small_beta*small_beta) )/ cite;

  return spe_heat;
}

double energy_cal(double & log_Z_up,double & log_Z_low,double & beta){
  double energy = 0;
  double cite = 8*pow(N_max,2) + 4*N_max;//系のスピンの数
  energy = -((log_Z_up - log_Z_low)/(2*small_beta))/cite;

  return energy;
}

double magnetization_cal(vector<double> & C_new,int & l){//自由端条件だと、必ず0になる。
  double magnetization = 0;

  // lは、拡大前のCTMの次元
  int k = q * l;//dimension of CTM
  vector<double> C2(k*k);// squared CTM
  vector<double> C4(k*k);// 4th power of CTM

  for (int i = 0; i < k; i++){
    for(int j=0; j<k; j++){
      double sum_C2 = 0;
      for (int theta = 0; theta < k; theta++) {
        sum_C2 += C_new[k*i + theta] * C_new[k*theta +j];
      }
      C2[k*i+j] = sum_C2;
    }
  }

  for (int i = 0; i < k; i++){
    for(int j=0; j<k; j++){
      double sum_C4 = 0;
      for (int theta = 0; theta < k; theta++) {
        sum_C4 += C2[k*i + theta] * C2[k*theta +j];
      }
      C4[k*i+j] = sum_C4;
    }
  }

  for (int sigma = 0; sigma < q; sigma++) {
    int z = 2*sigma - 1;
    for (int xi = 0; xi < l; xi++) {
        magnetization += z * C4[q*q*l*xi + q*l*sigma + q*xi + sigma];
    }
  }

  return magnetization;
}//中心付近のスピンの期待値を計算する

double energy_spin_cal(vector<double> & C_new,int & l){
  vector<double> U(q*q);//中心のスピン自由度以外を潰したテンソル

  // lは、拡大前のCTMの次元
  int k = q * l;//dimension of CTM
  vector<double> C2(k*k);// squared CTM
  vector<double> C3(k*k);// cubed CTM

  for (int i = 0; i < k; i++){
    for(int j=0; j<k; j++){
      double sum_C2 = 0;
      for (int theta = 0; theta < k; theta++) {
        sum_C2 += C_new[k*i + theta] * C_new[k*theta +j];
      }
      C2[k*i+j] = sum_C2;
    }
  }

  for (int i = 0; i < k; i++){
    for(int j=0; j < k; j++){
      double sum_C3 = 0;
      for (int theta = 0; theta < k; theta++) {
        sum_C3 += C2[k*i + theta] * C_new[k*theta +j];
      }
      C3[k*i+j] = sum_C3;
    }
  }

  for (int sigma_1 = 0; sigma_1 < q; sigma_1++) {
    for (int sigma_2 = 0; sigma_2 < q; sigma_2++) {
      double sum_U = 0;
      for (int xi = 0; xi < l; xi++) {
        for (int mu = 0; mu < l; mu++) {
          sum_U += C_new[q*q*l*xi + q*l*sigma_1 + q*mu + sigma_2] * C3[q*q*l*mu + q*l*sigma_2 + q*xi + sigma_1];
        }
      }
      U[q*sigma_1 + sigma_2] = sum_U;
    }
  }

  double uZ = 0;
  for (int sigma_1 = 0; sigma_1 < q; sigma_1++) {
    for (int sigma_2 = 0; sigma_2 < q; sigma_2++) {
        uZ += U[q*sigma_1 + sigma_2];
    }
  }

  double energy = 0;

  for (int sigma_1 = 0; sigma_1 < q; sigma_1++){
    for (int sigma_2 = 0; sigma_2 < q; sigma_2++){
      int z1 = 2*sigma_1 - 1;
      int z2 = 2*sigma_2 - 1;
      energy += -J * z1 * z2 * U[q*sigma_1 + sigma_2];
    }
  }

  energy = 2*energy / uZ;

  return energy;
}
