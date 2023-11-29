#ifndef ANA_FT
#define ANA_FT
#include <algorithm>
#include <sstream>
#include <vector>

// #include <stdio.h>
// #include <stdlib.h>

// #define folder_name "h0"//"stwmssnpn1"
// #define v0 1.0

//R,lo,ms,tauの順;
// double tau,Ms,R,Rbit=0,lo;
using std::cout;
using std::endl;
using std::pair;
using std::vector;

bool input(vector<vector<double>> &inpp,char *foldername_1) {
  char filename[128];
  snprintf(filename, 128, "%s/fais_R%.3f.dat",
           foldername_1, R);
  cout<<filename<<endl;
  std::ifstream file;
  file.open(filename);
  if (!file) {
    cout << "file is not found" << endl;
    return false;
  }
  std::string str;
  double inp_fl[4];
  vector<double> inp(5);
  while (getline(file, str)) {
    // if(str[0]=='#')continue;
    // cout<<str<<endl;
    std::stringstream ss;
    ss << str;
    ss >> inp[0] >> inp[1] >> inp[2] >> inp[3] >> inp[4];
    inpp.push_back(inp);
    // cout << inp[0] << endl;
  }
  file.close();
  return true;
}
void dft(vector<vector<double>> &inp, vector<vector<double>> &result) {
  double bun = 1. / inp.size();
  int size_inp = inp.size(), size_res = inp[0].size();
  vector<pair<double, double>> res_i(size_res);  // 0がsin、1がcos;
  vector<double> res_i_abs(size_res);
  double bit_q = 2 * M_PI / size_inp,
         bit_qans = 2 * M_PI / inp[size_inp - 1][0];
  pair<double, double> sincos;
  for (int k = 0; k < size_inp; k++) {
    for (auto &&ve : res_i) {
      ve.first = 0.;
      ve.second = 0.;
    }
    result[k][0] = bit_qans * k;
    // res_i[0].second = bit_qans * k;
    for (int j = 0; j < size_inp; j++) {
      sincos.first = sin(result[k][0] * inp[j][0]) * bun;
      sincos.second = cos(result[k][0] * inp[j][0]) * bun;
      for (int l = 1; l < size_res; l++) {
        res_i[l].first += inp[j][l] * sincos.first;
        res_i[l].second += inp[j][l] * sincos.second;
      }
    }
    for (int i = 1; i < size_res; i++) {
      result[k][i] = sqrt(res_i[i].first * res_i[i].first +
                          res_i[i].second * res_i[i].second);
    }
    // result[k]=res_i;
    // cout << k << endl;
  }
}
void dft_req(vector<vector<double>> &inp, vector<vector<double>> &result) {
  double bun = 1. / inp.size();
  int size_inp = inp.size(), size_res = inp[0].size();
  vector<pair<double, double>> res_i(size_res);  // 0がsin、1がcos;
  vector<double> res_i_abs(size_res);
  double ramda_min = inp[1][0] - inp[0][0],
         ramda_bit = inp[size_inp - 1][0] / size_inp, q;
  pair<double, double> sincos;
  for (int k = 0; k < size_inp; k++) {
    for (auto &&ve : res_i) {
      ve.first = 0.;
      ve.second = 0.;
    }
    result[k][0] = ramda_min + ramda_bit * k;
    q = 2 * M_PI / result[k][0];
    for (int j = 0; j < size_inp; j++) {
      sincos.first = sin(q * inp[j][0]) * bun;
      sincos.second = cos(q * inp[j][0]) * bun;
      for (int l = 1; l < size_res; l++) {
        res_i[l].first += inp[j][l] * sincos.first;
        res_i[l].second += inp[j][l] * sincos.second;
      }
    }
    for (int i = 1; i < size_res; i++) {
      result[k][i] = sqrt(res_i[i].first * res_i[i].first +
                          res_i[i].second * res_i[i].second);
    }
    // cout << k << endl;
  }
}
void output_ft(vector<vector<double>> result,char *folder_name_1) {
  char filename[128];
  snprintf(filename, 128,
           "%s/fais_qs_R%.3f.dat",
           folder_name_1, R);
  std::ofstream ofs;
  ofs.open(filename);
  ofs << "#波長　その他;";
  for (auto &&ve : result) {
    // ofs<<2*M_PI/ve[0]<<"\t";
    for (int i = 0; i < ve.size(); i++) {
      ofs << ve[i] << "\t";
    }
    ofs << endl;
  }
  ofs.close();
  
}
void output_ft_req(vector<vector<double>> result,char *foldername_1) {
  char filename[128];
  snprintf(filename, 128,"%s/fais_qs_r_R%.3f.dat", foldername_1, R);
  std::ofstream ofs;
  ofs.open(filename);
  ofs << "#波長　その他;";
  for (auto &&ve : result) {
    // ofs<<2*M_PI/ve[0]<<"\t";
    for (int i = 0; i < ve.size(); i++) {
      ofs << ve[i] << "\t";
    }
    ofs << endl;
  }
  ofs.close();
}
bool do_fts(char *folder_1){
  vector<vector<double>> inp;
  if (!input(inp,folder_1)) {
    return false;
  }
  // sort(inp.begin(), inp.end(),
  //      [](const vector<double> &alpha, const vector<double> &beta) {
  //        return alpha[0] < beta[0];
  //      });
  vector<vector<double>> result(
      inp.size(), vector<double>(inp[0].size()));
  dft(inp, result);
  cout << "No" << endl;
  output_ft(result,folder_1);
  dft_req(inp,result);
  output_ft_req(result,folder_1);
  return true;
}
#endif
