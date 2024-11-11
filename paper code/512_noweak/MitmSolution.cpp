#include "MitmSolution.h"

MitmSolution::MitmSolution(int Rounds) {
    DA.resize(Rounds + 1, vector<vector<vector<vector<int>>>>(5, vector<vector<vector<int>>>(5, vector<vector<int>>(64, vector<int>(3)))));
    DB.resize(Rounds, vector<vector<vector<vector<int>>>>(5, vector<vector<vector<int>>>(5, vector<vector<int>>(64, vector<int>(3)))));
    DC1_red.resize(Rounds, vector<vector<int>>(5, vector<int>(64)));
    DC1_blue.resize(Rounds, vector<vector<int>>(5, vector<int>(64)));
    DC12_red.resize(Rounds, vector<vector<int>>(5, vector<int>(64)));
    DC12_blue.resize(Rounds, vector<vector<int>>(5, vector<int>(64)));
    DP.resize(Rounds, vector<vector<vector<int>>>(5, vector<vector<int>>(64, vector<int>(3))));
    DP2.resize(Rounds, vector<vector<vector<int>>>(5, vector<vector<int>>(64, vector<int>(3))));
    DC2_red.resize(Rounds, vector<vector<vector<int>>>(5, vector<vector<int>>(5, vector<int>(64))));
    DC2_blue.resize(Rounds, vector<vector<vector<int>>>(5, vector<vector<int>>(5, vector<int>(64))));
    beta1.resize(Rounds, vector<vector<int>>(5, vector<int>(64)));
    beta2.resize(Rounds, vector<vector<vector<int>>>(5, vector<vector<int>>(5, vector<int>(64))));
    beta3.resize(Rounds, vector<vector<int>>(5, vector<int>(64)));
    /*dom.resize(2, vector<int>(64));*/
  //  dom1.resize(2, vector<int>(64));
    dom2.resize(5, vector<int>(64));
    obj.resize(3);
}
