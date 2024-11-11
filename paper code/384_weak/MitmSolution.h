#ifndef MITMSOLUTION_H
#define MITMSOLUTION_H
#pragma once
#include <vector>

using namespace std;

class MitmSolution {
public:
    int Rounds;
    int objective;
    vector<vector<vector<vector<vector<int>>>>> DA;          //DA对应每一轮的状态A
    vector<vector<vector<vector<int>>>> DP;                  //DP对应每一轮的状态C
    vector<vector<vector<int>>> DC1_red;                         //DC1对应A->C消耗的自由度
    vector<vector<vector<int>>> DC1_blue;
    vector<vector<vector<vector<int>>>> DP2;
    vector<vector<vector<int>>> DC12_red;
    vector<vector<vector<int>>> DC12_blue;
    vector<vector<vector<vector<vector<int>>>>> DB;
    vector<vector<vector<vector<int>>>> DC2_red;
    vector<vector<vector<vector<int>>>> DC2_blue;
    vector<int> obj;
   // vector<vector<int>> dom1;
    vector<vector<int>> dom2;
    vector<vector<vector<int>>> beta1;
    vector<vector<vector<vector<int>>>> beta2;
    vector<vector<vector<int>>> beta3;

    MitmSolution(int Rounds);

};

#endif