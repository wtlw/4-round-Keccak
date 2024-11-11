#include "MitmSHA3.h"

#define TESTOBJ 20
#define solfile "model1.sol"

Mitm::Mitm(int Rounds) {
    this->Rounds = Rounds;

    DA.resize(Rounds + 1, vector<vector<vector<vector<GRBVar>>>>(5, vector<vector<vector<GRBVar>>>(5, vector<vector<GRBVar>>(64, vector<GRBVar>(3)))));

    DB.resize(Rounds, vector<vector<vector<vector<GRBVar>>>>(5, vector<vector<vector<GRBVar>>>(5, vector<vector<GRBVar>>(64, vector<GRBVar>(3)))));

    DC1_red.resize(Rounds, vector<vector<GRBVar>>(5, vector<GRBVar>(64)));
    DC1_blue.resize(Rounds, vector<vector<GRBVar>>(5, vector<GRBVar>(64)));

    DC12_red.resize(Rounds, vector<vector<GRBVar>>(5, vector<GRBVar>(64)));
    DC12_blue.resize(Rounds, vector<vector<GRBVar>>(5, vector<GRBVar>(64)));

    DP.resize(Rounds, vector<vector<vector<GRBVar>>>(5, vector<vector<GRBVar>>(64, vector<GRBVar>(3))));

    DP2.resize(Rounds, vector<vector<vector<GRBVar>>>(5, vector<vector<GRBVar>>(64, vector<GRBVar>(3))));

    DC2_red.resize(Rounds, vector<vector<vector<GRBVar>>>(5, vector<vector<GRBVar>>(5, vector<GRBVar>(64))));
    DC2_blue.resize(Rounds, vector<vector<vector<GRBVar>>>(5, vector<vector<GRBVar>>(5, vector<GRBVar>(64))));


    dom1.resize(4, vector<GRBVar>(64));


    beta1.resize(Rounds, vector<vector<GRBVar>>(5, vector<GRBVar>(64)));
    beta2.resize(Rounds, vector<vector<vector<GRBVar>>>(5, vector<vector<GRBVar>>(5, vector<GRBVar>(64))));
    beta3.resize(Rounds, vector<vector<GRBVar>>(5, vector<GRBVar>(64)));

    obj.resize(3);
}

void Mitm::Value(GRBModel& model) {
    MitmFactory factory;

    for (int round = 0; round < Rounds + 1; round++)
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                for (int k = 0; k < 64; k++)
                    for (int l = 0; l < 3; l++) {
                        DA[round][i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "DA_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" + std::to_string(l));
                    }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                for (int k = 0; k < 64; k++)
                    for (int l = 0; l < 3; l++) {
                        DB[round][i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "DB_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" + std::to_string(l));
                    }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++)
                for (int l = 0; l < 3; l++) {
                    DP[round][i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "DP_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(l));
                }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++)
                for (int l = 0; l < 3; l++) {
                    DP2[round][i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "DP2_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(l));
                }

     
    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                for (int k = 0; k < 64; k++) {
                    DC2_red[round][i][j][k] = model.addVar(0.0, 2.0, 0.0, GRB_INTEGER, "DC2_red_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));
                }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                for (int k = 0; k < 64; k++) {
                    DC2_blue[round][i][j][k] = model.addVar(0.0, 2.0, 0.0, GRB_INTEGER, "DC2_blue_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));

                }
    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                DC1_red[round][i][k] = model.addVar(0.0, 2.0, 0.0, GRB_INTEGER, "DC1_red_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k));
            }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                DC1_blue[round][i][k] = model.addVar(0.0, 2.0, 0.0, GRB_INTEGER, "DC1_blue_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k));
            }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                DC12_red[round][i][k] = model.addVar(0.0, 2.0, 0.0, GRB_INTEGER, "DC12_red_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k));
            }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                DC12_blue[round][i][k] = model.addVar(0.0, 2.0, 0.0, GRB_INTEGER, "DC12_blue_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k));
            }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                beta1[round][i][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "beta1_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k));
            }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++)
                for (int j = 0; j < 5; j++) {
                    beta2[round][i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "beta2_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));
                }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                beta3[round][i][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "beta3_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k));
            }

    for (int i = 0; i < 4; i++)
        for (int k = 0; k < 64; k++) {
            dom1[i][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "dom1_" + std::to_string(i) + "_" + std::to_string(k));
        }


    vector<double> t1 = { 1.0, -1.0, 1.0, -1.0 };

 
    for (int k = 0; k < 64; k++) {
        for (int l = 0; l < 3; l++) {
            if (l != 1) {
                model.addConstr(DA[0][0][0][k][l] == DA[0][1][3][(k + 36) % 64][l], "");
                model.addConstr(DA[0][1][0][(k + 44) % 64][l] == DA[0][0][2][(k + 1) % 64][l], "");
                model.addConstr(DA[0][0][4][(k + 62) % 64][l] == DA[0][1][2][(k + 6) % 64][l], "");
                model.addConstr(DA[0][0][1][(k + 28) % 64][l] == DA[0][1][4][(k + 55) % 64][l], "");
            }
        }
    }


    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 64; k++) {
                if ((i == 0 && j == 0) || (i == 1 && j == 3) || (i == 0 && j == 4) || (i == 1 && j == 2) || (i == 1 && j == 0) || (i == 0 && j == 2) || (i == 0 && j == 1) || (i == 1 && j == 4)) {
                    model.addConstr(DA[0][i][j][k][1] == 1, "");
                    GRBLinExpr known = 0;
                    known += 1 * DA[0][i][j][k][0];
                    known += 1 * DA[0][i][j][k][2];
                    //model.addConstr(known == 1, "");
                    model.addConstr(known >= 1, "");
                }
                else if (i == 4 && j == 0) { 
                    model.addConstr(factory.linExprOf(t1, { DA[0][0][0][k][0], DA[0][1][0][k][0], DA[0][1][0][k][2], DA[0][0][0][k][2] }) <= 1, "");
                    model.addConstr(factory.linExprOf(t1, { DA[0][0][0][k][0], DA[0][1][0][k][0], DA[0][1][0][k][2], DA[0][0][0][k][2] }) >= -1, "");
                    model.addConstr(DA[0][i][j][k][1] == 1, "");
                    factory.Determine_ExistOne(DA[0][i][j][k][0], { DA[0][0][0][k][0], DA[0][1][0][k][0] }, model);
                    factory.Determine_ExistOne(DA[0][i][j][k][2], { DA[0][0][0][k][2], DA[0][1][0][k][2] }, model);
                }
                else if (i == 4 && j == 2) {
                    model.addConstr(factory.linExprOf(t1, { DA[0][0][2][k][0], DA[0][1][2][k][0], DA[0][1][2][k][2], DA[0][0][2][k][2] }) <= 1, "");
                    model.addConstr(factory.linExprOf(t1, { DA[0][0][2][k][0], DA[0][1][2][k][0], DA[0][1][2][k][2], DA[0][0][2][k][2] }) >= -1, "");
                    model.addConstr(DA[0][i][j][k][1] == 1, "");
                    factory.Determine_ExistOne(DA[0][i][j][k][0], { DA[0][0][2][k][0], DA[0][1][2][k][0] }, model);
                    factory.Determine_ExistOne(DA[0][i][j][k][2], { DA[0][0][2][k][2], DA[0][1][2][k][2] }, model);
                }
                else if (i == 4 && j == 4) {
                    model.addConstr(factory.linExprOf(t1, { DA[0][0][4][k][0], DA[0][1][4][k][0], DA[0][1][4][k][2], DA[0][0][4][k][2] }) <= 1, "");
                    model.addConstr(factory.linExprOf(t1, { DA[0][0][4][k][0], DA[0][1][4][k][0], DA[0][1][4][k][2], DA[0][0][4][k][2] }) >= -1, "");
                    model.addConstr(DA[0][i][j][k][1] == 1, "");
                    factory.Determine_ExistOne(DA[0][i][j][k][0], { DA[0][0][4][k][0], DA[0][1][4][k][0] }, model);
                    factory.Determine_ExistOne(DA[0][i][j][k][2], { DA[0][0][4][k][2], DA[0][1][4][k][2] }, model);
                }
                else {
                    model.addConstr(DA[0][i][j][k][0] == 1, "");
                    model.addConstr(DA[0][i][j][k][1] == 1, "");
                    model.addConstr(DA[0][i][j][k][2] == 1, "");
                }

            }
}


void Mitm::article_attack(GRBModel& model) {
    vector<int> DA_init(5 * 5 * 64 * 3, 5);
    vector<int> DA1_init(5 * 5 * 64 * 3, 5);


    std::ifstream file(solfile);
    if (!file.is_open()) {
        std::cerr << "Unable to open the file£¡" << std::endl;
    }

    std::string line;

    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "DA_0") {

            size_t pos = line.find('_', 4);

            int dim1_index = 0, dim2_index = 0, dim3_index = 0, dim4_index = 0, counter = 0;

            while (pos != std::string::npos) {
                std::string num_str;
                size_t next_pos = line.find('_', pos + 1);
                if (next_pos != std::string::npos) {
                    num_str = line.substr(pos + 1, next_pos - pos - 1);
                }
                else {
                    num_str = line.substr(pos + 1);
                }

                int number = std::stoi(num_str);
                counter++;
                if (counter == 1)
                    dim1_index = number;
                if (counter == 2)
                    dim2_index = number;
                if (counter == 3)
                    dim3_index = number;
                if (counter == 4)
                    dim4_index = number;

                pos = next_pos;
            }
            char last_char = line.back();
            DA_init[dim4_index + 3 * (dim3_index + 64 * (dim2_index + 5 * dim1_index))] = last_char - '0';
        }

        if (line.substr(0, 4) == "DA_1") {

            size_t pos = line.find('_', 4);

            int dim1_index = 0, dim2_index = 0, dim3_index = 0, dim4_index = 0, counter = 0;

            while (pos != std::string::npos) {
                std::string num_str;
                size_t next_pos = line.find('_', pos + 1);
                if (next_pos != std::string::npos) {
                    num_str = line.substr(pos + 1, next_pos - pos - 1);
                }
                else {
                    num_str = line.substr(pos + 1);
                }

                int number = std::stoi(num_str);
                counter++;
                if (counter == 1)
                    dim1_index = number;
                if (counter == 2)
                    dim2_index = number;
                if (counter == 3)
                    dim3_index = number;
                if (counter == 4)
                    dim4_index = number;

                pos = next_pos;
            }
            char last_char = line.back();
            DA1_init[dim4_index + 3 * (dim3_index + 64 * (dim2_index + 5 * dim1_index))] = last_char - '0';
        }
    }

    file.close();

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 64; k++)
                for (int l = 0; l < 3; l++) {
                    model.addConstr(DA[0][i][j][k][l] == DA_init[l + 3 * (k + 64 * (j + 5 * i))], "");
                    model.addConstr(DA[1][i][j][k][l] == DA1_init[l + 3 * (k + 64 * (j + 5 * i))], "");
                }
}

void Mitm::solve(GRBModel& model)
{
    vector<MitmSolution> solutions;
    string tikz;
    model.write("model.lp");
    model.optimize();
    model.write("model.sol");
    getAllFoundSolutions(model, solutions);

    tikz = generate_1(solutions[0]);


    ofstream fw;
    fw.open("pic\\tikzpic1.tex");

    if (!fw.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    fw << tikz;
    fw.close();

}


string Mitm::generate_1(MitmSolution& MitmSolution) {
    string output = "";
    // Header
    output += "\\documentclass{standalone}\n";
    output += "\\usepackage{tikz}\n";
    output += "\\usepackage{calc}\n";
    output += "\\usepackage{pgffor}\n";
    output += "\\usetikzlibrary{patterns}\n";
    output += "\\tikzset{base/.style = {draw=black, minimum width=0.02cm, minimum height=0.02cm, align=center, on chain},}\n";
    output += "\\begin{document}\n";
    output += "\\begin{tikzpicture}[scale = 0.45,every node/.style={scale=0.5}]\n";
    output += "\\makeatletter\n";

    const int dim1 = 5;
    const int dim2 = 5;
    const int dim3 = 64;
    const int dim4 = 3;
    double result = 0.0;
    std::ostringstream oss;
    std::string a;
    int orangenum = 0;

    vector<vector<vector<vector<double>>>> array(dim1, vector<vector<vector<double>>>(dim2, vector<vector<double>>(dim3, vector<double>(dim4, 0.0))));

    vector<vector<vector<double>>> array3D(dim1, vector<vector<double>>(dim2, vector<double>(dim3, 0.0)));

    const int arraySize = 576;
    string stringArray[arraySize];
    int consumered = 0, consumeblue = 0;


    for (int round = 0; round < MitmSolution.Rounds + 1; round++) {

        output += " \\node[align=center] at (" + std::to_string(-2) + "," + std::to_string(25 * (MitmSolution.Rounds - round) + 20.5) + ") [scale=1.5]{\\textbf{\\huge $A^{(" + std::to_string(round + 1) + ")}$}};\n";

        for (int k = 0; k < 64; k++)
        {
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                {
                    if (MitmSolution.DA[round][i][j][k][0] == 0 && MitmSolution.DA[round][i][j][k][1] == 0 && MitmSolution.DA[round][i][j][k][2] == 0)
                        output += "\\fill[color=" + WhiteColor + "] (" + std::to_string(6 * k + i) + "," + std::to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 18) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DA[round][i][j][k][0] == 0 && MitmSolution.DA[round][i][j][k][1] == 1 && MitmSolution.DA[round][i][j][k][2] == 1)
                        output += "\\fill[color=" + RedColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 18) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DA[round][i][j][k][0] == 1 && MitmSolution.DA[round][i][j][k][1] == 1 && MitmSolution.DA[round][i][j][k][2] == 0)
                        output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 18) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DA[round][i][j][k][0] == 0 && MitmSolution.DA[round][i][j][k][1] == 1 && MitmSolution.DA[round][i][j][k][2] == 0)
                        output += "\\fill[color=" + PurpleColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 18) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DA[round][i][j][k][0] == 1 && MitmSolution.DA[round][i][j][k][1] == 1 && MitmSolution.DA[round][i][j][k][2] == 1)
                        output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 18) + ") rectangle ++(1,1);\n";
                }
            output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 18) + ") grid ++(5,5);\n";
            output += " \\node[align=center] at (" + to_string(6 * k + 2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 17) + ")[scale=2] {{\\Large z=" + to_string(k) + "}};\n";
        }
        //DP


        if (round != MitmSolution.Rounds) {
            output += " \\node[align=center] at (" + to_string(-2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15.5) + ") [scale=1.5]{\\textbf{\\huge $C^{(" + to_string(round + 1) + ")}$}};\n";
            consumered = 0, consumeblue = 0;
            for (int k = 0; k < 64; k++)
            {
                for (int i = 0; i < 5; i++)
                {
                    if (MitmSolution.DP[round][i][k][0] == 0 && MitmSolution.DP[round][i][k][1] == 0 && MitmSolution.DP[round][i][k][2] == 0)
                        output += "\\fill[color=" + WhiteColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DP[round][i][k][0] == 0 && MitmSolution.DP[round][i][k][1] == 1 && MitmSolution.DP[round][i][k][2] == 1)
                        output += "\\fill[color=" + RedColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DP[round][i][k][0] == 1 && MitmSolution.DP[round][i][k][1] == 1 && MitmSolution.DP[round][i][k][2] == 0)
                        output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DP[round][i][k][0] == 0 && MitmSolution.DP[round][i][k][1] == 1 && MitmSolution.DP[round][i][k][2] == 0)
                        output += "\\fill[color=" + PurpleColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DP[round][i][k][0] == 1 && MitmSolution.DP[round][i][k][1] == 1 && MitmSolution.DP[round][i][k][2] == 1)
                        output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") rectangle ++(1,1);\n";
                }
                output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(5,1);\n";
                //output += " \\node[align=center] at ("+(6*k+2)+","+(21*(MitmSolution.Rounds-round)+13)+") {\\textbf{\\Large z="+k+"}};\n";
                for (int i = 0; i < 5; i++)
                {
                    //DC1
                    if (MitmSolution.DC1_red[round][i][k] == 2 && MitmSolution.beta1[round][i][k] == 0) {
                        consumered += 2;
                        output += "\\draw[line width=2pt, color=yellow](" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
                    }
                    if (MitmSolution.DC1_blue[round][i][k] == 2 && MitmSolution.beta1[round][i][k] == 0) {
                        consumeblue += 2;
                        output += "\\draw[line width=2pt, color=black](" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
                    }
                    if (MitmSolution.DC1_blue[round][i][k] == 2 && MitmSolution.beta1[round][i][k] == 1) {
                        consumeblue++;
                        consumered++;
                        orangenum++;
                        output += "\\draw[line width=2pt, color=orange](" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
                    }
                }
            }
            oss.str("");
            oss.clear();
            result = consumered / 2.0;
            oss << std::fixed << std::setprecision(1) << result;
            a = oss.str();
            output += "\\fill[color=" + RedColor + "] (" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ")  rectangle ++(1,1);\n";
            output += "\\draw(" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
            output += " \\node[align=center] at (" + to_string(6 * 63 + 4 + 2.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15.4) + ") {\\textbf{\\huge - " + a + "}};\n";
            oss.str("");
            oss.clear();
            result = consumeblue / 2.0;
            oss << std::fixed << std::setprecision(1) << result;
            a = oss.str();
            output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 16.5) + ")  rectangle ++(1,1);\n";
            output += "\\draw(" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 16.5) + ") grid ++(1,1);\n";
            output += " \\node[align=center] at (" + to_string(6 * 63 + 4 + 2.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + 16.9) + ") {\\textbf{\\huge - " + a + "}};\n";



            output += " \\node[align=center] at (" + to_string(-2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13.5) + ")[scale=1.5] {\\textbf{\\huge $D^{(" + to_string(round + 1) + ")}$}};\n";
            for (int k = 0; k < 64; k++)
            {
                for (int i = 0; i < 5; i++)
                {
                    if (MitmSolution.DP2[round][i][k][0] == 0 && MitmSolution.DP2[round][i][k][1] == 0 && MitmSolution.DP2[round][i][k][2] == 0)
                        output += "\\fill[color=" + WhiteColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DP2[round][i][k][0] == 0 && MitmSolution.DP2[round][i][k][1] == 1 && MitmSolution.DP2[round][i][k][2] == 1)
                        output += "\\fill[color=" + RedColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DP2[round][i][k][0] == 1 && MitmSolution.DP2[round][i][k][1] == 1 && MitmSolution.DP2[round][i][k][2] == 0)
                        output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DP2[round][i][k][0] == 0 && MitmSolution.DP2[round][i][k][1] == 1 && MitmSolution.DP2[round][i][k][2] == 0)
                        output += "\\fill[color=" + PurpleColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") rectangle ++(1,1);\n";
                    if (MitmSolution.DP2[round][i][k][0] == 1 && MitmSolution.DP2[round][i][k][1] == 1 && MitmSolution.DP2[round][i][k][2] == 1)
                        output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") rectangle ++(1,1);\n";
                }
                output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") grid ++(5,1);\n";
                //output += " \\node[align=center] at ("+(6*k+2)+","+(21*(MitmSolution.Rounds-round)+13)+") {\\textbf{\\Large z="+k+"}};\n";
                for (int i = 0; i < 5; i++)
                {
                    //DC12
                    if (MitmSolution.DC12_red[round][i][k] == 2 && MitmSolution.beta3[round][i][k] == 0) {
                        consumered += 2;
                        output += "\\draw[line width=2pt, color=yellow](" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") grid ++(1,1);\n";
                    }
                    if (MitmSolution.DC12_blue[round][i][k] == 2 && MitmSolution.beta3[round][i][k] == 0) {
                        consumeblue += 2;
                        output += "\\draw[line width=2pt, color=black](" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") grid ++(1,1);\n";
                    }
                    if (MitmSolution.DC12_blue[round][i][k] == 2 && MitmSolution.beta1[round][i][k] == 1) {
                        consumeblue++;
                        consumered++;
                        orangenum++;
                        output += "\\draw[line width=2pt, color=orange](" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") grid ++(1,1);\n";
                    }
                }
               
            }
        }
        // DB after theta
        if (round != MitmSolution.Rounds) {
            output += " \\node[align=center] at (" + to_string(-2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 9.5) + ")[scale=1.5] {\\textbf{\\huge ${\\Huge \\theta}^{(" + to_string(round + 1) + ")}$}};\n";
            for (int k = 0; k < 64; k++)
            {
                for (int i = 0; i < 5; i++)
                    for (int j = 0; j < 5; j++)
                    {
                        if (MitmSolution.DB[round][i][j][k][0] == 0 && MitmSolution.DB[round][i][j][k][1] == 0 && MitmSolution.DB[round][i][j][k][2] == 0)
                            output += "\\fill[color=" + WhiteColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DB[round][i][j][k][0] == 0 && MitmSolution.DB[round][i][j][k][1] == 1 && MitmSolution.DB[round][i][j][k][2] == 1)
                            output += "\\fill[color=" + RedColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DB[round][i][j][k][0] == 1 && MitmSolution.DB[round][i][j][k][1] == 1 && MitmSolution.DB[round][i][j][k][2] == 0)
                            output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DB[round][i][j][k][0] == 0 && MitmSolution.DB[round][i][j][k][1] == 1 && MitmSolution.DB[round][i][j][k][2] == 0)
                            output += "\\fill[color=" + PurpleColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DB[round][i][j][k][0] == 1 && MitmSolution.DB[round][i][j][k][1] == 1 && MitmSolution.DB[round][i][j][k][2] == 1)
                            output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                    }
            }
            for (int k = 0; k < 64; k++) {
                output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 7) + ") grid ++(5,5);\n";
                //output += " \\node[align=center] at ("+(6*k+2)+","+(25*(MitmSolution.Rounds-round)+6)+") {\\textbf{\\Large z="+k+"}};\n";
            }
            consumered = 0, consumeblue = 0;
            for (int k = 0; k < 64; k++)
            {
                for (int i = 0; i < 5; i++)
                    for (int j = 0; j < 5; j++)
                    {
                        //DC2
                        if (MitmSolution.DC2_red[round][i][j][k] == 2 && MitmSolution.beta2[round][i][j][k] == 0) {
                            consumered += 2;
                            output += "\\draw[line width=2pt, color=yellow](" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 7) + ") grid ++(1,1);\n";
                        }

                        if (MitmSolution.DC2_blue[round][i][j][k] == 2 && MitmSolution.beta2[round][i][j][k] == 0) {
                            consumeblue += 2;
                            output += "\\draw[line width=2pt, color=black](" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 7) + ") grid ++(1,1);\n";
                        }
                        if (MitmSolution.DC2_blue[round][i][j][k] == 2 && MitmSolution.beta2[round][i][j][k] == 1) {
                            consumeblue++;
                            consumered++;
                            orangenum++;
                            output += "\\draw[line width=2pt, color=orange](" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 7) + ") grid ++(1,1);\n";
                        }
                    }
            }
            oss.str("");
            oss.clear();
            result = consumered / 2.0;
            oss << std::fixed << std::setprecision(1) << result;
            a = oss.str();
            output += "\\fill[color=" + RedColor + "] (" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 9) + ")  rectangle ++(1,1);\n";
            output += "\\draw(" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 9) + ")  grid ++(1,1);\n";
            output += " \\node[align=center] at (" + to_string(6 * 63 + 4 + 2.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + 9.4) + ") {\\textbf{\\huge - " + a + "}};\n";
            oss.str("");
            oss.clear();
            result = consumeblue / 2.0;
            oss << std::fixed << std::setprecision(1) << result;
            a = oss.str();
            output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 10.5) + ")  rectangle ++(1,1);\n";
            output += "\\draw(" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 10.5) + ")  grid ++(1,1);\n";
            output += " \\node[align=center] at (" + to_string(6 * 63 + 4 + 2.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + 10.9) + ") {\\textbf{\\huge - " + a + "}};\n";

        }

        // DB
        if (round != MitmSolution.Rounds) {
            output += " \\node[align=center] at (" + to_string(-2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 2.5) + ")[scale=1.5] {\\textbf{\\huge $\\pi^{(" + to_string(round + 1) + ")}$}};\n";
            for (int k = 0; k < 64; k++)
            {
                for (int i = 0; i < 5; i++)
                    for (int j = 0; j < 5; j++)
                    {
                        if (MitmSolution.DB[round][i][j][k][0] == 0 && MitmSolution.DB[round][i][j][k][1] == 0 && MitmSolution.DB[round][i][j][k][2] == 0)
                            output += "\\fill[color=" + WhiteColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DB[round][i][j][k][0] == 0 && MitmSolution.DB[round][i][j][k][1] == 1 && MitmSolution.DB[round][i][j][k][2] == 1)
                            output += "\\fill[color=" + RedColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DB[round][i][j][k][0] == 1 && MitmSolution.DB[round][i][j][k][1] == 1 && MitmSolution.DB[round][i][j][k][2] == 0)
                            output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DB[round][i][j][k][0] == 0 && MitmSolution.DB[round][i][j][k][1] == 1 && MitmSolution.DB[round][i][j][k][2] == 0)
                            output += "\\fill[color=" + PurpleColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DB[round][i][j][k][0] == 1 && MitmSolution.DB[round][i][j][k][1] == 1 && MitmSolution.DB[round][i][j][k][2] == 1)
                            output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                    }
                output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round)) + ") grid ++(5,5);\n";
                //output += " \\node[align=center] at ("+(6*k+2)+","+(25*(MitmSolution.Rounds-round)-1)+") {\\textbf{\\Large z="+k+"}};\n";
            }
        }

        for (int k = 0; k < 64; k++)
        {
            if (MitmSolution.dom1[0][k] == 1)
            {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[0][2] + 64) % 64) + 0 + 0.5) + "," + to_string(18 + 2 + 0.5) + ") {\\textbf{\\Large $m$}};\n";
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[0][2] + 64) % 64) + 0 + 0.5) + "," + to_string(18 + 4 + 0.5) + ") {\\textbf{\\Large $m$}};\n";

                //output += "\\draw[line width=2pt, color=green]("+(6*((k-rho[3][0]+64)%64)+3)+","+(18+4)+") grid ++(1,1);\n";
                    //output += "\\draw[line width=2pt, color=green]("+(6*((k-rho[3][0]+64)%64)+3)+","+(18+1)+") grid ++(1,1);\n";
                    //output += "\\draw[line width=2pt, color=green]("+(6*((k-rho[0][2]+64)%64)+0)+","+(18+2)+") grid ++(1,1);\n";
                    //output += "\\draw[line width=2pt, color=green]("+(6*((k-rho[0][2]+64)%64)+0)+","+(18+4)+") grid ++(1,1);\n";
            }
            if (MitmSolution.dom1[1][k] == 1)
            {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[1][3] + 64) % 64) + 1 + 0.5) + "," + to_string(18 + 1 + 0.5) + ") {\\textbf{\\Large $m$}};\n";
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[1][3] + 64) % 64) + 1 + 0.5) + "," + to_string(18 + 3 + 0.5) + ") {\\textbf{\\Large $m$}};\n";

                //output += "\\draw[line width=2pt, color=green]("+(6*((k-rho[3][0]+64)%64)+3)+","+(18+4)+") grid ++(1,1);\n";
                    //output += "\\draw[line width=2pt, color=green]("+(6*((k-rho[3][0]+64)%64)+3)+","+(18+1)+") grid ++(1,1);\n";
                    //output += "\\draw[line width=2pt, color=green]("+(6*((k-rho[0][2]+64)%64)+0)+","+(18+2)+") grid ++(1,1);\n";
                    //output += "\\draw[line width=2pt, color=green]("+(6*((k-rho[0][2]+64)%64)+0)+","+(18+4)+") grid ++(1,1);\n";
            }
            if (MitmSolution.dom1[2][k] == 1)
            {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[3][0] + 64) % 64) + 3 + 0.5) + "," + to_string(18 + 4 + 0.5) + ") {\\textbf{\\Large $m$}};\n";
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[3][0] + 64) % 64) + 3 + 0.5) + "," + to_string(18 + 1 + 0.5) + ") {\\textbf{\\Large $m$}};\n";
            }
            if (MitmSolution.dom1[3][k] == 1)
            {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[4][1] + 64) % 64) + 4 + 0.5) + "," + to_string(18 + 3 + 0.5) + ") {\\textbf{\\Large $m$}};\n";
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[4][1] + 64) % 64) + 4 + 0.5) + "," + to_string(18 + 0 + 0.5) + ") {\\textbf{\\Large $m$}};\n";
            }
        }
        output += "\n";

    }

    int count = 0, count1 = 0;
    int round = -1;
    int allred = 0;
    int allblue = 0;
    if (round == -1) {

        output += " \\node[align=center] at (" + to_string(-2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 20.5) + ") [scale=1.5]{\\textbf{\\huge $A^{(" + to_string(round + 1) + ")}$}};\n";

        for (int k = 0; k < 64; k++)
        {
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                {
                    if (i == 0 || i == 1) {
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 0 && MitmSolution.DA[round + 1][i][j][k][1] == 0 && MitmSolution.DA[round + 1][i][j][k][2] == 0)
                            output += "\\fill[color=" + WhiteColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 18) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 0 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 1)
                            output += "\\fill[color=" + RedColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 18) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 1 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 0)
                            output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 18) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 0 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 0)
                            output += "\\fill[color=" + PurpleColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 18) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 1 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 1)
                            output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 18) + ") rectangle ++(1,1);\n";
                    }
                    else {
                        output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 18) + ") rectangle ++(1,1);\n";
                    }
                }
        }
        for (int k = 0; k < 64; k++)
        {
            output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 18) + ") grid ++(5,5);\n";
            output += " \\node[align=center] at (" + to_string(6 * k + 2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 17) + ")[scale=2] {{\\Large z=" + to_string(k) + "}};\n";
        }

        output += " \\node[align=center] at (" + to_string(-2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15.5) + ")[scale=1.5] {\\textbf{\\huge $C^{(" + to_string(round + 1) + ")}$}};\n";
        for (int k = 0; k < 64; k++)
        {
            for (int i = 0; i < 5; i++)
            {
                output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") rectangle ++(1,1);\n";
            }
            output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(5,1);\n";
        }
        for (int k = 0; k < 64; k++)
        {

            //DC1
            if (MitmSolution.DA[round + 1][0][0][k][0] == 0 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 1) {
                allred++;
                output += "\\draw[line width=2pt, color=yellow](" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 1 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 0) {
                allblue++;
                output += "\\draw[line width=2pt, color=black](" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
            }
            if (MitmSolution.DA[round + 1][1][0][k][0] == 0 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 1) {
                allred++;
                output += "\\draw[line width=2pt, color=yellow](" + to_string(6 * ((k - rho[1][1] + 64) % 64) + 1) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
            }
            if (MitmSolution.DA[round + 1][1][0][k][0] == 1 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 0) {
                allblue++;
                output += "\\draw[line width=2pt, color=black](" + to_string(6 * ((k - rho[1][1] + 64) % 64) + 1) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
            }
            if (MitmSolution.DA[round + 1][1][2][k][0] == 0 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 1) {
                allred++;
                output += "\\draw[line width=2pt, color=yellow](" + to_string(6 * ((k - rho[2][1] + 64) % 64) + 2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
            }
            if (MitmSolution.DA[round + 1][1][2][k][0] == 1 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 0) {
                allblue++;
                output += "\\draw[line width=2pt, color=black](" + to_string(6 * ((k - rho[2][1] + 64) % 64) + 2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
            }
            if (MitmSolution.DA[round + 1][0][1][k][0] == 0 && MitmSolution.DA[round + 1][0][1][k][1] == 1 && MitmSolution.DA[round + 1][0][1][k][2] == 1) {
                allred++;
                output += "\\draw[line width=2pt, color=yellow](" + to_string(6 * ((k - rho[3][0] + 64) % 64) + 3) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
            }
            if (MitmSolution.DA[round + 1][0][1][k][0] == 1 && MitmSolution.DA[round + 1][0][1][k][1] == 1 && MitmSolution.DA[round + 1][0][1][k][2] == 0) {
                allblue++;
                output += "\\draw[line width=2pt, color=black](" + to_string(6 * ((k - rho[3][0] + 64) % 64) + 3) + "," + to_string(25 * (MitmSolution.Rounds - round) + 15) + ") grid ++(1,1);\n";
            }

        }

        output += "\\fill[color=" + RedColor + "] (" + to_string(6 * 63 + 4 + 2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 21) + ")  rectangle ++(1,1);\n";
        output += "\\draw(" + to_string(6 * 63 + 4 + 2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 21) + ")  grid ++(1,1);\n";
        output += " \\node[align=center] at (" + to_string(6 * 63 + 4 + 4.2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 21.4) + ") {\\textbf{\\huge = " + to_string(allred * 2) + "}};\n";
        output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * 63 + 4 + 2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 19) + ")  rectangle ++(1,1);\n";
        output += "\\draw(" + to_string(6 * 63 + 4 + 2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 19) + ") grid ++(1,1);\n";
        output += " \\node[align=center] at (" + to_string(6 * 63 + 4 + 4.1) + "," + to_string(25 * (MitmSolution.Rounds - round) + 19.4) + ") {\\textbf{\\huge = " + to_string(allblue * 2) + "}};\n";
        output += " \\node[align=center] at (" + to_string(-2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13.5) + ") [scale=1.5] {\\textbf{\\huge $D^{(" + to_string(round + 1) + ")}$}};\n";
        for (int k = 0; k < 64; k++)
        {
            for (int i = 0; i < 5; i++)
            {
                output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") rectangle ++(1,1);\n";
            }
            output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 13) + ") grid ++(5,1);\n";
        }

        output += "\\fill[color=" + RedColor + "] (" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 16) + ")  rectangle ++(1,1);\n";
        output += "\\draw(" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 16) + ")  grid ++(1,1);\n";
        output += " \\node[align=center] at (" + to_string(6 * 63 + 4 + 2.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + 16.4) + ") {\\textbf{\\huge - " + to_string(allred) + "}};\n";
        output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 14.5) + ")  rectangle ++(1,1);\n";
        output += "\\draw(" + to_string(6 * 63 + 4 + 4) + "," + to_string(25 * (MitmSolution.Rounds - round) + 14.5) + ") grid ++(1,1);\n";
        output += " \\node[align=center] at (" + to_string(6 * 63 + 4 + 2.7) + "," + to_string(25 * (MitmSolution.Rounds - round) + 14.9) + ") {\\textbf{\\huge - " + to_string(allblue) + "}};\n";


        // DB after theta

        output += " \\node[align=center] at (" + to_string(-2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 9.5) + ")[scale=1.5] {{\\huge $\\theta^{(" + to_string(round + 1) + ")}$}};\n";
        for (int k = 0; k < 64; k++)
        {
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                {
                    if (i == 0 || i == 1) {
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 0 && MitmSolution.DA[round + 1][i][j][k][1] == 0 && MitmSolution.DA[round + 1][i][j][k][2] == 0)
                            output += "\\fill[color=" + WhiteColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 0 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 1) {
                            output += "\\fill[color=" + RedColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                            array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][0] = 0;
                            array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][1] = 1;
                            array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][2] = 1;
                        }
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 1 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 0) {
                            output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                            array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][0] = 1;
                            array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][1] = 1;
                            array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][2] = 0;
                        }
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 0 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 0)
                            output += "\\fill[color=" + PurpleColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 1 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 1) {
                            output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                            array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][0] = 1;
                            array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][1] = 1;
                            array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][2] = 1;
                        }
                    }
                    else {
                        output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * ((k - rho[(i + 3 * j) % 5][i] + 64) % 64) + ((i + 3 * j) % 5)) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - i) + 7) + ") rectangle ++(1,1);\n";
                        array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][0] = 1;
                        array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][1] = 1;
                        array[(i + 3 * j) % 5][i][(k - rho[(i + 3 * j) % 5][i] + 64) % 64][2] = 1;
                    }
                }


            int i = 0; int j = 0;
            if (MitmSolution.DA[round + 1][0][0][k][0] == 0 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 1 && MitmSolution.DA[round + 1][1][0][k][0] == 0 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 1 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 0 && MitmSolution.DA[round + 1][1][0][k][0] == 1 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 0 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 1 && MitmSolution.DA[round + 1][1][0][k][0] == 1 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64) + (((i + 1) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 1) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5][(k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 1 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 1 && MitmSolution.DA[round + 1][1][0][k][0] == 0 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64) + (((i + 2) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 2) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5][(k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i) % 5 + 3 * j) % 5][(i) % 5] + 64) % 64) + (((i) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5][(k - rho[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 1 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 0 && MitmSolution.DA[round + 1][1][0][k][0] == 1 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64) + (((i + 1) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 1) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5][(k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 1 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 1 && MitmSolution.DA[round + 1][1][0][k][0] == 1 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64) + (((i + 2) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 2) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5][(k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i) % 5 + 3 * j) % 5][(i) % 5] + 64) % 64) + (((i) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5][(k - rho[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5] + 64) % 64] = 1;
            }

            i = 0; j = 1;
            if (MitmSolution.DA[round + 1][0][1][k][0] == 0 && MitmSolution.DA[round + 1][0][1][k][1] == 1 && MitmSolution.DA[round + 1][0][1][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64) + (((i + 1) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 1) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5][(k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][1][k][0] == 1 && MitmSolution.DA[round + 1][0][1][k][1] == 1 && MitmSolution.DA[round + 1][0][1][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64) + (((i + 1) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 1) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5][(k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }

            i = 0; j = 2;
            if (MitmSolution.DA[round + 1][0][2][k][0] == 0 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 1 && MitmSolution.DA[round + 1][1][2][k][0] == 0 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 1 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 0 && MitmSolution.DA[round + 1][1][2][k][0] == 1 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 0 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 1 && MitmSolution.DA[round + 1][1][2][k][0] == 1 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64) + (((i + 1) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 1) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5][(k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 1 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 1 && MitmSolution.DA[round + 1][1][2][k][0] == 0 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64) + (((i + 2) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 2) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5][(k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i) % 5 + 3 * j) % 5][(i) % 5] + 64) % 64) + (((i) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5][(k - rho[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 1 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 0 && MitmSolution.DA[round + 1][1][2][k][0] == 1 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64) + (((i + 1) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 1) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5][(k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 1 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 1 && MitmSolution.DA[round + 1][1][2][k][0] == 1 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64) + (((i + 2) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 2) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5][(k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i) % 5 + 3 * j) % 5][(i) % 5] + 64) % 64) + (((i) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5][(k - rho[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5] + 64) % 64] = 1;
            }

            i = 0; j = 3;
            if (MitmSolution.DA[round + 1][1][3][k][0] == 0 && MitmSolution.DA[round + 1][1][3][k][1] == 1 && MitmSolution.DA[round + 1][1][3][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64) + (((i + 2) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 2) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5][(k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i) % 5 + 3 * j) % 5][(i) % 5] + 64) % 64) + (((i) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5][(k - rho[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][1][3][k][0] == 1 && MitmSolution.DA[round + 1][1][3][k][1] == 1 && MitmSolution.DA[round + 1][1][3][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64) + (((i + 2) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 2) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5][(k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i) % 5 + 3 * j) % 5][(i) % 5] + 64) % 64) + (((i) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5][(k - rho[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5] + 64) % 64] = 1;
            }

            i = 0; j = 4;
            if (MitmSolution.DA[round + 1][0][4][k][0] == 0 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 1 && MitmSolution.DA[round + 1][1][4][k][0] == 0 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 1 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 0 && MitmSolution.DA[round + 1][1][4][k][0] == 1 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 0 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 1 && MitmSolution.DA[round + 1][1][4][k][0] == 1 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64) + (((i + 1) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 1) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5][(k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 1 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 1 && MitmSolution.DA[round + 1][1][4][k][0] == 0 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64) + (((i + 2) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 2) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5][(k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i) % 5 + 3 * j) % 5][(i) % 5] + 64) % 64) + (((i) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5][(k - rho[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 1 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 0 && MitmSolution.DA[round + 1][1][4][k][0] == 1 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64) + (((i + 1) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 1) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5][(k - rho[((i + 1) % 5 + 3 * j) % 5][(i + 1) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64) + (((i + 4) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 4) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5][(k - rho[((i + 4) % 5 + 3 * j) % 5][(i + 4) % 5] + 64) % 64] = 1;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 1 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 1 && MitmSolution.DA[round + 1][1][4][k][0] == 1 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64) + (((i + 2) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i + 2) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count1++;
                array3D[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5][(k - rho[((i + 2) % 5 + 3 * j) % 5][(i + 2) % 5] + 64) % 64] = 1;
                output += "\\node[align=center] at (" + to_string(6 * ((k - rho[((i) % 5 + 3 * j) % 5][(i) % 5] + 64) % 64) + (((i) % 5 + 3 * j) % 5) + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - (i) % 5) + 7 + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count1++;
                array3D[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5][(k - rho[((i + 0) % 5 + 3 * j) % 5][(i + 0) % 5] + 64) % 64] = 1;
            }

        }
        for (int k = 0; k < 64; k++) {
            output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round) + 7) + ") grid ++(5,5);\n";
            //output += " \\node[align=center] at ("+(6*k+2)+","+(25*(MitmSolution.Rounds-round)+6)+") {\\textbf{\\Large z="+k+"}};\n";
        }


        // DB
        output += " \\node[align=center] at (" + to_string(-2) + "," + to_string(25 * (MitmSolution.Rounds - round) + 2.5) + ") [scale=1.5]{\\textbf{\\huge $\\pi^{(" + to_string(round + 1) + ")}$}};\n";
        for (int k = 0; k < 64; k++)
        {
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                {
                    if (i == 0 || i == 1) {
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 0 && MitmSolution.DA[round + 1][i][j][k][1] == 0 && MitmSolution.DA[round + 1][i][j][k][2] == 0)
                            output += "\\fill[color=" + WhiteColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 0 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 1)
                            output += "\\fill[color=" + RedColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 1 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 0)
                            output += "\\fill[color=" + BlueColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 0 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 0)
                            output += "\\fill[color=" + PurpleColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                        if (MitmSolution.DA[round + 1][i][j][k][0] == 1 && MitmSolution.DA[round + 1][i][j][k][1] == 1 && MitmSolution.DA[round + 1][i][j][k][2] == 1)
                            output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                    }
                    else {
                        output += "\\fill[color=" + GrayColor + "] (" + to_string(6 * k + i) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j)) + ") rectangle ++(1,1);\n";
                    }
                }
            output += "\\draw(" + to_string(6 * k) + "," + to_string(25 * (MitmSolution.Rounds - round)) + ") grid ++(5,5);\n";
            //output += " \\node[align=center] at ("+(6*k+2)+","+(25*(MitmSolution.Rounds-round)-1)+") {\\textbf{\\Large z="+k+"}};\n";

            int i = 0; int j = 0;
            if (MitmSolution.DA[round + 1][0][0][k][0] == 0 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 1 && MitmSolution.DA[round + 1][1][0][k][0] == 0 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 1 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 0 && MitmSolution.DA[round + 1][1][0][k][0] == 1 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 0 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 1 && MitmSolution.DA[round + 1][1][0][k][0] == 1 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 1) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 1 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 1 && MitmSolution.DA[round + 1][1][0][k][0] == 0 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + i % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 2) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 1 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 0 && MitmSolution.DA[round + 1][1][0][k][0] == 1 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 1) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][0][k][0] == 1 && MitmSolution.DA[round + 1][0][0][k][1] == 1 && MitmSolution.DA[round + 1][0][0][k][2] == 1 && MitmSolution.DA[round + 1][1][0][k][0] == 1 && MitmSolution.DA[round + 1][1][0][k][1] == 1 && MitmSolution.DA[round + 1][1][0][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * k + i % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 2) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
            }

            i = 0; j = 1;
            if (MitmSolution.DA[round + 1][0][1][k][0] == 0 && MitmSolution.DA[round + 1][0][1][k][1] == 1 && MitmSolution.DA[round + 1][0][1][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 1) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][1][k][0] == 1 && MitmSolution.DA[round + 1][0][1][k][1] == 1 && MitmSolution.DA[round + 1][0][1][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 1) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }

            i = 0; j = 2;
            if (MitmSolution.DA[round + 1][0][2][k][0] == 0 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 1 && MitmSolution.DA[round + 1][1][2][k][0] == 0 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 1 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 0 && MitmSolution.DA[round + 1][1][2][k][0] == 1 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 0 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 1 && MitmSolution.DA[round + 1][1][2][k][0] == 1 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 1) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 1 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 1 && MitmSolution.DA[round + 1][1][2][k][0] == 0 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + i % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 2) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 1 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 0 && MitmSolution.DA[round + 1][1][2][k][0] == 1 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 1) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][2][k][0] == 1 && MitmSolution.DA[round + 1][0][2][k][1] == 1 && MitmSolution.DA[round + 1][0][2][k][2] == 1 && MitmSolution.DA[round + 1][1][2][k][0] == 1 && MitmSolution.DA[round + 1][1][2][k][1] == 1 && MitmSolution.DA[round + 1][1][2][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * k + i % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 2) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
            }

            i = 0; j = 3;
            if (MitmSolution.DA[round + 1][1][3][k][0] == 0 && MitmSolution.DA[round + 1][1][3][k][1] == 1 && MitmSolution.DA[round + 1][1][3][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + i % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 2) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][1][3][k][0] == 1 && MitmSolution.DA[round + 1][1][3][k][1] == 1 && MitmSolution.DA[round + 1][1][3][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * k + i % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 2) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
            }

            i = 0; j = 4;
            if (MitmSolution.DA[round + 1][0][4][k][0] == 0 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 1 && MitmSolution.DA[round + 1][1][4][k][0] == 0 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 1 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 0 && MitmSolution.DA[round + 1][1][4][k][0] == 1 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 0 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 1 && MitmSolution.DA[round + 1][1][4][k][0] == 1 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 1) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 1 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 1 && MitmSolution.DA[round + 1][1][4][k][0] == 0 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + i % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 2) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 1 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 0 && MitmSolution.DA[round + 1][1][4][k][0] == 1 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 1) {
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 1) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 4) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
            }
            if (MitmSolution.DA[round + 1][0][4][k][0] == 1 && MitmSolution.DA[round + 1][0][4][k][1] == 1 && MitmSolution.DA[round + 1][0][4][k][2] == 1 && MitmSolution.DA[round + 1][1][4][k][0] == 1 && MitmSolution.DA[round + 1][1][4][k][1] == 1 && MitmSolution.DA[round + 1][1][4][k][2] == 0) {
                output += "\\node[align=center] at (" + to_string(6 * k + i % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 1}};\n";
                count++;
                output += "\\node[align=center] at (" + to_string(6 * k + (i + 2) % 5 + 0.5) + "," + to_string(25 * (MitmSolution.Rounds - round) + (4 - j) + 0.5) + ") {\\textbf{\\Large 0}};\n";
                count++;
            }
        }

    }

    int num = 0;

    for (int k = 0; k < dim3; ++k) {
        for (int i = 0; i < dim1; ++i) {
            for (int j = 0; j < dim2; ++j) {
                if (j == 0) {
                    if (i < 4) {
                        if ((array[i][j][k][0] == 1 && array[i][j][k][1] == 1 && array[i][j][k][2] == 0) || (array[i][j][k][0] == 0 && array[i][j][k][1] == 1 && array[i][j][k][2] == 1)) {
                            stringArray[2 * i + 9 * k] = "a" + to_string(num);
                            stringArray[2 * i + 9 * k + 1] = "a" + to_string(num);
                            num++;
                        }
                        else {
                            stringArray[2 * i + 9 * k] = "a" + to_string(num);
                            num++;
                            stringArray[2 * i + 9 * k + 1] = "a" + to_string(num);
                            num++;
                        }
                    }
                    else {
                        stringArray[2 * i + 9 * k] = "a" + to_string(num);
                        num++;

                    }

                }

            }
        }
    }

    num = 0;
    size_t pos1 = 0;

    const int rows = count;
    const int cols = 576 - allred - allblue;

    double** matrix = new double* [rows];
    for (int i = 0; i < rows; ++i) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = 0;
        }
    }

    for (int k = 0; k < dim3; ++k) {
        for (int i = 0; i < dim1; ++i) {
            for (int j = 0; j < dim2; ++j) {
                if (array3D[i][j][k] == 1) {
                    if ((i + 4) % 5 < 4) {
                        pos1 = stringArray[2 * ((i + 4) % 5) + 9 * k].find_first_not_of("0123456789", 1);
                        int postion1 = stoi(stringArray[2 * ((i + 4) % 5) + 9 * k].substr(1, pos1 - 1));
                        matrix[num][postion1] = 1;
                        pos1 = stringArray[2 * ((i + 4) % 5) + 9 * k + 1].find_first_not_of("0123456789", 1);
                        int postion2 = stoi(stringArray[2 * ((i + 4) % 5) + 9 * k + 1].substr(1, pos1 - 1));
                        matrix[num][postion2] = 1;
                    }
                    else {
                        pos1 = stringArray[2 * ((i + 4) % 5) + 9 * k].find_first_not_of("0123456789", 1);
                        int postion1 = stoi(stringArray[2 * ((i + 4) % 5) + 9 * k].substr(1, pos1 - 1));
                        matrix[num][postion1] = 1;
                    }
                    if ((i + 1) % 5 < 4) {
                        pos1 = stringArray[2 * ((i + 1) % 5) + 9 * ((k + 63) % 64)].find_first_not_of("0123456789", 1);
                        int postion1 = stoi(stringArray[2 * ((i + 1) % 5) + 9 * ((k + 63) % 64)].substr(1, pos1 - 1));
                        matrix[num][postion1] = 1;
                        pos1 = stringArray[2 * ((i + 1) % 5) + 9 * ((k + 63) % 64) + 1].find_first_not_of("0123456789", 1);
                        int postion2 = stoi(stringArray[2 * ((i + 1) % 5) + 9 * ((k + 63) % 64) + 1].substr(1, pos1 - 1));
                        matrix[num][postion2] = 1;
                    }
                    else {
                        pos1 = stringArray[2 * ((i + 1) % 5) + 9 * ((k + 63) % 64)].find_first_not_of("0123456789", 1);
                        int postion1 = stoi(stringArray[2 * ((i + 1) % 5) + 9 * ((k + 63) % 64)].substr(1, pos1 - 1));
                        matrix[num][postion1] = 1;
                    }
                    num++;
                }
            }
        }
    }

    // Footer
    output += "\\makeatother\n";
    output += "\\end{tikzpicture}\n";
    output += "\\end{document}\n";


    const string filename = "matrix.txt";

    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
    }


    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            outFile << matrix[i][j] << " ";
        }
        outFile << endl;
    }

    outFile.close();

    for (int i = 0; i < rows; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;

    return output;
}


void Mitm::getAllFoundSolutions(GRBModel& model, vector<MitmSolution>& solutions) {
    int solCount = model.get(GRB_IntAttr_SolCount);
    cout << solCount << endl;
    for (int solNb = 0; solNb < solCount; solNb++) {
        MitmSolution solution(Rounds);
        getSolution(model, solution);
        solutions.push_back(solution);
    }
}

void Mitm::getSolution(GRBModel& model, MitmSolution& MitmSolution) {
    for (int round = 0; round < Rounds; round++) {
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                for (int k = 0; k < 64; k++) {
                    MitmSolution.DC2_red[round][i][j][k] = static_cast<int>(DC2_red[round][i][j][k].get(GRB_DoubleAttr_Xn));
                    MitmSolution.DC2_blue[round][i][j][k] = static_cast<int>(DC2_blue[round][i][j][k].get(GRB_DoubleAttr_Xn));
                    for (int l = 0; l < 3; l++) {
                        MitmSolution.DA[round][i][j][k][l] = static_cast<int>(DA[round][i][j][k][l].get(GRB_DoubleAttr_Xn));
                        MitmSolution.DB[round][i][j][k][l] = static_cast<int>(DB[round][i][j][k][l].get(GRB_DoubleAttr_Xn));
                    }
                }
            }
        }
    }

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 64; k++)
                for (int l = 0; l < 3; l++) {
                    MitmSolution.DA[Rounds][i][j][k][l] = static_cast<int>(DA[Rounds][i][j][k][l].get(GRB_DoubleAttr_Xn));
                }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                MitmSolution.DC1_red[round][i][k] = static_cast<int>(DC1_red[round][i][k].get(GRB_DoubleAttr_Xn));
                MitmSolution.DC1_blue[round][i][k] = static_cast<int>(DC1_blue[round][i][k].get(GRB_DoubleAttr_Xn));
                MitmSolution.DC12_red[round][i][k] = static_cast<int>(DC12_red[round][i][k].get(GRB_DoubleAttr_Xn));
                MitmSolution.DC12_blue[round][i][k] = static_cast<int>(DC12_blue[round][i][k].get(GRB_DoubleAttr_Xn));

                for (int l = 0; l < 3; l++) {
                    MitmSolution.DP[round][i][k][l] = static_cast<int>(DP[round][i][k][l].get(GRB_DoubleAttr_Xn));
                    MitmSolution.DP2[round][i][k][l] = static_cast<int>(DP2[round][i][k][l].get(GRB_DoubleAttr_Xn));
                }
            }

    for (int i = 0; i < 4; i++)
        for (int k = 0; k < 64; k++)
            MitmSolution.dom1[i][k] = static_cast<int>(dom1[i][k].get(GRB_DoubleAttr_Xn));

    
    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                MitmSolution.beta1[round][i][k] = static_cast<int>(beta1[round][i][k].get(GRB_DoubleAttr_Xn));
            }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                MitmSolution.beta3[round][i][k] = static_cast<int>(beta3[round][i][k].get(GRB_DoubleAttr_Xn));
            }

    for (int round = 0; round < Rounds; round++) {
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                for (int k = 0; k < 64; k++) {
                    MitmSolution.beta2[round][i][j][k] = static_cast<int>(beta2[round][i][j][k].get(GRB_DoubleAttr_Xn));
                }
            }
        }
    }
    

    for (int i = 0; i < 3; i++)
        MitmSolution.obj[i] = static_cast<int>(obj[i].get(GRB_DoubleAttr_Xn));

    MitmSolution.Rounds = Rounds;

    MitmSolution.objective = static_cast<int>(model.get(GRB_DoubleAttr_PoolObjVal));
}

void Mitm::addconstraint(GRBModel& model) {
    MitmFactory factory;
    factory.addfivexor_red(DA, DP, DC1_red, DC1_blue, model);
    factory.addtwoxor_red(DP2, DP, DC12_red, DC12_blue, model);
    factory.addTheta_red(DA, DP2, DB, DC2_red, DC2_blue, model);
    factory.addSbox_nc(DB, DA, model);

    factory.addDoMnew(DA, dom1, model);

    vector<vector<GRBVar>> beta(4, vector<GRBVar>(64));

    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 64; ++k) {
            beta[i][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "beta_" + std::to_string(i) + "_" + std::to_string(k));
        }
    }

    // Assuming factory.betaConstraints is a method that takes DA and beta as parameters and applies constraints
    factory.betaConstraints(DA, beta, model);


    GRBLinExpr objective = 0;
    GRBLinExpr DoF_Lin = 0;
    GRBLinExpr DiF_Lin = 0;
    GRBLinExpr DoF_red = 0;
    GRBLinExpr DoF_blue = 0;
    //GRBLinExpr DoF_red_init = 0;
    //GRBLinExpr DoF_blue_init = 0;
    GRBLinExpr DoM = 0;
    GRBVar Obj = model.addVar(0.0, 256.0, 0.0, GRB_INTEGER, "Obj");
    //GRBVar DoF_red_init1 = model.addVar(0.0, 512.0, 0.0, GRB_INTEGER, "DoF_red_init");
    //GRBVar DoF_blue_init1 = model.addVar(0.0, 512.0, 0.0, GRB_INTEGER, "DoF_blue_init");

    objective += 1.0 * Obj;

    for (int k = 0; k < 64; k++) {
        for (int i = 0; i < 4; i++) {
            DoM += 1.0 * dom1[i][k];
        }
    }

   

    vector<double> t = { 1, 1, -2 };

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                model.addConstr(DC1_red[round][i][k] >= beta1[round][i][k], "");
                model.addConstr(DC1_blue[round][i][k] >= beta1[round][i][k], "");
                model.addConstr(factory.linExprOf(t, { DC1_red[round][i][k], DC1_blue[round][i][k], beta1[round][i][k] }) <= 2, "");

                model.addConstr(DC12_red[round][i][k] >= beta3[round][i][k], "");
                model.addConstr(DC12_blue[round][i][k] >= beta3[round][i][k], "");
                model.addConstr(factory.linExprOf(t, { DC12_red[round][i][k], DC12_blue[round][i][k], beta3[round][i][k] }) <= 2, "");


                for (int j = 0; j < 5; j++) {
                    model.addConstr(DC2_red[round][i][j][k] >= beta2[round][i][j][k], "");
                    model.addConstr(DC2_blue[round][i][j][k] >= beta2[round][i][j][k], "");
                    model.addConstr(factory.linExprOf(t, { DC2_red[round][i][j][k], DC2_blue[round][i][j][k], beta2[round][i][j][k] }) <= 2, "");
                }
            }

    obj[0] = model.addVar(0.0, 512.0, 0.0, GRB_INTEGER, "Obj1");
    obj[1] = model.addVar(0.0, 512.0, 0.0, GRB_INTEGER, "Obj2");

    for (int k = 0; k < 64; k++) {
        DoF_blue += 2.0 * DA[0][0][0][k][0];
        DoF_red += 2.0 * DA[0][0][0][k][2];
        DoF_blue += -2.0 * beta[0][k];
        DoF_red += -2.0 * beta[0][k];

        DoF_blue += 2.0 * DA[0][0][1][k][0];
        DoF_red += 2.0 * DA[0][0][1][k][2];
        DoF_blue += -2.0 * beta[1][k];
        DoF_red += -2.0 * beta[1][k];

        DoF_blue += 2.0 * DA[0][0][2][k][0];
        DoF_red += 2.0 * DA[0][0][2][k][2];
        DoF_blue += -2.0 * beta[2][k];
        DoF_red += -2.0 * beta[2][k];

        DoF_blue += 2.0 * DA[0][0][4][k][0];
        DoF_red += 2.0 * DA[0][0][4][k][2];
        DoF_blue += -2.0 * beta[3][k];
        DoF_red += -2.0 * beta[3][k];

        /*DoF_blue_init += 2.0 * DA[0][0][0][k][0];
        DoF_red_init += 2.0 * DA[0][0][0][k][2];
        DoF_blue_init += -2.0 * beta[0][k];
        DoF_red_init += -2.0 * beta[0][k];

        DoF_blue_init += 2.0 * DA[0][0][1][k][0];
        DoF_red_init += 2.0 * DA[0][0][1][k][2];
        DoF_blue_init += -2.0 * beta[1][k];
        DoF_red_init += -2.0 * beta[1][k];

        DoF_blue_init += 2.0 * DA[0][0][2][k][0];
        DoF_red_init += 2.0 * DA[0][0][2][k][2];
        DoF_blue_init += -2.0 * beta[2][k];
        DoF_red_init += -2.0 * beta[2][k];

        DoF_blue_init += 2.0 * DA[0][0][4][k][0];
        DoF_red_init += 2.0 * DA[0][0][4][k][2];
        DoF_blue_init += -2.0 * beta[3][k];
        DoF_red_init += -2.0 * beta[3][k];*/
    }

    for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 64; k++) {
                DoF_red += -1.0 * DC1_red[round][i][k];
                DoF_red += 1.0 * beta1[round][i][k];
                DoF_blue += -1.0 * DC1_blue[round][i][k];
                DoF_blue += 1.0 * beta1[round][i][k];
                DoF_red += -1.0 * DC12_red[round][i][k];
                DoF_red += 1.0 * beta3[round][i][k];
                DoF_blue += -1.0 * DC12_blue[round][i][k];
                DoF_blue += 1.0 * beta3[round][i][k];
                for (int j = 0; j < 5; j++) {
                    DoF_red += -1.0 * DC2_red[round][i][j][k];
                    DoF_red += 1.0 * beta2[round][i][j][k];
                    DoF_blue += -1.0 * DC2_blue[round][i][j][k];
                    DoF_blue += 1.0 * beta2[round][i][j][k];
                }
            }

    model.addConstr(DoF_red >= 2, "");
    model.addConstr(DoF_blue >= 2, "");
    DoF_Lin += DoF_red;
    DoF_Lin += DoF_blue;
    DoF_Lin -= 2.0 * DoM;
    model.addConstr(DoF_Lin >= 4, "");
    DiF_Lin += DoF_blue;
    DiF_Lin -= DoF_red;
    model.addConstr(DiF_Lin <= 6, "");
    model.addConstr(DiF_Lin >= -6, "");
    //  model.addConstr(2.0 * objective <= DoF_blue, "");
    //   model.addConstr(2.0 * objective <= DoF_red, "");
    


    obj[2] = model.addVar(0.0, 256.0, 0.0, GRB_INTEGER, "Obj3");

    //for (int i = 0; i < 2; i++) 
    //    for (int k = 0; k < 64; k++) {
    //        DoM += 1.0*dom1[i][k];
    //    }

    model.addConstr(DoM >= 1, "");
    model.addConstr(DoF_red == obj[0], "");
    model.addConstr(DoF_blue == obj[1], "");
    model.addConstr(DoM == obj[2], "");
    model.addConstr(objective <= DoM, "");
    //model.addConstr(DoF_red_init1 == DoF_red_init, "");
    //model.addConstr(DoF_blue_init1 == DoF_blue_init, "");
    model.addConstr(objective, GRB_LESS_EQUAL, TESTOBJ, "");
    model.setObjective(objective, GRB_MAXIMIZE);
}

void MitmFactory::addfivexor_red(vector<vector<vector<vector<vector<GRBVar>>>>> DA,
    vector<vector<vector<vector<GRBVar>>>> DP,
    vector<vector<vector<GRBVar>>> DC1_red, vector<vector<vector<GRBVar>>> DC1_blue, GRBModel& model) {

    std::vector<std::vector<std::vector<std::vector<GRBVar>>>> DA_Allone(DP.size(), std::vector<std::vector<std::vector<GRBVar>>>(5, std::vector<std::vector<GRBVar>>(64, std::vector<GRBVar>(3))));

    for (int round = 0; round < DP.size(); ++round) {
        for (int i = 0; i < 5; ++i) {
            for (int k = 0; k < 64; ++k) {
                for (int l = 0; l < 3; ++l) {
                    DA_Allone[round][i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "DA_Allone_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(l));

                    if (round == 0) {
                        if (i == 0) {
                            if (l != 1)
                                Determine_AllOne(DA_Allone[round][i][k][l], { DA[round][i][0][k][l], DA[round][i][1][k][l], DA[round][i][2][k][l], DA[round][i][4][k][l] }, model);
                            else
                                model.addConstr(DA_Allone[round][i][k][l] == 1, "");
                        }
                        else if (i == 1) {
                            if (l != 1)
                                Determine_AllOne(DA_Allone[round][i][k][l], { DA[round][i][0][k][l], DA[round][i][2][k][l], DA[round][i][3][k][l], DA[round][i][4][k][l] }, model);
                            else
                                model.addConstr(DA_Allone[round][i][k][l] == 1, "");
                        }
                        else if (i == 4) {
                            if (l != 1)
                                Determine_AllOne(DA_Allone[round][i][k][l], { DA[round][i][0][k][l], DA[round][i][2][k][l], DA[round][i][4][k][l] }, model);
                            else
                                model.addConstr(DA_Allone[round][i][k][l] == 1, "");
                        }
                        else
                            model.addConstr(DA_Allone[round][i][k][l] == 1, "");
                    }
                    else {
                        Determine_AllOne(DA_Allone[round][i][k][l], { DA[round][i][0][k][l], DA[round][i][1][k][l], DA[round][i][2][k][l], DA[round][i][3][k][l], DA[round][i][4][k][l] }, model);
                    }
                }

                model.addConstr(DA_Allone[round][i][k][1] == DP[round][i][k][1], "");
                model.addConstr(DA_Allone[round][i][k][1] >= DP[round][i][k][0], "");
                model.addConstr(DA_Allone[round][i][k][0] <= DP[round][i][k][0], "");
                model.addConstr(DA_Allone[round][i][k][2] <= DP[round][i][k][2], "");
                model.addConstr(DA_Allone[round][i][k][1] >= DP[round][i][k][2], "");

                
                vector<double> t1 = { 1.0, 2.0, -2.0 };
                model.addConstr(linExprOf(t1, { DC1_red[round][i][k], DA_Allone[round][i][k][0], DP[round][i][k][0] }) == 0, "");
                model.addConstr(linExprOf(t1, { DC1_blue[round][i][k], DA_Allone[round][i][k][2], DP[round][i][k][2] }) == 0, "");
                

            }
        }
    }
}

void MitmFactory::addtwoxor_red(vector<vector<vector<vector<GRBVar>>>> DP2,
    vector<vector<vector<vector<GRBVar>>>> DP,
    vector<vector<vector<GRBVar>>> DC12_red, vector<vector<vector<GRBVar>>> DC12_blue, GRBModel& model) {

    std::vector<std::vector<std::vector<std::vector<GRBVar>>>> DP_Allone(DP.size(), std::vector<std::vector<std::vector<GRBVar>>>(5, std::vector<std::vector<GRBVar>>(64, std::vector<GRBVar>(3))));


    for (int round = 0; round < DP.size(); ++round) {
        for (int i = 0; i < 5; ++i) {
            for (int k = 0; k < 64; ++k) {
                for (int l = 0; l < 3; ++l) {
                    DP_Allone[round][i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "DP_Allone_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(l));
                    Determine_AllOne(DP_Allone[round][i][k][l], { DP[round][(i + 4) % 5][k][l], DP[round][(i + 1) % 5][(k + 63) % 64][l] }, model);
                }

                model.addConstr(DP_Allone[round][i][k][1] == DP2[round][i][k][1], "");
                model.addConstr(DP_Allone[round][i][k][1] >= DP2[round][i][k][0], "");
                model.addConstr(DP_Allone[round][i][k][0] <= DP2[round][i][k][0], "");
                model.addConstr(DP_Allone[round][i][k][2] <= DP2[round][i][k][2], "");
                model.addConstr(DP_Allone[round][i][k][1] >= DP2[round][i][k][2], "");


                vector<double> t1 = { 1.0, 2.0, -2.0 };
                model.addConstr(linExprOf(t1, { DC12_red[round][i][k], DP_Allone[round][i][k][0], DP2[round][i][k][0] }) == 0, "");
                model.addConstr(linExprOf(t1, { DC12_blue[round][i][k], DP_Allone[round][i][k][2], DP2[round][i][k][2] }) == 0, "");
                

            }
        }
    }

}


void MitmFactory::addTheta_red(vector<vector<vector<vector<vector<GRBVar>>>>> DA,
    vector<vector<vector<vector<GRBVar>>>> DP2,
    vector<vector<vector<vector<vector<GRBVar>>>>> DB,
    vector<vector<vector<vector<GRBVar>>>> DC2_red, vector<vector<vector<vector<GRBVar>>>> DC2_blue, GRBModel& model) {

    std::vector<std::vector<std::vector<std::vector<std::vector<GRBVar>>>>> DP2_Allone(DP2.size(), std::vector<std::vector<std::vector<std::vector<GRBVar>>>>(5, std::vector<std::vector<std::vector<GRBVar>>>(5, std::vector<std::vector<GRBVar>>(64, std::vector<GRBVar>(3)))));

    for (int round = 0; round < DP2.size(); ++round) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                for (int k = 0; k < 64; ++k) {
                    for (int l = 0; l < 3; ++l) {
                        DP2_Allone[round][i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "DP2_Allone_" + std::to_string(round) + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" + std::to_string(l));
                        Determine_AllOne(DP2_Allone[round][i][j][k][l], { DA[round][i][j][k][l], DP2[round][i][k][l] }, model);
                    }

                    model.addConstr(DP2_Allone[round][i][j][k][1] == DB[round][j][(2 * i + 3 * j) % 5][(k + rho[i][j]) % 64][1], "");
                    model.addConstr(DP2_Allone[round][i][j][k][1] >= DB[round][j][(2 * i + 3 * j) % 5][(k + rho[i][j]) % 64][0], "");
                    model.addConstr(DP2_Allone[round][i][j][k][0] <= DB[round][j][(2 * i + 3 * j) % 5][(k + rho[i][j]) % 64][0], "");
                    model.addConstr(DP2_Allone[round][i][j][k][1] >= DB[round][j][(2 * i + 3 * j) % 5][(k + rho[i][j]) % 64][2], "");
                    model.addConstr(DP2_Allone[round][i][j][k][2] <= DB[round][j][(2 * i + 3 * j) % 5][(k + rho[i][j]) % 64][2], "");

                   
                    vector<double> t1 = { 1.0, 2.0, -2.0 };
                    model.addConstr(linExprOf(t1, { DC2_red[round][i][j][k], DP2_Allone[round][i][j][k][0], DB[round][j][(2 * i + 3 * j) % 5][(k + rho[i][j]) % 64][0] }) == 0, "");
                    model.addConstr(linExprOf(t1, { DC2_blue[round][i][j][k], DP2_Allone[round][i][j][k][2], DB[round][j][(2 * i + 3 * j) % 5][(k + rho[i][j]) % 64][2] }) == 0, "");
                    
                }
            }
        }
    }
}


void MitmFactory::addSbox_nc(vector<vector<vector<vector<vector<GRBVar>>>>> DB,
    vector<vector<vector<vector<vector<GRBVar>>>>> DA, GRBModel& model) {

    vector<vector<double>> t = {
    {1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1},
            {0, 1, -1, 1, 0, 0, 1, 0, 0, -1, -1, 1},
            {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1},
            {0, -1, 0, 0, 0, -1, 0, 0, -1, 0, 1, 0},
            {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0},
            {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0},
            {0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -1, 0},
            {0, 1, -1, 0, 1, -1, 0, 1, -1, 0, -1, 1},
            {0, -1, 0, -1, 0, 0, -1, 0, 0, 0, 1, 0},
            {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1},
            {-1, 1, 0, -1, 1, 0, -1, 1, 0, 1, -1, 0},
            {0, -1, 0, 0, -1, 0, -1, 1, -1, 0, 1, 0},
            {0, -1, 0, -1, 1, -1, 0, -1, 0, 0, 1, 0},
            {0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0},
            {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0},
            {0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0}
    };

    double c[19] = { 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, -2, 0, 0, -2, -2, 0, 0, 0, 0 };

    for (size_t round = 0; round < DB.size(); ++round) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                for (int k = 0; k < 64; ++k) {
                    for (int q = 0; q < 19; ++q) {
                        model.addConstr(linExprOf(t[q],
                            { DB[round][i][j][k][0], DB[round][i][j][k][1], DB[round][i][j][k][2],
                            DB[round][(i + 1) % 5][j][k][0], DB[round][(i + 1) % 5][j][k][1], DB[round][(i + 1) % 5][j][k][2],
                            DB[round][(i + 2) % 5][j][k][0], DB[round][(i + 2) % 5][j][k][1], DB[round][(i + 2) % 5][j][k][2],
                            DA[round + 1][i][j][k][0], DA[round + 1][i][j][k][1], DA[round + 1][i][j][k][2] }
                        ) >= c[q], "");
                    }
                }
            }
        }
    }
}


void MitmFactory::addDoMnew(vector<vector<vector<vector<vector<GRBVar>>>>> DA,
    vector<vector<GRBVar>> dom1, GRBModel& model) {

    const int rows1 = 4;
    const int cols = 64;

    std::vector<std::vector<GRBVar>> Nowhite1(rows1, std::vector<GRBVar>(cols));
    std::vector<std::vector<GRBVar>> Existred1(rows1, std::vector<GRBVar>(cols));
    std::vector<std::vector<GRBVar>> Existblue1(rows1, std::vector<GRBVar>(cols));


    for (int i = 0; i < rows1; i++) {
        for (int k = 0; k < cols; k++) {
            Nowhite1[i][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "Nowhite1_" + std::to_string(i) + "_" + std::to_string(k));
            Existred1[i][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "Existred1_" + std::to_string(i) + "_" + std::to_string(k));
            Existblue1[i][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "Existblue1_" + std::to_string(i) + "_" + std::to_string(k));
        }
    }

    int r = DA.size() - 1;

    for (int k = 0; k < 64; k++) {
        Determine_AllOne(Nowhite1[0][k], { DA[r][0][2][(k - rho[0][2] + 64) % 64][1], DA[r][0][0][(k - rho[0][2] + 64) % 64][1] }, model);
        Determine_AllOne(Nowhite1[1][k], { DA[r][1][3][(k - rho[1][3] + 64) % 64][1], DA[r][1][1][(k - rho[1][3] + 64) % 64][1] }, model);
        Determine_AllOne(Nowhite1[2][k], { DA[r][3][3][(k - rho[3][0] + 64) % 64][1], DA[r][3][0][(k - rho[3][0] + 64) % 64][1] }, model);
        Determine_AllOne(Nowhite1[3][k], { DA[r][4][4][(k - rho[4][1] + 64) % 64][1], DA[r][4][1][(k - rho[4][1] + 64) % 64][1] }, model);
        
        Determine_Existzero(Existred1[0][k], { DA[r][0][2][(k - rho[0][2] + 64) % 64][0], DA[r][0][0][(k - rho[0][2] + 64) % 64][0] }, model);
        Determine_Existzero(Existred1[1][k], { DA[r][1][3][(k - rho[1][3] + 64) % 64][0], DA[r][1][1][(k - rho[1][3] + 64) % 64][0] }, model);
        Determine_Existzero(Existred1[2][k], { DA[r][3][3][(k - rho[3][0] + 64) % 64][0], DA[r][3][0][(k - rho[3][0] + 64) % 64][0] }, model);
        Determine_Existzero(Existred1[3][k], { DA[r][4][4][(k - rho[4][1] + 64) % 64][0], DA[r][4][1][(k - rho[4][1] + 64) % 64][0] }, model);

        Determine_Existzero(Existblue1[0][k], { DA[r][0][2][(k - rho[0][2] + 64) % 64][2], DA[r][0][0][(k - rho[0][2] + 64) % 64][2] }, model);
        Determine_Existzero(Existblue1[1][k], { DA[r][1][3][(k - rho[1][3] + 64) % 64][2], DA[r][1][1][(k - rho[1][3] + 64) % 64][2] }, model);
        Determine_Existzero(Existblue1[2][k], { DA[r][3][3][(k - rho[3][0] + 64) % 64][2], DA[r][3][0][(k - rho[3][0] + 64) % 64][2] }, model);
        Determine_Existzero(Existblue1[3][k], { DA[r][4][4][(k - rho[4][1] + 64) % 64][2], DA[r][4][1][(k - rho[4][1] + 64) % 64][2] }, model);


        Determine_AllOne(dom1[0][k], { Nowhite1[0][k], Existred1[0][k], Existblue1[0][k] }, model);
        Determine_AllOne(dom1[1][k], { Nowhite1[1][k], Existred1[1][k], Existblue1[1][k] }, model);
        Determine_AllOne(dom1[2][k], { Nowhite1[2][k], Existred1[2][k], Existblue1[2][k] }, model);
        Determine_AllOne(dom1[3][k], { Nowhite1[3][k], Existred1[3][k], Existblue1[3][k] }, model);
    }
}


void MitmFactory::betaConstraints(vector<vector<vector<vector<vector<GRBVar>>>>> DA,
    vector<vector<GRBVar>> beta, GRBModel& model) {
    vector<double> t = { 1, 1, -2 };

    for (int k = 0; k < 64; k++) {
        model.addConstr(DA[0][0][0][k][0] >= beta[0][k], "");
        model.addConstr(DA[0][0][0][k][2] >= beta[0][k], "");
        model.addConstr(linExprOf(t, { DA[0][0][0][k][0], DA[0][0][0][k][2], beta[0][k] }) <= 1, "");

        model.addConstr(DA[0][0][1][k][0] >= beta[1][k], "");
        model.addConstr(DA[0][0][1][k][2] >= beta[1][k], "");
        model.addConstr(linExprOf(t, { DA[0][0][1][k][0], DA[0][0][1][k][2], beta[1][k] }) <= 1, "");

        model.addConstr(DA[0][0][2][k][0] >= beta[2][k], "");
        model.addConstr(DA[0][0][2][k][2] >= beta[2][k], "");
        model.addConstr(linExprOf(t, { DA[0][0][2][k][0], DA[0][0][2][k][2], beta[2][k] }) <= 1, "");

        model.addConstr(DA[0][0][4][k][0] >= beta[3][k], "");
        model.addConstr(DA[0][0][4][k][2] >= beta[3][k], "");
        model.addConstr(linExprOf(t, { DA[0][0][4][k][0], DA[0][0][4][k][2], beta[3][k] }) <= 1, "");
    }
}

void MitmFactory::Determine_Allzero(GRBVar mainVar, initializer_list<GRBVar> vars, GRBModel& model) {
    GRBLinExpr expr = 0;
    expr += 1.0 * mainVar;
    for (auto& var : vars) {
        expr += 1.0 * var;
        model.addConstr(linExprOf({ mainVar,  var }) <= 1, "");
    }

    model.addConstr(expr >= 1, "");
}

void MitmFactory::Determine_AllOne(GRBVar mainVar, initializer_list<GRBVar> vars, GRBModel& model) {
    int size = vars.size();
    GRBLinExpr exprm = 0;
    exprm += 1.0 * mainVar;

    GRBLinExpr exprp = 0;
    exprp += -1.0 * size * mainVar;

    for (auto& var : vars) {
        exprm += -1.0 * var;
        exprp += 1.0 * var;
    }

    model.addConstr(exprm >= 1 - size, "");
    model.addConstr(exprp >= 0, "");
}

void MitmFactory::Determine_ExistOne(GRBVar mainVar, initializer_list<GRBVar> vars, GRBModel& model) {
    GRBLinExpr expr1 = 0;
    expr1 += 1.0 * mainVar;
    GRBLinExpr exprm = 0;
    exprm += vars.size() * mainVar;
    for (auto& var : vars) {
        expr1 += -1.0 * var;
        exprm += -1.0 * var;
    }

    model.addConstr(exprm >= 0, "");
    model.addConstr(expr1 <= 0, "");
}

void MitmFactory::Determine_Existzero(GRBVar mainVar, initializer_list<GRBVar> vars, GRBModel& model) {
    int size = vars.size();
    GRBLinExpr expr = 0;
    expr += 1.0 * mainVar;
    for (auto& var : vars) {
        expr += 1.0 * var;
        model.addConstr(linExprOf({ mainVar, var }) >= 1, "");
    }

    model.addConstr(expr <= size, "");
}

GRBLinExpr MitmFactory::linExprOf(vector<double> coeffs, initializer_list<GRBVar> vars) {
    int size = vars.size();
    GRBLinExpr ofVars = 0;
    GRBVar* myArray = new GRBVar[size];
    copy(vars.begin(), vars.end(), myArray);
    ofVars.addTerms(coeffs.data(), myArray, size);
    delete[] myArray;
    return ofVars;
}

GRBLinExpr MitmFactory::linExprOf(double constant, initializer_list<GRBVar> vars) {
    GRBLinExpr ofVars = linExprOf(vars);
    ofVars += constant;
    return ofVars;
}

GRBLinExpr MitmFactory::linExprOf(initializer_list<GRBVar> vars) {
    GRBLinExpr expr = 0;
    for (GRBVar var : vars)
        expr += 1.0 * var;
    return expr;
}

MitmFactory::MitmFactory() {

}