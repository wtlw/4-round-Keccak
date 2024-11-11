#include "MitmSHA3.h"



int main(int argc, char* argv[]) {
    GRBEnv* env =  new GRBEnv();                        //创建环境
    int Rounds = 2;                                     //设置搜索轮数
   // env->set(GRB_IntParam_OutputFlag, 0);
   // env->set(GRB_DoubleParam_TimeLimit, 100.0);
    env->set(GRB_IntParam_Threads, 8);

    GRBModel model = GRBModel(*env);                    //创建新的MILP模型

    Mitm mitm(Rounds);                                  //创建对应的数组

    mitm.Value(model);                                  //给每一个数组赋值，并给出初始变量间的关系
    mitm.article_attack(model);
    mitm.addconstraint(model);                          //建立变量传播模型
    mitm.solve(model);                                  //求解模型，并给出tikz代码

    delete env;
    
    return 0;
}



