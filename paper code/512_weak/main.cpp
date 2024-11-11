#include "MitmSHA3.h"



int main(int argc, char* argv[]) {
    GRBEnv* env =  new GRBEnv();                        
    int Rounds = 2;                                     
   // env->set(GRB_IntParam_OutputFlag, 0);
   // env->set(GRB_DoubleParam_TimeLimit, 100.0);
   // env->set(GRB_IntParam_Threads, 8);

    GRBModel model = GRBModel(*env);                    

    Mitm mitm(Rounds);                                  

    mitm.Value(model);                                  
    mitm.article_attack(model);
    mitm.addconstraint(model);                          
    mitm.solve(model);                                  

    delete env;
    
    return 0;
}



