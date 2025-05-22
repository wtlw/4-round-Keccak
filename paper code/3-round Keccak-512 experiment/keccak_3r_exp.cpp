//
//  main.cpp
//  3r experiments on keccak-512
//  


#include <iostream>
#include <cstdlib>
#include <ctime>
#include <map>
#include <cmath>
#include <vector>

using namespace std;

typedef unsigned char UINT8;
typedef unsigned long int UINT32;
typedef unsigned long long int UINT64;


#define random(x) (rand())%x;
#define nrRounds 3
UINT64 KeccakRoundConstants[nrRounds];//these are constant,
#define nrLanes 25
unsigned int KeccakRhoOffsets[nrLanes];//these are constant,

#define index(x, y) (((x)%5)+5*((y)%5))
#define ROL64(a, offset) ((offset != 0) ? ((((UINT64)a) << offset) ^ (((UINT64)a) >> (64-offset))) : a)
#define ROR64(a, offset) ((offset != 0) ? ((((UINT64)a) >> offset) ^ (((UINT64)a) << (64-offset))) : a)
#define TAKE_BIT(x, pos) (((x) >> (pos)) & 0x1)


void KeccakPermutationOnWords(UINT64 *state);
void theta(UINT64 *A);
void rho(UINT64 *A);
void rho_inv(UINT64 *A);
void pi(UINT64 *A);
void chi(UINT64 *A);
void chi_inv(UINT64 *A);
void iota(UINT64 *A, unsigned int indexRound);


void KeccakPermutationOnWords(UINT64 state[], int round)
{
    unsigned int i;


    for(i=0; i<round; i++) {
        theta(state);
        rho(state);
        pi(state);
        chi(state);
        iota(state, i);
    }
}


void theta(UINT64 *A)
{
    unsigned int x, y;
    UINT64 C[5], D[5];//C are the Xors of the five bits in every column. D are the Xors of the ten bits in right-behind column and right column

    for(x=0; x<5; x++) {
        C[x] = 0;
        for(y=0; y<5; y++)
            C[x] ^= A[index(x, y)];
    }
    for(x=0; x<5; x++)
        D[x] = ROL64(C[(x+1)%5], 1) ^ C[(x+4)%5];
    for(x=0; x<5; x++)
        for(y=0; y<5; y++)
            A[index(x, y)] ^= D[x];
}

void rho(UINT64 *A)
{
    unsigned int x, y;

    for(x=0; x<5; x++)
        for(y=0; y<5; y++)
            A[index(x, y)] = ROL64(A[index(x, y)], KeccakRhoOffsets[index(x, y)]);
}

void rho_inv1(UINT64 *A)
{
    unsigned int x;

    for(x=0; x<5; x++)
        A[index(x, 0)] = ROR64(A[index(x, 0)], KeccakRhoOffsets[index(x, x)]);
}

void rho_inv2(UINT64 *A)
{
    unsigned int x;

    for(x=0; x<5; x++)
        A[index(x+2, 1)] = ROR64(A[index(x+2, 1)], KeccakRhoOffsets[index(x, x+2)]);
}

void pi(UINT64 *A)
{
    unsigned int x, y;
    UINT64 tempA[25];

    for(x=0; x<5; x++)
        for(y=0; y<5; y++)
            tempA[index(x, y)] = A[index(x, y)];
    for(x=0; x<5; x++)
        for(y=0; y<5; y++)
            A[index(0*x+1*y, 2*x+3*y)] = tempA[index(x, y)];
}

void pi_inv(UINT64 *A)
{
    unsigned int x;
    UINT64 tempA[25];

    for(x=0; x<5; x++)
        tempA[index(x, 0)] = A[index(x, 0)];
    for(x=0; x<5; x++)
        A[index(x, x)] = tempA[index(x, 0)];
}

void chi(UINT64 *A)
{
    unsigned int x, y;
    UINT64 C[5];

    for(y=0; y<5; y++) {
        for(x=0; x<5; x++)
            C[x] = A[index(x, y)] ^ ((~A[index(x+1, y)]) & A[index(x+2, y)]);
        for(x=0; x<5; x++)
            A[index(x, y)] = C[x];
    }
}

void chi_inv(UINT64 *A)
{
    unsigned int x;
    UINT64 C[5];
    
    for(x=0; x<5; x++)
        C[x] = A[index(x, 0)] ^ ((~A[index(x+1, 0)]) & (A[index(x+2, 0)] ^ ((~A[index(x+3, 0)]) & A[index(x+4, 0)])));
    for(x=0; x<5; x++)
        A[index(x, 0)] = C[x];
    
}

void iota(UINT64 *A, unsigned int indexRound)
{
    A[index(0, 0)] ^= KeccakRoundConstants[indexRound];
}



int LFSR86540(UINT8 *LFSR)
{
    int result = ((*LFSR) & 0x01) != 0;
    if (((*LFSR) & 0x80) != 0)
        // Primitive polynomial over GF(2): x^8+x^6+x^5+x^4+1
        (*LFSR) = ((*LFSR) << 1) ^ 0x71;
    else
        (*LFSR) <<= 1;
    return result;
}

void KeccakInitializeRoundConstants()
{
    UINT8 LFSRstate = 0x01;
    unsigned int i, j, bitPosition;

    for(i=0; i<nrRounds; i++) {
        KeccakRoundConstants[i] = 0;
        for(j=0; j<7; j++) {
            bitPosition = (1<<j)-1; //2^j-1
            if (LFSR86540(&LFSRstate))
                KeccakRoundConstants[i] ^= ((UINT64)1<<bitPosition);
        }
    }
}

void KeccakInitializeRhoOffsets()
{
    unsigned int x, y, t, newX, newY;

    KeccakRhoOffsets[index(0, 0)] = 0;
    x = 1;
    y = 0;
    for(t=0; t<24; t++) {
        KeccakRhoOffsets[index(x, y)] = ((t+1)*(t+2)/2) % 64;
        newX = (0*x+1*y) % 5;
        newY = (2*x+3*y) % 5;
        x = newX;
        y = newY;
    }
    //for (x=0;x<5;x++)
        //for (y=0;y<5;y++)
            //cout << KeccakRhoOffsets[index(x, y)] << endl;
}

void KeccakInitialize()
{
    KeccakInitializeRoundConstants();
    KeccakInitializeRhoOffsets();
}

void displaystate(UINT64 *state, int lanes)
{
    unsigned int i;
    for(i=0;i<lanes;i++)
    {
        //printf("%08x ",(unsigned int)(state[i]));
        printf("%016llx ",(state[i]));
        if((i+1)%5==0) printf("\n");
    }
    printf("\n");
}


int main(int argc, const char * argv[])
{
    srand((unsigned)time(NULL));
    
    KeccakInitialize();
    
    unsigned int indexb[6][2]={0};
    unsigned int indexr[6][2]={0};
    unsigned int index1[24][2]={0};
    unsigned int index0[24][2]={0};
    
    indexr[0][0]=1; indexr[0][1]=25;
    indexr[1][0]=0; indexr[1][1]=26;
    indexr[2][0]=2; indexr[2][1]=47;
    indexr[3][0]=0; indexr[3][1]=53;
    indexr[4][0]=3; indexr[4][1]=62;
    indexr[5][0]=2; indexr[5][1]=63;
    

    indexb[0][0]=3; indexb[0][1]=16;
    indexb[1][0]=1; indexb[1][1]=27;
    indexb[2][0]=3; indexb[2][1]=27;
    indexb[3][0]=1; indexb[3][1]=43;
    indexb[4][0]=0; indexb[4][1]=55;
    indexb[5][0]=2; indexb[5][1]=57;
    
    
    index0[0][0]=8; index0[0][1]=0; //+5 +6
    index0[1][0]=8; index0[1][1]=6;
    index0[2][0]=9; index0[2][1]=6;
    index0[3][0]=6; index0[3][1]=9;
    index0[4][0]=6; index0[4][1]=11;
    index0[5][0]=14; index0[5][1]=14;
    index0[6][0]=11; index0[6][1]=15;
    index0[7][0]=11; index0[7][1]=17;
    
    index0[8][0]=7; index0[8][1]=20;//+5 +38
    index0[9][0]=7; index0[9][1]=22;
    index0[10][0]=9; index0[10][1]=24;
    index0[11][0]=12; index0[11][1]=26;
    index0[12][0]=12; index0[12][1]=28;
    index0[13][0]=13; index0[13][1]=28;
    index0[14][0]=14; index0[14][1]=32;
    index0[15][0]=9; index0[15][1]=35;

    index0[16][0]=7; index0[16][1]=38;
    index0[17][0]=13; index0[17][1]=38;//+5 +38
    index0[18][0]=14; index0[18][1]=43;
    index0[19][0]=12; index0[19][1]=44;
    index0[20][0]=13; index0[20][1]=44;
    index0[21][0]=6; index0[21][1]=46;
    index0[22][0]=11; index0[22][1]=52;
    index0[23][0]=8; index0[23][1]=54;

    
    index1[0][0]=4; index1[0][1]=0; //-20 -41
    index1[1][0]=1; index1[1][1]=4;
    index1[2][0]=0; index1[2][1]=5;
    index1[3][0]=0; index1[3][1]=7;
    index1[4][0]=20; index1[4][1]=8;
    index1[5][0]=2; index1[5][1]=9;
    index1[6][0]=20; index1[7][1]=10;
    index1[7][0]=24; index1[7][1]=12;
    
    index1[8][0]=2; index1[8][1]=20;//-20 -55
    index1[9][0]=0; index1[9][1]=23;
    index1[10][0]=20; index1[10][1]=26;
    index1[11][0]=22; index1[11][1]=29;
    index1[12][0]=4; index1[12][1]=35;
    index1[13][0]=24; index1[13][1]=39;
    index1[14][0]=24; index1[14][1]=41;
    index1[15][0]=21; index1[15][1]=43;

    index1[16][0]=22; index1[8][1]=47;//-20 -55
    index1[17][0]=1; index1[9][1]=52;
    index1[18][0]=21; index1[10][1]=53;
    index1[19][0]=2; index1[11][1]=55;
    index1[20][0]=22; index1[12][1]=58;
    index1[21][0]=21; index1[13][1]=59;
    index1[22][0]=1; index1[14][1]=62;
    index1[23][0]=4; index1[15][1]=62;
    
    
    UINT64 InitialState[25]={0};
    UINT64 TempState[25]={0};
    UINT64 FinalState[10]={0};
    UINT64 ThetaState1[5]={0};
    UINT64 ThetaState2[5]={0};

    
    UINT64 rightkey=0;
    
    //Init the 1600-bit state with 0
    for(UINT64 i=0;i<25;i++){
        InitialState[i]=0;
    }
    //Randomly chosse the gray bits A_{0,0,z}=A_{0,1,z}, A_{2,0,z}=A_{2,1,z},
    for(UINT64 i=0;i<64;i++){
        UINT64 temp=random(2);
        if(temp){
            InitialState[0] |= (UINT64(1)<<i);
            InitialState[5] |= (UINT64(1)<<i);
        }
        temp=random(2);
        if(temp){
            InitialState[1] |= (UINT64(1)<<i);
            InitialState[6] |= (UINT64(1)<<i);
        }
        temp=random(2);
        if(temp){
            InitialState[2] |= (UINT64(1)<<i);
            InitialState[7] |= (UINT64(1)<<i);
        }
        temp=random(2);
        if(temp){
            InitialState[3] |= (UINT64(1)<<i);
            InitialState[8] |= (UINT64(1)<<i);
        }
    }
    
    //Randomly set red bits and blue bits in M2
    for(UINT64 i=0;i<6;i++){
        UINT64 temp=random(2);
        if(temp){
            InitialState[indexr[i][0]] |= (UINT64(1)<<indexr[i][1]);
            InitialState[indexr[i][0]+5] |= (UINT64(1)<<indexr[i][1]);
        }
        else {
            InitialState[indexr[i][0]] &= ROL64(~UINT64(1),indexr[i][1]);
            InitialState[indexr[i][0]+5] &= ROL64(~UINT64(1),indexr[i][1]);
        }
        //generate the right key
        rightkey = (rightkey<<1) ^ temp;
    }
    
    for(UINT64 i=0;i<6;i++){
        UINT64 temp=random(2);
        if(temp){
            InitialState[indexb[i][0]] |= (UINT64(1)<<indexb[i][1]);
            InitialState[indexb[i][0]+5] |= (UINT64(1)<<indexb[i][1]);
        }
        else {
            InitialState[indexb[i][0]] &= ROL64(~UINT64(1),indexb[i][1]);
            InitialState[indexb[i][0]+5] &= ROL64(~UINT64(1),indexb[i][1]);
        }
        //generate the right key
        rightkey = (rightkey<<1) ^ temp;
    }
    
    //Set two bits padding
    InitialState[8] |= (UINT64(1)<<63);
    InitialState[3] |= (UINT64(1)<<63);
    InitialState[8] |= (UINT64(1)<<62);
    InitialState[3] |= (UINT64(1)<<62);
    
    //Set 64 bit conditions satified with M1
    for(UINT64 i=0;i<24;i++)
    {
        //set conditions=1
        InitialState[index1[i][0]] |= (UINT64(1)<<index1[i][1]);
        
        //set conditions=0
        InitialState[index0[i][0]] &= ROL64(~UINT64(1),index0[i][1]);
    }

    //displaystate(InitialState,25);
    
    //Computer the hash value
    for(UINT64 i=0;i<25;i++){
        TempState[i]=InitialState[i];
    }
    KeccakPermutationOnWords(TempState, 3);
    //displaystate(TempState,25);
    for(UINT64 i=0;i<10;i++){
        FinalState[i]=TempState[i];
    }
    
    //Inverse the Sbox
    for(UINT64 i=0;i<5;i++){
        ThetaState1[i]=TempState[i];
        ThetaState2[i]=TempState[i+5];
    }
    iota(ThetaState1, 2);
    chi_inv(ThetaState1);
    chi_inv(ThetaState2);
    rho_inv1(ThetaState1);
    rho_inv2(ThetaState2);
   
    map<UINT64, vector<UINT64>> TableU;
    //Generate table U
    for(UINT64 j=0;j<(UINT64(1)<<6);j++){
        for(UINT64 k=0;k<25;k++){
            TempState[k]=InitialState[k];
        }
        //Set blue=0
        for(UINT64 k=0;k<6;k++){
            TempState[indexb[k][0]] &= ROL64(~UINT64(1),indexb[k][1]);
            TempState[indexb[k][0]+5] &= ROL64(~UINT64(1),indexb[k][1]);
        }
        //Traverse the red bits
        for(UINT64 k=0;k<6;k++) {
            UINT64 temp1=(j>>(5-k))&1;
            if(temp1){
                TempState[indexr[k][0]] |= (UINT64(1)<<indexr[k][1]);
                TempState[indexr[k][0]+5] |= (UINT64(1)<<indexr[k][1]);
            }
            else {
                TempState[indexr[k][0]] &= ROL64(~UINT64(1),indexr[k][1]);
                TempState[indexr[k][0]+5] &= ROL64(~UINT64(1),indexr[k][1]);
            }
                    
        }

        //Compute A^(2)
        KeccakPermutationOnWords(TempState, 2);

        
        UINT64 Matchpoint[8] = {0};
        UINT64 MatchRandG = 0;
        //Compute f'm
        Matchpoint[0] = TAKE_BIT(TempState[12],31) ^ TAKE_BIT(TempState[22],31) ^ TAKE_BIT(ThetaState1[2],31) ^ TAKE_BIT(ThetaState2[2],31);
        MatchRandG = Matchpoint[0] & 0x1;
    
        Matchpoint[1] = TAKE_BIT(TempState[0],33) ^ TAKE_BIT(TempState[10],31) ^ TAKE_BIT(ThetaState1[0],33) ^ TAKE_BIT(ThetaState2[0],33);
        MatchRandG = (MatchRandG << 1) ^ (Matchpoint[1] & 0x1);
    
        Matchpoint[2] = TAKE_BIT(TempState[6],36) ^ TAKE_BIT(TempState[16],36) ^ TAKE_BIT(ThetaState1[1],36) ^ TAKE_BIT(ThetaState2[1],36);
        MatchRandG = (MatchRandG << 1) ^ (Matchpoint[2] & 0x1);
    
        Matchpoint[3] = TAKE_BIT(TempState[12],37) ^TAKE_BIT(TempState[22],37) ^ TAKE_BIT(ThetaState1[2],37) ^ TAKE_BIT(ThetaState2[2],37);
        MatchRandG = (MatchRandG << 1) ^ (Matchpoint[3] & 0x1);
        
        Matchpoint[4] = TAKE_BIT(TempState[6],42) ^ TAKE_BIT(TempState[16],42) ^ TAKE_BIT(ThetaState1[1],42) ^ TAKE_BIT(ThetaState2[1],42);
        MatchRandG = (MatchRandG << 1) ^ (Matchpoint[4] & 0x1);
    
        Matchpoint[5] = TAKE_BIT(TempState[6],47) ^ TAKE_BIT(TempState[16],47) ^ TAKE_BIT(ThetaState1[1],47) ^ TAKE_BIT(ThetaState2[1],47);      
        MatchRandG = (MatchRandG << 1) ^ (Matchpoint[5] & 0x1);
    
        Matchpoint[6] = TAKE_BIT(TempState[6],55) ^ TAKE_BIT(TempState[16],55) ^ TAKE_BIT(ThetaState1[1],55) ^ TAKE_BIT(ThetaState2[1],55);     
        MatchRandG = (MatchRandG << 1) ^ (Matchpoint[6] & 0x1);
    
        Matchpoint[7] = TAKE_BIT(TempState[9],60) ^ TAKE_BIT(TempState[24],60) ^ TAKE_BIT(ThetaState1[1],60) ^ TAKE_BIT(ThetaState2[1],60);  
        MatchRandG = (MatchRandG << 1) ^ (Matchpoint[7] & 0x1) ;

        
        //cout << MatchRandG << endl;
        //Store the value for 8 red bits
        if(TableU.find(MatchRandG) != TableU.end())
        {
            vector<UINT64> ttmp =TableU[MatchRandG];
            ttmp.push_back(j);
            TableU[MatchRandG] = ttmp;
        } else {
            vector<UINT64> ttmp;
            ttmp.push_back(j);
            TableU[MatchRandG] = ttmp;
        }
        
    }
    
    UINT64 MatchNum=0;
    //f'''m
    UINT64 Conste[8] = {0};
    //Fix Red, traverse the blue bits
    for(UINT64 j=0;j<(UINT64(1)<<6);j++){
        for(UINT64 k=0;k<25;k++){
            TempState[k]=InitialState[k];
        }
        
        for(UINT64 k=0;k<6;k++)
        {
            UINT64 temp1=(j>>(5-k))&1;
            if(temp1){
                TempState[indexb[k][0]] |= (UINT64(1)<<indexb[k][1]);
                TempState[indexb[k][0]+5] |= (UINT64(1)<<indexb[k][1]);
            }
            else{
                TempState[indexb[k][0]] &= ROL64(~UINT64(1),indexb[k][1]);
                TempState[indexb[k][0]+5] &= ROL64(~UINT64(1),indexb[k][1]);
            }
                    
        }
        KeccakPermutationOnWords(TempState, 2);
        
        UINT64 Matchpoint[8] = {0};
        UINT64 MatchB = 0;
        //Compute f''m
        Matchpoint[0] = TAKE_BIT(TempState[12],31) ^ TAKE_BIT(TempState[22],31) ^ TAKE_BIT(ThetaState1[2],31) ^ TAKE_BIT(ThetaState2[2],31);
        //Compute f'''m
        if (j==0)
           Conste[0]=Matchpoint[0];
        //f''m^f'''m
        MatchB = (Matchpoint[0]^Conste[0]) & 0x1;
    
        Matchpoint[1] = TAKE_BIT(TempState[0],33) ^ TAKE_BIT(TempState[10],31) ^ TAKE_BIT(ThetaState1[0],33) ^ TAKE_BIT(ThetaState2[0],33);
        if (j==0)
            Conste[1]=Matchpoint[1];
        MatchB = (MatchB << 1) ^ ((Matchpoint[1]^Conste[1]) & 0x1);
    
        Matchpoint[2] = TAKE_BIT(TempState[6],36) ^ TAKE_BIT(TempState[16],36) ^ TAKE_BIT(ThetaState1[1],36) ^ TAKE_BIT(ThetaState2[1],36);
        if (j==0)
            Conste[2]=Matchpoint[2];
        MatchB = (MatchB << 1) ^ ((Matchpoint[2]^Conste[2]) & 0x1);
    
        Matchpoint[3] = TAKE_BIT(TempState[12],37) ^TAKE_BIT(TempState[22],37) ^ TAKE_BIT(ThetaState1[2],37) ^ TAKE_BIT(ThetaState2[2],37);
        if (j==0)
            Conste[3]=Matchpoint[3];
        MatchB = (MatchB << 1) ^ ((Matchpoint[3]^Conste[3]) & 0x1);
        
        Matchpoint[4] = TAKE_BIT(TempState[6],42) ^ TAKE_BIT(TempState[16],42) ^ TAKE_BIT(ThetaState1[1],42) ^ TAKE_BIT(ThetaState2[1],42);
        if (j==0)
            Conste[4]=Matchpoint[4];
        MatchB = (MatchB << 1) ^ ((Matchpoint[4]^Conste[4]) & 0x1);
    
        Matchpoint[5] = TAKE_BIT(TempState[6],47) ^ TAKE_BIT(TempState[16],47) ^ TAKE_BIT(ThetaState1[1],47) ^ TAKE_BIT(ThetaState2[1],47); 
        if (j==0)
            Conste[5]=Matchpoint[5];
        MatchB = (MatchB << 1) ^ ((Matchpoint[5]^Conste[5]) & 0x1);
        
        Matchpoint[6] = TAKE_BIT(TempState[6],55) ^ TAKE_BIT(TempState[16],55) ^ TAKE_BIT(ThetaState1[1],55) ^ TAKE_BIT(ThetaState2[1],55);
        if (j==0)
            Conste[6]=Matchpoint[6];
        MatchB = (MatchB << 1) ^ ((Matchpoint[6]^Conste[6]) & 0x1);
        
        Matchpoint[7] = TAKE_BIT(TempState[9],60) ^ TAKE_BIT(TempState[24],60) ^ TAKE_BIT(ThetaState1[1],60) ^ TAKE_BIT(ThetaState2[1],60);  
        if (j==0)
            Conste[7]=Matchpoint[7];
        MatchB = (MatchB << 1) ^ ((Matchpoint[7]^Conste[7]) & 0x1);
        
        
        if(TableU.find(MatchB) != TableU.end()) {
            //cout << "Find the Match!\n";
            vector<UINT64> ttmp = TableU[MatchB];
            MatchNum += ttmp.size();
            for (UINT64 k=0; k<ttmp.size(); k++) {
                if(((ttmp[k]<<6)^j)==rightkey)
                    cout << "The original values for red and blue cells remain!" << endl;
            }
        }
    }
    
    cout << "In total, " << MatchNum << " matches are found!" << endl;
    cout << "In total, 2^" << log(double(MatchNum)) / log(2.0) << " matches are found!" << endl;
    return 0;

}


