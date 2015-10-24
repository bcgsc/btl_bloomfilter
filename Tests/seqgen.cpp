#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

const int alphNum = 4;
const char iTb[alphNum] = {'A', 'C', 'G', 'T'};

void makeRead(const size_t rNum, const size_t rLen) {
    srand(time(NULL));
    size_t bTc[alphNum] = {0, 0, 0, 0};
    string myrSeq;
    myrSeq.resize(rLen);
    ofstream rFile("reads.fa");

    for (size_t j=0; j < rNum; j++) {
        for (int i=0; i< rLen; i++) {
            int ranInd = rand() % alphNum;
            myrSeq[i] = iTb[ranInd];
            ++bTc[ranInd];
        }
        rFile << ">" << j << "\n" << myrSeq << "\n";
    }
    
    for (int i=0; i< alphNum; i++)
        cerr << bTc[i] << "\n";
    rFile.close();
}


void makeGenome(const size_t seqLen) {
    srand(time(NULL));
    size_t bTc[alphNum] = {0, 0, 0, 0};
    string mygSeq;
    mygSeq.resize(seqLen);
    for (int i=0; i< seqLen; i++) {
        int ranInd = rand() % alphNum;
        mygSeq[i] = iTb[ranInd];
        ++bTc[ranInd];
    }
    ofstream gFile("genome.fa");
    gFile << ">1\n" << mygSeq << "\n";
    for (int i=0; i< alphNum; i++)
        cerr << bTc[i] << "\n";
    gFile.close();
}

int main(int argc, const char* argv[]) {
    //makeGenome(3000000000+32-1);
    makeRead(1000000000, 32);
    return 0;
}