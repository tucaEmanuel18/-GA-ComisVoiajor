// ComisVoiajor.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <list>
#include <chrono>
#include <fstream>
#include <algorithm>
using namespace std;
ifstream fin("date.in.txt");
#define INFINIT 999999999

const float mutateProb = 0.2;
const float crossProb = 0.2;
const int popSize = 100;
const int nrVertex = 3;
vector<vector<int>> pop(popSize * 2 + 2);
vector<float> popEval(popSize * 2 + 2);
list<int> unitPerm({0,1,2,3,4,5});
float distanceMatrix[nrVertex + 1][nrVertex + 1];

// Generarea populatiei initiale
vector<int> cromInit()
{
    int nrVertexLocal = nrVertex;
    vector<int> crom;
    for (; nrVertexLocal > 1; --nrVertexLocal)
    {
        crom.push_back(rand() % nrVertexLocal);
    }
    return crom;
}
void popInit()
{ 
    for (int ii = 0; ii < popSize; ++ii)
    {
        pop.push_back(cromInit());
    }
}


// Evaluarea populatiei
vector<int> cromReal(const vector<int>& crom)
{
    vector<int> real;
    list<int> unusedVertex = unitPerm;
    auto it = unusedVertex.begin();
    for (auto ii : crom)
    {
        advance(it, ii);
        
        real.push_back(*it);
        unusedVertex.erase(it);

        it = unusedVertex.begin();
    }
    real.push_back(*it);
    return real;
}
float evalCrom(const vector<int>& crom)
{
    vector<int> realCrom = cromReal(crom);
    float sum = 0;
    int currentVertex = realCrom[0];
    for (int nextVertex = 1; nextVertex < realCrom.size(); ++nextVertex)
    {
        sum += distanceMatrix[currentVertex][nextVertex];
        currentVertex = nextVertex;
    }
    sum += distanceMatrix[currentVertex][realCrom[0]];

    return sum;
}
float evalFitnessCrom(const vector<int>& crom)
{
    float result = evalCrom(crom);
    if (result == 0)
    {
        result += 0.001;
    }
    return 1 / result;
}
void eval()
{
    popEval.clear();
    for (auto ii : pop)
    {
        popEval.push_back(evalFitnessCrom(ii));
    }
}

// Cei mai buni cromozomi
vector<int> bestTwoCrom()
{
    //Salvam cei mai buni 2 cromozomi pentru elitism
    float maxF0 = float(-1 * INFINIT);
    float maxF1 = float(-1 * INFINIT);
    vector<int> pos = { 0, 0 };
    for (unsigned int i = 0; i < popEval.size(); ++i)
    {
        if (popEval[i] > maxF0)
        {
            maxF1 = maxF0;
            pos[1] = pos[0];
            maxF0 = popEval[i];
            pos[0] = i;
        }
        else
            if (popEval[i] > maxF1)
            {
                maxF1 = popEval[i];
                pos[1] = i;
            }
    }
    return pos;
}

int bestCrom()
{
    float maxF = float(-1 * INFINIT);
    int pos = 0;
    for (unsigned int i = 0; i < popEval.size(); ++i)
    {
        if (popEval[i] > maxF)
        {
            maxF = popEval[i];
            pos = i;
        }
    }
    return pos;
}
// Mutatia
void mutate()
{
    //printf("\nMUTATE:");
    float p;
    int nrCrom = 0;
    for (auto crom : pop)
    {
        int nrLocus = 0;
        for (auto locus : crom)
        {

            //char old = locus;
            p = ((float)rand() / (RAND_MAX));
            if (p < mutateProb)
            {
                pop[nrCrom][nrLocus] = (pop[nrCrom][nrLocus] + 1) % (nrVertex - nrLocus);
            }


            nrLocus++;
        }
        nrCrom++;
    }
}

//Crossover
void Cross(int index1, int index2) {
    //Declare the children
    std::vector<int> child1;
    std::vector<int> child2;
    child1.reserve(nrVertex);
    child2.reserve(nrVertex);

    //Swap la punctDeSwap
    int punctDeSwap = std::rand() % (nrVertex - 1); //random index from 0 to cityNumber-2

    std::vector<int> parent1 = pop[index1];
    std::vector<int> parent2 = pop[index2];

    for (int ii = 0; ii <= punctDeSwap; ++ii) {
        child1.push_back(parent1[ii]);
        child2.push_back(parent2[ii]);
    }
    for (int jj = punctDeSwap + 1; jj < nrVertex; ++jj) {
        child1.push_back(parent2[jj]);
        child2.push_back(parent1[jj]);
    }
    //Replace the parents with their children in the population
    pop.push_back(child1);
    pop.push_back(child2);
}


void crossOver() {
    std::vector<std::pair<float, int>> crossList;
    crossList.reserve(popSize + 2);
    for (int i = 0; i < pop.size(); ++i) {
        //For each chromosome, generate random probability between 0 and 1
        double rand = (double(std::rand()) / double(RAND_MAX));
        //Push crossList the pair that consists of probability and chromosome index
        crossList.push_back(std::pair<float, int>(rand, i));
    }
    //Sort the list after probability
    std::sort(crossList.begin(), crossList.end());

    //Crossover 2 by 2 till the line
    int j = 0;
    while (crossList[j].first < crossProb) {
        if (crossList[j + 1].first >= crossProb) break;
        else {
            Cross(crossList[j].second, crossList[j + 1].second);
        }
        j = j + 2; //Next pair
    }

}

void select() {
    vector<vector<int>> newPop(popSize * 2 + 2);
    //Global Eval vector for evaluation of every chromosome`s fitness
    vector<float> cumulatedSelectionProbability(popSize * 2 + 2);
    vector<float> selectionProbability(popSize * 2 + 2);
   
    //Total fitness
    float T = 0;
    for (int ii = 0; ii < popSize; ++ii) {
        T += popEval[ii];
    }

    //Individual selection probability
    for (int ii = 0; ii < popSize; ++ii) {
        selectionProbability[ii] = popEval[ii] / T; // Current fitness impartit la total fitness
    }

    //Accumulated selection probability
    cumulatedSelectionProbability[0] = 0;
    for (int ii = 0; ii < popSize - 1; ++ii) {
        cumulatedSelectionProbability[ii + 1] = cumulatedSelectionProbability[ii] + selectionProbability[ii];
    }
    cumulatedSelectionProbability[popSize - 1] = 1;
      
    //Selection procedure
    double ran;
    for (int ii = 0; ii < popSize; ++ii) {
        //Generate a random number between 0 and 1
        ran = (double(std::rand()) / double(RAND_MAX));

        for (int jj = 0; jj < popSize - 1; ++jj) {

            //The jj chromosome corresponds to the random sector in roata norocului
            if (ran >= cumulatedSelectionProbability[jj] && ran <= cumulatedSelectionProbability[jj + 1]) {
                //Push the jj chromosome in our population
                newPop.push_back(pop[jj]);
                break; //Go to next ii
            }
        }
    }
    pop.clear();
    pop.reserve(popSize * 2 + 2);
    pop = newPop;
}

// Printari
void printRealCrom(vector<int> crom)
{
    vector<int> real = cromReal(crom);
    for (auto ii : real)
    {
        printf("%d", ii);
    }
}
void printCrom(vector<int> crom)
{
    for (auto ii : crom)
    {
        printf("%d", ii);
    }
}
void printPop()
{
    printf("\nPopulation:");
    int nrCrom = 0;
    for (auto ii : pop)
    {
        printf("\nCrom_%d: ", nrCrom++);
        printCrom(ii);
    }
}
//Algoritmul Genetic
void algGenetic(int nrGeneratii)
{
	int t = 0;
	// Generam populatia initiala aleator
	popInit();
    
    // evaluam populatia initiala (avem nevoie pentru elitism)
	eval();

	// Elitism (alegem cei mai buni 2 cromozomi)
	vector<int> positionBest = bestTwoCrom();
	vector<vector<int>> bestTwo;
	bestTwo.push_back(pop[positionBest[0]]);
	bestTwo.push_back(pop[positionBest[1]]);

    while (t < nrGeneratii)
    {
        // aplicam mutatia peste populatie
        mutate();

        // adaugam in populatie pe cei doi indivizi buni (Elitism)
        pop.push_back(bestTwo[0]);
        pop.push_back(bestTwo[1]);

        // CrossOVER
        crossOver();
        // evaluam populatia obtinuta
         eval();

        // Elitism (alegem cei mai buni 2 cromozomi)
        positionBest = bestTwoCrom();
        bestTwo[0] = pop[positionBest[0]];
        bestTwo[1] = pop[positionBest[1]];

        // SELECT
        select();
        t++;

        eval();
        int posBest = bestCrom();

        printf("\nGeneratia %d\nBestCrom:", t);
        printRealCrom(pop[posBest]);
        printf("\nResult = %f", evalCrom(pop[posBest]));

    }

    // Gasirea celui mai bun individ din ultima generatie.
	eval();
    int posBest = bestCrom();

    printRealCrom(pop[posBest]);
    printf("\nResult = %f", evalCrom(pop[posBest]));
  

	//printf("\n##################   FINAL    #######\n");
	//printCromValue(pop[pos]);
}



int main()
{
    srand(unsigned int(time(NULL))); 

    ifstream infile;
    infile.open("date.in");
 
    for (int ii = 0; ii < nrVertex; ++ii)
    {
        for (int jj = 0; jj < nrVertex; ++jj)
        {
            int value;
            infile >> value;
           distanceMatrix[ii][jj] = value;
        }
    }
    algGenetic(15);

    infile.close();

}
