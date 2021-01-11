// ComisVoiajor.cpp : This file contains the 'main' function. Program execution begins and ends there.
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
#include <ratio>
using namespace std;

#define INFINIT 999999999
ofstream outfile("br17.out");
const float defaultMutateProb = 0.01;
float mutateProb = defaultMutateProb;
float oldBestEval;
int blockCounter = 0;
const float crossProb = 0.7;
const int popSize = 250;
const int nrVertex = 17;
vector<vector<int>> pop;
vector<float> popEval;
list<int> unitPerm;
void unitPermInit()
{
    for (int i = 0; i < nrVertex; i++)
    {
        unitPerm.push_back(i);
    }
}
int distanceMatrix[nrVertex + 1][nrVertex + 1];

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
    for (int i = 1; i < realCrom.size(); ++i)
    {
        int nextVertex = realCrom[i];
        sum += distanceMatrix[currentVertex][nextVertex];
        currentVertex = nextVertex;
    }
  
    sum += distanceMatrix[currentVertex][realCrom[0]];

    return sum;
}
float evalFitnessCrom(float cost, float min, float max)
{
    float result = (max - cost) / (max - min) + 0.01;
    return pow(result, 3);
}
void eval()
{
    popEval.clear();
    vector<float> popCost;
    float min, max;
    min = INFINIT;
    max = -1 * INFINIT;
    for(auto ii : pop)
    {
        float actualCost = evalCrom(ii);
        if(actualCost > max)
        {
            max = actualCost;
        }
        if(actualCost < min)
        {
            min = actualCost;
        }
        popCost.push_back(actualCost);
    }
    for (auto ii : popCost)
    {
        popEval.push_back(evalFitnessCrom(ii, min, max));
    }
}

// Cei mai buni cromozomi
vector<int> bestThreeCrom()
{
    //Salvam cei mai buni 2 cromozomi pentru elitism
    float maxF0 = float(-1 * INFINIT);
    float maxF1 = float(-1 * INFINIT);
    float maxF2 = float(-1 * INFINIT);
    vector<int> pos = { 0, 0, 0 };
    for (unsigned int i = 0; i < popEval.size(); ++i)
    {
        if (popEval[i] > maxF0)
        {
            maxF2 = maxF1;
            pos[2] = pos[1];
            maxF1 = maxF0;
            pos[1] = pos[0];
            maxF0 = popEval[i];
            pos[0] = i;
        }
        else
            if (popEval[i] > maxF1)
            {
                maxF2 = maxF1;
                pos[2] = pos[1];
                maxF1 = popEval[i];
                pos[1] = i;
            }
            else
                if(popEval[i] > maxF2)
                {
                    maxF2 = popEval[i];
                    pos[2] = i;
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
void mutate2()
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


void Permute(vector<int>& crom, int nrLocus, int newNrLocus)
{
    int aux = crom[nrLocus];
    crom[nrLocus] = crom[newNrLocus];
    crom[newNrLocus] = aux;
}
void mutate1()
{
    int nrCrom = 0;
    for (auto crom : pop)
    {
        vector<int> newCrom = cromReal(crom);
        int nrLocus = 0;
        for (auto locus : newCrom)
        {
            //char old = locus;
           int p = ((float)rand() / (RAND_MAX));
            if (p < mutateProb)
            {
                int j = rand() % nrVertex;
                Permute(newCrom, nrLocus, j);
            }

            nrLocus++;
        }
        vector<int> leftVertex;
        for(int i = 0; i < nrVertex; ++i)
        {
            leftVertex.push_back(i);
        }
        vector<int> transCrom;
        for(auto ii : newCrom)
        {
            if(transCrom.size() < nrVertex - 1)
            {
                auto findPos = find(leftVertex.begin(), leftVertex.end(), ii); 
                transCrom.push_back(distance(leftVertex.begin(), findPos));
                leftVertex.erase(findPos);
            }
        }
        pop[nrCrom] = transCrom;
        nrCrom++;
    }
}

//Crossover
/*void Cross2(int index1, int index2) {
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
    for (int jj = punctDeSwap + 1; jj < nrVertex - 1; ++jj) {
        child1.push_back(parent2[jj]);
        child2.push_back(parent1[jj]);
    }
    //Replace the parents with their children in the population
    pop.push_back(child1);
    pop.push_back(child2);
}*/

void Cross(int index1, int index2) {
    //Declare the children
    std::vector<int> child;
    child.reserve(nrVertex);

    // declare aux data structure
    bool isValid1, isValid2;
    int nextVertex;

    //create set of leftVertext
    vector<int> leftVertex;
    for(int i = 0; i < nrVertex; ++i)
    {
        leftVertex.push_back(i);
    }


    // obtain the real path
    std::vector<int> parent1 = cromReal(pop[index1]);
    std::vector<int> parent2 = cromReal(pop[index2]);

    // obtain the start vertex with rand
    int vertex = rand() % nrVertex;

    child.push_back(vertex);
    auto erasePos = find(leftVertex.begin(), leftVertex.end(), vertex);
    leftVertex.erase(erasePos);


    while(leftVertex.size() > 0)
    {
        // obtain nextVertex1
        auto pos1 = find(parent1.begin(), parent1.end(), vertex);
        pos1++;
        if(pos1 == parent1.end())
            pos1 = parent1.begin();
        int nextVertex1 = *pos1;

        //obtain nextVertex2
        auto pos2 = find(parent2.begin(), parent2.end(), vertex);
        pos2++;
        if(pos2 == parent2.end())
            pos2 = parent2.begin();
        int nextVertex2 = *pos2;

        // obtain the nextVertexes validity
        isValid1 = true;
        isValid2 = true;
        for(auto ii : child)
        {
            if(ii == nextVertex1)
            {
                isValid1 = false;
            }
            if(ii == nextVertex2)
            {
                isValid2 = false;
            }
        }

        // obtain nextVertex
        if(isValid1)
        {
            if(isValid2)
            {
                if(distanceMatrix[vertex][nextVertex1] <= distanceMatrix[vertex][nextVertex2])
                {
                    nextVertex = nextVertex1;
                }
                else
                {
                    nextVertex = nextVertex2;
                }
            }
            else
            {
                nextVertex = nextVertex1;
            }
        }
        else
        {
            if(isValid2)
            {
                nextVertex = nextVertex2;
            }
            else
            {
                int randLeftVertex = rand()  % leftVertex.size();
                nextVertex = leftVertex[randLeftVertex];
            }
        }
        // add next vertex to child Path
        child.push_back(nextVertex);
        auto erasePos = find(leftVertex.begin(), leftVertex.end(), nextVertex);
        leftVertex.erase(erasePos);


        vertex = nextVertex;
    }

    // translate child path to representation
    leftVertex.clear(); 
    for(int i = 0; i < nrVertex; ++i)
    {
        leftVertex.push_back(i);
    }
    vector<int> newChild;
    for(auto ii : child)
    {
        if(newChild.size() < nrVertex - 1)
        {
            auto findPos = find(leftVertex.begin(), leftVertex.end(), ii); 
            newChild.push_back(distance(leftVertex.begin(), findPos));
            leftVertex.erase(findPos);
        }
    }
    //Add children in the population
    pop.push_back(newChild);
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
        else 
        {
            Cross(crossList[j].second, crossList[j + 1].second);
        }
        j = j + 2; //Next pair
    }

}

void select() {
    vector<vector<int>> newPop;
    //Global Eval vector for evaluation of every chromosome`s fitness
    vector<float> cumulatedSelectionProbability(popSize * 2 + 2);
    vector<float> selectionProbability(popSize * 2 + 2);
   
    //Total fitness
    float T = 0;
    for (int ii = 0; ii < pop.size(); ++ii) {
        T += popEval[ii];
    }

    //Individual selection probability
    for (int ii = 0; ii < pop.size(); ++ii) {
        selectionProbability[ii] = popEval[ii] / T; // Current fitness impartit la total fitness
    }

    //Accumulated selection probability
    cumulatedSelectionProbability[0] = 0;
    for (int ii = 0; ii < pop.size() - 1; ++ii) {
        cumulatedSelectionProbability[ii + 1] = cumulatedSelectionProbability[ii] + selectionProbability[ii];
    }
    cumulatedSelectionProbability[pop.size() - 1] = 1;
      
    //Selection procedure
    double ran;
    for (int ii = 0; ii < popSize; ++ii) 
    {
        //Generate a random number between 0 and 1
        ran = (double(std::rand()) / double(RAND_MAX));

        for (int jj = 0; jj < pop.size() - 1; ++jj) 
        {

            //The jj chromosome corresponds to the random sector in roata norocului
            if (ran >= cumulatedSelectionProbability[jj] && ran <= cumulatedSelectionProbability[jj + 1]) 
            {
                //Push the jj chromosome in our population
               /* if(evalCrom(pop[ii]) >= evalCrom(pop[jj]))
                {
                    newPop.push_back(pop[jj]);
                    
                }
                else
                {
                    newPop.push_back(pop[ii]);
                }
                break; //Go to next ii  */
                newPop.push_back(pop[jj]);
                break;  
            }
        }
    }
    pop.clear();
    pop = newPop;
}

// Printari
void printRealCrom(vector<int> crom)
{
    vector<int> real = cromReal(crom);
    for (auto ii : real)
    {
        printf("%d ", ii);
    }
}

void printRealPop()
{
    for(auto ii : pop)
    {
        printRealCrom(ii);
        printf("\n");
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
void printEval()
{
    int nrCrom = 0;
    for (auto ii : popEval)
    {
        printf("\n EvalCrom_%d = %f", nrCrom++, ii);
    }
}
//Algoritmul Genetic
void algGenetic(int nrGeneratii)
{
	int t = 0;
	// Generam populatia initiala aleator
	popInit();
    //printf("\nInitial Population\n");
    //printPop();


   // evaluam populatia initiala (avem nevoie pentru elitism)
   // printf("\nEval:\n");
	eval();
  //  printEval();

	// Elitism (alegem cei mai buni 2 cromozomi)
	vector<int> positionBest = bestThreeCrom();
	vector<vector<int>> bestThree;
	bestThree.push_back(pop[positionBest[0]]);
	bestThree.push_back(pop[positionBest[1]]);
    bestThree.push_back(pop[positionBest[2]]);
    while (t < nrGeneratii)
    {
        // aplicam mutatia peste populatie
        int mutateVar = rand() % 2;
        if(mutateVar)
        {
            mutate1();
        }
        else
        {
            mutate2();
        }
        

        // adaugam in populatie pe cei doi indivizi buni (Elitism)
        pop.push_back(bestThree[0]);
        pop.push_back(bestThree[1]);
        pop.push_back(bestThree[2]);

        // CrossOVER
        crossOver();
        // evaluam populatia obtinuta
         eval();

        // Elitism (alegem cei mai buni 2 cromozomi)
        positionBest = bestThreeCrom();
        bestThree[0] = pop[positionBest[0]];
        bestThree[1] = pop[positionBest[1]];
        bestThree[2] = pop[positionBest[2]];

        // SELECT
        select();
        t++;

        eval();
        int posBest = bestCrom();
        float evalBest = evalCrom(pop[posBest]);
        if(t == 1)
        {
            oldBestEval = evalBest;
        }
        else
        {
            if(oldBestEval == evalBest)
            {
                blockCounter++;
            }
            else
            {
                blockCounter = 0;
                oldBestEval = evalBest;
            }

            if(blockCounter == 0)
            {
                mutateProb = defaultMutateProb;
            }
            else
            {
                if(blockCounter % 5 == 0 && mutateProb < 0.75)
                {
                    mutateProb *= 1.2;
                }   
            }
        }




        /*printf("\nGeneratia %d\n", t);
        printf("\nResult = %f", evalCrom(pop[posBest]));
        printf("\nMutateProb = %f", mutateProb);
*/
    }

    // Gasirea celui mai bun individ din ultima generatie.
	eval();
    int posBest = bestCrom();
    //printf("\n");
   // printRealCrom(pop[posBest]);
    float evalResult = evalCrom(pop[posBest]);
   //printf("\nResult = %f", evalResult);
   outfile << evalResult;
  

	//printf("\n##################   FINAL    #######\n");
	//printCromValue(pop[pos]);
}



int main(int argc, char* argv[])
{

    if(argc != 2)
    {
        printf("Please get the name of test file!\n");
        exit(0);
    }
    else
    {
        srand(time(NULL)); 
        unitPermInit();
        //printf("Name of testFile: %s\n", argv[1]);
        ifstream infile;
        infile.open(argv[1]);
        int nr;
        for(int i = 0; i < nrVertex; i++)
        {
            for(int j = 0; j < nrVertex; j++)
            {
                infile >> nr;
                //printf("%d\n", nr);
                distanceMatrix[i][j] = nr;
            }       
        }
        /*for(int i = 0; i < nrVertex; i++)
        {
            for(int j = 0; j < nrVertex; j++)
            {
                printf("%d ", distanceMatrix[i][j]);
            }
            printf("\n");
        }*/
       
       for(int i = 0; i < 2; i++)
       {
            auto start = std::chrono::high_resolution_clock::now(); 
            algGenetic(5000);
            auto stop = std::chrono::high_resolution_clock::now(); 
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
            outfile << '/' << duration.count() << endl;
       }


        infile.close();
    }

 

}