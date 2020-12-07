//Comentariu
#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>

const float crossProb = 0.2;
const int populationSize = 100;
const int cityNumber = 8;
std::vector<std::vector<int>> globalPopulation(populationSize);

void generateCandidate(std::vector<int>& candidateParam) {
	//Reserve the necessary capacity for the vector
	candidateParam.reserve(cityNumber);
	for (int i = 0; i < cityNumber; ++i) {
		//random number between 0 and cityNumber - i
		int rand = std::rand() % (cityNumber - i);
		candidateParam.push_back(rand);
	}
}

void Initialization() {
	for (auto it = globalPopulation.begin(); it != globalPopulation.end(); ++it) {
		//Generate candidate one by one
		generateCandidate(*it);
	}
}

//Cross function
void Cross(int index1, int index2) {
	//Declare the children
	std::vector<int> child1;
	std::vector<int> child2;
	child1.reserve(cityNumber);
	child2.reserve(cityNumber);

	//Swap la punctDeSwap
	int punctDeSwap = std::rand() % (cityNumber-1); //random index from 0 to cityNumber-2

	std::vector<int> parent1 = globalPopulation[index1];
	std::vector<int> parent2 = globalPopulation[index2];
	
	for (int ii = 0; ii <= punctDeSwap; ++ii) {
		child1.push_back(parent1[ii]);
		child2.push_back(parent2[ii]);
	}
	for (int jj = punctDeSwap+1; jj < cityNumber; ++jj) {
		child1.push_back(parent2[jj]);
		child2.push_back(parent1[jj]);
	}
	//Replace the parents with their children in the population
	globalPopulation.push_back(child1);
	globalPopulation.push_back(child2);
}


void crossOver() {
	std::vector<std::pair<float, int>> crossList;
	crossList.reserve(populationSize);
	for (int i = 0; i < globalPopulation.size(); ++i) {
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

void Selection() {
	//Global Eval vector for evaluation of every chromosome`s fitness
	std::vector<float> cumulatedSelectionProbability(populationSize);
	std::vector<float> selectionProbability(populationSize);

	//Total fitness
	float T = 0;
	for (int ii = 0; ii < populationSize; ++ii) {
		T += evalPop[ii];
	}

	//Individual selection probability
	for (int ii = 0; ii < populationSize; ++ii) {
		selectionProbability[ii] = evalPop[ii] / T; // Current fitness impartit la total fitness
	}

	//Accumulated selection probability
	cumulatedSelectionProbability[0] = 0;
	for (int ii = 0; ii < populationSize - 1; ++ii) {
		cumulatedSelectionProbability[ii + 1] = cumulatedSelectionProbability[ii] + selectionProbability[ii];
	}
	cumulatedSelectionProbability[populationSize - 1] = 1;

	//Selection procedure
	double ran;
	for (int ii = 0; ii < populationSize; ++ii) {
		//Generate a random number between 0 and 1
		ran = (double(std::rand()) / double(RAND_MAX));

		for (int jj = 0; jj < populationSize - 1; ++jj) {
			
			//The jj chromosome corresponds to the random sector in roata norocului
			if (ran >= cumulatedSelectionProbability[jj] && ran <= cumulatedSelectionProbability[jj + 1]) {
					//Push the jj chromosome in our population
					globalPopulation[ii] = globalPopulation[jj];
					break; //Go to next ii
			}
		}
	}
}

int main() {
	std::srand(std::time(nullptr));
	Initialization();
	crossOver();
	Selection();
	return 0;
}