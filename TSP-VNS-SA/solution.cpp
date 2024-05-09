#include <bits/stdc++.h>
#include "utils.hpp"
#include <fstream>
using namespace std;
int maxIterations = 100;
int seed = 1714104515;
int best_known=0;
int vns_ans=0;
int sa_ans=0;
double calculate_cost(vector<vector<double>> *dist,vector<int>&soln) {
    double sum = 0;
    for(int i = 0; i < soln.size() - 1; i++) {
        sum += dist[0][soln[i] - 1][soln[i + 1] - 1];
    }
    sum += dist[0][soln[soln.size() - 1] - 1][soln[0] - 1];
    return sum;
}
vector<int> generate_random_solution(int size) {
    vector<int> solution(size);
    for(int i = 0; i < size; i++) {
        solution[i] = i + 1;
    }
    srand(seed);
    printf("\ntime %d\n",time(0));
    random_shuffle(solution.begin(), solution.end());
    return solution;
}

void reverseTour(vector<int>& tour, int i, int j) {
    while (i < j) {
        swap(tour[i], tour[j]);
        i++;
        j--;
    }
}

void shake(vector<int>& tour, int k) {
    int choice = rand() % 2;
    if (choice == 0) {
        // 2-opt neighborhood change
        int i = rand() % tour.size();
        int j;
        do {
            j = rand() % tour.size();
        } while (abs(i - j) < 2);
        reverseTour(tour, i, j);
    } else {
        // 3-opt neighborhood change
        int i, j, l;
        do {
            i = rand() % tour.size();
            j = rand() % tour.size();
            l = rand() % tour.size();
        } while (abs(i - j) < 2 || abs(j - l) < 2 || abs(i - l) < 2);
        reverseTour(tour, i, j);
        reverseTour(tour, j, l);
    }
}

void twoOptLocalSearch(vector<vector<double>> *dist, vector<int>& tour) {
    bool improved = true;

    while (improved) {
        improved = false;
        for (int i = 0; i < tour.size() - 2; i++) {
            for (int j = i + 2; j < tour.size(); j++) {
                // Check if reversing the segment between indices i and j improves the tour
                double dist_before = (*dist)[tour[i] - 1][tour[i + 1] - 1] + (*dist)[tour[j] - 1][tour[(j + 1) % tour.size()] - 1];
                double dist_after = (*dist)[tour[i] - 1][tour[j] - 1] + (*dist)[tour[i + 1] - 1][tour[(j + 1) % tour.size()] - 1];

                if (dist_after < dist_before) {
                    // Reverse the segment
                    reverseTour(tour, i + 1, j);
                    improved = true;
                }
            }
        }
    }
}

vector<int> variable_neighborhood_search(const char * problem,int k_max){
    best_known = 0;
    vector<vector<double>> *dist =get_matrix(problem,&best_known);
    int n_cities = (*dist).size();
    cout<<"Best_known cost: "<<best_known<<endl;
    vector<int>initial_solution = generate_random_solution(n_cities);
    double current_cost = calculate_cost(dist,initial_solution);
    cout<<"Initial solution: ";
    for(int i = 0; i < initial_solution.size(); i++) {
        cout << initial_solution[i] << ' ';
    }
    cout << '\n';
    cout<<"Initial cost: "<<current_cost<<endl;
    int k = 1;
    int iteration = 0;
    while (iteration < maxIterations) {
        while (k <= k_max) {
            vector<int> currentTour = initial_solution;

            shake(currentTour, k);

            twoOptLocalSearch(dist, currentTour);
            double currentLength = calculate_cost(dist, currentTour);
            double bestLength = calculate_cost(dist, initial_solution);

            if (currentLength < bestLength) {
                initial_solution = currentTour;
                k = 1;
            } else {
                k++;
            }
        }
        iteration++;
    }
    vns_ans=calculate_cost(dist,initial_solution);
    cout<<"Final cost: "<<vns_ans<<endl;
    return initial_solution;
}



vector<int> alternate_solution(vector<int> solution) {
    vector<int> new_solution = solution;
    int i = rand() % solution.size();
    int j = rand() % solution.size();
    swap(new_solution[i], new_solution[j]);
    return new_solution;
}
float acceptance_probability(float prev, float curr, float temperature) {
    if (curr < prev) {
        return 1.0;
    }
    return exp((prev - curr) / temperature);
}

vector<int> simulated_annealing(const char * problem) {
    best_known = 0;
    vector<vector<double>> *dist =get_matrix(problem,&best_known);
    int n_cities = (*dist).size();
    vector<int> solution = generate_random_solution(n_cities);
    double current_cost = calculate_cost(dist,solution);
    cout<<"Initial solution: ";
    for(int i = 0; i < solution.size(); i++) {
        cout << solution[i] << ' ';
    }
    cout << '\n';
    cout<<"Initial cost: "<<current_cost<<endl;
    float temperature = 100000;
    float cooling_rate = 0.99995;
    float stopping_temperature = 0.000001;
    while(temperature > stopping_temperature) {
        vector<int> new_solution = alternate_solution(solution);
        float new_cost = calculate_cost(dist,new_solution);
        if(new_cost < current_cost || acceptance_probability(current_cost, new_cost, temperature) > ((double)rand() / RAND_MAX)) {
            solution = new_solution;
            current_cost = new_cost;
        }
        temperature *= cooling_rate;
    }
    sa_ans=calculate_cost(dist,solution);
    cout<<"Final cost: "<<sa_ans<<endl;
    return solution;
}
int main(int argc, const char ** argv) {
    std::ofstream vns_output("vns_output.txt");
    std::ofstream sa_output("sa_output.txt");
    sa_output << argv[1] << '\n';
    vns_output << argv[1] << '\n';
    int k_max=5000;
    if(argc == 3) {
        seed = atoi(argv[2]);
        std::cout << "seed" << seed << std::endl;
    }
    cout<<"--------------------------------Variable Neighborhood Search---=-----------------------------------\n";
    vector<int> solution = variable_neighborhood_search(argv[1],k_max);
    for (const auto& i : solution) {
        cout << i << ' ';
        vns_output << i << ' ';
    }
    cout << '\n';
    cout<<"-------------------------------------------Simulated Annealing-------------------------------------\n";
    solution = simulated_annealing(argv[1]);
    for (const auto& i : solution) {
        cout << i << ' ';
        sa_output << i << ' ';
    }
    cout << '\n';
    vns_output.close();
    sa_output.close();

    cout<<"----------------Stats----------------\n";
    cout<<"Best Known Solution: "<<best_known<<endl;
    cout<<"VNS Solution: "<<vns_ans<<endl;
    cout<<"SA Solution: "<<sa_ans<<endl;
    cout<<"Better to use "<<(vns_ans<sa_ans?"VNS":"SA") <<" for this "<<endl;
    return 0;
}