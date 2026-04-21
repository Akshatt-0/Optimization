#include <ilcplex/ilocplex.h>
#include <iostream>
#include <vector>
#include <cmath>

ILOSTLBEGIN

const int n = 7;              // number of nodes
const int p = 4;              // number of hubs to open
const int r = 2;              // number of interdictions allowed

// Coordinates for nodes
double coords[n][2] = {
    {0, 0}, {2, 1}, {5, 2}, {6, 5}, {8, 3}, {1, 7}, {3, 6}
};

// Demand matrix
double demand[n][n] = {
  { 0, 5, 8, 6, 4, 7, 9 },
  { 6, 0, 3, 7, 5, 8, 4 },
  { 4, 7, 0, 6, 9, 5, 3 },
  { 9, 2, 6, 0, 4, 7, 6 },
  { 7, 8, 5, 3, 0, 4, 2 },
  { 3, 5, 7, 6, 5, 0, 8 },
  { 6, 4, 5, 7, 8, 2, 0 }
};

// Euclidean distance function
double distance(int i, int j) {
    return std::sqrt(std::pow(coords[i][0] - coords[j][0], 2) + std::pow(coords[i][1] - coords[j][1], 2));
}

int main() {
    IloEnv env;
    try {
        IloModel model(env);

        IloBoolVarArray z(env, n);  // 1 if hub at node k is interdicted
        for (int i = 0; i < n; ++i) {
            z[i] = IloBoolVar(env, 0, 1);
        }

        IloBoolVarArray y(env, n);  // 1 if node i is selected as hub
        for (int i = 0; i < n; ++i) {
            y[i] = IloBoolVar(env, 0, 1);
        }

        // Assignment variables: x[i][j][k][m] is fraction from i to j via hubs k, m
        IloArray<IloArray<IloArray<IloNumVarArray>>> x(env, n);
        for (int i = 0; i < n; ++i) {
            x[i] = IloArray<IloArray<IloNumVarArray>>(env, n);
            for (int j = 0; j < n; ++j) {
                x[i][j] = IloArray<IloNumVarArray>(env, n);
                for (int k = 0; k < n; ++k) {
                    x[i][j][k] = IloNumVarArray(env, n);
                    for (int m = 0; m < n; ++m) {
                        x[i][j][k][m] = IloNumVar(env, 0.0, 1.0);
                    }
                }
            }
        }

        // Objective: maximize minimum total cost under interdiction
        IloNumVar totalCost(env);
        IloExpr costExpr(env);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    for (int m = 0; m < n; ++m) {
                        double c = distance(i, k) + 0.5 * distance(k, m) + distance(m, j);
                        costExpr += demand[i][j] * c * x[i][j][k][m];
                    }
                }
            }
        }
        model.add(totalCost == costExpr);
        model.add(IloMaximize(env, totalCost));

        // (1) Exactly p hubs are opened
        IloExpr hubOpen(env);
        for (int k = 0; k < n; ++k) hubOpen += y[k];
        model.add(hubOpen == p);
        hubOpen.end();

        // (2) Exactly r hubs are interdicted
        IloExpr interdictCount(env);
        for (int k = 0; k < n; ++k) interdictCount += z[k];
        model.add(interdictCount == r);
        interdictCount.end();

        // (3) A hub can only be interdicted if it is opened
        for (int k = 0; k < n; ++k) {
            model.add(z[k] <= y[k]);
        }

        // CAC2: x[i][j][q][s] <= 2 - z_k - z_m for hubs k,m
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                for (int k = 0; k < n; ++k) {
                    for (int m = 0; m < n; ++m) {
                        if (k == m) continue;
                        double d_ijkm = distance(i,k) + 0.5 * distance(k,m) + distance(m,j);

                        IloExpr lhs(env);
                        for (int q = 0; q < n; ++q) {
                            for (int s = 0; s < n; ++s) {
                                if (q == s) continue;
                                double d_ijqs = distance(i,q) + 0.5 * distance(q,s) + distance(s,j);
                                if (d_ijqs > d_ijkm) {
                                    lhs += x[i][j][q][s];
                                }
                            }
                        }

                        model.add(lhs <= 2 - z[k] - z[m]);
                        lhs.end();
                    }
                }
            }
        }

        // Each demand pair must be fully assigned
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                IloExpr sum(env);
                for (int k = 0; k < n; ++k) {
                    for (int m = 0; m < n; ++m) {
                        sum += x[i][j][k][m];
                    }
                }
                model.add(sum == 1);
            }
        }

        // ∑ x_{ijkm} + ∑ x_{ijmk} ≤ z_k  ∀ i,j ∈ N; k ∈ N
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    IloExpr sum(env);
                    for (int m = 0; m < n; ++m) {
                        sum += x[i][j][k][m]; // k as first hub
                        if (m != k) {
                            sum += x[i][j][m][k]; // k as second hub (excluding duplicate)
                        }
                    }
                    model.add(sum <= z[k]);
                }
            }
        }

        // Solve
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        if (cplex.solve()) {
            std::cout << "Maximized total cost under interdiction: " << cplex.getValue(totalCost) << "\n";
            std::cout << "Hubs: ";
            for (int i = 0; i < n; ++i) if (cplex.getValue(y[i]) > 0.5) std::cout << i << " ";
            std::cout << "\nInterdicted hubs: ";
            for (int i = 0; i < n; ++i) if (cplex.getValue(z[i]) > 0.5) std::cout << i << " ";
            std::cout << "\n";
        } else {
            std::cout << "No solution found.\n";
        }
    } catch (IloException& e) {
        std::cerr << "Concert exception: " << e << "\n";
    } catch (...) {
        std::cerr << "Unknown exception caught.\n";
    }
    env.end();
    return 0;
}
