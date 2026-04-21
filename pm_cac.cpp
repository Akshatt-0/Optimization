#include <cmath>
#include <vector>
#include <iostream>

#include "utils.h"
#include "cc.h"

int main() {
    IloEnv env;
    try {
        // Define problem data
        std::vector<Point> nodes = {
            {0, 0}, {2, 1}, {1, 3}, {5, 0}, {6, 2}, {5, 4}
        };

        int n = nodes.size();
        int p = 3;  // number of facilities to open

        // Model
        IloModel model(env);

        // Decision variables
        IloArray<IloNumVarArray> x(env, n);  // x[i][j] = 1 if client i assigned to facility j
        for (int i = 0; i < n; ++i) {
            x[i] = IloNumVarArray(env, n);
            for (int j = 0; j < n; ++j) {
                x[i][j] = IloNumVar(env, 0, 1, ILOINT);
            }
        }

        IloNumVarArray y(env, n, 0, 1, ILOINT);  // y[j] = 1 if facility is open at j

        // Objective: minimize total distance
        IloExpr totalCost(env);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                totalCost += distance(nodes[i], nodes[j]) * x[i][j];
            }
        }
        model.add(IloMinimize(env, totalCost));
        totalCost.end();

        // Constraint 1: each client assigned to exactly one facility
        for (int i = 0; i < n; ++i) {
            IloExpr assign(env);
            for (int j = 0; j < n; ++j) {
                assign += x[i][j];
            }
            model.add(assign == 1);
            assign.end();
        }

        // Constraint 2: assignment only to open facilities
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                model.add(x[i][j] <= y[j]);
            }
        }

        // Constraint 3: open exactly p facilities
        IloExpr facilityCount(env);
        for (int j = 0; j < n; ++j) {
            facilityCount += y[j];
        }
        model.add(facilityCount == p);
        facilityCount.end();

        // Constraint 4: CAC(CC)
        addClosestAssignmentConstraint(env, model, nodes, xij);

        // Solve
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());  // silence solver

        if (cplex.solve()) {
            std::cout << "Total cost: " << cplex.getObjValue() << "\n";
            for (int j = 0; j < n; ++j) {
                if (cplex.getValue(y[j]) > 0.5) {
                    std::cout << "Facility opened at node " << j << "\n";
                }
            }
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (cplex.getValue(x[i][j]) > 0.5) {
                        std::cout << "Client " << i << " assigned to facility " << j << "\n";
                    }
                }
            }
        } else {
            std::cout << "No solution found.\n";
        }
    } catch (IloException& e) {
        std::cerr << "Concert exception: " << e << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception caught." << std::endl;
    }

    env.end();
    return 0;
}
