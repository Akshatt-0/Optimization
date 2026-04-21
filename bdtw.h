#ifndef BDTW_H
#define BDTW_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <cmath>
#include "utils.h"

// CAC type: BDTW
void addClosestAssignmentConstraint(IloEnv env, IloModel& model, const std::vector<Point>& nodes,
                                    IloArray<IloNumVarArray>& x, IloNumVarArray& y) {
    int n = nodes.size();

    double M = 100;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            IloExpr lhs(env);
            double dij = distance(nodes[i], nodes[j]);
            lhs = (M - dij) * y[j];
                for (int a = 0; a < n; ++a) {
                    lhs += distance(nodes[i], nodes[a]) * x[i][a];
                }
                model.add(lhs <= M);
                lhs.end();
            }
        }
}
#endif // BDTW