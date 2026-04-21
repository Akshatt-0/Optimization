#ifndef RR_H
#define RR_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <cmath>
#include "utils.h"

// CAC type: RR
void addClosestAssignmentConstraint(IloEnv env, IloModel& model, const std::vector<Point>& nodes,
                                    IloArray<IloNumVarArray>& x, IloNumVarArray& y) {
    int n = nodes.size();

    for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                IloExpr lhs(env);
                double dij = distance(nodes[i], nodes[j]);
                for (int a = 0; a < n; ++a) {
                    if (distance(nodes[i], nodes[a]) <= dij) {
                        lhs += x[i][a];
                    }
                }
                model.add(lhs >= y[j]);
                lhs.end();
            }
    }
}
#endif // RR