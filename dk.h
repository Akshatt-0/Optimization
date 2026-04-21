#ifndef DK_H
#define DK_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <cmath>
#include "utils.h"

// CAC type: DK
void addClosestAssignmentConstraint(IloEnv env, IloModel& model, const std::vector<Point>& nodes,
                                    IloArray<IloNumVarArray>& x, IloNumVarArray& y) {
    int n = nodes.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double dij = distance(nodes[i], nodes[j]);
            for (int a = 0; a < n; ++a) {
                if (distance(nodes[i], nodes[a]) < dij) {
                    IloExpr lhs(env);
                    lhs = x[i][j] + y[a];
                    model.add(lhs <= 1);
                    lhs.end();
                }
            }
        }
    }
}
#endif // DK