#include <ilcplex/ilocplex.h>
#include <vector>
#include <cmath>


struct Point {
    double x, y;
    Point(double _x, double _y) : x(_x), y(_y) {}
};

double euclideanDistance(const Point& a, const Point& b) {
    return std::sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

int main() {
    IloEnv env;
    try {
        IloModel model(env);
        int n = 6;        // total nodes
        int p = 2;        // total hubs
        double alpha = 0.75;  // inter-hub discount factor

        // Sample coordinates
        std::vector<Point> nodes = {
            {0, 0}, {2, 0}, {4, 0},
            {1, 1}, {3, 1}, {5, 1}
        };

        // Variables
        IloArray<IloNumVarArray> x(env, n);  // assignment vars
        for (int i = 0; i < n; ++i)
            x[i] = IloNumVarArray(env, n, 0, 1, ILOINT);  // x[i][k]

        IloNumVarArray y(env, n, 0, 1, ILOINT);  // y[k]

        // Objective function
        IloExpr totalCost(env);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    for (int l = 0; l < n; ++l) {
                        double d1 = euclideanDistance(nodes[i], nodes[k]);
                        double d2 = euclideanDistance(nodes[k], nodes[l]) * alpha;
                        double d3 = euclideanDistance(nodes[l], nodes[j]);

                        totalCost += (d1 + d2 + d3) * x[i][k] * x[j][l];
                    }
                }
            }
        }
        model.add(IloMinimize(env, totalCost));
        totalCost.end();

        // Constraint 1: assign each node to 1 hub
        for (int i = 0; i < n; ++i) {
            IloExpr assign(env);
            for (int k = 0; k < n; ++k)
                assign += x[i][k];
            model.add(assign == 1);
            assign.end();
        }

        // Constraint 2: node assigned to hub only if hub is open
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k)
                model.add(x[i][k] <= y[k]);

        // Constraint 3: open exactly p hubs
        IloExpr hubCount(env);
        for (int k = 0; k < n; ++k)
            hubCount += y[k];
        model.add(hubCount == p);
        hubCount.end();

        // Solve
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());  // mute output
        if (cplex.solve()) {
            std::cout << "Total cost: " << cplex.getObjValue() << "\n";
            for (int k = 0; k < n; ++k)
                if (cplex.getValue(y[k]) > 0.5)
                    std::cout << "Hub at node " << k << "\n";

            for (int i = 0; i < n; ++i)
                for (int k = 0; k < n; ++k)
                    if (cplex.getValue(x[i][k]) > 0.5)
                        std::cout << "Node " << i << " assigned to hub " << k << "\n";
        } else {
            std::cout << "No solution found.\n";
        }

    } catch (IloException& e) {
        std::cerr << "Concert exception: " << e << "\n";
    } catch (...) {
        std::cerr << "Unknown error\n";
    }
    env.end();
    return 0;
}
