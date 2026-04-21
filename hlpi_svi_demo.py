
# hlpi_svi_demo.py
# A tiny, brute-force demo of the Supervalid Inequality (SVI) cutting-plane idea for a toy HLPI instance.
from itertools import combinations

# Problem size and parameters
n = 5
nodes = list(range(n))
p = 3
r = 1

coords = {
    0: (0,0),
    1: (1,0),
    2: (0,1),
    3: (1,1),
    4: (2,1),
}
def euclid(a,b):
    ax, ay = coords[a]
    bx, by = coords[b]
    return ((ax-bx)**2 + (ay-by)**2)**0.5
d = [[euclid(i,j) for j in nodes] for i in nodes]

W = [[0, 3, 2, 0, 1],
     [3, 0, 2, 1, 1],
     [2, 2, 0, 2, 0],
     [0, 1, 2, 0, 4],
     [1, 1, 0, 4, 0]]

delta = 1.0
alpha = 0.6
gamma = 1.0
w1 = 1.0
w2 = 1.0

def route_cost_given_hubs(open_hubs, d, W, delta, alpha, gamma):
    hubs = list(open_hubs)
    total = 0.0
    for i in range(n):
        for j in range(n):
            if i == j or W[i][j] == 0:
                continue
            best = float('inf')
            for k in hubs:
                for m in hubs:
                    c = delta*d[i][k] + alpha*d[k][m] + gamma*d[m][j]
                    if c < best:
                        best = c
            total += W[i][j] * best
    return total

def solve_MP_forbidden_pairs(forbidden_pairs):
    best_cost = float('inf')
    best_set = None
    for hubset in combinations(nodes, p):
        violates = any(set(pair).issubset(hubset) for pair in forbidden_pairs)
        if violates:
            continue
        pre = route_cost_given_hubs(hubset, d, W, delta, alpha, gamma)
        if pre < best_cost:
            best_cost = pre
            best_set = tuple(sorted(hubset))
    return best_cost, best_set

def solve_SP(open_hubs):
    hubs = list(open_hubs)
    worst_cost = -1.0
    worst_interdict = None
    worst_Q = None
    for h in hubs:
        surviving = tuple(sorted(x for x in hubs if x != h))
        post = route_cost_given_hubs(surviving, d, W, delta, alpha, gamma)
        if post > worst_cost:
            worst_cost = post
            worst_interdict = h
            worst_Q = surviving
    return worst_cost, worst_interdict, worst_Q

def main():
    forbidden_pairs = []
    incumbent = None
    iteration = 0
    while True:
        pre_cost, hubset = solve_MP_forbidden_pairs(forbidden_pairs)
        if hubset is None:
            break
        post_cost, interdict, Q = solve_SP(hubset)
        obj = w1*pre_cost + w2*post_cost
        print(f"[Iter {iteration}] hubs={hubset}, pre={pre_cost:.4f}, interdict={interdict}, Q={Q}, post={post_cost:.4f}, obj={obj:.4f}")
        if (incumbent is None) or (obj < incumbent[0]):
            incumbent = (obj, hubset, interdict, Q, pre_cost, post_cost)
        if Q not in forbidden_pairs:
            forbidden_pairs.append(Q)
        iteration += 1
    if incumbent:
        obj, hubset, interdict, Q, pre_cost, post_cost = incumbent
        print("\\nBest: obj={:.4f}, hubs={}, interdict={}, Q={}, pre={:.4f}, post={:.4f}".format(
            obj, hubset, interdict, Q, pre_cost, post_cost))
    else:
        print("No feasible solution (unexpected).")

if __name__ == "__main__":
    main()
