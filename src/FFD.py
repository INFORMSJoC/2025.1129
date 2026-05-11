from gurobipy import GRB
import gurobipy as gp

def model_FFD(m: int, L: int):
    model = gp.Model("FFD")
    model.Params.LogToConsole = 0

    s = model.addVars(m, lb=0.001, ub=1.0, name="s")
    p = model.addVars(L, lb=0.0, ub=1.0, name="p")
    x_FFD = model.addVars(m, L, vtype=GRB.BINARY, name="x_FFD")
    x_OPT = model.addVars(m, L, vtype=GRB.BINARY, name="x_OPT")

    # Normalization constraints
    model.addConstr(s[0] == 1.0)
    model.addConstrs(p[j] >= p[j + 1] for j in range(L - 1))
    model.addConstrs(s[i] >= s[i + 1] for i in range(m - 1))

    # Constraints
    model.addConstrs(gp.quicksum(x_OPT[i, j] for i in range(m)) == 1 for j in range(L))

    opt_load = model.addVars(m, lb=0.0, ub=1.0*L)
    model.addConstrs(opt_load[i]*s[i] == gp.quicksum(p[j] * x_OPT[i, j] for j in range(L)) for i in range(m))
    C_opt = model.addVar(lb=0.9, ub=1.1)
    model.addGenConstrMax(C_opt, [opt_load[i] for i in range(m)])
    model.addConstr(C_opt == 1.0)

    y_FFD = model.addVars(L, vtype=GRB.BINARY, name="y_FFD")
    model.addConstr(gp.quicksum(y_FFD[j] for j in range(L)) >= 1)
    model.addConstrs(y_FFD[j] + gp.quicksum(x_FFD[i, j] for i in range(m)) == 1 for j in range(L))

    gamma = model.addVar(lb=1.0, ub=7.0, name="C_FFD")
    model.setObjective(gamma, GRB.MAXIMIZE)

    # FFD
    for i in range(m):
        for j in range(L):
            model.addConstr(
                p[j] + gp.quicksum(x_FFD[i, k]*p[k] for k in range(j)) >= gamma * s[i] - 100*gp.quicksum(x_FFD[h, j] for h in range(i, m))
            )

    return model

if __name__ == "__main__":
    # Reproduce the first column of Table 3 of the paper.
    for m, L in [
        (2,4),
        (3,7),
        (4,9),
    ]:
        model = model_FFD(m=m, L=L)
        model.optimize()
        C_FFD = model.getVarByName("C_FFD").X
        lemma = 1 + (m-1)/(L+1)
        assert C_FFD > lemma, f"Expected C_FFD > {lemma:.4f}, but got C_FFD={C_FFD:.4f}"
        print(f"m={m}, L={L}, C_FFD={C_FFD:.4f}")
