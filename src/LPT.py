from gurobipy import GRB
import gurobipy as gp

def model_LPT(m: int, L: int):
    model = gp.Model("LPT")
    model.Params.LogToConsole = 0

    s = model.addVars(m, lb=0.001, ub=1.0, name="s")
    p = model.addVars(L, lb=0.0, ub=1.0, name="p")
    x_LPT = model.addVars(m, L, vtype=GRB.BINARY, name="x_LPT")
    x_OPT = model.addVars(m, L, vtype=GRB.BINARY, name="x_OPT")

    # Gurobi does not allow to multiply three variables together,
    # so we introduce temporary variables to store the product x[i, j] * p[j].
    y_LPT = {}
    for i in range(m):
        for j in range(L):
            y_LPT[i, j] = model.addVar(lb=0.0, ub=1.0)
            model.addConstr(y_LPT[i, j] == p[j] * x_LPT[i, j])

    # Normalization constraints
    model.addConstr(s[0] == 1.0)
    model.addConstrs(p[j] >= p[j + 1] for j in range(L - 1))
    model.addConstrs(s[i] >= s[i + 1] for i in range(m - 1))

    # Constraints
    model.addConstrs(gp.quicksum(x_LPT[i, j] for i in range(m)) == 1 for j in range(L))
    model.addConstrs(gp.quicksum(x_OPT[i, j] for i in range(m)) == 1 for j in range(L))

    opt_load = model.addVars(m, lb=0.0, ub=1.0*L)
    model.addConstrs(opt_load[i]*s[i] == gp.quicksum(p[j] * x_OPT[i, j] for j in range(L)) for i in range(m))
    C_opt = model.addVar(lb=0.9, ub=1.1)
    model.addGenConstrMax(C_opt, [opt_load[i] for i in range(m)])
    model.addConstr(C_opt == 1.0)

    lpt_load = model.addVars(m, lb=0.0, ub=1.0*L)
    model.addConstrs(lpt_load[i]*s[i] == gp.quicksum(y_LPT[i, j] for j in range(L)) for i in range(m))
    C_lpt = model.addVar(lb=1.0, ub=2.0, name="C_lpt")
    model.addGenConstrMax(C_lpt, [lpt_load[i] for i in range(m)])
    model.setObjective(C_lpt, GRB.MAXIMIZE)


    # Add the LPT constraints
    B = 100.0
    for i in range(m):
        for itilde in range(m):
            if i == itilde:
                continue
            for j in range(L):
                term1 = s[itilde] * (p[j] + gp.quicksum(y_LPT[i,k] for k in range(j)))
                term2 = s[i] * (p[j] + gp.quicksum(y_LPT[itilde,k] for k in range(j)))
                model.addConstr(term1 <= term2 + B * (1 - x_LPT[i,j]))

    return model

if __name__ == "__main__":
    # Reproduce Table 1 of the paper.
    for m, L in [
        (2,3),
        (3,5),
        (4,6),
        # (5,8),
        # (6,10),
        # (7,12),
    ]:
        model = model_LPT(m=m, L=L)
        model.optimize()
        C_lpt = model.getVarByName("C_lpt").X
        lemma = 1 + (m-1)/(L+1)
        assert C_lpt > lemma, f"Expected C_lpt > {lemma:.4f}, but got C_lpt={C_lpt:.4f}"
        print(f"m={m}, L={L}, C_LPT={C_lpt:.4f}")