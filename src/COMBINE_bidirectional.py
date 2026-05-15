from gurobipy import GRB
import gurobipy as gp

def model_bCOMBINE(m: int, L: int):
    model = gp.Model("bidirectionalCOMBINE")
    model.Params.LogToConsole = 0

    s = model.addVars(m, lb=0.001, ub=1.0, name="s")
    p = model.addVars(L, lb=0.0, ub=1.0, name="p")
    x_FFD = model.addVars(m, L, vtype=GRB.BINARY, name="x_FFD")
    x_rFFD = model.addVars(m, L, vtype=GRB.BINARY, name="x_rFFD")
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
    model.addConstrs(gp.quicksum(x_OPT[i, j] for i in range(m)) == 1 for j in range(L))
    model.addConstrs(gp.quicksum(x_LPT[i, j] for i in range(m)) == 1 for j in range(L))

    opt_load = model.addVars(m, lb=0.0, ub=1.0*L)
    model.addConstrs(opt_load[i]*s[i] == gp.quicksum(p[j] * x_OPT[i, j] for j in range(L)) for i in range(m))
    C_opt = model.addVar(lb=0.9, ub=1.1)
    model.addGenConstrMax(C_opt, [opt_load[i] for i in range(m)])
    model.addConstr(C_opt == 1.0)

    y_FFD = model.addVars(L, vtype=GRB.BINARY, name="y_FFD")
    model.addConstr(gp.quicksum(y_FFD[j] for j in range(L)) >= 1)
    model.addConstrs(y_FFD[j] + gp.quicksum(x_FFD[i, j] for i in range(m)) == 1 for j in range(L))

    y_rFFD = model.addVars(L, vtype=GRB.BINARY, name="y_rFFD")
    model.addConstr(gp.quicksum(y_rFFD[j] for j in range(L)) >= 1)
    model.addConstrs(y_rFFD[j] + gp.quicksum(x_rFFD[i, j] for i in range(m)) == 1 for j in range(L))

    gamma_FFD = model.addVar(lb=1.0, ub=7.0, name="C_FFD")
    gamma_rFFD = model.addVar(lb=1.0, ub=7.0, name="C_rFFD")

    # "Forward"-FFD
    for i in range(m):
        for j in range(L):
            model.addConstr(
                p[j] + gp.quicksum(x_FFD[i, k]*p[k] for k in range(j)) >= gamma_FFD * s[i] - 100*gp.quicksum(x_FFD[h, j] for h in range(i, m))
            )

    
    # reversed FFD
    for i in range(m):
        for j in range(L):
            model.addConstr(
                p[j] + gp.quicksum(x_rFFD[i, k]*p[k] for k in range(j)) >= gamma_rFFD * s[i] - 100*gp.quicksum(x_rFFD[h, j] for h in range(i+1))
            )


    # LPT
    lpt_load = model.addVars(m, lb=0.0, ub=1.0*L)
    model.addConstrs(lpt_load[i]*s[i] == gp.quicksum(y_LPT[i, j] for j in range(L)) for i in range(m))
    C_lpt = model.addVar(lb=1.0, ub=2.0, name="C_lpt")
    model.addGenConstrMax(C_lpt, [lpt_load[i] for i in range(m)])

    B = 100.0
    for i in range(m):
        for itilde in range(m):
            if i == itilde:
                continue
            for j in range(L):
                term1 = s[itilde] * (p[j] + gp.quicksum(y_LPT[i,k] for k in range(j)))
                term2 = s[i] * (p[j] + gp.quicksum(y_LPT[itilde,k] for k in range(j)))
                model.addConstr(term1 <= term2 + B * (1 - x_LPT[i,j]))


    C = model.addVar(lb=1.0, ub=7.0, name="C")
    model.addConstr(C <= C_lpt)
    model.addConstr(C <= gamma_FFD)
    model.addConstr(C <= gamma_rFFD)
    model.setObjective(C, GRB.MAXIMIZE)

    return model

if __name__ == "__main__":
    # Reproduce the fifth column of Table 3 of the paper.
    for m, L in [
        (2,7),
        (3,8),
    ]:
        model = model_bCOMBINE(m=m, L=L)
        model.optimize()
        C = model.getVarByName("C").X
        lemma = 1 + (m-1)/(L+1)
        assert C > lemma, f"Expected C > {lemma:.4f}, but got C={C:.4f}"
        print(f"m={m} and L={L}: {C:.4f}")
