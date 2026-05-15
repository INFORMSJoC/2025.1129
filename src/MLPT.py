from gurobipy import GRB
import gurobipy as gp

def model_MLPT(m: int, L: int, h: int):
    model = gp.Model("MLPT")
    model.Params.LogToConsole = 0

    s = model.addVars(m, lb=0.001, ub=1.0, name="s")
    p = model.addVars(L, lb=0.0, ub=1.0, name="p")
    x_MLPT = model.addVars(m, L, vtype=GRB.BINARY, name="x_MLPT")
    x_OPT = model.addVars(m, L, vtype=GRB.BINARY, name="x_OPT")

    # Gurobi does not allow to multiply three variables together,
    # so we introduce temporary variables to store the product x[i, j] * p[j].
    y_MLPT = {}
    for i in range(m):
        for j in range(L):
            y_MLPT[i, j] = model.addVar(lb=0.0, ub=1.0)
            model.addConstr(y_MLPT[i, j] == p[j] * x_MLPT[i, j])

    # Normalization constraints
    model.addConstr(s[0] == 1.0)
    model.addConstrs(p[j] >= p[j + 1] for j in range(L - 1))
    model.addConstrs(s[i] >= s[i + 1] for i in range(m - 1))

    # Constraints
    model.addConstrs(gp.quicksum(x_MLPT[i, j] for i in range(m)) == 1 for j in range(L))
    model.addConstrs(gp.quicksum(x_OPT[i, j] for i in range(m)) == 1 for j in range(L))

    opt_load = model.addVars(m, lb=0.0, ub=1.0*L)
    model.addConstrs(opt_load[i]*s[i] == gp.quicksum(p[j] * x_OPT[i, j] for j in range(L)) for i in range(m))
    C_opt = model.addVar(lb=0.9, ub=1.1)
    model.addGenConstrMax(C_opt, [opt_load[i] for i in range(m)])
    model.addConstr(C_opt == 1.0)

    Mlpt_load = model.addVars(m, lb=0.0, ub=1.0*L)
    model.addConstrs(Mlpt_load[i]*s[i] == gp.quicksum(y_MLPT[i, j] for j in range(L)) for i in range(m))
    C_Mlpt = model.addVar(lb=1.0, ub=7.0, name="C_Mlpt")
    model.addGenConstrMax(C_Mlpt, [Mlpt_load[i] for i in range(m)])
    model.setObjective(C_Mlpt, GRB.MAXIMIZE)

    # As jobs (0,...,h-1) are sorted optimally and C^=1, we known that the completion times of the first h jobs in the HEU solution are at most 1
    model.addConstrs(gp.quicksum(x_MLPT[i,j]*p[j] for j in range(h)) <= s[i] for i in range(m))
    
    # Finally we the makespan is at most the sum of all the remaining jobs plus the partial completion times
    remaining_jobs_total = gp.quicksum(p[j] for j in range(h, L))
    model.addConstrs(C_Mlpt*s[i] <= gp.quicksum(x_MLPT[i,j]*p[j] for j in range(h))+remaining_jobs_total for i in range(m))

    # Add the LPT constraints
    B = 100.0
    for i in range(m):
        for itilde in range(m):
            if i == itilde:
                continue
            for j in range(h, L):
                term1 = s[itilde] * (p[j] + gp.quicksum(y_MLPT[i,k] for k in range(j)))
                term2 = s[i] * (p[j] + gp.quicksum(y_MLPT[itilde,k] for k in range(j)))
                model.addConstr(term1 <= term2 + B * (1 - x_MLPT[i,j]))


    return model

if __name__ == "__main__":
    # Reproduce Table 2 of the paper.
    for m, L, h in [
        (2,4,3),
        (2,5,4),
        (2,6,5),
        (2,7,6),
        (2,8,7),
        (3,5,4),
        (3,6,5),
        (3,7,6),
        (3,8,7),
        (4,7,5),
        (4,8,6),
        # (4,8,7)
    ]:
        model = model_MLPT(m=m, L=L, h=h)
        model.optimize()
        C_Mlpt = model.getVarByName("C_Mlpt").X
        lemma = 1 + (m-1)/(L+1)
        assert C_Mlpt > lemma, f"Expected C_Mlpt > {lemma:.4f}, but got C_Mlpt={C_Mlpt:.4f}"
        print(f"m={m}, L={L}, C_MLPT={C_Mlpt:.4f}")