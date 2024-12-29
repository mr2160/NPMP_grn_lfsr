import numpy as np 

def solve_model(T,state):
    D1, C, notD1, notC, D1andC, notD1andC, Q1, Q1_, D2, C2, notD2, notC2, D2andC2, notD2andC2, Q2, Q2_ = state
    dD1 = -D1*0
    dC = -C*0
    dnotD1 = -notD1*0.7+40*(1)/(1+((D1/10)**4))
    dnotC = -notC*0.7+40*(1)/(1+((C/10)**4))
    dD1andC = -D1andC*0.7+40*(((C/10)**4))/(1+((notD1/10)**4)+((C/10)**4)+((notD1/10)**4)*((C/10)**4))
    dnotD1andC = -notD1andC*0.7+40*(((C/10)**4))/(1+((D1/10)**4)+((C/10)**4)+((D1/10)**4)*((C/10)**4))
    dQ1 = -Q1*0.7+40*(1)/(1+((notD1andC/10)**4)+((Q1_/10)**4)+((notD1andC/10)**4)*((Q1_/10)**4))
    dQ1_ = -Q1_*0.7+40*(1)/(1+((D1andC/10)**4)+((Q1/10)**4)+((D1andC/10)**4)*((Q1/10)**4))
    dD2 = -D2*0.7+40*(((Q1/10)**4))/(1+((Q1/10)**4))
    dC2 = -C2*0.7+40*(((notC/10)**4))/(1+((notC/10)**4))
    dnotD2 = -notD2*0.7+40*(1)/(1+((D2/10)**4))
    dnotC2 = -notC2*0.7+40*(1)/(1+((C2/10)**4))
    dD2andC2 = -D2andC2*0.7+40*(((C2/10)**4))/(1+((notD2/10)**4)+((C2/10)**4)+((notD2/10)**4)*((C2/10)**4))
    dnotD2andC2 = -notD2andC2*0.7+40*(((C2/10)**4))/(1+((D2/10)**4)+((C2/10)**4)+((D2/10)**4)*((C2/10)**4))
    dQ2 = -Q2*0.7+40*(1)/(1+((notD2andC2/10)**4)+((Q2_/10)**4)+((notD2andC2/10)**4)*((Q2_/10)**4))
    dQ2_ = -Q2_*0.7+40*(1)/(1+((D2andC2/10)**4)+((Q2/10)**4)+((D2andC2/10)**4)*((Q2/10)**4))
    return np.array([dD1, dC, dnotD1, dnotC, dD1andC, dnotD1andC, dQ1, dQ1_, dD2, dC2, dnotD2, dnotC2, dD2andC2, dnotD2andC2, dQ2, dQ2_])

def solve_model_steady(state):
    return solve_model(0, state)
