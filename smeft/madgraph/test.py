import sympy as sp

def solve_equation(equation):
    # 'equation' must be a string with python-like math and equal 0. 
    # x = sp.symbols('x')
    # return sp.solve(eval(equation), x)
    return sp.solveset(equation, 'x', sp.S.Reals)

print solve_equation('(1+0.000365*x+0.000001*x^2) - 11.493111')