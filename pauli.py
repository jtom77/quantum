import sympy as sp

theta = sp.symbols('theta', real=True)

X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.Matrix([[1, 0], [0, -1]])
H = 1/sp.sqrt(2)*sp.Matrix([[1, 1], [1, -1]])
S = sp.Matrix([[1, 0], [0, sp.I]])
T = sp.Matrix([[1, 0], [0, sp.exp(sp.I*sp.pi / 4)]])
RX = sp.Matrix([[sp.cos(theta/2), -sp.I * sp.sin(theta/2)], 
                [-sp.I * sp.sin(theta/2), sp.cos(theta/2)]])
RY = sp.Matrix([[sp.cos(theta/2), -sp.sin(theta/2)], 
                [sp.sin(theta/2), sp.cos(theta/2)]])
RZ = sp.Matrix([[sp.exp(-sp.I*theta/2), 0], 
                [0, sp.exp(sp.I*theta/2)]])

def Rx(phi):
    return RX.subs(theta, phi)

def Ry(phi):
    return RY.subs(theta, phi)

def Rz(phi):
    return RZ.subs(theta, phi)

ket0 = sp.Matrix([1, 0])
ket1 = sp.Matrix([0, -1])
ketp = 1/sp.sqrt(2)* sp.Matrix([1, 1])
ketm = 1/sp.sqrt(2)* sp.Matrix([1, -1])


def normalize_state(state):
    rho = -sp.arg(state[0]) + sp.arg(state[1])
    return 1/state.norm() * sp.Matrix([sp.Abs(state[0]), sp.exp(sp.I * rho)*sp.Abs(state[1])])

def bloch_coords(state, numeric=False):
    psi = normalize_state(state)
    x = 2 * sp.acos(psi[0])
    y = sp.arg(psi[1])
    if numeric:
        x = x.evalf()
        y = y.evalf()
    return (x, y)
