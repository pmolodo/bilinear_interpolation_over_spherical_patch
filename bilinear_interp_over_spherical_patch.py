"""
Deriving the formula for integrating a bilinear interpolation over a spherical patch in spherical coordinates
"""

import math

import sympy

from sympy import cos, sin


def eq_equal(equation1, equation2):
    return sympy.simplify(equation1.expand()) == sympy.simplify(equation2.expand())


def linear_interp(variable, var_min, var_max, val_min, val_max):
    return (val_max - val_min) / (var_max - var_min) * (variable - var_min) + val_min


def bilinear_interp(
    var1,
    var1_min,
    var1_max,
    var2,
    var2_min,
    var2_max,
    val_minmin,
    val_minmax,
    val_maxmin,
    val_maxmax,
):
    interp_var1_var2min = linear_interp(
        var1, var1_min, var1_max, val_minmin, val_maxmin
    )
    interp_var1_var2max = linear_interp(
        var1, var1_min, var1_max, val_minmax, val_maxmax
    )
    return sympy.simplify(
        linear_interp(
            var2, var2_min, var2_max, interp_var1_var2min, interp_var1_var2max
        )
    )


# define our spherical coordinates such that:
#  theta goes from 0 to pi  ("lattitude", except 0 is a pole, not the equator)
#  phi goes from 0 to 2*pi  ("longitude")


def spherical_surface_integral(
    theta, phi, func_theta_phi, theta_min, theta_max, phi_min, phi_max
):
    inner_theta_integral = sympy.simplify(
        sympy.integrate(
            func_theta_phi * sympy.sin(theta), (theta, theta_min, theta_max)
        )
    )
    outer_phi_integral = sympy.simplify(
        sympy.integrate(inner_theta_integral, (phi, phi_min, phi_max))
    )
    return outer_phi_integral


def bilinear_interp_over_spherical_patch(
    theta,
    phi,
    theta_min,
    theta_max,
    phi_min,
    phi_max,
    val_thetamin_phimin,
    val_thetamax_phimin,
    val_thetamin_phimax,
    val_thetamax_phimax,
):
    bilinear_func = bilinear_interp(
        theta,
        theta_min,
        theta_max,
        phi,
        phi_min,
        phi_max,
        val_thetamin_phimin,
        val_thetamin_phimax,
        val_thetamax_phimin,
        val_thetamax_phimax,
    )
    return spherical_surface_integral(
        theta, phi, bilinear_func, theta_min, theta_max, phi_min, phi_max
    )


theta, phi = sympy.symbols("theta phi")

# theta min/max = A/B
#   phi min/max = S/T
A, B, S, T = sympy.symbols("A B S T")

theta_min = A
theta_max = B
phi_min = S
phi_max = T

# interpolated corner values
E, F, G, H = sympy.symbols("E F G H")

val_thetamin_phimin = E
val_thetamax_phimin = F
val_thetamin_phimax = G
val_thetamax_phimax = H


# My formulation of the bilinear function
bilinear_func_ABSTEFGH = (
    (E - F - G + H) * theta * phi
    + ((G - H) * S + (F - E) * T) * theta
    + ((F - H) * A + (G - E) * B) * phi
    + H * A * S
    - F * A * T
    - G * B * S
    + E * B * T
) / ((B - A) * (T - S))


# Validate that my forumulation of the bilinear function is correct
bilinear_func_calc = bilinear_interp(
    theta,
    theta_min,
    theta_max,
    phi,
    phi_min,
    phi_max,
    val_thetamin_phimin,
    val_thetamin_phimax,
    val_thetamax_phimin,
    val_thetamax_phimax,
)
assert eq_equal(bilinear_func_ABSTEFGH, bilinear_func_calc)

# Now define in terms of coefficients P, Q, R and constants K, L

P_eq = E - F - G + H
Q_eq = (G - H) * S + (F - E) * T
R_eq = (F - H) * A + (G - E) * B
K_eq = H * A * S - F * A * T - G * B * S + E * B * T
L_eq = (B - A) * (T - S)

P, Q, R, K, L = sympy.symbols("P Q R K L")

PQRKL_expansion = dict(P=P_eq, Q=Q_eq, R=R_eq, K=K_eq, L=L_eq)

bilinear_func_PQRKL = (P * theta * phi + Q * theta + R * phi + K) / L

# double check that bilinear_func_PQRKL is the same as bilinear_func1
assert bilinear_func_PQRKL.subs(PQRKL_expansion) == bilinear_func_ABSTEFGH

spherical_surf_formula_ABSTEFGH = bilinear_interp_over_spherical_patch(
    theta,
    phi,
    theta_min=A,
    theta_max=B,
    phi_min=S,
    phi_max=T,
    val_thetamin_phimin=E,
    val_thetamax_phimin=F,
    val_thetamin_phimax=G,
    val_thetamax_phimax=H,
)

spherical_surf_formula_ABSTPQRKL = spherical_surface_integral(
    theta,
    phi,
    bilinear_func_PQRKL,
    theta_min=A,
    theta_max=B,
    phi_min=S,
    phi_max=T,
)

# Another validation of spherical_surf_formula_ABSTPQRKL
assert eq_equal(
    spherical_surf_formula_ABSTEFGH,
    spherical_surf_formula_ABSTPQRKL.subs(PQRKL_expansion),
)


def gamma(x):
    return sin(x) - x * cos(x)


C, D = sympy.symbols("C D")

C_eq = gamma(B) - gamma(A)
D_eq = cos(A) - cos(B)

CD_expansion = dict(C=C_eq, D=D_eq)

spherical_surf_formula_CDSTPQRKL = (
    (P * C + R * D) / 2 * (T**2 - S**2) + (Q * C + K * D) * (T - S)
) / L

assert eq_equal(
    spherical_surf_formula_CDSTPQRKL.subs(CD_expansion),
    spherical_surf_formula_ABSTPQRKL,
)


def calc_bilinear_interp_over_spherical_patch(
    theta_min,
    theta_max,
    phi_min,
    phi_max,
    val_thetamin_phimin,
    val_thetamax_phimin,
    val_thetamin_phimax,
    val_thetamax_phimax,
    use_sympy=False,
):
    if use_sympy:
        sin = sympy.sin
        cos = sympy.cos
    else:
        sin = math.sin
        cos = math.cos

    A = theta_min
    B = theta_max
    S = phi_min
    T = phi_max
    E = val_thetamin_phimin
    F = val_thetamax_phimin
    G = val_thetamin_phimax
    H = val_thetamax_phimax

    sinA = sin(A)
    cosA = cos(A)

    sinB = sin(B)
    cosB = cos(B)

    def gamma(x, sinx, cosx):
        return sinx - x * cosx

    gammaA = sinA - A * cosA
    gammaB = sinB - B * cosB

    C = gammaB - gammaA
    D = cosA - cosB

    P = E - F - G + H
    Q = (G - H) * S + (F - E) * T
    R = (F - H) * A + (G - E) * B
    K = H * A * S - F * A * T - G * B * S + E * B * T
    L = (B - A) * (T - S)

    return ((P * C + R * D) / 2 * (T**2 - S**2) + (Q * C + K * D) * (T - S)) / L


spherical_surf_formula_CDSTPQRKL_func = calc_bilinear_interp_over_spherical_patch(
    theta_min=A,
    theta_max=B,
    phi_min=S,
    phi_max=T,
    val_thetamin_phimin=E,
    val_thetamax_phimin=F,
    val_thetamin_phimax=G,
    val_thetamax_phimax=H,
    use_sympy=True,
)

assert eq_equal(
    spherical_surf_formula_CDSTPQRKL_func,
    spherical_surf_formula_ABSTEFGH,
)


# Note - it seems that our choice of spherical coordinate is affecting our final
# result, which seems wrong... ?

# ie, I've done this assuming 0 <= theta <= pi
