import numpy as np
from src.pump_performance import CentrifugalPumpSpecification


def flow_coefficient(n_d, Q_d, H_d):
    nq = n_d * np.sqrt(Q_d) / H_d ** (3 / 4)
    return 0.003 * nq - 0.0013


def pump_design(
    fluid,
    n_d,
    Q_d,
    H_d,
    head_coefficient_psi,
    flow_coefficient_phi,
    beta2,
    t,
    z,
    dh1,
    i,
):
    nq = n_d * np.sqrt(Q_d) / H_d ** (3 / 4)
    g = 9.81
    omega_d = n_d / 60 * 2 * np.pi
    d2 = 2 / omega_d * np.sqrt(g * H_d / head_coefficient_psi)
    u2 = omega_d * d2 / 2
    b2 = Q_d / (flow_coefficient_phi * u2 * (np.pi * d2 - t * z))
    A2 = np.pi * d2 * b2 - t * z * b2
    c2m = Q_d / A2
    w_theta_2 = c2m / np.tan(beta2)

    ns = 51.55 * nq
    Kd = 1e-13 * ns ** 3 - 7e-9 * ns ** 2 + 1e-4 * ns + 0.2661
    d1 = Kd * d2
    b1 = (d1 - dh1) / 2
    A1 = np.pi * ((dh1 + d1) / 2) * b1 - t * z * b1
    c1m = Q_d / A1
    beta1 = np.arctan(c1m / w_theta_2) + i

    u1 = omega_d * d1 / 2
    alpha1 = np.deg2rad(90)
    c_theta_1 = c1m / np.tan(alpha1)
    w_theta_1 = u1 - c_theta_1
    w1 = np.sqrt(c1m ** 2 + w_theta_1 ** 2)
    sigma = 1 - np.sin(beta2) / (z ** 0.7)
    w_theta_2 = c2m / np.tan(beta2)
    w2 = np.sqrt(c2m ** 2 + w_theta_2 ** 2)
    c_theta_2_slip = sigma * (u2 - w_theta_2)

    diameter_size_check = d1 < d2
    blade_height_size_check = b2 < b1

    Kw = (w1 - w2) / w1
    Ks = (u2 - u1) / c_theta_2_slip

    pump_design_result = CentrifugalPumpSpecification(
        fluid=fluid,
        rho=1000,
        impeller_inner_diameter=0.05,
        impeller_outer_diameter=0.240,
        inlet_passage_width=0.0175,
        outlet_passage_width=0.0065,
        blade_clearance=0.05e-3,
        blade_thickness=0.005,
        number_of_impeller_blades=5,
        roughness=0.045e-3,
        beta_1_deg=36,
        beta_2_deg=23,
        alpha_1_deg=90,
        volute_throat_area=np.pi * 0.03 ** 2,
        volute_height=1e-3,
        diffuser_area_ratio=2,
    )

    print("Specific speed: {}".format(nq))
    print("Outer diameter: {}".format(d2))
    print("Inner diameter: {}".format(d1))
    print("Beta 1: {}".format(np.rad2deg(beta1)))
    print("diameter size check {}".format(diameter_size_check))
    print("blade height size check {}".format(blade_height_size_check))
    print("Kw factor below 0.25: {}".format(Kw))
    print("Ks factor below 0.9: {}".format(Ks))
    print()


n_d = 1000
Q_d = 40 / 1000
H_d = 10

pump_design(
    fluid="INCOMP::Water",
    n_d=n_d,
    Q_d=Q_d,
    H_d=H_d,
    head_coefficient_psi=0.5,
    flow_coefficient_phi=flow_coefficient(n_d=n_d, Q_d=Q_d, H_d=H_d),
    beta2=10,
    t=8e-3,
    z=6,
    dh1=0.06,
    i=0,
)
