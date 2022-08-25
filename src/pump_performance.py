import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import root_scalar, minimize_scalar


class CentrifugalPumpSpecification:
    def __init__(
        self,
        fluid: str,
        rho: float,
        impeller_inner_diameter: float,
        impeller_outer_diameter: float,
        inlet_passage_width: float,
        outlet_passage_width: float,
        blade_clearance: float,
        blade_thickness: float,
        number_of_impeller_blades: int,
        roughness: float,
        beta_1_deg: float,
        beta_2_deg: float,
        volute_throat_area: float,
        volute_height: float,
        diffuser_area_ratio: float,
        alpha_1_deg: float = 90,
    ):
        self.fluid = fluid
        self.rho = rho
        self.d1 = impeller_inner_diameter
        self.d2 = impeller_outer_diameter
        self.b1 = inlet_passage_width
        self.b2 = outlet_passage_width
        self.s = blade_clearance
        self.t = blade_thickness
        self.z = number_of_impeller_blades
        self.epsilon = roughness
        self.beta1 = np.radians(beta_1_deg)
        self.beta2 = np.radians(beta_2_deg)
        self.alpha1 = np.radians(alpha_1_deg)
        self.A4 = volute_throat_area
        self.A5 = self.A4 * diffuser_area_ratio
        self.b3 = volute_height


class CentrifugalPump:
    def __init__(self, data: CentrifugalPumpSpecification):
        self.fluid = data.fluid
        self.g = 9.81
        self.rho = data.rho
        self.d1 = data.d1
        self.d2 = data.d2
        self.b1 = data.b1
        self.b2 = data.b2
        self.b3 = data.b3
        self.s = data.s
        self.t = data.t
        self.z = data.z
        self.epsilon = data.epsilon
        self.beta1 = data.beta1
        self.beta2 = data.beta2
        self.alpha1 = data.alpha1
        self.r1 = self.d1 / 2
        self.r2 = self.d2 / 2
        self.A1 = np.pi * self.d1 * self.b1
        self.A2 = np.pi * self.d2 * self.b2
        self.A1r = self.A1 - self.t * self.z * self.b1
        self.A2r = self.A2 - self.t * self.z * self.b2
        self.A4 = data.A4
        self.A5 = data.A5

    def _friction_colebrook(self, f, dh, Re):
        res = 1 / np.sqrt(f) + 2 * np.log10(
            (self.epsilon / dh) / 3.7 + 2.51 / (Re * np.sqrt(f))
        )
        return res

    def _friction_factor(self, a1, a2, t, v_avg):
        P = 101325
        dh = 2 * (a2 * self.b2 + a1 * self.b1) / (a1 + a2 + self.b1 + self.b1)
        v = PropsSI("V", "T", t, "P", P, self.fluid)
        Re = self.rho * v_avg * dh / v
        f0 = 0.25 * (np.log((self.epsilon / dh) / 3.7 + 5.74 / Re ** 0.9)) ** (-2)
        f = root_scalar(
            method="bisect",
            f=self._friction_colebrook,
            x0=f0,
            args=(dh, Re),
            bracket=[0.0001, 300],
        ).root
        return f

    def _pump_friction_loss(self, v1: float, v2: float, t: float):
        L = (self.d2 - self.d1) / (2 * np.sin(self.beta2))
        a1 = np.pi * self.d1 / self.z * np.sin(self.beta1)
        a2 = np.pi * self.d2 / self.z * np.sin(self.beta2)
        dh = 2 * (a2 * self.b2 + a1 * self.b1) / (a1 + a2 + self.b1 + self.b1)
        v_avg = (v1 + v2) / 2
        f = self._friction_factor(a1=a1, a2=a2, t=t, v_avg=v_avg)
        Hl_f = f * L * v_avg ** 2 / (dh * 2 * self.g)

        return Hl_f

    def _pump_incidence_loss(self, beta1, beta1_fl, w1):
        lmda = np.sin(beta1_fl) / np.sin(2 * beta1_fl - beta1) - np.sqrt(
            (np.sin(beta1_fl) ** 2 - np.sin(beta1) * np.sin(2 * beta1_fl - beta1))
            / (np.sin(2 * beta1_fl - beta1)) ** 2
        )
        Hl_inc = (
            w1 ** 2
            / (2 * self.g)
            * (1 / lmda) ** 2
            * (1 - lmda * np.sin(beta1_fl) / np.sin(beta1)) ** 2
        )
        return Hl_inc

    def _pump_leak_loss_open(self, Hth):
        a2 = np.pi * self.d2 / self.z * np.sin(self.beta2)
        nq = 20
        Hl_leak_ratio = (2.5 * self.s / self.d2) / (
            np.sqrt(self.b2 / self.d2 * (1 - self.d1 / self.d2))
            * (self.t / a2)
            * nq ** 0.1
            * np.sin(self.beta2) ** 1.2
            * np.sin(self.beta1) ** 0.4
        )
        return Hl_leak_ratio * Hth

    def _pump_leak_loss_closed(self, Q, Hth, omega, w1: float, w2: float, t: float):
        Hl_stat = Hth - (omega / 2) ** 2 * (self.d2 ** 2 - self.d1 ** 2) / (8 * self.g)
        a1 = np.pi * self.d1 / self.z * np.sin(self.beta1)
        a2 = np.pi * self.d2 / self.z * np.sin(self.beta2)
        w_avg = (w1 + w2) / 2
        f = self._friction_factor(a1=a1, a2=a2, t=t, v_avg=w_avg)
        L_gap = (self.d2 - self.d1) / 2
        v_gap = np.sqrt(2 * self.g * Hl_stat / (f * L_gap / self.s + 0.5 + 0.5))
        A_gap = np.pi * self.d2 * self.s
        Q_leak = v_gap * A_gap
        eta_vol = Q / (Q + Q_leak)
        return eta_vol

    def _pump_diffusion_loss(self, w1):
        Hl_d = 0.25 * w1 ** 2 / (2 * self.g)
        return Hl_d

    def _volute_friction_loss(self, Q, t):
        P = 101325
        dh_volute = 4 * self.A4 / (4 * self.b3)
        d3 = self.d2 + 2 * self.b3
        c4m = Q / self.A4
        v = PropsSI("V", "T", t, "P", P, self.fluid)
        Re = self.rho * c4m * dh_volute / v
        f0 = 0.25 * (np.log((self.epsilon / dh_volute) / 3.7 + 5.74 / Re ** 0.9)) ** (
            -2
        )
        f = root_scalar(
            method="bisect",
            f=self._friction_colebrook,
            x0=f0,
            args=(dh_volute, Re),
            bracket=[0.0001, 100],
        ).root

        Hl_f_v = f * (4 * d3 * self.b3 * c4m ** 2) / (dh_volute ** 2 * 2 * self.g)
        return Hl_f_v

    def _volute_radial_loss(self, c2r):
        Hl_radial_v = c2r ** 2 / (2 * self.g)
        return Hl_radial_v

    def _volute_head_loss(self, Q, c_theta_2):
        c4m = Q / self.A4
        d4 = self.d2 + self.b3 + np.sqrt(self.A4 / np.pi) / 2
        c_theta_4 = c_theta_2 * self.d2 / d4
        Hl_h_v = 0.8 * (c_theta_4 ** 2 - c4m ** 2) / (2 * self.g)
        return Hl_h_v

    def _volute_diffuser_loss(self, Q):
        cp = 2 * ((self.A4 / self.A5) - (self.A4 / self.A5) ** 2)
        delta_diff = 1 - cp
        c4m = Q / self.A4
        Hl_diff_v = delta_diff * c4m ** 2 / (2 * self.g)
        return Hl_diff_v

    def _power_disk_friction(self, t, u2):
        P = 101325
        v = PropsSI("V", "T", t, "P", P, self.fluid)
        k = 7.3e-4 * (2 * v * 1e6 / (u2 * self.d2)) ** (1 / 8)
        P_disk = k * self.rho * u2 ** 3 * self.d2 * (self.d2 + 5 * self.s)
        return P_disk

    def pump_head(self, Q, n, T):
        omega = n / 60 * 2 * np.pi

        c1m = Q / self.A1r
        c2m = Q / self.A2r

        u1 = omega * self.r1
        u2 = omega * self.r2

        c_theta_1 = c1m / np.tan(self.alpha1)
        c1 = np.sqrt(c1m ** 2 + c_theta_1 ** 2)
        w_theta_1 = u1 - c_theta_1
        w1 = np.sqrt(c1m ** 2 + w_theta_1 ** 2)
        beta1_fl = np.arctan(c1m / w_theta_1)
        incidence = self.beta1 - beta1_fl

        sigma = 1 - np.sin(self.beta2) / (self.z ** 0.7)

        w_theta_2 = c2m / np.tan(self.beta2)
        w2 = np.sqrt(c2m ** 2 + w_theta_2 ** 2)
        c_theta_2_slip = sigma * (u2 - w_theta_2)
        w_theta_2_slip = u2 - c_theta_2_slip
        beta2_fl = np.arctan(c2m / w_theta_2_slip)

        Hth = 1 / self.g * (u2 * c_theta_2_slip - u1 * c_theta_1)
        Hl_f = self._pump_friction_loss(v1=w1, v2=w2, t=T)
        Hl_inc = self._pump_incidence_loss(beta1=self.beta1, beta1_fl=beta1_fl, w1=w1)
        Hl_leak = self._pump_leak_loss_open(Hth=Hth)
        Hl_d = self._pump_diffusion_loss(w1=w1)
        Hl_f_v = self._volute_friction_loss(Q=Q, t=T)
        Hl_r_v = self._volute_radial_loss(c2r=c2m)
        Hl_h_v = self._volute_head_loss(Q=Q, c_theta_2=c_theta_2_slip)
        Hl_diff_v = self._volute_diffuser_loss(Q=Q)
        H_out = (
            Hth - Hl_f - Hl_inc - Hl_leak - Hl_d - Hl_f_v - Hl_r_v - Hl_h_v - Hl_diff_v
        )

        P_disk = self._power_disk_friction(t=T, u2=u2)
        P_hyd = Hth * self.g * self.rho * Q
        P_out = H_out * self.g * self.rho * Q
        P_in = P_hyd + P_disk
        eta = P_out / P_in

        return [
            H_out,
            np.rad2deg(self.beta1 - beta1_fl),
            Hth,
            Hl_f,
            Hl_inc,
            Hl_d,
            Hl_f_v,
            Hl_r_v,
            Hl_h_v,
            Hl_diff_v,
            P_in,
            P_out,
            eta,
        ]
