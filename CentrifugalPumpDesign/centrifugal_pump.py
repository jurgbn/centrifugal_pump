# AUTOGENERATED! DO NOT EDIT! File to edit: ..\00_core.ipynb.

# %% auto 0
__all__ = ['CentrifugalPumpSpecification', 'view_pump_design_parameters', 'data_for_20sim', 'PumpDesign', 'CentrifugalPump']

# %% ..\00_core.ipynb 4
from CoolProp.CoolProp import PropsSI
from scipy.optimize import root_scalar
from collections import namedtuple
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.interpolate import griddata
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt


# %% ..\00_core.ipynb 9
@dataclass
class CentrifugalPumpSpecification:
        fluid: str
        specific_speed: float
        impeller_inner_diameter: float
        impeller_outer_diameter: float
        inlet_blade_height: float
        outlet_blade_height: float
        impeller_inlet_area: float
        impeller_outlet_area: float
        blade_clearance: float
        blade_thickness: float
        number_of_impeller_blades: int
        roughness: float
        blade_inlet_angle_rad: float
        blade_outlet_angle_rad: float
        volute_throat_area: float
        volute_height: float
        diffuser_area_ratio: float
        inlet_flow_angle_rad: float
        open_impeller: bool = True
        diameter_size_check: bool = False
        blade_height_size_check: bool = False
        Kw: float = 0.0
        Ks: float = 0.0
        Kw_check: bool = False
        Ks_check: bool = False


# %% ..\00_core.ipynb 10
def view_pump_design_parameters(pump_design: CentrifugalPumpSpecification):
    print('Inlet blade angle: {}'.format(round(np.rad2deg(pump_design.blade_inlet_angle_rad)), 2))
    print('Outlet blade angle: {}'.format(round(np.rad2deg(pump_design.blade_outlet_angle_rad)), 2))
    print('Impeller inner diameter: {} mm'.format(round(pump_design.impeller_inner_diameter, 3) * 1000))
    print('Impeller outer diameter {} mm'.format(round(pump_design.impeller_outer_diameter,3) * 1000))
    print('Inlet blade height {} mm'.format(round(pump_design.inlet_blade_height, 3) * 1000))
    print('Outlet blade height {} mm'.format(round(pump_design.outlet_blade_height, 3) * 1000))
    print('Impeller inlet area {} m2'.format(round(pump_design.impeller_inlet_area, 4)))
    print('Impeller outlet area {} m2'.format(round(pump_design.impeller_outlet_area, 4)))
    print('Volute throat {} m2'.format(round(pump_design.volute_throat_area, 5)))


def data_for_20sim(pump_design: CentrifugalPumpSpecification):
        print('string fluid = \'{}\';'.format(pump_design.fluid))
        print('real impeller_inner_diameter = {};'.format(pump_design.impeller_inner_diameter))
        print('real impeller_outer_diameter = {};'.format(pump_design.impeller_outer_diameter))
        print('real impeller_inlet_passage_witdh = {};'.format(pump_design.inlet_blade_height))
        print('real impeller_outlet_passage_width = {};'.format(pump_design.outlet_blade_height))
        print('real blade_clearance = {};'.format(pump_design.blade_clearance))
        print('real blade_thickness = {};'.format(pump_design.blade_thickness))
        print('real number_of_impeller_blades = {};'.format(pump_design.number_of_impeller_blades))
        print('real impeller_material_roughness = {};'.format(pump_design.roughness))
        print('real beta_1_deg = {};'.format(np.rad2deg(pump_design.blade_inlet_angle_rad)))
        print('real beta_2_deg = {};'.format(np.rad2deg(pump_design.blade_outlet_angle_rad)))
        print('real alpha_1_deg = {};'.format(np.rad2deg(pump_design.inlet_flow_angle_rad)))
        print('real casting_volute_throat_area = {};'.format(pump_design.volute_throat_area))
        print('real casting_volute_height = {};'.format(pump_design.volute_height))
        print('real casting_diffusor_area_ratio = {};'.format(pump_design.diffuser_area_ratio))
        print('real specific_speed = {};'.format(pump_design.specific_speed))
        if pump_design.open_impeller:
            open_impeller_str = 'true'
        else:
            open_impeller_str = 'false'
        print('boolean open_impeller = {};'.format(open_impeller_str))


class PumpDesign:
    def __init__(self,
                 fluid: str,
                 n_d: float,
                 Q_d: float,
                 H_d: float,
                 head_coefficient: float = 0.5,
                 flow_coefficient_phi: float = 0,
                 epsilon: float=0.045e-3):
        self.fluid = fluid
        self.n_d = n_d
        self.Q_d = Q_d
        self.H_d = H_d
        self.nq = n_d * np.sqrt(Q_d) / (H_d ** (3 / 4))
        if flow_coefficient_phi == 0:
            self.flow_coefficient_phi = 3.7838e-3 * self.nq - 0.02919
        else:
            self.flow_coefficient_phi=flow_coefficient_phi
        self.epsilon = epsilon
        self.head_coefficient_psi = head_coefficient

    def pump_design(self,
                    blade_outlet_angle_deg,
                    number_of_blades,
                    inlet_impeller_hub_diameter,
                    incidence,
                    plot_figure = False,
                    flow_inlet_angle_deg = 90,
                    blade_thickness = 8e-3,
                    volute_height=5e-3,
                    diffusor_area_ratio=2,
                    blade_clearance = 0.05e-3) -> CentrifugalPumpSpecification:

        t = blade_thickness
        z = number_of_blades
        dh1 = inlet_impeller_hub_diameter
        i = incidence
        b3 = volute_height
        AR = diffusor_area_ratio
        s = blade_clearance
        beta2 = np.deg2rad(blade_outlet_angle_deg)
        alpha1 = np.deg2rad(flow_inlet_angle_deg)

        g = 9.81
        omega_d = self.n_d / 60 * 2 * np.pi
        d2 = 2 / omega_d * np.sqrt(g * self.H_d / self.head_coefficient_psi)
        u2 = omega_d * d2 / 2

        b2 = self.Q_d / (self.flow_coefficient_phi * u2 * (np.pi * d2 - t * z))
        A2 = np.pi * d2 * b2 - t * z * b2
        c2m = self.Q_d / A2


        ns = 51.55 * self.nq

        Kd = 1e-13 * ns ** 3 - 7e-9 * ns ** 2 + 1e-4 * ns + 0.2661
        d1 = Kd * d2
        b1 = (d1 - dh1) / 2
        A1 = np.pi * d1 * b1 - t * z * b1
        c1m = self.Q_d / A1

        u1 = omega_d * d1 / 2
        c_theta_1 = c1m / np.tan(alpha1)
        w_theta_1 = u1 - c_theta_1
        beta1 = np.arctan(c1m / w_theta_1) + i
        w1 = np.sqrt(c1m ** 2 + w_theta_1 ** 2)
        sigma = 1 - np.sin(beta2) / (z ** 0.7)
        w_theta_2 = c2m / np.tan(beta2)
        w2 = np.sqrt(c2m ** 2 + w_theta_2 ** 2)
        c_theta_2 = u2 - w_theta_2
        c_theta_2_slip = sigma * (u2 - w_theta_2)

        c2m_check = w_theta_2 * np.tan(beta2)
        diameter_size_check = d1 < d2
        blade_height_size_check = b2 < b1

        Kw = (w1 - w2) / w1
        Ks = (u2 - u1) / c_theta_2

        L0 = 1/360 * c_theta_2_slip * d2/2
        A4_check = self.Q_d / c_theta_2_slip * 1.25
        A4 = 0

        theta_list = []
        d_n_list = []
        d_1_list = []
        d_2_list = []
        for theta in range(1, 360, 1):
            coeff_1 = np.pi * L0 * 360**2 /(16 * self.Q_d * theta)
            coeff_2 = -(2 * np.pi * L0 * 360**2 * d2 / (16 * self.Q_d * theta) + 1 / 4)
            coeff_3 = (d2**2 * np.pi * L0 * 360**2 / (16 * self.Q_d * theta) - d2/4)
            """coeff_1 = np.pi * c_theta_2_slip * d2/2
            coeff_2 = - theta / (2 * np.pi) * self.Q_d
            coeff_3 = - d2 / 2 * theta / (2 * np.pi) * self.Q_d"""
            d_n_root = np.roots([coeff_1, coeff_2, coeff_3])
            d_n = d_n_root[0]
            c_theta_n = L0 * 360 / (d2/2 + (d_n - d2) / 4)
            d_v = (d_n - d2) / 2
            A4 = np.pi * (d_v * 1.25 /2)**2
            theta_list.append(theta/(360/(2*np.pi)))
            d_n_list.append(max(d_n/2 * 1.25, d2/2 + b3))
            d_2_list.append(d2/2)
            d_1_list.append(d1/2)


        if plot_figure:
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.plot(theta_list, d_n_list)
            ax.plot(theta_list, d_1_list)
            ax.plot(theta_list, d_2_list)
            # ax.set_rmax(0.25)
            ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
            ax.grid(True)

            ax.set_title("A line plot on a polar axis", va='bottom')
            plt.show()

        pump_design_parameters = CentrifugalPumpSpecification(
            fluid=self.fluid,
            specific_speed = self.nq,
            impeller_inner_diameter = d1,
            impeller_outer_diameter = d2,
            inlet_blade_height = b1,
            outlet_blade_height = b2,
            impeller_inlet_area= A1,
            impeller_outlet_area = A2,
            volute_height = b3,
            blade_clearance = s,
            blade_thickness = t,
            number_of_impeller_blades = z,
            roughness = self.epsilon,
            blade_inlet_angle_rad= beta1,
            blade_outlet_angle_rad = beta2,
            inlet_flow_angle_rad = alpha1,
            volute_throat_area = A4,
            diffuser_area_ratio=AR,
            diameter_size_check=diameter_size_check,
            blade_height_size_check=blade_height_size_check,
            Kw=Kw,
            Ks=Ks,
            Kw_check=Kw<0.25,
            Ks_check=Ks<0.9)

        return pump_design_parameters





# %% ..\00_core.ipynb 22
class CentrifugalPump:
    def __init__(self, fluid: str, temperature: float, data: CentrifugalPumpSpecification):
        self.open_impeller = data.open_impeller
        self.g = 9.81
        self.P = 101325
        self.T = temperature
        self.fluid = fluid
        self.nq = data.specific_speed
        self.d1 = data.impeller_inner_diameter
        self.d2 = data.impeller_outer_diameter
        self.b1 = data.inlet_blade_height
        self.b2 = data.outlet_blade_height
        self.b3 = data.volute_height
        self.s = data.blade_clearance
        self.t = data.blade_thickness
        self.z = data.number_of_impeller_blades
        self.epsilon = data.roughness
        self.beta1 = data.blade_inlet_angle_rad
        self.beta2 = data.blade_outlet_angle_rad
        self.alpha1 = data.inlet_flow_angle_rad
        self.r1 = self.d1 / 2
        self.r2 = self.d2 / 2
        self.A1 = np.pi * self.d1 * self.b1
        self.A2 = np.pi * self.d2 * self.b2
        self.A1r = self.A1 - self.t * self.z * self.b1
        self.A2r = self.A2 - self.t * self.z * self.b2
        self.A4 = data.volute_throat_area
        self.A5 = self.A4 * data.diffuser_area_ratio

    def _friction_colebrook(self, f, dh, Re):
        res = 1 / np.sqrt(f) + 2 * np.log10(
            (self.epsilon / dh) / 3.7 + 2.51 / (Re * np.sqrt(f))
        )
        return res

    def _friction_factor(self, a1, a2, temperature, v_avg):
        dh = 2 * (a2 * self.b2 + a1 * self.b1) / (a1 + a2 + self.b1 + self.b1)
        v = PropsSI("V", "T", temperature, "P", self.P, self.fluid)
        rho = PropsSI('D', 'T', temperature, 'P', self.P, self.fluid)
        Re = rho * v_avg * dh / v
        f0 = 0.25 * (np.log((self.epsilon / dh) / 3.7 + 5.74 / Re ** 0.9)) ** (-2)
        f = root_scalar(
            method="bisect",
            f=self._friction_colebrook,
            x0=f0,
            args=(dh, Re),
            bracket=[0.0001, 300],
        ).root
        return f

    def _pump_friction_loss(self, v1: float, v2: float, temperature: float):
        L = (self.d2 - self.d1) / (2 * np.sin(self.beta2))
        a1 = np.pi * self.d1 / self.z * np.sin(self.beta1)
        a2 = np.pi * self.d2 / self.z * np.sin(self.beta2)
        dh = 2 * (a2 * self.b2 + a1 * self.b1) / (a1 + a2 + self.b1 + self.b1)
        v_avg = (v1 + v2) / 2
        f = self._friction_factor(a1=a1, a2=a2, temperature=temperature, v_avg=v_avg)
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
        Hl_leak_ratio = (2.5 * self.s / self.d2) / (
            np.sqrt(self.b2 / self.d2 * (1 - self.d1 / self.d2))
            * (self.t / a2)
            * self.nq ** 0.1
            * np.sin(self.beta2) ** 1.2
            * np.sin(self.beta1) ** 0.4
        )
        return Hl_leak_ratio * Hth

    def _pump_leak_loss_closed(self, Q, Hth, omega, w1: float, w2: float, temperature: float):
        Hl_stat = Hth - (omega / 2) ** 2 * (self.d2 ** 2 - self.d1 ** 2) / (8 * self.g)
        a1 = np.pi * self.d1 / self.z * np.sin(self.beta1)
        a2 = np.pi * self.d2 / self.z * np.sin(self.beta2)
        w_avg = (w1 + w2) / 2
        f = self._friction_factor(a1=a1, a2=a2, temperature=temperature, v_avg=w_avg)
        L_gap = (self.d2 - self.d1) / 2
        v_gap = np.sqrt(2 * self.g * Hl_stat / (f * L_gap / self.s + 0.5 + 0.5))
        A_gap = np.pi * self.d2 * self.s
        Q_leak = v_gap * A_gap
        eta_vol = Q / (Q + Q_leak)
        return eta_vol

    def _pump_diffusion_loss(self, w1):
        Hl_d = 0.25 * w1 ** 2 / (2 * self.g)
        return Hl_d

    def _volute_friction_loss(self, Q, temperature):
        dh_volute = 4 * self.A4 / (4 * self.b3)
        d3 = self.d2 + 2 * self.b3
        c4m = Q / self.A4
        v = PropsSI("V", "T", temperature, "P", self.P, self.fluid)
        rho = PropsSI('D', 'T', temperature, 'P', self.P, self.fluid)
        Re = rho * c4m * dh_volute / v
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
        d4 = self.d2 + self.b3 + np.sqrt(self.A4 / np.pi) * 2
        c_theta_4 = c_theta_2 * self.d2 / d4
        Hl_h_v = 0.8 * (c_theta_4 ** 2 - c4m ** 2) / (2 * self.g)
        return Hl_h_v

    def _volute_diffuser_loss(self, Q):
        cp = 2 * ((self.A4 / self.A5) - (self.A4 / self.A5) ** 2)
        delta_diff = 1 - cp
        c4m = Q / self.A4
        Hl_diff_v = delta_diff * c4m ** 2 / (2 * self.g)
        return Hl_diff_v

    def _power_disk_friction(self, temperature, u2):
        v = PropsSI("V", "T", temperature, "P", self.P, self.fluid)
        rho = PropsSI('D', 'T', temperature, 'P', self.P, self.fluid)
        k = 7.3e-4 * (2 * v * 1e6 / (u2 * self.d2)) ** (1 / 8)
        P_disk = k * rho * u2 ** 3 * self.d2 * (self.d2 + 5 * self.s)
        return P_disk

    def pump_head(self, Q, n, temperature):
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
        Hl_f = self._pump_friction_loss(v1=w1, v2=w2, temperature=temperature)
        Hl_inc = self._pump_incidence_loss(beta1=self.beta1, beta1_fl=beta1_fl, w1=w1)
        if self.open_impeller:
            Hl_leak = self._pump_leak_loss_open(Hth=Hth)
        else:
            Hl_leak = self._pump_leak_loss_closed(Q=Q, Hth=Hth, omega=omega, w1=w1, w2=w2, temperature=temperature)
        Hl_d = self._pump_diffusion_loss(w1=w1)
        Hl_f_v = self._volute_friction_loss(Q=Q, temperature=temperature)
        Hl_r_v = self._volute_radial_loss(c2r=c2m)
        Hl_h_v = self._volute_head_loss(Q=Q, c_theta_2=c_theta_2_slip)
        Hl_diff_v = self._volute_diffuser_loss(Q=Q)
        """H_out = (
            Hth - Hl_f - Hl_inc - Hl_leak - Hl_d - Hl_f_v - Hl_r_v - Hl_h_v - Hl_diff_v
        )"""
        H_out = (
            Hth - Hl_f - Hl_inc - Hl_leak - Hl_d - Hl_f_v - Hl_r_v- Hl_diff_v
        )

        rho = PropsSI('D', 'T', temperature, 'P', self.P, self.fluid)

        P_disk = self._power_disk_friction(temperature=temperature, u2=u2)
        P_hyd = Hth * self.g * rho * Q
        P_out = H_out * self.g * rho * Q
        P_in = P_hyd + P_disk
        eta = P_out / P_in

        PumpPerformance = namedtuple('PumpPerformance', ['head',
                                                         'head_incidence',
                                                         'head_theoretical',
                                                         'pump_friction_loss',
                                                         'pump_incidence_loss',
                                                         'pump_diffusion_loss',
                                                         'volute_friction_loss',
                                                         'volute_radial_loss',
                                                         'volute_head_loss',
                                                         'volute_diffuser_loss',
                                                         'power_in',
                                                         'power_out',
                                                         'pump_efficiency'])

        return PumpPerformance(H_out,
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
                               eta)

    def create_performance_figures(self, min_Q, max_Q, dQ, min_n, max_n, dn):
        n = []
        Q = []
        H = []
        Hth = []
        inc_loss = []
        inc = []
        Hl_f = []
        Hl_d = []
        Hl_f_v = []
        Hl_r_v = []
        Hl_h_v = []
        Hl_diff_v = []
        power_in = []
        power_out = []
        eta = []

        for speed in np.arange(min_n, max_n, dn):
            for flow in np.arange(min_Q, max_Q, dQ):

                pump_performance = self.pump_head(Q=flow, n=speed, temperature=self.T)
                if pump_performance.head < 0:
                    break
                n.append(speed)
                Q.append(flow * 3600)
                H.append(pump_performance.head)
                Hth.append(pump_performance.head_theoretical)
                inc_loss.append(pump_performance.pump_incidence_loss)
                inc.append(pump_performance.head_incidence)
                Hl_f.append(pump_performance.pump_friction_loss)
                Hl_d.append(pump_performance.pump_diffusion_loss)
                Hl_f_v.append(pump_performance.volute_friction_loss)
                Hl_r_v.append(pump_performance.volute_radial_loss)
                Hl_h_v.append(pump_performance.volute_head_loss)
                Hl_diff_v.append(pump_performance.volute_diffuser_loss)
                power_in.append(pump_performance.power_in / 1000)
                power_out.append(pump_performance.power_out / 1000)
                eta.append(pump_performance.pump_efficiency)
        pump_dict = {'Speed':n, 'Flow [l/min]':Q, 'Head [m]':H, 'Head theoretical':Hth, 'Incidence loss':inc_loss, 'Incidence':inc, 'Head loss friction':Hl_f, 'Head loss volute friction':Hl_f_v, 'Head loss volute radial':Hl_r_v, 'Head loss volute':Hl_h_v, 'Head loss volute diffusor': Hl_diff_v, 'Power in':power_in, 'Power out': power_out, 'Pump efficiency':eta}
        pump_df = pd.DataFrame(data=pump_dict)

        H_min = 0
        H_max = max(H)
        dh = 1
        grid_x, grid_y = np.mgrid[min_Q * 3600:max_Q * 3600:0.0001 * 3600, H_min:H_max:0.1]

        points = np.column_stack((Q, H))
        values = eta

        grid_z0 = griddata(points=points, values=values, xi=(grid_x, grid_y), method='nearest')
        grid_z1 = griddata(points=points, values=values, xi=(grid_x, grid_y), method='linear')
        grid_z2 = griddata(points=points, values=values, xi=(grid_x, grid_y), method='cubic')

        # print(grid_z2)

        self.fig_flow_head_th = px.line(pump_df, x="Flow [l/min]", y='Head theoretical', color="Speed", title="Pump theoretical head")
        self.fig_flow_head_th.show()

        self.fig_flow_head = px.line(pump_df, x="Flow [l/min]", y="Head [m]", color="Speed", title="Pump head")
        self.fig_flow_head.show()

        self.fig_flow_power_in = px.line(pump_df, x="Flow [l/min]", y="Power in", color="Speed", title="Power in")
        self.fig_flow_power_in.show()

        self.fig_flow_power_out = px.line(pump_df, x="Flow [l/min]", y="Power out", color="Speed", title="Power out")
        self.fig_flow_power_out.show()

        self.fig_flow_eta =go.Figure(data = go.Contour(z=grid_z2.T, x=np.arange(min_Q * 3600, max_Q * 3600, dQ * 3600).T, y=np.arange(H_min, H_max, 0.1)))
        self.fig_flow_eta.show()

        self.fig_flow_incloss = px.line(pump_df, x="Flow [l/min]", y='Incidence loss', color="Speed", title="Incidence loss")
        self.fig_flow_incloss.show()

        self.fig_flow_friction = px.line(pump_df, x="Flow [l/min]", y='Head loss friction', color="Speed", title="Friction loss")
        self.fig_flow_friction.show()

        self.fig_flow_volute_friction = px.line(pump_df, x="Flow [l/min]", y='Head loss volute friction', color="Speed", title="Volute friction loss")
        self.fig_flow_volute_friction.show()

        self.fig_flow_volute_radial= px.line(pump_df, x="Flow [l/min]", y='Head loss volute radial', color="Speed", title="Volute radial loss")
        self.fig_flow_volute_radial.show()

        self.fig_flow_volute_head= px.line(pump_df, x="Flow [l/min]", y='Head loss volute', color="Speed", title="Volute head loss")
        self.fig_flow_volute_head.show()

        self.fig_flow_volute_diffusor= px.line(pump_df, x="Flow [l/min]", y='Head loss volute diffusor', color="Speed", title="Volute diffusor loss")
        self.fig_flow_volute_diffusor.show()


