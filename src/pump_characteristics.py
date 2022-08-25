import numpy as np
from matplotlib import pyplot as plt
from src.pump_performance import CentrifugalPumpSpecification, CentrifugalPump

""" Constants"""

g = 9.81

""" Fluid data """

rho = 1000
t = 300


""" Pump charateristics """

pump_data = CentrifugalPumpSpecification(
    fluid="INCOMP::Water",
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
pump = CentrifugalPump(data=pump_data)


""" Pump operation parameters """
plt.figure(0)
plt.figure(1)
plt.figure(2)
plt.figure(3)
plt.figure(4)
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
legend = []

for n in np.arange(500, 3000, 500):
    legend.append(n)
    Q = []
    H = []
    Hth = []
    inc_loss = []
    inc = []
    power_in = []
    power_out = []
    eta = []
    for flow in np.arange(0.0001, 0.1, 0.0001):
        Q.append(flow * 1000 * 60)
        H_list = pump.pump_head(Q=flow, n=n, T=t)
        H.append(H_list[0])
        Hth.append(H_list[2])
        inc_loss.append(H_list[4])
        inc.append(H_list[1])
        power_in.append(H_list[10] / 1000)
        power_out.append(H_list[11] / 1000)
        eta.append(H_list[12])
    plt.figure(0)
    plt.plot(Q, H, label=str(n) + " rpm")
    plt.legend()

    plt.figure(1)
    color = "tab:red"
    ax1.plot(Q, inc_loss, color=color)
    color = "tab:blue"
    ax2.plot(Q, inc, color=color)

    plt.figure(2)
    plt.plot(Q, power_in, label=str(n) + " rpm")
    plt.legend()

    plt.figure(3)
    plt.plot(Q, eta, label=str(n) + " rpm")
    plt.legend()

    plt.figure(4)
    plt.plot(Q, Hth, label=str(n) + " rpm")
    plt.legend()


plt.figure(1)
color = "tab:red"
ax1.set_xlabel("Flow Q (l/min)")
ax1.set_ylabel("Incidence head loss [m]", color=color)
ax1.tick_params(axis="y", labelcolor=color)
ax1.set_ylim([0, 10])
color = "tab:blue"
ax2.set_ylabel("Incidence [deg]", color=color)
ax2.tick_params(axis="y", labelcolor=color)
fig.tight_layout()

plt.figure(0)
plt.ylim([0, None])
plt.xlim([0, 3000])
plt.xlabel("Flow rate [l/m]")
plt.ylabel("Pump head [m]")

plt.figure(2)
plt.ylim([0, None])
plt.xlim([0, 3000])
plt.xlabel("Flow rate [l/m]")
plt.ylabel("Pump power [kW]")

plt.figure(3)
plt.ylim([0, 1])
plt.xlim([0, 3000])
plt.xlabel("Flow rate [l/m]")
plt.ylabel("Pump efficiency")

plt.figure(4)
plt.ylim([0, None])
plt.xlim([0, 3000])
plt.xlabel("Flow rate [l/m]")
plt.ylabel("Pump head th [m]")


plt.show()
