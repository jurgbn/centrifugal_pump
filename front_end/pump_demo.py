from CentrifugalPumpDesign.centrifugal_pump import PumpDesign, CentrifugalPump
import streamlit as st
import numpy as np

st.title("Initial pump design tool")
# Create a pump object
st.subheader("Pump design point")
design_speed_rpm = st.slider("Design Speed [rpm]", 100, 3600, 1750)
design_head_m = st.slider("Design Head [m]", 1, 100, 10)
design_flow_m3_per_hour = st.slider("Design Flow [m3/h]", 1.0, 2e2, 90.0)
pump_design = PumpDesign(
    fluid="water",
    n_d=design_speed_rpm,
    Q_d=design_flow_m3_per_hour / 3600,
    H_d=design_head_m,
)

st.subheader("Pump design parameters")
blade_outlet_angle_deg = st.slider("Blade Outlet Angle [deg]", 1, 90, 45)
number_of_blades = st.slider("Number of Blades", 1, 10, 4)
inlet_impeller_hub_diameter_m = st.slider(
    "Inlet Impeller Hub Diameter [m]", 0.01, 1.0, 0.05
)
pump_design_result = pump_design.pump_design(
    blade_outlet_angle_deg=blade_outlet_angle_deg,
    number_of_blades=number_of_blades,
    inlet_impeller_hub_diameter=inlet_impeller_hub_diameter_m,
    incidence=np.deg2rad(0),
)
st.write(
    "Impeller outer diameter {}[m]".format(
        round(pump_design_result.impeller_outer_diameter, 2)
    )
)
st.write(
    "Impeller inlet angle {}[deg]".format(
        round(np.rad2deg(pump_design_result.blade_inlet_angle_rad), 2)
    )
)
st.write(
    "Kw check {} < 0.25 is {}".format(
        round(pump_design_result.Kw, 3), pump_design_result.Kw_check
    )
)
st.write(
    "Ks check {} < 0.9 is {}".format(
        round(pump_design_result.Ks, 3), pump_design_result.Ks_check
    )
)

st.subheader("Pump performance")

pump_performance = CentrifugalPump(
    fluid="INCOMP::Water", temperature=300, data=pump_design_result
)

test_flow_m3_per_hour = st.slider("Test Flow [m3/h]", 1.0, 2e2, 90.0)
test_speed_rpm = st.slider("Test Speed [rpm]", 100, 3600, 1750)
pump_performance_result = pump_performance.pump_head(
    Q=test_flow_m3_per_hour / 3600, n=test_speed_rpm, temperature=300
)
st.write(pump_performance_result)
