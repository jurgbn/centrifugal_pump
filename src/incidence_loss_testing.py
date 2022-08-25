import numpy as np
from matplotlib import pyplot as plt


def pump_incidence_loss(beta1, beta1_fl, w1):
    g = 9.81
    print(np.rad2deg(beta1 - beta1_fl))
    print(np.sin(2 * beta1_fl - beta1))
    lmda = np.sin(beta1_fl) / np.sin(2 * beta1_fl - beta1) - np.sqrt(
        (np.sin(beta1_fl) ** 2 - np.sin(beta1) * np.sin(2 * beta1_fl - beta1))
        / (np.sin(2 * beta1_fl - beta1)) ** 2
    )
    print(lmda)
    print("\n")
    Hl_inc = (
        w1 ** 2
        / (2 * g)
        * (1 / lmda) ** 2
        * (1 - lmda * np.sin(beta1_fl) / np.sin(beta1)) ** 2
    )
    return [Hl_inc, lmda]


beta1 = 36

lmda = []
Hl = []
inc = []

for incidence in np.arange(-30, 30, 1):
    beta1_fl = np.deg2rad(beta1 - incidence)
    res = pump_incidence_loss(beta1=np.deg2rad(beta1), beta1_fl=beta1_fl, w1=1)
    inc.append(incidence)
    lmda.append(res[1])
    Hl.append(res[0])

plt.figure(0)
plt.plot(inc, lmda)

plt.show
