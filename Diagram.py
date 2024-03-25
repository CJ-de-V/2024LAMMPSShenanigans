# Simple Diagram drawer for the scenarios I need exact representations for
import numpy as np
from matplotlib import pyplot as plt

lb = 1
Np = 33
Nt = 16
angle = 80 * (np.pi / 180)
lb_x = lb * np.cos(angle)
lb_y = lb * np.sin(angle)
tco = np.asarray([[0, x * lb] for x in range(Nt + 1)])  # Tether coordinates
pco = np.asarray(
    [[((-Np + 1) / 2 + x) * lb_x, Nt * lb + lb_y * np.abs((Np - 1) / 2 - x)] for x in range(Np)])  # PolymerCoordinates
maxy = np.max([np.max(tco[:, 1]), np.max(pco[:, 1])])  # max y value
maxx = np.max([np.max(tco[:, 0]), np.max(pco[:, 0])]) # Max x value
wco = []
if angle != 0:
    wco = np.asarray([[-maxy / np.tan(angle), maxy], [0, 0], [maxy / np.tan(angle), maxy]]) # Cone
else:
    wco = np.asarray([[-maxx, 0], [0, 0], [maxx, 0]])  # Cone

# tether and polymer
plt.plot(tco[:, 0], tco[:, 1], marker='o')
plt.plot(pco[:, 0], pco[:, 1], marker='o')
plt.plot(wco[:, 0], wco[:, 1], color='black')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.savefig('figure.png')
plt.show()
