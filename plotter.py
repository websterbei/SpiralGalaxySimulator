import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

fname = sys.argv[1]

xs = []
ys = []
zs = []

plotEnd = 100000;

with open(fname, 'r') as f:
    for line in f.readlines():
        tmp = line.rstrip("\n").split(" ")
        xs.append(float(tmp[0]))
        ys.append(float(tmp[1]))
        zs.append(float(tmp[2]))

fig = plt.figure()
ax = fig.add_subplot(2, 2, 1, projection='3d')
ax.plot(xs[:plotEnd],ys[:plotEnd],zs[:plotEnd])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

ax = fig.add_subplot(2, 2, 2)
ax.plot(xs[:plotEnd],ys[:plotEnd])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax = fig.add_subplot(2, 2, 3)
ax.plot(xs[:plotEnd],zs[:plotEnd])
ax.set_xlabel("x")
ax.set_ylabel("z")
ax = fig.add_subplot(2, 2, 4)
ax.plot(ys[:plotEnd],zs[:plotEnd])
ax.set_xlabel("y")
ax.set_ylabel("z")

plt.savefig(fname.rstrip("txt")+"jpg")
