import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

fname = sys.argv[1]

xs = []
ys = []
zs = []

with open(fname, 'r') as f:
    for line in f.readlines():
        tmp = line.rstrip("\n").split(" ")
        length = len(tmp)
        count = 0
        while(count<length):
            xs.append(float(tmp[count+0]))
            ys.append(float(tmp[count+1]))
            zs.append(float(tmp[count+2]))
            count=count+3

plotEnd = min(100000, len(xs));

fig = plt.figure()
ax = fig.add_subplot(2, 2, 1, projection='3d')
ax.scatter(xs[:plotEnd],ys[:plotEnd],zs[:plotEnd], s=1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

ax = fig.add_subplot(2, 2, 2)
ax.scatter(xs[:plotEnd],ys[:plotEnd], s=1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax = fig.add_subplot(2, 2, 3)
ax.scatter(xs[:plotEnd],zs[:plotEnd], s=1)
ax.set_xlabel("x")
ax.set_ylabel("z")
ax = fig.add_subplot(2, 2, 4)
ax.scatter(ys[:plotEnd],zs[:plotEnd], s=1)
ax.set_xlabel("y")
ax.set_ylabel("z")

plt.savefig(fname.rstrip("txt")+"png")
