import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

n = int(sys.argv[1])

def plotter(fname):
    xs = []
    ys = []
    zs = []
    with open(fname, 'r') as f:
        for line in f.readlines():
            tmp = line.rstrip("\n").split(" ")
            xs.append(float(tmp[0]))
            ys.append(float(tmp[1]))
            zs.append(float(tmp[2]))

    plotEnd = len(xs);

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
    ax.axis([-100, 100, -100, 100])
    ax = fig.add_subplot(2, 2, 3)
    ax.scatter(xs[:plotEnd],zs[:plotEnd], s=1)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.axis([-100, 100, -100, 100])
    ax = fig.add_subplot(2, 2, 4)
    ax.scatter(ys[:plotEnd],zs[:plotEnd], s=1)
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.axis([-100, 100, -100, 100])

    plt.savefig(fname.rstrip("txt"))

for i in range(n+1):
    fname = "./pics/multiSim"+str(i);
    plotter(fname)
