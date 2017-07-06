import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
from PIL import Image

n = int(sys.argv[1])
scale = 200000;
width = 1000
height = 1000

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
    ax.scatter(xs[:10:plotEnd],ys[:10:plotEnd],zs[:10:plotEnd], s=1)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax = fig.add_subplot(2, 2, 2)
    ax.scatter(xs[:plotEnd],ys[:plotEnd], s=1)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.axis([-scale, scale, -scale, scale])
    ax = fig.add_subplot(2, 2, 3)
    ax.scatter(xs[:plotEnd],zs[:plotEnd], s=1)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.axis([-scale, scale, -scale, scale])
    ax = fig.add_subplot(2, 2, 4)
    ax.scatter(ys[:plotEnd],zs[:plotEnd], s=1)
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.axis([-scale, scale, -scale, scale])

    plt.savefig(fname.rstrip("txt"))

    img = densityMap(xs, ys)
    im = Image.new('L', (width+1, height+1))
    im.putdata(img.flatten().tolist())
    im.save(fname.rstrip("txt")+"Map"+".jpg")

def densityMap(x, y):
    maxX = max(x)
    maxY = max(y)
    map = np.zeros((width+1, height+1))
    for i in range(len(x)):
        xcor = int(float(x[i])/maxX*width/2)+width/2
        ycor = int(float(y[i])/maxY*height/2)+height/2
        map[xcor][ycor]+=500
    return map

for i in range(n+1):
    fname = "./pics/multiSim"+str(i);
    plotter(fname)
