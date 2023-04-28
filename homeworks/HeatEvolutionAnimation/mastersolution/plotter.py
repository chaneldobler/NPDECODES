import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

M = 20 # number of frames

data = np.loadtxt("solution.txt")
axis_length = int(np.sqrt(data.shape[1]))
data = data.reshape((data.shape[0], axis_length, axis_length))

step = int(data.shape[0] / M)
def animate(i):
    plt.clf()
    plt.imshow(data[i*step], vmin=0.)
    plt.colorbar()
    print(f"frame {i}", end="\r")

fig = plt.figure()
ani = animation.FuncAnimation(fig, animate, frames= M)
ani.save("animation.mp4")