import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np

sol = np.loadtxt("solution.txt")

if sol.shape[0] == 441:
    plt.imshow(np.reshape(sol, (21, 21)), vmin=-1, vmax=1)
    plt.colorbar()
    plt.savefig("solution.png")
else:
    sol = sol.reshape(sol.shape[0], 21, 21)
    fig = plt.figure()

    def animate(i):
        plt.clf()
        plt.imshow(sol[i, :, :], vmin=-1, vmax=1)
        plt.colorbar()
    ani = anim.FuncAnimation(fig, animate, frames=200, interval=10)
    ani.save("output.mp4")