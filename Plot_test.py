import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import animation

from matplotlib.animation import FuncAnimation, PillowWriter
import scipy.interpolate


def show_3d_plot():
    # 3D
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    x, y, z = np.meshgrid(np.arange(-0.8, 1, 0.2),
                          np.arange(-0.8, 1, 0.2),
                          np.arange(-0.8, 1, 0.8))

    u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
    w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
         np.sin(np.pi * z))

    ax.quiver(x, y, z, u, v, w, length=0.1, color='black')

    plt.show()


def animated_2d_plot():
    X, Y = np.mgrid[:2 * np.pi:10j, :2 * np.pi:5j]
    U = np.cos(X)
    V = np.sin(Y)

    fig, ax = plt.subplots(1, 1)
    Q = ax.quiver(X, Y, U, V, pivot='mid', color='r', units='inches')

    ax.set_xlim(-1, 7)
    ax.set_ylim(-1, 7)

    def update_quiver(num, Q, X, Y):
        """updates the horizontal and vertical vector components by a
        fixed increment on each frame
        """

        U = np.cos(X + num * 0.1)
        V = np.sin(Y + num * 0.1)

        Q.set_UVC(U, V)

        return Q,

    # you need to set blit=False, or the first set of arrows never gets
    # cleared on subsequent frames
    anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q, X, Y),
                                   interval=50, blit=False)
    fig.tight_layout()
    plt.show()


def test_animation():
    x1 = np.arange(0, -0.2, -0.002)
    y1 = np.arange(0, -0.2, -0.002)
    x2 = np.arange(3.9, 3.7, -0.002)
    y2 = np.arange(0, 1, 0.01)
    x3 = np.arange(0, 1.8, 0.018)
    y3 = np.array(x3 ** 2)

    fig, ax = plt.subplots()

    def animate(i):
        ax.clear()
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)
        line, = ax.plot(x1[0:i], y1[0:i], color='blue', lw=1)
        line2, = ax.plot(x2[0:i], y2[0:i], color='red', lw=1)
        line3, = ax.plot(x3[0:i], y3[0:i], color='purple', lw=1)
        point1, = ax.plot(x1[i], y1[i], marker='.', color='blue')
        point2, = ax.plot(x2[i], y2[i], marker='.', color='red')
        point3, = ax.plot(x3[i], y3[i], marker='.', color='purple')
        return line, line2, line3, point1, point2, point3,

    ani = FuncAnimation(fig, animate, interval=40, blit=True, repeat=True, frames=100)
    ani.save("TLI.gif", dpi=300, writer=PillowWriter(fps=25))


def get_plot_slices(x, y, u, v, num_slices=100, slice_step=10,movement_factor =0.1):
    xx = x[0::slice_step, 0::slice_step]
    yy = y[0::slice_step, 0::slice_step]
    uu = u[0::slice_step, 0::slice_step]
    vv = v[0::slice_step, 0::slice_step]

    result = []
    result.append((xx, yy, uu, vv))
    for i in range(num_slices):
        # Reshape u, v for interpolation
        u_known_reshaped = u.reshape(-1)
        v_known_reshaped = v.reshape(-1)
        points = np.column_stack((x.ravel(), y.ravel()))

        x_new = xx + uu*movement_factor
        y_new = yy + vv*movement_factor
        new_points = np.column_stack((x_new.ravel(), y_new.ravel()))
        # Interpolate new x, y to find corresponding u, v
        u_new = scipy.interpolate.griddata(points, u_known_reshaped, new_points, method='nearest')
        v_new = scipy.interpolate.griddata(points, v_known_reshaped, new_points, method='nearest')

        # Reshape interpolated u, v
        u_new = u_new.reshape(x_new.shape)
        v_new = v_new.reshape(y_new.shape)

        xx = x_new
        yy = y_new
        uu = u_new
        vv = v_new
        result.append((xx, yy, uu, vv))

    return result

def main():
    # full plot of all known forces
    x, y = np.meshgrid(np.linspace(-5, 5, 100), np.linspace(-5, 5, 100))

    u = -y / np.sqrt(x ** 2 + y ** 2)
    v = x / np.sqrt(x ** 2 + y ** 2)

    # full plot
    fig, ax = plt.subplots()
    q = ax.quiver(x, y, u, v)
    plt.show()
    data = get_plot_slices(x, y, u, v)

    def animate(i, plot_data):
        ax.clear()
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        x, y, u, v = plot_data[i]
        q = ax.quiver(x, y, u, v)

        return q,

    ani = FuncAnimation(fig, animate, interval=40, blit=True, repeat=True, frames=100,fargs=[data])
    ani.save("TLI.gif", dpi=300, writer=PillowWriter(fps=25))

    print("")

    exit()



if __name__ == '__main__':
    main()
