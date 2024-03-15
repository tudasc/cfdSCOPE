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

    x, y, z = np.meshgrid(np.linspace(-5, 5, 100), np.linspace(-5, 5, 100), np.linspace(-5, 5, 100))

    u = -y / np.sqrt(x ** 2 + y ** 2)
    v = x / np.sqrt(x ** 2 + y ** 2)
    w = np.zeros(shape=u.shape)

    FRAMES = 50

    data = get_plot_slices_3d(x, y, z, u, v, w,num_slices=FRAMES)



    def animate(i, plot_data):
        ax.clear()
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_zlim(-5, 5)
        x, y, z, u, v, w = plot_data[i]
        q = ax.quiver(x, y, z, u, v, w)

        return q,

    print("calculated interpolation, plotting data...")

    ani = FuncAnimation(fig, animate, interval=40, blit=True, repeat=True, frames=FRAMES, fargs=[data])
    ani.save("TLI.gif", dpi=300, writer=PillowWriter(fps=25))
    print("3D plot")
    exit()


def get_plot_slices_3d(x, y, z, u, v, w, num_slices=100, slice_step=10, movement_factor=0.1):
    # Reshape u, v for interpolation
    u_known_reshaped = u.reshape(-1)
    v_known_reshaped = v.reshape(-1)
    w_known_reshaped = w.reshape(-1)
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))

    xx = x[0::slice_step, 0::slice_step, 0::slice_step]
    yy = y[0::slice_step, 0::slice_step, 0::slice_step]
    zz = z[0::slice_step, 0::slice_step, 0::slice_step]
    uu = u[0::slice_step, 0::slice_step, 0::slice_step]
    vv = v[0::slice_step, 0::slice_step, 0::slice_step]
    ww = v[0::slice_step, 0::slice_step, 0::slice_step]

    result = []
    result.append((xx, yy, zz, uu, vv, ww))
    for i in range(num_slices):
        x_new = xx + uu * movement_factor
        y_new = yy + vv * movement_factor
        z_new = zz + ww * movement_factor
        new_points = np.column_stack((x_new.ravel(), y_new.ravel(), z_new.ravel()))
        # Interpolate new x, y,z to find corresponding u, v,w
        u_new = scipy.interpolate.griddata(points, u_known_reshaped, new_points, method='nearest')
        v_new = scipy.interpolate.griddata(points, v_known_reshaped, new_points, method='nearest')
        w_new = scipy.interpolate.griddata(points, w_known_reshaped, new_points, method='nearest')

        # Reshape interpolated u, v
        u_new = u_new.reshape(x_new.shape)
        v_new = v_new.reshape(y_new.shape)
        w_new = w_new.reshape(z_new.shape)

        xx = x_new
        yy = y_new
        zz = z_new
        uu = u_new
        vv = v_new
        ww = w_new
        result.append((xx, yy, zz, uu, vv, ww))

    return result


def get_plot_slices(x, y, u, v, num_slices=100, slice_step=10, movement_factor=0.1):
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

        x_new = xx + uu * movement_factor
        y_new = yy + vv * movement_factor
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
    show_3d_plot()

    exit()
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
        q = ax.quiver(x, y, u, v, np.linalg.norm((v, u)))

        return q,

    ani = FuncAnimation(fig, animate, interval=40, blit=True, repeat=True, frames=100, fargs=[data])
    ani.save("TLI.gif", dpi=300, writer=PillowWriter(fps=25))

    print("")

    exit()


if __name__ == '__main__':
    main()
