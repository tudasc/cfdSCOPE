import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
import argparse
from tqdm import tqdm
import scipy.interpolate

# to show a progress bar during rendering
progress_bar = None


# for progressing the progress bar
def progress_callback_func(current_frame: int, total_frames: int):
    progress_bar.update(1)


def main():
    parser = argparse.ArgumentParser("Visualizes CFD data")
    parser.add_argument("input", type=str, action="store", help="name of input file. "
                                                                "If SAMPLE is given as filename, it will generate sample data to visualize and does NOT read the SAMPLE file (if present).")
    parser.add_argument("--nth_point", type=int, default=1, action="store",
                        help="distance of grid points to visualize with arrows (100 grid points and nth_point=20 ==> 5 arrows visualized)",
                        required=False)
    parser.add_argument("--output", type=str, default="visu.gif", action="store", help="name of output file",
                        required=False)
    parser.add_argument("--frames", type=int, default=100, action="store", help="number of frames to animate",
                        required=False)
    twod = parser.add_argument_group('arguments to use if a 2D plot is desired')
    twod.add_argument("--two_d", action="store_true", help="Plot only 2D", required=False)
    twod.add_argument("--slice_dim", default='z', choices=['x', 'y', 'z'], action="store",
                      help="which dimension to slice in case of 2D plot", required=False)
    twod.add_argument("--slice_point", type=int, default=0, action="store",
                      help="where to slice dimension (given in grid Points)", required=False)

    ARGS = parser.parse_args()

    if ARGS.input == "SAMPLE":
        x, y, z, u, v, w = get_sample_3d_input_data()
    else:
        x, y, z, u, v, w = read_3d_input_data(ARGS.input)

    if ARGS.two_d:
        x, y, u, v = get_2d_from_3d(x, y, z, u, v, w, ARGS.slice_dim, ARGS.slice_point)
        plot_2d(x, y, u, v, ARGS.frames, ARGS.nth_point, ARGS.output)
    else:
        plot_3d(x, y, z, u, v, w, ARGS.frames, ARGS.nth_point, ARGS.output)

    print("Visualization written to %s" % ARGS.output)


def get_sample_3d_input_data():
    x, y, z = np.meshgrid(np.linspace(-5, 5, 100), np.linspace(-5, 5, 100), np.linspace(-5, 5, 100))

    u = -y / np.sqrt(x ** 2 + y ** 2)
    v = x / np.sqrt(x ** 2 + y ** 2)
    w = np.zeros(shape=u.shape)
    return x, y, z, u, v, w


def read_3d_input_data(filename):
    # parse header
    f = open(filename)
    header = f.readline().strip().split(',')
    header = [int(h) for h in header]  # as int
    assert len(header) == 4

    num_x = header[0]
    num_y = header[1]
    num_z = header[2]
    resolution = header[3]

    data = np.genfromtxt(filename, delimiter=',', skip_header=1)

    # pressure would be the 0th value, but we don't visualize it
    # the pressure array would contain one nan for each line in the file that needs to be removed
    # ( this happens due to the trailing comma)
    # p = data[:, 0:-1:4]
    # p = p.reshape(num_x, num_y, num_z)
    # p is not needed but this would be the correct code to read it in
    u = data[:, 1::4]  # every 4th value is u
    # Reshape into the correct shape
    u = u.reshape(num_x, num_y, num_z)
    v = data[:, 2::4]
    v = v.reshape(num_x, num_y, num_z)
    w = data[:, 3::4]
    w = w.reshape(num_x, num_y, num_z)

    x, y, z = np.meshgrid(np.linspace(0, resolution * num_x, num_x),
                          np.linspace(0, resolution * num_y, num_y),
                          np.linspace(0, resolution * num_z, num_z))

    return x, y, z, u, v, w


# interpolates the arrow movement for animation
# moves the arrow according to the forces and interpolates the new force at the result position
# does not change or "calculate" any forces, it is just used to visualize the existing force field with movement
def interpolate_movement_3d(x, y, z, u, v, w, num_slices=100, slice_step=10, movement_factor=0.1):
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
    print("Interpolating movement:")
    for i in tqdm(range(num_slices)):
        x_new = xx + uu * movement_factor
        y_new = yy + vv * movement_factor
        z_new = zz + ww * movement_factor
        new_points = np.column_stack((x_new.ravel(), y_new.ravel(), z_new.ravel()))
        # Interpolate new x,y,z to find corresponding u,v,w
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


# generate a 3D animation
def plot_3d(x, y, z, u, v, w, frames, nth_point, outfile):
    # 3D
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    data = interpolate_movement_3d(x, y, z, u, v, w, num_slices=frames, slice_step=nth_point)

    def animate(i, plot_data):
        ax.clear()
        x, y, z, u, v, w = plot_data[i]
        q = ax.quiver(x, y, z, u, v, w)

        return q,

    print("Render Animation:")
    global progress_bar
    progress_bar = tqdm(total=frames)

    ani = FuncAnimation(fig, animate, interval=40, blit=True, repeat=True, frames=frames, fargs=[data])
    ani.save(outfile, dpi=300, writer=PillowWriter(fps=25), progress_callback=progress_callback_func)

    progress_bar.close()
    progress_bar = None


# get a 2D slice from the 3D data
def get_2d_from_3d(x, y, z, u, v, w, slice_dim, slice_point):
    assert slice_dim in ["x", "y", "z"]

    if slice_dim == "x":
        return x[slice_point, :, :], y[slice_point, :, :], v[slice_point, :, :], w[slice_point, :, :]
    if slice_dim == "y":
        return x[:, slice_point, :], y[:, slice_point, :], u[:, slice_point, :], w[:, slice_point, :]
    if slice_dim == "z":
        return x[:, :, slice_point], y[:, :, slice_point], u[:, :, slice_point], v[:, :, slice_point]


# same as above but for 2D
def interpolate_movement_2d(x, y, u, v, num_slices=100, slice_step=10, movement_factor=0.1):
    xx = x[0::slice_step, 0::slice_step]
    yy = y[0::slice_step, 0::slice_step]
    uu = u[0::slice_step, 0::slice_step]
    vv = v[0::slice_step, 0::slice_step]

    result = []
    result.append((xx, yy, uu, vv))
    print("Interpolating movement:")
    for i in tqdm(range(num_slices)):
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


def plot_2d(x, y, u, v, frames, nth_point, outfile):
    data = interpolate_movement_2d(x, y, u, v, num_slices=frames, slice_step=nth_point)
    fig, ax = plt.subplots()

    def animate(i, plot_data):
        ax.clear()
        x, y, u, v = plot_data[i]
        q = ax.quiver(x, y, u, v, np.linalg.norm((v, u)))

        return q,

    print("Render Animation:")
    global progress_bar
    progress_bar = tqdm(total=frames)

    ani = FuncAnimation(fig, animate, interval=40, blit=True, repeat=True, frames=frames, fargs=[data])
    ani.save(outfile, dpi=300, writer=PillowWriter(fps=25), progress_callback=progress_callback_func)

    progress_bar.close()
    progress_bar = None


if __name__ == '__main__':
    main()
