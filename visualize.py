import os
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
import argparse
from tqdm import tqdm
import scipy.interpolate

# quiver_size_args = {'length': 0.5, 'arrow_length_ratio': 0.3}
quiver_size_args = {'length': 1.5, 'normalize': True}

# to show a progress bar during rendering
progress_bar = None


# for progressing the progress bar
def progress_callback_func(current_frame: int, total_frames: int):
    progress_bar.update(1)


def main():
    parser = argparse.ArgumentParser("Visualizes CFD data")
    parser.add_argument("input", type=str, action="store", help="name of input file. ")
    parser.add_argument("--static", action="store_true", help="dont animate the plot.")
    #parser.add_argument("--nth_point", type=int, default=1, action="store",
    #                   help="distance of grid points to visualize with arrows (100 grid points and nth_point=20 ==> 5 arrows visualized)",
    #                   required=False)
    parser.add_argument("--num_particles", type=int, action="store", help="number of particles to trace", default=1000)
    parser.add_argument("--output", type=str, default="visu.gif", action="store", help="name of output file",
                        required=False)
    parser.add_argument("--frames_per_step", type=int, default=4, action="store", help="number of frames per time step to animate",
                        required=False)
    #parser.add_argument("--fps", type=int, default=25, action="store", help="frames per second",
    #                    required=False)
    parser.add_argument("--time_step", type=float, action="store", help="Time step size in seconds", required=True)
    twod = parser.add_argument_group('arguments to use if a 2D plot is desired')
    twod.add_argument("--two_d", action="store_true", help="Plot only 2D", required=False)
    twod.add_argument("--slice_dim", default='z', choices=['x', 'y', 'z'], action="store",
                      help="which dimension to slice in case of 2D plot", required=False)
    twod.add_argument("--slice_point", type=int, default=0, action="store",
                      help="where to slice dimension (given in grid Points)", required=False)

    ARGS = parser.parse_args()

    # x, y, z, u, v, w = read_3d_input_data(ARGS.input)
    data = read_input_files(ARGS.input)

    if ARGS.two_d:
        x, y, u, v = get_2d_from_3d(x, y, z, u, v, w, ARGS.slice_dim, ARGS.slice_point)
        if ARGS.static:
            plot_2d_static(x, y, u, v, ARGS.nth_point, ARGS.output)
        else:
            plot_2d(x, y, u, v, ARGS.frames_per_step, ARGS.nth_point, ARGS.output)
    else:
        if ARGS.static:
            plot_3d_static(x, y, z, u, v, w, ARGS.nth_point, ARGS.output)
        else:
            plot_3d(data, ARGS.time_step, ARGS.frames_per_step, ARGS.num_particles, ARGS.output)

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

    xx = []
    yy = []
    zz = []
    p = []
    u = []
    v = []
    w = []
    for z in range(num_z):
        for y in range(num_y):
            for x in range(num_x):
                xx.append(x * resolution)
                yy.append(y * resolution)
                zz.append(z * resolution)

                p.append(data[z * num_y + y, x * 4 + 0])
                u.append(data[z * num_y + y, x * 4 + 1])
                v.append(data[z * num_y + y, x * 4 + 2])
                w.append(data[z * num_y + y, x * 4 + 3])

    s = (num_x, num_y, num_z)

    return (np.array(xx).reshape(s), np.array(yy).reshape(s), np.array(zz).reshape(s),
            np.array(u).reshape(s), np.array(v).reshape(s), np.array(w).reshape(s))


def read_input_files(dir):
    pattern = r"fields_(\d+)\.txt"

    result = {}
    for fname in os.listdir(dir):
        match = re.match(pattern, fname)
        if match:
            index = int(match.group(1))  # Extract the index as an integer
            data = read_3d_input_data(os.path.join(dir, fname))
            assert index not in result  # not read same file twice
            result[index] = data

    # Sort the list of tuples based on keys
    dict_items_sorted = sorted(list(result.items()), key=lambda x: x[0])

    return [elem[1] for elem in dict_items_sorted]

# interpolates the arrow movement for animation
# moves the arrow according to the forces and interpolates the new force at the result position
# does not change or "calculate" any forces, it is just used to visualize the existing force field with movement
def interpolate_movement_3d(data, xx, yy, zz, x_bound, y_bound, z_bound, time_step, num_slices=100):
    result = []
    xx=xx.ravel()
    yy = yy.ravel()
    zz = zz.ravel()
    result.append((xx, yy, zz, np.zeros_like(xx)))  # starting position

    # don't show arrows moving out of bounds

    dt = time_step / num_slices


    print("Interpolating movement:")
    progress_bar = tqdm(total=num_slices * len(data))
    target_shape = xx.shape
    assert xx.shape == yy.shape == zz.shape
    for x, y, z, u, v, w in data:
        # interpolation of particle forces
        u_interp = scipy.interpolate.LinearNDInterpolator((x.ravel(), y.ravel(), z.ravel()), u.ravel())
        v_interp = scipy.interpolate.LinearNDInterpolator((x.ravel(), y.ravel(), z.ravel()), v.ravel())
        w_interp = scipy.interpolate.LinearNDInterpolator((x.ravel(), y.ravel(), z.ravel()), w.ravel())

        for i in range(num_slices):
            # uu = u_interp((xx.ravel(), yy.ravel(), zz.ravel())).reshape(target_shape)
            # vv = v_interp((xx.ravel(), yy.ravel(), zz.ravel())).reshape(target_shape)
            # ww = w_interp((xx.ravel(), yy.ravel(), zz.ravel())).reshape(target_shape)
            uu = u_interp((xx, yy, zz))
            vv = v_interp((xx, yy, zz))
            ww = w_interp((xx, yy, zz))
            # move particles in plot
            xx = xx + uu * dt
            yy = yy + vv * dt
            zz = zz + ww * dt

            # stop particles at field bounds when they "hit the wall"
            xx[xx > x_bound] = x_bound
            yy[yy > y_bound] = y_bound
            zz[zz > z_bound] = z_bound

            vals = np.append([uu], np.append([vv], [ww], axis=0), axis=0)
            speed = np.linalg.norm(vals, axis=0)

            result.append((xx, yy, zz, speed))  # resulting position
            progress_bar.update(1)

    return result
def generate_animation_data_3d(data, nth_point, num_slices=100, slice_step=10, movement_factor=0.1):

    result = []

    for x, y, z, u, v, w in data:
        xx = x[0::nth_point, 0::nth_point, 0::nth_point]
        yy = y[0::nth_point, 0::nth_point, 0::nth_point]
        zz = z[0::nth_point, 0::nth_point, 0::nth_point]

        uu = u[0::nth_point, 0::nth_point, 0::nth_point]
        vv = v[0::nth_point, 0::nth_point, 0::nth_point]
        ww = w[0::nth_point, 0::nth_point, 0::nth_point]

        for i in range(num_slices):
            result.append((xx, yy, zz, uu, vv, ww))

    return result


#
def plot_3d_static(x, y, z, u, v, w, nth_point, outfile):  #
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.set_xlim([0, x.max()])
    ax.set_ylim([0, y.max()])
    ax.set_zlim([0, z.max()])

    q = ax.quiver(x, y, z, u, v, w, **quiver_size_args)

    plt.savefig(outfile, dpi=300)


# generate a 3D animation
def plot_3d(data, time_step, frames_per_step, num_particles, outfile):
    # 3D
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    x = data[0][0]
    y = data[0][1]
    z = data[0][2]

    # don't show arrows moving out of bounds
    xlim = x.max()
    ylim = y.max()
    zlim = z.max()

    # take the sample of points to visualize
    #xx = x[0::nth_point, 0::nth_point, 0::nth_point]
    #yy = y[0::nth_point, 0::nth_point, 0::nth_point]
    #zz = z[0::nth_point, 0::nth_point, 0::nth_point]
    xx = np.random.rand(num_particles) * xlim
    yy = np.random.rand(num_particles) * ylim
    zz = np.random.rand(num_particles) * zlim

    anim_data = interpolate_movement_3d(data, xx, yy, zz,xlim,ylim,zlim, time_step, num_slices=frames_per_step)

    #data = generate_animation_data_3d(data, nth_point, num_slices=frames, slice_step=nth_point)


    def animate(i, plot_data, xlim, ylim, zlim):
        ax.clear()
        #x, y, z = plot_data[i]
        x, y, z, speed = plot_data[i]

        ax.set_xlim([0, xlim])
        ax.set_ylim([0, ylim])
        ax.set_zlim([0, zlim])

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        s = ax.scatter(x, y, z, c=speed, cmap='plasma', marker='.')
        return s,

        #q = ax.quiver(x, y, z, u, v, w, **quiver_size_args)
        #return q,

    print("Render Animation:")
    global progress_bar
    progress_bar = tqdm(total=len(anim_data))

    time_per_anim_frame = time_step / frames_per_step
    time_per_anim_frame_ms = (time_step / frames_per_step) * 1000
    anim_fps = 1.0 / time_per_anim_frame

    print(f"Animation interval: {time_per_anim_frame_ms}")
    print(f"FPS: {anim_fps}")

    ani = FuncAnimation(fig, animate, interval=time_per_anim_frame_ms, blit=True, repeat=True, frames=len(anim_data),
                        fargs=[anim_data, xlim, ylim, zlim])
    ani.save(outfile, dpi=300, writer=PillowWriter(fps=anim_fps), progress_callback=progress_callback_func)

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


def plot_2d_static(x, y, u, v, nth_point, outfile):  #
    fig = plt.figure()
    ax = fig.add_subplot()

    ax.set_xlim([0, x.max()])
    ax.set_ylim([0, y.max()])

    q = ax.quiver(x, y, u, v, **quiver_size_args)
    plt.savefig(outfile, dpi=300)


def plot_2d(x, y, u, v, frames, nth_point, outfile):
    data = interpolate_movement_2d(x, y, u, v, num_slices=frames, slice_step=nth_point)
    fig, ax = plt.subplots()

    xlim = x.max()
    ylim = y.max()

    def animate(i, plot_data, xlim, ylim):
        ax.clear()
        x, y, u, v = plot_data[i]

        ax.set_xlim([0, xlim])
        ax.set_ylim([0, ylim])

        s = ax.scatter(x, y, z, marker='.')

        return s,

    print("Render Animation:")
    global progress_bar
    progress_bar = tqdm(total=frames)

    ani = FuncAnimation(fig, animate, interval=40, blit=True, repeat=True, frames=frames, fargs=[data, xlim, ylim])
    ani.save(outfile, dpi=300, writer=PillowWriter(fps=25), progress_callback=progress_callback_func)

    progress_bar.close()
    progress_bar = None


if __name__ == '__main__':
    main()
