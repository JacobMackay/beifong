import os
import numpy as np
import mitsuba

# Set the desired mitsuba variant
# mitsuba.set_variant('scalar_rgb')
# mitsuba.set_variant('scalar_mono')
mitsuba.set_variant('scalar_spectral')
# mitsuba.set_variant('gpu_spectral')
# mitsuba.set_variant('gpu_rgb')
# mitsuba.set_variant('packet_spectral')

from mitsuba.core import Bitmap, Struct, Thread
from mitsuba.core.xml import load_file

from matplotlib import pyplot as plt

# Animation bits
from mitsuba.core.xml import load_dict
from mitsuba.core import Vector3f, Transform4f
from mitsuba.core import ScalarTransform4f
import matplotlib.animation as animation
from matplotlib import cm
# matplotlib.use('Agg')

# ---------------------------------------------------------
# Other settings
# N_BINS = 500
N_BINS = 400
DR = 0.1
# resx = 25
# resy = 1600
resx = 600
resy = 600
spp = 50
power = 1000

# =============================================================================
# Load up dicts
targ_size = Transform4f.scale([100e-3, 50e-3, 1])
# targ_size = Transform4f.scale([500e-3, 500e-3, 1])
# targ_size = Transform4f.scale([200e-3, 200e-3, 1])
targ_loc = Transform4f.look_at(origin=[2, 0, 20e-2],
                               target=[0, 0, 20e-2],
                               up=[0, 0, 1])
targ_transform = targ_loc*targ_size
gnd_size = Transform4f.scale([20, 20, 1])
gnd_loc = Transform4f.look_at(origin=[0, 0, 0],
                              target=[0, 0, 1],
                              up=[0, 1, 0])
gnd_transform = gnd_loc*gnd_size
emit_size = Transform4f.scale([20e-3, 50e-3, 1])
# emit_size = Transform4f.scale([6.5e-3, 6.5e-3, 1])
# emit_loc = Transform4f.look_at(origin=[20e-2, 0, 20e-2],
#                                target=[1, 0, 20e-2],
#                                up=[0, 0, 1])
# emit_transform = emit_loc*emit_size
emit_transform = emit_size

# base_transform = Transform4f.look_at(origin=[20e-2, 0, 20e-2],
#                                      target=[1, 0, 20e-2],
#                                      up=[0, 0, 1])
# base_transform = Transform4f.look_at(origin=[0, 0, 20e-2],
#                                      target=[1, 0, 20e-2],
#                                      up=[0, 0, 1])

base_transform = Transform4f.look_at(origin=[0, 0, 1],
                                     target=[1, 0, 1],
                                     up=[0, 0, 1])

# car_loc = Transform4f.look_at(origin=[-2, 0, 1],
#                                target=[0, 0, 1],
#                                up=[0, 0, 1])
car_loc = Transform4f.look_at(origin=[-10, 0, 100e-2],
                               target=[0, 1, 100e-2],
                               up=[0, 0, 1])

# car_size = Transform4f.scale([0.0005, 0.0005, 0.0005])
car_size = Transform4f.scale([0.001, 0.001, 0.001])
# car_size = Transform4f.scale([0.0015, 0.0015, 0.0015])
# car_size = Transform4f.scale([0.0025, 0.0025, 0.0025])
car_transform = car_loc*car_size
# car_transform = car_loc

bsdfs = load_dict({
    'type': 'twosided',
    'id': 'material',
    'bsdf': {
        # 'type': 'diffuse',
        'type': 'plastic',
        'specular_reflectance': {
        # 'reflectance': {
            'type': 'spectrum',
            'value': [(8292683, 0.5), (8717949, 0.5)],
        },
    },
})

bsdft = load_dict({
    'type': 'twosided',
    'id': 'material',
    'bsdf': {
        # 'type': 'diffuse',
        'type': 'conductor',
        # 'reflectance': {
        'specular_reflectance': {
            'type': 'spectrum',
            'value': [(8292683, 1), (8717949, 1)],
        },
    },
})

bsdfc = load_dict({
    'type': 'twosided',
    'id': 'material',
    'bsdf': {
        # 'type': 'diffuse',
        # 'reflectance': {
        'type': 'roughconductor',
        # 'type': 'plastic',
        'specular_reflectance': {
            'type': 'spectrum',
            'value': [(8292683, 1), (8717949, 1)],
        },
    },
})

targ = load_dict({
    'type': 'rectangle',
    'id': 'target',
    'to_world': targ_transform,
    'bsdf': bsdfs,
})
gnd = load_dict({
    'type': 'rectangle',
    'id': 'gnd',
    'to_world': gnd_transform,
    'bsdf': bsdfs,
})

# car = load_dict({
#     'type': 'ply',
#     'filename': './bus_ply.ply',
#     # 'filename': './Car-body.ply',
#     # 'filename': './Motorbike_ply.ply',
#     'to_world': car_transform,
#     'bsdf': bsdfc,
# })

car = load_dict({
    'type': 'obj',
    'filename': './Bus.obj',
    # 'filename': './Car-body.ply',
    # 'filename': './Motorbike_ply.ply',
    'to_world': car_transform,
    'bsdf': bsdfc,
})

ints = load_dict({
    'type': 'range',
    'integrator': {
        'type': 'pathlength',
    },
    'dr': DR,
    'bins': N_BINS,
})
# emit = load_dict({
#     'type': 'spot',
#     'cutoff_angle': 25,
#     'beam_width': 20,
#     'intensity': {
#         'type': 'spectrum',
#         'value': [(8292683, 1), (8717949, 1)],
#     }
# })
emit = load_dict({
    'type': 'rectangle',
    'to_world': emit_transform,
    'emitter': {
        'type': 'wignertransmitter',
        'radiance': {
            'type': 'spectrum',
            'value': [(8292683, power), (8717949, power)],
        },
    },
})
film = load_dict({
    'type': 'hdrfilm',
    'rfilter': {'type': 'box'},
    'width': resx,
    'height': resy,
})
sampler = load_dict({
    'type': 'independent',
    'sample_count': spp,
    'id': 'sampler',
})
sen = load_dict({
    'type': 'perspective',
    'near_clip': 0.100000,
    'far_clip': 100.0,
    'fov_axis': 'x',
    'fov': 5.0,
    'sampler': sampler,
    'film': film,
})
# =============================================================================

# =============================================================================
# Do a render to test

# lorigin = Vector3f(0, 0, 0)
# lorigin = Vector3f(20e-2, 0, 20e-2)
# lorigin = Vector3f(0, 0, 20e-2)
lorigin = Vector3f(0, 0, 1)
# boresight = Vector3f(0, -1, 0)
# boresight = Vector3f(-1, 0, 20e-2)
boresight = Vector3f(-1, 0, 1)
# boresight = Vector3f(0, 1, 20e-2)

# (origin=[20e-2, 0, 20e-2],
#                                target=[1, 0, 20e-2],
#                                up=[0, 0, 1])
i = 10
STEP_SIZE = 0
# rotation_cur = Transform4f.rotate(Vector3f(0, 0, 1), i*STEP_SIZE)
rotation_cur = Transform4f.rotate(Vector3f(0, 0, 1), i*STEP_SIZE)
new_boresight_c = rotation_cur.transform_vector(boresight)
new_up_c = rotation_cur.transform_vector(Vector3f(0, 0, 1))
to_world_cur = Transform4f.look_at(lorigin, new_boresight_c, new_up_c)

sen = load_dict({
    'type': 'perspective',
    'near_clip': 0.100000,
    'far_clip': 100.0,
    'fov_axis': 'y',
    'fov': 75.0,
    # 'fov_axis': 'x',
    # 'fov': 45.0,
    'sampler': sampler,
    'film': film,
    'to_world': to_world_cur,
})

emit = load_dict({
    'type': 'rectangle',
    'emitter': {
        'type': 'wignertransmitter',
        'radiance': {
            'type': 'spectrum',
            'value': [(8292683, power), (8717949, power)],
        },
    },
    'to_world': to_world_cur*emit_transform,
})

# emit = load_dict({
#     'type': 'spot',
#     'cutoff_angle': 25,
#     'beam_width': 20,
#     'intensity': {
#         'type': 'spectrum',
#         'value': [(8292683, 1), (8717949, 1)],
#     },
#     'to_world': to_world_cur,
# })

scene = load_dict({
    'type': 'scene',
    'integrator': ints,
    'sensor': sen,
    'emitter': emit,
    'so': targ,
    's1': gnd,
    's2': car,
})

# scene.integrator().render(scene, scene.sensors()[0])
scene.integrator().render(scene, sen)

bmp2 = film.bitmap(raw=True)
bmp22_np = np.array(bmp2)

# After rendering, the rendered data is stored in the film
film = scene.sensors()[0].film()

# Write out rendering as high dynamic range OpenEXR file
# film.set_destination_file('./' + name + '.exr')
# film.develop()

bmp2 = film.bitmap(raw=True)
bmp22_np = np.array(bmp2)
print(bmp22_np.shape)
# print(bmp22_np.tolist())

fig1 = plt.figure(figsize=(14, 7))
axes1 = fig1.add_subplot()

axes1.imshow(bmp22_np[:,:,0:3], interpolation='nearest')

fig1.show()

# Section for images ----------------------------------------------------------
bmp_ss = bmp22_np[:, :, 0:5]
bmpout = Bitmap(bmp_ss, Bitmap.PixelFormat.XYZAW)
bmpout.write('./transients/' + 'bmp_ss.exr')
bmpout.convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True).write('./transients/' + 'bmp_ss' + '.jpg')

exit()
# =============================================================================

# ------------------------------------------------------------
# Set up the animation
STEP_SIZE = 1
N_FRAMES = 360//STEP_SIZE + 1  # // is the floor operator
npixels = resx*resy
N_TIMES = N_FRAMES
# N_SS_CHAN = 3
# N_TRA_CHAN = 3

lesum = np.zeros([N_TIMES, N_BINS])

# Dunno
# fig1, ax1 = plt.subplots()
# ims1 = []
# FPS = 1
# ------------------------------------------------------------

# =============================================================================
# Begin the rotations

# lorigin = Vector3f(0, 0, 0)
# lorigin = Vector3f(20e-2, 0, 20e-2)
# lorigin = Vector3f(0, 0, 20e-2)
# # boresight = Vector3f(0, -1, 0)
# # boresight = Vector3f(1, 0, 20e-2)
# # boresight = Vector3f(0, 1, 20e-2)
# boresight = Vector3f(0, 1, 20e-2)

lorigin = Vector3f(0, 0, 1)
# boresight = Vector3f(0, -1, 0)
# boresight = Vector3f(-1, 0, 20e-2)
boresight = Vector3f(-1, 0, 1)

# (origin=[20e-2, 0, 20e-2],
#                                target=[1, 0, 20e-2],
#                                up=[0, 0, 1])
for i in range(1, N_FRAMES):
    # rotation_cur = Transform4f.rotate(Vector3f(0, 0, 1), i*STEP_SIZE)
    rotation_cur = Transform4f.rotate(Vector3f(0, 0, 1), i*STEP_SIZE)
    new_boresight_c = rotation_cur.transform_vector(boresight)
    new_up_c = rotation_cur.transform_vector(Vector3f(0, 0, 1))
    to_world_cur = Transform4f.look_at(lorigin, new_boresight_c, new_up_c)

    sen = load_dict({
        'type': 'perspective',
        'near_clip': 0.100000,
        'far_clip': 100.0,
        'fov_axis': 'y',
        'fov': 75.0,
        'sampler': sampler,
        'film': film,
        'to_world': to_world_cur,
    })

    emit = load_dict({
        'type': 'rectangle',
        'emitter': {
            'type': 'wignertransmitter',
            'radiance': {
                'type': 'spectrum',
                'value': [(8292683, power), (8717949, power)],
            },
        },
        'to_world': to_world_cur*emit_transform,
    })

    # emit = load_dict({
    #     'type': 'spot',
    #     'cutoff_angle': 25,
    #     'beam_width': 20,
    #     'intensity': {
    #         'type': 'spectrum',
    #         'value': [(8292683, 1000), (8717949, 1000)],
    #     },
    #     'to_world': to_world_cur,
    # })

    scene = load_dict({
        'type': 'scene',
        'integrator': ints,
        'sensor': sen,
        'emitter': emit,
        'so': targ,
        's1': gnd,
        's2': car,
    })

    # scene.integrator().render(scene, scene.sensors()[0])
    scene.integrator().render(scene, sen)

    bmp2 = film.bitmap(raw=True)
    bmp22_np = np.array(bmp2)

    # Section for steady state-------------------------------
    # bmp_ss = bmp22_np[:, :, 0:3]
    # bmpout = Bitmap(bmp_ss, Bitmap.PixelFormat.XYZ)
    # bmpout.convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True)
    # -------------------------------------------------------

    # Section for radar-------------------------------
    for j in range(0, N_BINS):
        idx_lo = 5 + j + 2*j
        idx_hi = 5 + j+3 + 2*j
        # bmp_tra = sum(bmp22_np[0, 0, idx_lo:idx_hi])
        # bmp_tra = bmp22_np[0, 0, idx_lo]
        lesum[i, j] = (sum(sum(sum(bmp22_np[:, :, idx_lo:idx_hi]))))
        # rng_scan[i, j] = np.abs(bmp_tra)
        # rng_scan[i, j] = bmp_tra

# =====================================================================
# Rendering done

# ---------------------------------------------------------
# Outputs
lesum = 10*np.log10(np.abs(lesum) + np.finfo(float).eps)
rs = np.arange(0, DR*N_BINS, DR)
rs = rs/2

# thetas = np.deg2rad(np.arange(0, N_FRAMES*STEP_SIZE, STEP_SIZE))
thetas = np.deg2rad(np.arange(0+180, N_FRAMES*STEP_SIZE+180, STEP_SIZE))
rr, tt = np.meshgrid(rs, thetas)

# minc = lesum[(lesum > -140)].min()
minc = lesum[(lesum > -200)].min()
maxc = lesum.max()

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
pol = ax.contourf(tt, rr, lesum, cmap=cm.coolwarm, vmin=minc, vmax=maxc)
fig.colorbar(pol)

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(rr, np.rad2deg(tt), lesum
#                        , cmap=cm.coolwarm, vmin=minc, vmax=maxc)
# fig.colorbar(surf)


# xLabel = ax.set_xlabel('Range [m]')
# yLabel = ax.set_ylabel('Angle [deg]')
# xLabel = ax.set_xlabel('X [m]')
# yLabel = ax.set_ylabel('Y [m]')
# zLabel = ax.set_zlabel('Return Value [dBW]')
title = ax.set_title('Signal Return Ground + Car + Plate')
plt.show()
