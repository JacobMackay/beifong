"""Testing for getting a radar sweep."""
# import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from matplotlib import cm
import mitsuba
mitsuba.set_variant('scalar_spectral')
# mitsuba.set_variant('scalar_rgb')
# mitsuba.set_variant('packet_rgb')
from mitsuba.core.xml import load_dict, load_string
from mitsuba.core import Vector3f, Transform4f


targ_size = Transform4f.scale([1, 1, 1])
# targ_loc = Transform4f.look_at(origin=[0, -4, 0],
#                                target=[0, 0, 0],
#                                up=[0, 0, 1])
targ_loc = Transform4f.look_at(origin=[0, 4, 0],
                               target=[0, 0, 0],
                               up=[0, 0, 1])
targ_transform = targ_loc*targ_size
gnd_size = Transform4f.scale([20, 20, 1])
gnd_loc = Transform4f.look_at(origin=[0, 0, -0.5],
                              target=[0, 0, 0.5],
                              up=[0, 1, 0])
gnd_transform = gnd_loc*gnd_size
sen_size = Transform4f.scale([0.5, 0.5, 1])
sen_loc = Transform4f.look_at(origin=[0, -1, 0],
                              target=[0, 0, 0],
                              up=[0, 0, 1])
sen_transform = sen_loc*sen_size
txa_size = Transform4f.scale([0.5, 0.5, 1])
txa_loc = Transform4f.look_at(origin=[0, -1, 0],
                              target=[0, 0, 0],
                              up=[0, 0, 1])
txa_transform = txa_loc*txa_size

SPP = 1000

BINS = 50
DR = 0.2

max_range = BINS*DR + DR
min_range = DR
# remember that the camera also samples 4 wavelengths


# bsdfs_d = {
#     "type": "twosided",
#     "id": "material",
#     "bsdf": {
#         "type": "diffuse",
#         "reflectance": {
#             "type": "spectrum",
#             "value": 1,
#         },
#     },
# }
# targ_d = {
#     "type": "rectangle",
#     "id": "target",
#     "to_world": targ_transform,
#     "bsdf": bsdfs_d,
# }
# gnd_d = {
#     "type": "rectangle",
#     "id": "gnd",
#     "to_world": gnd_transform,
#     "bsdf": bsdfs_d,
# }
# ints_d = {
#     "type": "time",
#     "integrator": {
#         "type": "pathtime",
#     }
# }
# emit_s_d = {
#     "type": "spot",
#     "cutoff_angle": 25,
#     "beam_width": 20,
#     "intensity": {
#         "type": "spectrum",
#         "value": 100,
#     }
# }
# film_d = {
#     "type": "hdrfilm",
#     "rfilter": {"type": "box"},
#     "width": 1,
#     "height": 1,
# }
# sampler_d = {
#     "type": "independent",
#     "sample_count": SPP,
# }
# # lolsen = load_dict({
# #     "type": "fluxmeter",
# #     "id": "rx",
# #     # "type": "radiancemeter",
# #     "sampler": sampler2,
# #     "film": film2,
# #     # "shape": {
# #     #     "type": "rectangle",
# #     #     "to_world": sen_transform,
# #     # },
# # })
# # print(lolsen)
# # sen2 = load_dict({
# #     "type": "rectangle",
# #     "id": "receiveAntenna",
# #     "to_world": sen_transform,
# #     "sensor": [],
# # })
# # print(sen2)
# # sen2["sensor"] = lolsen
# # print(sen2)
# sen_d = {
#     "type": "perspective",
#     "near_clip": 0.100000,
#     "far_clip": 100.0,
#     "fov_axis": "x",
#     "fov": 45.0,
#     "sampler": sampler_d,
#     "film": film_d,
# }
# scene_d = {
#     "type": "scene",
#     "integrator": ints_d,
#     "sensor": sen_d,
#     "emitter": emit_s_d,
#     "so": targ_d,
#     "s1": gnd_d,
# }
#
# film = load_dict(film_d)
# sen = load_dict(sen_d)
# scene = load_dict(scene_d)

bsdfs = load_dict({
    "type": "twosided",
    "id": "material",
    "bsdf": {
        "type": "diffuse",
        "reflectance": {
            "type": "spectrum",
            "value": 1,
        },
    },
})
targ = load_dict({
    "type": "rectangle",
    "id": "target",
    "to_world": targ_transform,
    "bsdf": bsdfs,
})
gnd = load_dict({
    "type": "rectangle",
    "id": "gnd",
    "to_world": gnd_transform,
    "bsdf": bsdfs,
})
ints = load_dict({
    "type": "range",
    "integrator": {
        "type": "pathlength",
    },
    "dr": DR,
    "bins": BINS,
})
emit_s = load_dict({
    "type": "spot",
    "cutoff_angle": 25,
    "beam_width": 20,
    "intensity": {
        "type": "spectrum",
        # "value": 10e-3*0.5e-6/SPP,   # W/sr * t
        "value": 1000,   # W/sr * t
    }
})
emit_r = load_dict({
    "type": "rectangle",
    "id": "txa",
    "to_world": txa_transform,
    "emitter": {
        "type": "area",
        "radiance": {
            "type": "spectrum",
            # "value": 10e-3*0.5e-6/SPP,   # W/sr * t
            "value": 1000,   # W/sr * t
        },
    },
})
film = load_dict({
    "type": "hdrfilm",
    "rfilter": {"type": "box"},
    "width": 1,
    "height": 1,
})
sampler = load_dict({
    "type": "independent",
    "sample_count": SPP,
})
# sen_r = load_dict({
#     "type": "rectangle",
#     "id": "rxa",
#     "to_world": sen_transform,
#     "sensor": {
#         "type": "fluxmeter",
#         # "type": "irradiancemeter",
#         "sampler": sampler,
#         "film": film,
#     },
# })
# <matrix value="0 -0.53 0 -1.79 0.92 0 0 8.03 0 0 0.53 0 0 0 0 1"/>
# sen_r_str = ("<shape version='2.0.0' type='rectangle' id='txa'> "
#              "<transform name='to_world'> "
#              "<scale x='0.05' y='0.05'/> "
#              "<lookat origin='0, 0, 0' "
#              "target='0, -1, 0' "
#              "up='0, 0, 1'/> "
#              "</transform> "
#              "<emitter type='area'> "
#              "<spectrum name='radiance' value='1.0'/> "
#              "</emitter> "
#              "</shape>")
# sen_r = load_string(sen_r_str)
# print(sen_r)
sen_s = load_dict({
    "type": "perspective",
    "near_clip": min_range,
    "far_clip": max_range,
    "fov_axis": "x",
    "fov": 45,
    "sampler": sampler,
    "film": film,
})
scene = load_dict({
    "type": "scene",
    "integrator": ints,
    "sensor": sen_s,
    # "sensor": sen_r,
    # "emitter": emit_s,
    "emitter": emit_r,
    "so": targ,
    "s1": gnd,
})

# Next TODO:
# Get 'shape' emit_ster
# Get 'shape fluxmeter'
# Figure out units and return values...
# Get angular sampling working for bsdfs
# Start with tx and rx beamwidths, but leave the environment alone
# Paper: Simple compare what we've done on this trial to simulated reults
# Renderer that accounts for beamwidth
# Include doppler: each object would also have to have velocity properties.
# Lets not deal with this for now.
# Make intergrator that includes range

# To get a paper ready by october...
# I've made a radar sim based off mitsuba2
# it handles beam shape*** is this an interesting algorithm?
# i've added an algorithm for long range
# Note it doesn't do multipath (yet ;) )
# ground modelled as rough plastic/dielectric
# reflector modelled as conductor
# need an ior/wavelength profile for materials...
# unless I do diffuse. Conductor has an optin of none which is perfectly
# specular
# Probably use rough/smooth conductors for metallic objects, and rough
# dielectric for the rest
# Choose one of the txrx..even tho it beamforms..
# Beam forming radar ray tracer
# Simple wdf only on tx and for multipath.
# Specify tx as shapes in locs. Add a phase offset
# Radar beamforming simulation with the wdf.
# Updated experiment plan: reflector at small height, move toward/away and get
#  fade/glint

# implement a narrowband variant: Specify wavelength ranges
# don't need to, this is being handled! Now wait.

# Beamforming and multipath for groundbased imaging radar simulations.
# So I need wdf + multipath
# Compare diffuse tx array, angular tx array and wdf tx array.
# Visualise for comparison a radar 'image' and then experiments.

STEP_SIZE = 5
N_FRAMES = 360//STEP_SIZE + 1

npx, npy = film.crop_size()
# npx, npy = film_d["width"], film_d["height"]
npixels = npx*npy
N_BINS = BINS
N_TIMES = N_FRAMES
N_SS_CHAN = 3
N_TRA_CHAN = 3

rng_scan = np.empty([N_TIMES, N_BINS])

FPS = 1

lorigin = Vector3f(0, 0, 0)
boresight = Vector3f(0, -1, 0)
# for i in range(1, N_FRAMES):
for i in range(1, N_FRAMES):
    rotation_cur = Transform4f.rotate(Vector3f(0, 0, 1), i*STEP_SIZE)
    new_boresight_c = rotation_cur.transform_vector(boresight)
    new_up_c = rotation_cur.transform_vector(Vector3f(0, 0, 1))
    to_world_cur = Transform4f.look_at(lorigin, new_boresight_c, new_up_c)

    # sen2 = load_dict({
    #     "type": "rectangle",
    #     "id": "receiveAntenna",
    #     "to_world": to_world_cur*sen_size,
    #     "sensor": lolsen,
    # })

    # sen.update({"to_world": to_world_cur})
    sen = load_dict({
        "type": "perspective",
        "near_clip": min_range,
        "far_clip": max_range,
        "fov_axis": "x",
        "fov": 45,
        "sampler": sampler,
        "film": film,
        "to_world": to_world_cur,
    })

    # emit_s = load_dict({
    #     "type": "spot",
    #     "intensity": {
    #         "type": "spectrum",
    #         # "value": 10e-3*0.5e-6/SPP,
    #         "value": 100,
    #     },
    #     "cutoff_angle": 25,
    #     "beam_width": 20,
    #     "to_world": to_world_cur,
    # })
    emit_r = load_dict({
        "type": "rectangle",
        "id": "txa",
        "to_world": to_world_cur*txa_size,
        "emitter": {
            "type": "area",
            "radiance": {
                "type": "spectrum",
                # "value": 10e-3*0.5e-6/SPP,   # W/sr * t
                "value": 100,   # W/sr * t
            },
        },
    })
    # emit_s.update({"to_world": to_world_cur})
    # emit_s2["to_world"] = to_world_cur
    # sen = load_dict(sen_d)
    # scene = load_dict(scene)

    scene = load_dict({
        "type": "scene",
        "integrator": ints,
        "sensor": sen,
        # "emitter": emit_s,
        "emitter": emit_r,
        "so": targ,
        "s1": gnd,
    })

    # Could make alternate render interface which returns values rather than
    # storing them on a film...nahhhh
    scene.integrator().render(scene, sen)

    bmp2 = film.bitmap(raw=True)
    bmp22_np = np.array(bmp2)

    for j in range(0, N_BINS):
        idx_lo = 5 + j + 2*j
        idx_hi = 5 + j+3 + 2*j
        # bmp_tra = sum(bmp22_np[0, 0, idx_lo:idx_hi])
        bmp_tra = bmp22_np[0, 0, idx_lo]
        # rng_scan[i, j] = np.abs(bmp_tra)
        rng_scan[i, j] = bmp_tra

# dt = 0.5e-9
dr = DR
# rng_scan = 10*np.log10(rng_scan/dt + 100*np.finfo(float).eps)
# print(rng_scan)
rng_scan = 10*np.log10(abs(rng_scan)/SPP + 100*np.finfo(float).eps)
# ts = np.arange(0, dt*50, dt)
# rs = ts*3e8
rs = np.arange(0, dr*BINS, dr)
thetas = np.deg2rad(np.arange(0, N_FRAMES*STEP_SIZE, STEP_SIZE))
rr, tt = np.meshgrid(rs, thetas)

minc = rng_scan[(rng_scan > -140)].min()
maxc = rng_scan.max()



# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# ax.contourf(tt, rr, rng_scan)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(rr, np.rad2deg(tt), rng_scan
                       , cmap=cm.coolwarm, vmin=minc, vmax=maxc)
fig.colorbar(surf)


xLabel = ax.set_xlabel("Range [m]")
yLabel = ax.set_ylabel("Angle [deg]")
# xLabel = ax.set_xlabel("X [m]")
# yLabel = ax.set_ylabel("Y [m]")
zLabel = ax.set_zlabel("Return Value [dBW]")
title = ax.set_title("Signal Return Diffuse Ground + Plate")
plt.show()

# Power back should be pretty trivial.
# I sample a bunch of rays with starting time (pulse). This should also give us
# pulse width effects. The length of the pulse*P_tx
# maybe i should consider energy conservation instead. pulse time*peak power =E
# So then each ray has an energy associated with it. Integrate over a
# pixel(space), and hemisphere (direction),
# So total energy is peak power*transmit time. Each ray leaves with total
# energy/total samples. We send these through the scene and collect energy at
# the end. The power in each bin is the energy collected in that bin/bin time.

# Still not sure. Should be correct if transmit energy
