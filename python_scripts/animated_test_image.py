"""File for testing how to create an animated image using mitsuba/beifong."""
# import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.animation as animation
matplotlib.use("Agg")
# plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
import mitsuba
mitsuba.set_variant('scalar_spectral')
# import enoki as ek
# import cv2
# from mitsuba.core.xml import load_file, load_dict, load_string
from mitsuba.core.xml import load_file, load_dict
from mitsuba.core import Vector3f, Transform4f
from mitsuba.core import Bitmap, Struct, ScalarTransform4f

# Set the desired mitsuba variant
# mitsuba.set_variant('scalar_rgb')
# mitsuba.set_variant('scalar_mono')
# mitsuba.set_variant('scalar_spectral')
# mitsuba.set_variant('gpu_spectral')
# mitsuba.set_variant('gpu_rgb')
# mitsuba.set_variant('packet_spectral')

# from mitsuba.core import Bitmap, Struct, Thread, *
# from mitsuba.core import *
# from mitsuba.render import *
# from mitsuba.core.xml import *
# from mitsuba.python.util import traverse

# from matplotlib import image as img

# I want to be able to read from a file, AND add or update elements

# Absolute or relative path to the XML file
FILENAME = './animated_trans_image.xml'
NAME = 'animated_trans_image'

targ_size = Transform4f.scale([0.1, 0.1, 1])
targ_loc = Transform4f.look_at(origin=[0, -4, 0],
                               target=[0, 0, 0],
                               up=[0, 0, 1])
targ_transform = targ_loc*targ_size
gnd_size = Transform4f.scale([20, 20, 1])
gnd_loc = Transform4f.look_at(origin=[0, 0, -0.5],
                              target=[0, 0, 0.5],
                              up=[0, 1, 0])
gnd_transform = gnd_loc*gnd_size

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
    "type": "time",
    "integrator": {
        "type": "pathtime",
    }
})
emit2 = load_dict({
    "type": "spot",
    "cutoff_angle": 25,
    "beam_width": 20,
    "intensity": {
        "type": "spectrum",
        "value": 100,
    }
})
film2 = load_dict({
    "type": "hdrfilm",
    "rfilter": {"type": "box"},
    "width": 400,
    "height": 400,
})
sampler2 = load_dict({
    "type": "independent",
    "sample_count": 4,
    "id": "sampler2",
})
sen2 = load_dict({
    "type": "perspective",
    "near_clip": 0.100000,
    "far_clip": 100.0,
    "fov_axis": "x",
    "fov": 45.0,
    "sampler": sampler2,
    "film": film2,
})


STEP_SIZE = 5
N_FRAMES = 360//STEP_SIZE  # // is the floor operator

npx, npy = film2.crop_size()
npixels = npx*npy
N_BINS = 50
N_TIMES = N_FRAMES
N_SS_CHAN = 3
N_TRA_CHAN = 3
# bmp_tra_np = np.empty([npx, npy, N_TRA_CHAN, N_BINS, N_TIMES])
# bmp_ss_np = np.empty([npx, npy, N_SS_CHAN, N_TIMES])
# Initialise frame storage

# print(bmp_tra_np.shape, bmp_ss_np.shape)

# These arrays are too big...

fig1, ax1 = plt.subplots()
ims1 = []
# fig2, ax2 = plt.subplots()
# ims2 = []
FPS = 1

lorigin = Vector3f(0, 0, 0)
boresight = Vector3f(0, -1, 0)
for i in range(1, N_FRAMES):
    rotation_cur = Transform4f.rotate(Vector3f(0, 0, 1), i*STEP_SIZE)
    new_boresight_c = rotation_cur.transform_vector(boresight)
    new_up_c = rotation_cur.transform_vector(Vector3f(0, 0, 1))
    to_world_cur = Transform4f.look_at(lorigin, new_boresight_c, new_up_c)

    sen2 = load_dict({
        "type": "perspective",
        "near_clip": 0.100000,
        "far_clip": 100.0,
        "fov_axis": "x",
        "fov": 45.0,
        "to_world": to_world_cur,
        "sampler": sampler2,
        "film": film2,
    })

    emit2 = load_dict({
        "type": "spot",
        "intensity": {
            "type": "spectrum",
            "value": 100,
        },
        "cutoff_angle": 25,
        "beam_width": 20,
        "to_world": to_world_cur,
    })

    scene2 = load_dict({
        "type": "scene",
        "integrator": ints,
        "sensor": sen2,
        "emitter": emit2,
        "so": targ,
        "s1": gnd,
    })

    scene2.integrator().render(scene2, sen2)

    bmp2 = film2.bitmap(raw=True)
    bmp22_np = np.array(bmp2)
    bmp_ss = bmp22_np[:, :, 0:3]
    bmpout = Bitmap(bmp_ss, Bitmap.PixelFormat.XYZ)
    bmpout.convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True)
    # bmpout = Bitmap(bmp_ss, Bitmap.PixelFormat.XYZAW)
    # bmpout.write('./animated/frame_ss_0' + str(i) +'_' + name + '.exr')
    # bmpout.convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8,
    #     srgb_gamma=True).write('./animated/frame_ss_0' + str(i) +'_' + name +
    #     '.jpg')

    # Append new frame
    # bmp_ss_np[:, :, :, i] = np.array(bmpout)
    im1 = ax1.imshow(np.array(bmpout), animated=True)
    ims1.append([im1])

    for j in range(0, N_BINS):
        idx_lo = 5 + j + 2*j
        idx_hi = 5 + j+3 + 2*j
        bmp_tra = bmp22_np[:, :, idx_lo:idx_hi]
        bmpout2 = Bitmap(bmp_tra, Bitmap.PixelFormat.XYZ)
        bmpout2.convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8,
                        srgb_gamma=True)
        # im2 = ax2.imshow(np.array(bmpout2), animated=True)
        # ims2.append([im2])
        # print(bmp_tra.shape)
        # print(bmp_tra)
        # bmp_tra = Bitmap(bmp_tra, Bitmap.PixelFormat.Y)
        # bmp_tra = Bitmap(bmp_tra, Bitmap.PixelFormat.XYZ)
        # bmp_tra.convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8,
        #     srgb_gamma=True).write('./animated/transients/frame_ss_0' +str(i)
        #     + '_' + 'tra_0' + str(j) + '.jpg')

ani = animation.ArtistAnimation(fig1, ims1, interval=60//FPS)
ani.save('./animated/ss.gif', writer='imagemagick')
# ani2 = animation.ArtistAnimation(fig2, ims2, interval=60//FPS)
# ani2.save('./animated/tra.gif', writer='imagemagick')

# Section for radar -----------------------------------------------------------
# bmp_tra_np = np.empty([npixels, nbins, ntimes])
#
# # render scene at a short burst of times and put result in corresponding bins
# bmp22_np_re = np.reshape(bmp22_np, [npixels, nbins*3+5, ntimes])
# spp = scene.sensors()[0].sampler().sample_count()
# for i in range(0,npixels):
#     for k in range(0,ntimes):
#         for j in range(0,nbins):
#             idx_lo = 5 + j + 2*j
#             idx_hi = 5 + j+3 + 2*j
#             bmp_tra = sum(bmp22_np_re[i,idx_lo:idx_hi,k]);
#             bmp_tra_np[i,j,k] = bmp_tra/spp + np.finfo(float).eps
#
# bmp_tra_np = 10*np.log10(bmp_tra_np)
#
# dt = 0.5e-9;
# ts = np.arange(0,dt*50,dt)
# rs = ts*3e8
#
#
# plt.plot(rs,bmp_tra_np[0,:,0])
# plt.xlabel("Range [m]")
# plt.ylabel("Return Value [dB]")
# plt.title('Signal Return')
# plt.show()
