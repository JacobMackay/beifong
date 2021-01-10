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

# Absolute or relative path to the XML file
filename = './trans_image.xml'
# filename = './trans_image_rad.xml'
name = 'trans_image'

# Add the scene directory to the FileResolver's search path
Thread.thread().file_resolver().append(os.path.dirname(filename))

# Load the actual scene
scene = load_file(filename)

print(scene)

nbins = 100

# Call the scene's integrator to render the loaded scene
scene.integrator().render(scene, scene.sensors()[0])

# After rendering, the rendered data is stored in the film
film = scene.sensors()[0].film()


# Write out rendering as high dynamic range OpenEXR file
# film.set_destination_file('./' + name + '.exr')
# film.develop()

bmp2 = film.bitmap(raw=True)
bmp22_np = np.array(bmp2)
print(bmp22_np.shape)
# print(bmp22_np.tolist())

npx, npy = film.crop_size()
npixels = npx*npy
ntimes = 1

# Section for images ----------------------------------------------------------
bmp_ss = bmp22_np[:, :, 0:5]
bmpout = Bitmap(bmp_ss, Bitmap.PixelFormat.XYZAW)
bmpout.write('./transients/' + 'bmp_ss.exr')
bmpout.convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True).write('./transients/' + 'bmp_ss' + '.jpg')

# for i in range(0,nbins):
#     idx_lo = 5 + i + 2*i
#     idx_hi = 5 + i+3 + 2*i
#     print(idx_lo, idx_hi)
#     bmp_tra = bmp22_np[:, :, idx_lo:idx_hi]
#     bmp_tra = Bitmap(bmp_tra, Bitmap.PixelFormat.XYZ)
#     bmp_tra.convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True).write('./transients/' + 'bmp_tra_' + str(i) + '.jpg')

lesum = np.zeros([nbins])

for i in range(0,nbins):
    idx_lo = 5 + i + 2*i
    idx_hi = 5 + i+3 + 2*i
    print(idx_lo, idx_hi)
    bmp_tra = bmp22_np[:, :, idx_lo:idx_hi]
    # print(sum(bmp22_np[i, idx_lo:idx_hi, k]))
    # print(sum(sum(sum(bmp22_np[:, :, idx_lo:idx_hi]))))
    lesum[i] = (sum(sum(sum(bmp22_np[:, :, idx_lo:idx_hi]))))
    # bmp_tra = Bitmap(bmp_tra, Bitmap.PixelFormat.XYZ)
    # bmp_tra.convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True).write('./transients/' + 'bmp_tra_' + str(i) + '.jpg')

lesum = 10*np.log10(np.abs(lesum + np.finfo(float).eps))

dr = 0.1 / 2
rs = np.arange(0, dr*nbins, dr)

plt.plot(rs,lesum)
plt.xlabel("Range [m]")
plt.ylabel("Return Value [dB]")
plt.title('Signal Return')
plt.show()
