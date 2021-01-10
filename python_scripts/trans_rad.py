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
# filename = './trans_rad.xml'
filename = './trans_image_rad.xml'
name = 'trans_rad'

# Add the scene directory to the FileResolver's search path
Thread.thread().file_resolver().append(os.path.dirname(filename))

# Load the actual scene
scene = load_file(filename)

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
nbins = 50
ntimes = 1

# Section for radar ------------------------------------------------------------
bmp_tra_np = np.empty([npixels, nbins, ntimes])

# render scene at a short burst of times and put result in corresponding bins
bmp22_np_re = np.reshape(bmp22_np, [npixels, nbins*3+5, ntimes])
spp = scene.sensors()[0].sampler().sample_count()
for i in range(0,npixels):
    for k in range(0,ntimes):
        for j in range(0,nbins):
            idx_lo = 5 + j + 2*j
            idx_hi = 5 + j+3 + 2*j
            bmp_tra = sum(bmp22_np_re[i,idx_lo:idx_hi,k])
            bmp_tra_np[i,j,k] = bmp_tra/spp + np.finfo(float).eps

bmp_tra_np = 10*np.log10(bmp_tra_np)

dt = 0.5e-9
ts = np.arange(0,dt*50,dt)
rs = ts*3e8


plt.plot(rs,bmp_tra_np[0,:,0])
plt.xlabel("Range [m]")
plt.ylabel("Return Value [dB]")
plt.title('Signal Return')
plt.show()
