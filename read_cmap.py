import numpy as np
from bitarray import bitarray
import matplotlib.pyplot as plt
from libContactMapper import contact_mapper
import time
cm = contact_mapper()

# save_path ='/home/soliareofastora/genomics_data/database/cmap/6bz7.bin'
save_path = '/home/soliareofastora/genomics_data/database/cmap/1_468_6fqx.3.A_60806a71dafdd895ba140273.bin'



start = time.time()
cpp_map = cm.load_cmap(save_path)
print((time.time() - start))


start = time.time()
loaded_bits = bitarray(endian="little")
with open(save_path, 'rb') as f:
    loaded_bits.fromfile(f)

n = int((1+np.sqrt(8*len(loaded_bits)+1))/2)
py_map = np.zeros((n, n), dtype=np.bool)
np.fill_diagonal(py_map, 1)
for k in range(int((n*(n-1) )/ 2)):
    if loaded_bits[k]:
        i = int(n - 2 - int(np.sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5))
        j = int(k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2)
        py_map[i][j] = 1
        py_map[j][i] = 1
print((time.time() - start))


plt.imshow(py_map, vmin=0, vmax=1)
plt.title("py")
plt.show()

plt.imshow(cpp_map, vmin=0, vmax=1)
plt.title("cpp")
plt.show()

diff = py_map ^ cpp_map
if np.sum(diff) > 0:
    print("difference:", np.sum(diff))
    plt.imshow(diff, vmin=0, vmax=1)
    plt.title("difference")
    plt.show()





