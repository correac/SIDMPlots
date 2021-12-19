import swiftsimio as sw
from velociraptor import load
import unyt
import numpy as np

folder = "/scratch/tangosidm/Tango_sims/L025N752/DMONLY/SigmaVelDep20Isotropic"
#folder = "/Users/camila/SimulationData/mahti/L025N752/DMONLY/SigmaVelDep20Isotropic"

for snap in range(37):
    filename = folder+"/snapshot_%04i.hdf5"%snap
    print(filename)
    # Now load the snapshot with this mask
    data = sw.load(filename)
    h = data.dark_matter.sidm_search_radius.value
    select = h != 0.03
    print(np.max(h[select]))
#vr_file = folder+"/subhalo_0036.properties"

# mask = sw.mask(filename)
# # The full metadata object is available from within the mask
# size = unyt.unyt_array([1, 1, 1], 'Mpc')
#
# boxsize = mask.metadata.boxsize
# # load_region is a 3x2 list [[left, right], [bottom, top], [front, back]]
# #load_region = [[0.0 * b, 0.5 * b] for b in size]
#
# catalogue = load(vr_file)
# x = catalogue.positions.xcminpot[0].to("Mpc").value
# y = catalogue.positions.ycminpot[0].to("Mpc").value
# z = catalogue.positions.zcminpot[0].to("Mpc").value
# origin = unyt.unyt_array([24.7, 24.8, 24.9], 'Mpc')
#
# region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(size, origin)]
#
# # Constrain the mask
# mask.constrain_spatial(region)
#
# # Now load the snapshot with this mask
# data = sw.load(filename, mask=mask)
#
# print(np.max(data.dark_matter.coordinates))
