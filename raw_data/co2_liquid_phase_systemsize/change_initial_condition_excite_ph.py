import numpy as np
import sys
import xml.etree.ElementTree as ET
import json as js
from collections.abc import Mapping, Sequence

class ScientificNotationEncoder(js.JSONEncoder):
    def iterencode(self, o, _one_shot=False):
        if isinstance(o, float):
            return "{:e}".format(o)
        elif isinstance(o, Mapping):
            return "{{{}}}".format(', '.join('"{}" : {}'.format(str(ok), self.iterencode(ov))
                                             for ok, ov in o.items()))
        elif isinstance(o, Sequence) and not isinstance(o, str):
            return "[{}]".format(', '.join(map(self.iterencode, o)))
        return ', '.join(super().iterencode(o, _one_shot))


def operate_one_file(xml_path, store_folder="./", idx_displacement=0, r_displacement=0.0):
    # Open original file
    tree = ET.parse(xml_path)
    # Read this file
    root = tree.getroot()
    a = root.find("system")
    suba = a.find("beads")
    # load mass and momenta
    momenta = suba.find("p")
    mass = suba.find("m")
    positions = suba.find("q")
    momenta_value = momenta.text
    mass_value = mass.text
    positions_value = positions.text
    momenta_value = js.loads(momenta_value)
    mass_value = js.loads(mass_value)
    positions_value = js.loads(positions_value)

    # change the positions array, displace it
    positions_value[idx_displacement] = r_displacement
    positions.text = js.dumps(positions_value, cls=ScientificNotationEncoder)

    # Write back to file
    pure_path = xml_path.split('/')[-1].split('_')[-1]
    tree.write(store_folder+"/init_" + pure_path)



if __name__ == "__main__":

    path, n_increase_molecule, n_increase_ph = sys.argv[1], sys.argv[2], sys.argv[3]
    #print("working under %s with n_increase_molecule = %s and n_increase_ph = %s" %(path, n_increase_molecule, n_increase_ph))
    n_increase_molecule = int(n_increase_molecule)
    n_increase_ph = int(n_increase_ph)
    n_co2 = 36 * 36 * n_increase_molecule
    idx_excite = 3 * (3 * n_co2 + 12 * n_increase_ph - 1) + 1  # excite the y-axis of the l=12 cavity mode
    #print("will excite No.%d idx" %idx_excite)
    r_displacement = 100.0 * 1.8897259886 * n_increase_molecule**0.5
    #print("displaced to %.5E" %r_displacement)
    #print("will displace the photon coordinate to %.3E Angstrom" %(r_displacement/1.8897259886))

    #print("Dealing with %s" %path)
    operate_one_file(xml_path=path, idx_displacement=idx_excite, r_displacement=r_displacement)