from GRS.targets.target import *
from GRS.targets.config_tool import *

class FromStructure(Target):
    def __init__(self,t,fingerprints,typ):
        super().__init__(t,fingerprints,typ)
        self.target_type = typ
        self.fingerprints = fingerprints
        self.t = t
        s = self.get_structure()
        self.set_structure(s)

    def get_structure(self):
        if self.target_typ == 'file':
            s = read(self.t)
        elif self.target_typ == 'ase':
            s = self.t
        elif self.target_typ == 'json':
            fsats = fsnap_atoms()
            fsats.read_json(self.t)
            s = fsats.get_ASE()
        else:
            raise TypeError("type %s not implemented" % self.target_typ)
        return s

    def set_structure(self,s):
        struct_elems = sorted(list(set([elemi for elemi in s.symbols]))) 
        fingerprint_space_elems = self.fingerprints.fingerprint_settings['elements']
        all_contained = all([selem in fingerprint_space_elems for selem in struct_elems])
        assert all_contained, "One or more of the elements in your structure is not in yaur fingerprint space. Rebuild your fingerprint space to include these new elements before continuing"
        self.target_structure = s

    def set_target_fingerprints(self):


