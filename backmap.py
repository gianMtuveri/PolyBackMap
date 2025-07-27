import json
import numpy as np
from itertools import product
import string
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation as R

def align_fragment_to_beads(fragment_coords, src_vec, tgt_vec):
    src_vec = src_vec / np.linalg.norm(src_vec)
    tgt_vec = tgt_vec / np.linalg.norm(tgt_vec)
    rot, _ = R.align_vectors([tgt_vec], [src_vec])
    return rot.apply(fragment_coords)

def load_mapping(mapping_file):
    with open(mapping_file) as f:
        return json.load(f)

def center_fragment_on_bead(frag_mol, cg_pos, cg_vec=None):
    conf = frag_mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(frag_mol.GetNumAtoms())])
    com = coords.mean(axis=0)
    centered = coords - com
    if cg_vec is not None:
        aligned = align_fragment_to_beads(centered, np.array([1.0, 0.0, 0.0]), cg_vec)
    else:
        aligned = centered
    translated = aligned + cg_pos
    for i in range(frag_mol.GetNumAtoms()):
        conf.SetAtomPosition(i, translated[i])
    return frag_mol

def write_gro(atom_list, box, output_file):
    with open(output_file, 'w') as f:
        f.write("Backmapped structure\n")
        f.write(f"{len(atom_list)}\n")
        for i, atom in enumerate(atom_list):
            f.write(f"{atom['resid']:<5}{atom['resname']:<5}{atom['atom_name']:<5}{i+1:5d}"
                    f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}\n")
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")

def write_itp(atom_list, bond_list, output_file):
    with open(output_file, 'w') as f:
        f.write("[ atoms ]\n")
        f.write("; nr type resnr resid atom cgnr charge mass\n")
        for i, atom in enumerate(atom_list):
            f.write(f"{i+1:5d} {atom['atom_name']:>4} {atom['resid']:4d} {atom['resname']:>4} "
                    f"{atom['atom_name']:>4} {i+1:5d} 0.000 12.011\n")
        f.write("\n[ bonds ]\n")
        f.write("; ai aj funct length force_constant\n")
        for bond in bond_list:
            f.write(f"{bond[0]:5d} {bond[1]:5d} 1 0.15 1000\n")

def generate_chain_ids(n):
    upper = list(string.ascii_uppercase)
    lower = list(string.ascii_lowercase)
    all_ids = upper + lower
    if n > len(all_ids):
        raise ValueError(f"Too many chains: {n} requested, but only {len(all_ids)} unique single-letter chain IDs available.")
    return all_ids[:n]

def write_pdb_with_ter_and_conect(atom_list, bond_list, output_file):
    with open(output_file, 'w') as f:
        prev_chain = None
        for i, atom in enumerate(atom_list):
            atom_index = i + 1
            chain_id = atom.get("chain_id", "A")[:2]
            if prev_chain is not None and chain_id != prev_chain:
                f.write(f"TER   {atom_index:5d}      {atom['resname']:<4} {atom['resid']:4d}\n")
            f.write(
                f"ATOM  {atom_index:5d} {atom['atom_name']:<4} {atom['resname']:<3} {chain_id:<2}{atom['resid']:4d}   "
                f"{atom['x']*10:8.3f}{atom['y']*10:8.3f}{atom['z']*10:8.3f}  1.00  0.00\n")
            prev_chain = chain_id
        f.write(f"TER   {len(atom_list)+1:5d}      {atom_list[-1]['resname']:<4} {atom_list[-1]['resid']:4d}\n")
        for bond in bond_list:
            ai, aj = bond
            f.write(f"CONECT{ai:5d}{aj:5d}\n")
        f.write("END\n")

def backmap_to_gromacs(cg_gro_file, mapping_file, output_gro, output_itp, output_pdb):
    u = mda.Universe(cg_gro_file)
    cg_atoms = u.atoms
    mapping = load_mapping(mapping_file)
    atom_list = []
    bond_list = []
    current_resid = 1
    atom_offset = 0
    atoms_in_chain = 0
    atoms_per_chain = 218
    num_chains = 46
    chain_ids = generate_chain_ids(num_chains)
    current_chain = 0

    for bead in cg_atoms:
        bead_type = bead.name
        cg_pos = bead.position
        if bead_type not in mapping:
            raise ValueError(f"Bead {bead_type} not in mapping.")

        info = mapping[bead_type]
        frag = Chem.MolFromMol2File(info["fragment_file"], removeHs=False)
        if frag is None:
            raise ValueError(f"Could not read fragment {info['fragment_file']}")

        AllChem.EmbedMolecule(frag, randomSeed=42)

        next_idx = bead.index + 1
        if next_idx < len(cg_atoms):
            cg_vec = u.atoms[next_idx].position - bead.position
        else:
            cg_vec = np.array([1.0, 0.0, 0.0])

        resname = info["resname"]
        if resname.upper() == "PLA":
            cg_vec = -cg_vec
            atom_names = info["atom_names"]
            idx_o1 = atom_names.index("C1")
            idx_c2 = atom_names.index("C2")
            conf = frag.GetConformer()
            pos_o1 = np.array(conf.GetAtomPosition(idx_o1))
            pos_c2 = np.array(conf.GetAtomPosition(idx_c2))
            frag_vec = pos_c2 - pos_o1
            frag = center_fragment_on_bead(frag, cg_pos, cg_vec=frag_vec)
        else:
            frag = center_fragment_on_bead(frag, cg_pos, cg_vec)

        conf = frag.GetConformer()
        for j, atom_name in enumerate(info["atom_names"]):
            pos = conf.GetAtomPosition(j)
            atom_list.append({
                "resid": current_resid,
                "resname": info["resname"],
                "atom_name": atom_name,
                "x": pos.x / 10.0,
                "y": pos.y / 10.0,
                "z": pos.z / 10.0,
                "chain_id": chain_ids[current_chain],
            })

        for bond in frag.GetBonds():
            ai = bond.GetBeginAtomIdx() + 1 + atom_offset
            aj = bond.GetEndAtomIdx() + 1 + atom_offset
            bond_list.append((ai, aj))

        atom_offset += frag.GetNumAtoms()
        current_resid += 1
        atoms_in_chain += frag.GetNumAtoms()
        if atoms_in_chain >= atoms_per_chain:
            current_chain += 1
            atoms_in_chain = 0

    box = u.dimensions[:3] / 10.0
    write_gro(atom_list, box, output_gro)
    write_itp(atom_list, bond_list, output_itp)
    write_pdb_with_ter_and_conect(atom_list, bond_list, output_pdb)
    print(f" Wrote backmapped structure to {output_gro}")
    print(f" Wrote topology (atoms + bonds) to {output_itp}")
    print(f" Wrote PDB with TER and CONECT to {output_pdb}")

if __name__ == "__main__":
    backmap_to_gromacs(
        "micelle/step5.2_production_pull_final_now.gro",
        "mapping_noh.json",
        "micelle/backmapped_v6_noh.gro",
        "micelle/topology_v6_noh.itp",
        "micelle/backmapped_v6_noh.pdb"
    )
