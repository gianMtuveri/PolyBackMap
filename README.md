# PolyBackMap: Backmapping Coarse-Grained (CG) Polymer Structures to Atomistic Resolution

This repository provides a Python script to backmap coarse-grained (CG) structures (in GROMACS `.gro` format) into atomistic models, using user-defined mapping rules and fragment templates. The process reconstructs atomistic coordinates and topology from a CG trajectory, outputting GROMACS and PDB-compatible files.

## Features

- **Backmapping from CG to atomistic** using a mapping JSON file and fragment templates.
- **Supports GROMACS `.gro` input and output.**
- **Generates GROMACS `.itp` topology** and **PDB files** with correct `TER` and `CONECT` records.
- **Handles multiple chains** with custom chain identifiers.
- **Alignment and placement** of atomistic fragments onto CG beads, maintaining vector orientation.
- **Compatible with MDAnalysis and RDKit** for structural manipulation.

## Installation

Install the required dependencies (Python 3.7+ recommended):

```bash
pip install numpy mdanalysis rdkit scipy
```

> **Note:**  
> RDKit is best installed via [conda](https://anaconda.org/conda-forge/rdkit) due to binary dependencies:
> ```bash
> conda install -c conda-forge rdkit
> ```

## Usage

1. **Prepare Inputs:**
   - A CG `.gro` file (e.g., `micelle/step5.2_production_pull_final_now.gro`)
   - A mapping file in JSON format (`mapping_noh.json`) that specifies bead types, atomistic fragments (Mol2 format), and atom names.
   - Atomistic fragment files (Mol2 format) referenced in your mapping file.

2. **Run the script:**

```bash
python backmap.py
```

By default, the script will look for the above files and produce:
- Atomistic `.gro` file (e.g., `micelle/backmapped_v6_noh.gro`)
- GROMACS topology `.itp` file (`micelle/topology_v6_noh.itp`)
- PDB file with chain and connectivity info (`micelle/backmapped_v6_noh.pdb`)

You can modify the paths and filenames in the `if __name__ == "__main__":` section of `backmap.py`.

## Mapping File Format

The mapping file (`mapping_noh.json`) should be a JSON dictionary mapping CG bead names to fragment and atomistic information. For example:

```json
{
  "BEAD1": {
    "fragment_file": "fragments/bead1.mol2",
    "resname": "RES",
    "atom_names": ["C1", "C2", "O1", ...]
  },
  "BEAD2": {
    ...
  }
}
```

- `fragment_file`: Path to the atomistic fragment in Mol2 format.
- `resname`: Residue name for output files.
- `atom_names`: List of atom names, matching the order in the fragment.

## Output Files

- **GROMACS `.gro` file**: Atomistic coordinates, for use in simulations or visualization.
- **GROMACS `.itp` file**: Simple topology with atoms and bonds.
- **PDB file**: For visualization, with chain IDs, `TER`, and `CONECT` records.

## Customization

- Change chain length and number of chains by editing `atoms_per_chain` and `num_chains` in `backmap.py`.
- Adjust mapping rules and fragment templates in your JSON and Mol2 files.

## Example Command Line

To use custom input/output paths, edit the following section in `backmap.py`:

```python
if __name__ == "__main__":
    backmap_to_gromacs(
        "micelle/step5.2_production_pull_final_now.gro",
        "mapping_noh.json",
        "micelle/backmapped_v6_noh.gro",
        "micelle/topology_v6_noh.itp",
        "micelle/backmapped_v6_noh.pdb"
    )
```

## Notes

- The script assumes all fragments and mapping entries are correct and consistent.
- The atom mass and charge in `.itp` are placeholders and may need to be adapted for your force field.
- For large systems, memory and computation time can increase.

## License

[MIT License](LICENSE)
