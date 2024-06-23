import numpy as np
import MDAnalysis as mda
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue

def adjust_orientation(atoms, phosphorus_atom_name='P'):
    """Adjust the orientation of a molecule so that the phosphorus atom is always facing upwards."""
    phosphorus_atom = next((atom for atom in atoms if atom.name.strip().startswith(phosphorus_atom_name)), None)
    if phosphorus_atom:
        print(f"Phosphorus atom found: {phosphorus_atom.name} at {phosphorus_atom.coord}")
        if phosphorus_atom.coord[2] < 0:
            for atom in atoms:
                atom.coord[2] *= -1
            print("Molecule flipped to make phosphorus atom face upwards.")
    else:
        print("Phosphorus atom not found in the molecule.")
    return atoms

def load_molecule(filename):
    """Load a molecule from a PDB file, center it at the origin, and adjust its orientation."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('Lipid', filename)
    atoms = [atom.copy() for atom in structure.get_atoms()]
    
    if not atoms:
        raise ValueError(f"No atoms found in the file: {filename}")
    
    geom_center = np.mean([atom.coord for atom in atoms], axis=0)
    print(f"Centering molecule. Geometric center: {geom_center}")

    # Center the molecule at the origin
    for atom in atoms:
        atom.coord -= geom_center
    print(f"After centering: {[atom.coord for atom in atoms]}")

    # Adjust orientation so that the phosphorus atom faces upwards
    atoms = adjust_orientation(atoms)
    print(f"After orientation adjustment: {[atom.coord for atom in atoms]}")

    return atoms

def rotate_around_z(atoms, angle):
    """Rotate atoms around the z-axis by the given angle."""
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle),  np.cos(angle), 0],
        [0, 0, 1]
    ])
    for atom in atoms:
        atom.coord = np.dot(rotation_matrix, atom.coord)
    return atoms

def calculate_lipid_numbers(total_lipids, ratios, ntypes):
    """Calculate the number of each type of lipid based on the provided ratios."""
    if ntypes == 1:
        return [total_lipids]
    
    ratio_values = [int(r) for r in ratios.split(':')]
    if len(ratio_values) != ntypes:
        raise ValueError("The number of ratios provided does not match the number of lipid types.")
    ratio_sum = sum(ratio_values)
    lipid_numbers = [int(total_lipids * r / ratio_sum) for r in ratio_values]
    
    # Ensure the total number of lipids matches the total_lipids
    while sum(lipid_numbers) < total_lipids:
        for i in range(len(lipid_numbers)):
            if sum(lipid_numbers) < total_lipids:
                lipid_numbers[i] += 1
    
    return lipid_numbers

def generate_membrane(params):
    struct = Structure.Structure("Membrane")
    model = Model.Model(0)
    struct.add(model)
    chain = Chain.Chain('A')
    model.add(chain)
    
    boxlen = float(params['boxlen'])
    nlipidtot = int(params['nlipidtot'])
    nside = int(np.ceil(np.sqrt(nlipidtot)))
    stepval = boxlen / nside
    random_rotation = params.get('random_rotation', 'False').lower() == 'true'
    
    # Calculate the number of each type of lipid
    ntypes = int(params['ntype'])
    lipid_numbers = calculate_lipid_numbers(nlipidtot, params.get('lipid_ratios', '1:1'), ntypes)
    
    molecules = []
    resnames = []
    for i in range(ntypes):
        try:
            atoms = load_molecule(params[f'filename_{i+1}'])
            molecules.append(atoms)
            resnames.append(params[f'resname_{i+1}'])
        except ValueError as e:
            print(e)
            continue  # Skip this molecule or handle it differently

    if not molecules:
        print("No valid molecules loaded, aborting membrane generation.")
        return None
    
    resid = 1
    atom_counter = 1
    for layer in range(2):
        z_shift = float(params['refz1']) if layer == 0 else float(params['refz2'])
        placed_lipids = 0
        lipid_counts = [0] * len(molecules)
        for i in range(nside):
            for j in range(nside):
                if placed_lipids >= nlipidtot:
                    break
                
                # Select the lipid type based on the remaining counts
                mol_type = np.random.choice([idx for idx, count in enumerate(lipid_counts) if count < lipid_numbers[idx]])
                lipid_counts[mol_type] += 1
                atoms = [atom.copy() for atom in molecules[mol_type]]

                x, y = i * stepval, j * stepval
                z = z_shift

                if layer == 0 and random_rotation:
                    # Randomly rotate the molecule around the z-axis for the first layer
                    angle = np.random.uniform(0, 2 * np.pi)
                    atoms = rotate_around_z(atoms, angle)
                elif layer == 1:
                    # Flip the molecule along the z-axis for the second layer
                    for atom in atoms:
                        atom.coord[2] *= -1
                    z = -z_shift  # Adjust z position for the second layer

                print(f"Placing molecule at ({x}, {y}, {z}) in layer {layer}, residue ID: {resid}")

                residue = Residue.Residue((' ', resid, ' '), resnames[mol_type], ' ')
                chain.add(residue)
                for atom in atoms:
                    unique_atom_id = atom.name + str(atom_counter)
                    atom.id = unique_atom_id
                    atom.coord += np.array([x, y, z])
                    print(f"Atom {atom.name} final coordinates: {atom.coord}")
                    residue.add(atom)
                    atom_counter += 1

                resid += 1
                placed_lipids += 1

    return struct

def write_pdb_file(structure, filename):
    """Write the PDB file with correct formatting."""
    atom_serial_number = 1
    with open(filename, 'w') as f:
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        line = (
                            "ATOM  {:>5} {:<4} {:>3} {:>1}{:>4}    "
                            "{:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}\n"
                        ).format(
                            atom_serial_number,          # atom serial number
                            atom.name.strip(),           # atom name
                            residue.resname.strip(),     # residue name
                            chain.id,                    # chain identifier
                            int(residue.id[1]),          # residue sequence number
                            atom.coord[0],               # x coordinate
                            atom.coord[1],               # y coordinate
                            atom.coord[2],               # z coordinate
                            atom.occupancy,              # occupancy
                            atom.bfactor,                # temperature factor
                            atom.element.strip()         # element symbol
                        )
                        f.write(line)
                        atom_serial_number += 1
        f.write("END\n")

def fix_pdb_format(filename):
    """Fix the PDB format by removing extra spaces before x coordinate."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    with open(filename, 'w') as f:
        for line in lines:
            if line.startswith("ATOM"):
                line = line[:30] + line[31:]  # Remove the extra space before x coordinate
            f.write(line)

def generate_topology(pdb_filename, top_filename):
    """Generate a topology file from a PDB file using MDAnalysis."""
    u = mda.Universe(pdb_filename)
    residue_counts = {}

    for res in u.residues:
        resname = res.resname
        if resname in residue_counts:
            residue_counts[resname] += 1
        else:
            residue_counts[resname] = 1
    
    with open(top_filename, 'w') as f:
        for resname, count in residue_counts.items():
            f.write(f"{resname} {count}\n")

def read_config(filename="input.txt"):
    """Read configuration from a text file, skipping comment lines."""
    params = {}
    with open(filename, 'r') as file:
        for line in file:
            stripped_line = line.strip()
            # Skip empty lines and comment lines
            if stripped_line and not stripped_line.startswith("#"):
                key, value = stripped_line.split(':', 1)  # Split only at the first colon
                params[key.strip()] = value.strip()
    return params

def main():
    params = read_config()
    structure = generate_membrane(params)
    if structure is not None:
        write_pdb_file(structure, params['outfilename'])
        fix_pdb_format(params['outfilename'])
        generate_topology(params['outfilename'], 'topol.top')

if __name__ == "__main__":
    main()
