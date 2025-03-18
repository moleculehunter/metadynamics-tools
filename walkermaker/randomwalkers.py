import numpy as np
import MDAnalysis as mda
from scipy.spatial.transform import Rotation as R

def load_structure(gro_file):
    """Load a GROMACS .gro file."""
    return mda.Universe(gro_file)

def get_random_plane_atoms(universe, dna_selection):
    """Randomly selects four atoms defining a reference plane."""
    dna_atoms = universe.select_atoms(dna_selection)
    if len(dna_atoms) < 4:
        raise ValueError("Not enough atoms in selection to define a plane.")
    return dna_atoms[np.random.choice(len(dna_atoms), 4, replace=False)]

def compute_plane_normal(atoms):
    """Computes the normal vector of a plane defined by four atoms."""
    p1, p2, p3, _ = atoms.positions
    normal = np.cross(p2 - p1, p3 - p1)
    return normal / np.linalg.norm(normal)

def is_collision_free(new_position, dna_atoms, min_distance=2.0):
    """Checks if the new ligand position collides with DNA atoms."""
    for atom in dna_atoms.positions:
        if np.linalg.norm(new_position - atom) < min_distance:
            return False
    return True

def displace_and_rotate_ligand(universe, ligand_selection, dna_atoms, normal, min_distance=5.0, max_distance=15.0, max_angle=360.0):
    """Moves one ligand to a higher distance in random directions relative to the chosen plane and applies random rotation, ensuring no collision with DNA atoms."""
    ligands = universe.select_atoms(ligand_selection)
    if len(ligands) < 2:
        raise ValueError("Less than two ligands found in selection.")
    
    ligand = ligands[:len(ligands)//2]  # Select only one ligand to move
    displacement = np.random.uniform(min_distance, max_distance)
    
    # Generate a random move vector in a random direction
    while True:
        move_vector = np.random.randn(3)
        move_vector -= normal * np.dot(move_vector, normal)  # Ensure movement stays relative to the plane
        move_vector = move_vector / np.linalg.norm(move_vector) * displacement
        new_position = ligand.center_of_mass() + move_vector
        if is_collision_free(new_position, dna_atoms):
            break
    
    # Apply random rotation
    rotation_axis = np.random.rand(3) - 0.5  # Random rotation axis
    rotation_axis /= np.linalg.norm(rotation_axis)
    angle = np.random.uniform(0, max_angle)
    
    rot = R.from_rotvec(angle * rotation_axis)
    rotated_positions = rot.apply(ligand.positions - ligand.center_of_mass()) + ligand.center_of_mass()
    
    ligand.positions = rotated_positions + move_vector

def save_structure(universe, output_file):
    """Save modified structure to a new .gro file."""
    universe.atoms.write(output_file)

def generate_walkers(input_gro, output_prefix, dna_sel='resname DG', ligand_sel='resname A51', num_walkers=40, min_distance=5.0, max_distance=30.0, max_angle=360.0):
    """Generates multiple walkers with ligands spread at higher distances, in random directions above randomly chosen planes, avoiding collisions with DNA atoms."""
    u = load_structure(input_gro)
    
    for i in range(num_walkers):
        walker = u.copy()
        dna_atoms = get_random_plane_atoms(walker, dna_sel)
        normal = compute_plane_normal(dna_atoms)
        displace_and_rotate_ligand(walker, ligand_sel, dna_atoms, normal, min_distance, max_distance, max_angle)
        save_structure(walker, f"{output_prefix}_walker{i+1}.gro")
    
    print(f"Generated {num_walkers} walkers with random ligand positioning and rotation in random planes, avoiding collisions with DNA.")

# Example Usage
generate_walkers("A51zero.gro", "walker")

