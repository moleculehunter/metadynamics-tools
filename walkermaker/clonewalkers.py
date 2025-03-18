import numpy as np
import shutil
import glob

def read_gro(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    return lines

def write_gro(filename, lines):
    with open(filename, 'w') as f:
        f.writelines(lines)

def extract_coordinates(lines, resname, index=1):
    coords = []
    atom_indices = []
    count = 0
    current_residue = None
    
    for i, line in enumerate(lines[2:-1]):  # Skip header and final line
        if len(line) > 20:
            res = line[5:10].strip()
            res_num = int(line[0:5])  # Get residue number
            
            if res == resname:
                if current_residue is None or current_residue != res_num:
                    count += 1
                    current_residue = res_num
                
                if count == index:
                    coords.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
                    atom_indices.append(i + 2)  # Store line indices for modification
    
    return np.array(coords), atom_indices

def modify_coordinates(lines, resname, target_origin, index=1):
    new_lines = lines.copy()
    ligand_coords, atom_indices = extract_coordinates(lines, resname, index=index)
    
    if len(ligand_coords) == 0:
        raise ValueError(f"Ligand {resname} not found in the file.")
    
    # Compute ligand centroid
    ligand_centroid = np.mean(ligand_coords, axis=0)
    
    # Compute shift vector
    shift_vector = target_origin - ligand_centroid
    
    for i, atom_index in enumerate(atom_indices):
        if len(new_lines[atom_index]) > 20:
            x, y, z = ligand_coords[i] + shift_vector
            new_lines[atom_index] = f"{new_lines[atom_index][:20]}{x:8.3f}{y:8.3f}{z:8.3f}{new_lines[atom_index][44:]}"
    
    return new_lines

def process_multiple_files(source_prefix, source_ligand, source_index, target_file, target_ligand, target_index): 
    source_files = sorted(glob.glob(f"{source_prefix}*.gro"))  # Get all matching source files
    
    for i, source_file in enumerate(source_files, start=0):
        output_file = f"{target_ligand}dn{i}.gro"  # Change name as needed (up or dn for output file)
        align_ligand(source_file, source_ligand, target_file, target_ligand, output_file, source_index, target_index)
    
    print("All files have been processed.")

def align_ligand(source_file, source_ligand, target_file, target_ligand, output_file, source_index=1, target_index=1):
    # Read files
    source_lines = read_gro(source_file)
    target_lines = read_gro(target_file)
    
    # Get centroid of the selected reference ligand
    source_coords, _ = extract_coordinates(source_lines, source_ligand, index=source_index)
    
    if source_coords.size == 0:
        raise ValueError(f"Selected ligand {source_ligand} not found in file {source_file}")
    
    # Compute centroid as target position
    target_origin = np.mean(source_coords, axis=0)
    
    # Modify target file for selected ligand
    modified_target_lines = modify_coordinates(target_lines, target_ligand, target_origin, index=target_index)
    
    # Save modified file
    write_gro(output_file, modified_target_lines)
    print(f"Saved modified file as {output_file}")

# Example usage - user can select reference file, reference ligand, target file, and target ligand
process_multiple_files("A51dn", "A51", 2, "C40zero.gro", "C40", 2)

