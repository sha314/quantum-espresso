import numpy as np
import pickle
import matplotlib.pyplot as plt
import getpass
username = getpass.getuser()

#File inputs and naming output file
signature = "nb3s4-run7"
folder_path = "/home/{}/quantum-espresso/".format(username) + signature + "/"
bands_path = folder_path + "bands.dat"
fermi_path = folder_path + "dos.dat"
nscf_path = folder_path + "input/nscf.in"
branches_path = folder_path + "input/bands.in"

output_dir = "/home/{}/materials-project/".format(username)
filename = output_dir + "QE-{}.pkl".format(signature)


#Getting kpoints and ebands from bands.dat
kpoints = []
ebands = []
try:
    with open(bands_path, 'r') as f:
        #Skip the first line of file
        next(f)
        current_k = None
        current_energies = []
        content = f.read()
        if not content:
            print("File bands.dat is empty")
        else:
            lines = content.splitlines()
            for line in lines:
                #Any line with 3 strings is kpoint 
                #Any line with strings not 0 or 3 will be energies
                values = line.strip().split()
                if not values:
                    continue
                elif len(values) == 3:
                    if current_k is not None:
                        kpoints.append(current_k)
                        ebands.append(current_energies)
                        current_energies = []
                    current_k = [float(v) for v in values]
                else:
                    current_energies.extend([float(v) for v in values])
        
            if current_k is not None:
                kpoints.append(current_k)
                ebands.append(current_energies)
except FileNotFoundError:
    print("Incorrect file path for bands.dat file")

#dos with energies and e_fermi from dos.dat
dos = []
try:
    with open(fermi_path, 'r') as f:
        content = f.read()
        if not content:
            print("File bands.dat is empty")
        else:
            lines = content.splitlines()
            #Taking fermi energy from first line
            first_line = lines[0]
            val = first_line.strip().split()
            e_fermi = float(val[-2])

            for line in lines:
                #For lines with 3 strings taking only first 2 (E, dos)
                values = line.strip().split()
                if len(values) == 3:
                    dos.append([float(values[0]), float(values[1])])
                else:
                    continue
except FileNotFoundError:
    print("Incorrect file path for dos.dat")
    

#Getting branches from bands.in file
branches = []
gamma = 'Γ'

try:
    with open(branches_path, 'r') as f:
        content = f.read()
        if not content:
            print("File bands.in is empty")
        else:
            content_lines = content.splitlines()
            #Reading from bottom of file, taking lines with length 6 until break
            lines = reversed(content_lines)
            entries = []
            for line in lines:
                if len(line.strip().split()) == 6:
                    entries.append(line)
                else:
                    break

            #Looping through paths and stripping info
            prev_line = None
            index_counter = 0  

            for line in reversed(entries):
                curr_line = line.strip().split()
                
                if curr_line[-1] == gamma:
                    curr_line[-1] = '\\Gamma'

                if prev_line is not None:
                    segment_length = int(prev_line[-3])  
                    if segment_length != 0:

                        branch_dict = {
                            'name': f"{prev_line[-1]}-{curr_line[-1]}",
                            'start_index': index_counter,
                            'end_index': index_counter + segment_length - 1
                        }
                        branches.append(branch_dict)
                        index_counter = branch_dict['end_index'] + 1
                    else:
                        None
                prev_line = curr_line
except FileNotFoundError:
    print("Incorrect file path for bands.in file")

#Getting cell params, atomic species, atomic positions & mesh size from nscf.in file

cell_param = []
atomic_species = []
atomic_positions = []

#Helper function to count number of lines with same string length 
def get_data_length(i, lines):
    start = i + 1
    while start < len(lines) and not lines[start].strip():
        start += 1

    if start >= len(lines):
        return [], start

    expected_len = len(lines[start].strip().split())
    data = []
    j = start

    while j < len(lines):
        split_line = lines[j].strip().split()

        # if not split_line:  
        #     j += 1
        #     continue

        if len(split_line) == expected_len:
            data.append(split_line)
            j += 1
        else:
            break  

    return data, j

try:
    with open(nscf_path, 'r') as f:
        content = f.read()
        if not content:
            print("File nscf.in is empty")
        else:
            lines = content.splitlines()
#Getting mesh size from bottom of the file
            rev_lines = reversed(lines)
            vals = next(rev_lines).strip().split()
            mesh = vals[:3]

#Parsing through file looking for keywords then taking following lines with same format
            i = 0
            while i < len(lines):
                if "CELL_PARAMETERS" in lines[i]:
                    data, j = get_data_length(i, lines)
                    cell_param.append(lines[i].strip().split())

                    for line in data:
                        cell_param.append(line)
                    
                    i=j

                elif "ATOMIC_SPECIES" in lines[i]:
                    data, j = get_data_length(i, lines)
                    for line in data:
                        atomic_species.append(line)
                    i=j
                
                elif "ATOMIC_POSITIONS" in lines[i]:
                    atomic_positions.append(lines[i].strip().split())
                    data, j = get_data_length(i, lines)
                    for line in data:
                        atomic_positions.append(line)
                    i=j

                else:
                    i += 1
except FileNotFoundError:
    print("Incorrect file path for nscf.in file")

ebands = np.array(ebands).T


#Combining all data into one dictionary
data = {"kpoints": kpoints, 
        "bands": ebands, 
        "e_fermi": e_fermi, 
        "dos": dos, 
        "mesh_size": mesh, 
        "branches":branches,
        "cell_parameters": cell_param,
        "atomic_species": atomic_species,
        "atomic_positions": atomic_positions}


#Saving pickle file
with open(filename, 'wb') as f:
    pickle.dump(data, f)