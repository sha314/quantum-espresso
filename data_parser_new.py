import numpy as np
import pickle
import matplotlib.pyplot as plt
import getpass

username = getpass.getuser()

#File inputs and naming output file
signature = "nb3s4-run8"
folder_path = "/home/{}/quantum-espresso/".format(username) + signature + "/"
bands_path = folder_path + "bands.dat"
fermi_path = folder_path + "dos.dat"
nscf_path = folder_path + "input/nscf.in"
branches_path = folder_path + "input/bands.in"

output_dir = "/home/{}/materials-project/".format(username)
filename = output_dir + "QE-{}.pkl".format(signature)

## Trying to make data parser that turns repeated steps into their own function
## Edit find_keywords function to add any other keyword classifications
## Run the parse_file on your file then do any necessary post-processing

def find_keywords(prev,line):
    #Assumes that keywords will be "&...", "#..." or follow empty space or "/"
    #If need to update formatting of keywords, change the if statements
    if any("&" in s for s in line):
        return True
    if any("#" in s for s in line):
        return True
    elif any("/" in s for s in prev):
        return True
    elif len(prev) == 0 and len(line) != 0:
        return True
    else:
        return False

#Helper function to count number of lines with same string length 
def get_data_length(i, lines):
    start = i + 1
    while start < len(lines) and not lines[start].strip():
        start += 1

    if start >= len(lines):
        return [], start

    data = []
    j = start

    while j < len(lines):
        split_line = lines[j].strip().split()

        # if not split_line:  
        #     j += 1
        #     continue

        if len(split_line) != 0:
            data.append(split_line)
            j += 1
        else:
            break  

    return data, j



def read_keywords(lines, keywords, include_keyword):
    data_blocks = {}
    i = 0
    #Will turn desired data into list of lists

    formatted_keywords = set(k.lower() for k in keywords)

    while i < len(lines):
        line = lines[i].strip()
        vals = line.split()

        if find_keywords(lines[i-1].strip().split() if i > 0 else [], vals):
            formatted_vals = [v.lower().lstrip("&") for v in vals]
            found_keyword = formatted_keywords & set(formatted_vals)

            if found_keyword:
                keyword = list(found_keyword)[0]
                data, j = get_data_length(i, lines)
                if include_keyword:
                    # Store the keyword line as first entry
                    data = [vals] + data

                data_blocks[keyword] = data
                i = j
                continue 
        i += 1

    return data_blocks



def parse_file(file_path, keywords_list, include_keyword=False):
    """
    Method will parse the file inputted looking for all keywords that are in the
    keywords list. If include_keyword is true it will include the keyword line
    as the first entry in the data. 
    """
    with open(file_path, "r") as f:
        lines = f.readlines()

        file_data = read_keywords(lines, keywords_list, include_keyword)
    
    return file_data

#Retriving data from the 4 files
nscf_data = parse_file(nscf_path, ["CELL_PARAMETERS", 
                                   "ATOMIC_SPECIES", 
                                   "ATOMIC_POSITIONS", 
                                   "K_POINTS"])
bands_data = parse_file(bands_path, ["plot"]) 
branch_data = parse_file(branches_path, ["K_POINTS"]) 
dos_data = parse_file(fermi_path, ["E"], True) 


#Post data processing

#Splitting bands_data into kpoints and bands
bands_data_clean = {"kpoints":[], "bands":[]}
for line in bands_data["plot"]:
    if len(line) == 3:
        bands_data_clean["kpoints"].append(line)
    elif len(line) == 10:
        bands_data_clean["bands"].append(line)
    else: 
        continue

bands_data_clean["bands"] = np.array(bands_data_clean["bands"]).T
#Renaming kpoints from nscf as mesh size
nscf_data["mesh_size"] = nscf_data.pop("k_points")

#Cleaning up branch_data
branch_list = branch_data["k_points"][1:]
branch_data_clean = {"branches":[]}
gamma = 'Γ'
#Looping through paths and stripping info
prev_line = None
index_counter = 0  

for line in branch_list:
    if line[-1] == gamma:
        line[-1] = '\\Gamma'

    if prev_line is not None:
        segment_length = int(prev_line[-3])  
        if segment_length != 0:

            branch_dict = {
                'name': f"{prev_line[-1]}-{line[-1]}",
                        'start_index': index_counter,
                        'end_index': index_counter + segment_length - 1
                    }
                    #print(branch_dict)
            branch_data_clean['branches'].append(branch_dict)
            index_counter = branch_dict['end_index'] + 1
        else:
                None
    prev_line = line


#Cleaning dos & getting E-fermi#  

dos_data_clean = {"dos":[], "e_fermi": None}
for i, line in enumerate(dos_data["e"]):
    if i == 0:
        dos_data_clean["e_fermi"] = float(line[-2])
    else:
        dos_data_clean["dos"].append(line[:2])

#print(dos_data_clean["e_fermi"])

#Combine Data Together
data = bands_data_clean|nscf_data|dos_data_clean|branch_data_clean


print(data["bands"].shape)


#Saving pickle file
# with open(filename, 'wb') as f:
#     pickle.dump(data, f)
#     print(f"Pickled data saved to: {filename}")

