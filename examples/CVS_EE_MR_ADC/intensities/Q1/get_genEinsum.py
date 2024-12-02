def extract_einsum_lines(file_path):

    with open(file_path, 'r') as file:

        lines = file.readlines()

    equations = []
    capture = False
    #current_section = []
 
    for line in lines:

        if "genEinsum equations" in line:

            current_section = []
            capture = True
            continue

        elif capture and line.strip().startswith("------"):

            capture = False

            if current_section:
                equations.append(current_section)

            continue

        elif capture and line.strip().startswith("TY"):

            current_section.append(line.strip())
            #equations.append(line.strip())

    return equations

##use
file_path = 'TY_q1_ccaa.dat'
extracted_lines = extract_einsum_lines(file_path)

print('einsum equations from ' + file_path) 

#for line in extracted_lines:
#    print(line)

for section in extracted_lines:
    print('')
    for line in section:
        print (line)
