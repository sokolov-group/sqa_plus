import subprocess

amplitudes_list = ['t2_aaee', 't2_caee', 't2_ccaa', 't2_ccea', 't2_ccee']

processes = []

for amplitude in amplitudes_list:
    with open('V1__master_part2.py', 'r') as file:
        content = file.read()

    amp_select = content.replace('selected_amplitude = None', "selected_amplitude = '{}'".format(amplitude))

    input_file = 'V1__{}.py'.format(amplitude)
    with open(input_file, 'w') as file:
        file.write(amp_select)

    output_file = 'V1__{}.dat'.format(amplitude)
    with open(output_file, 'w') as file:
        process = subprocess.Popen(['python', input_file], stdout=file)
        processes.append(process)

    print("Python file for amplitude='{}' has been started, output file will be saved to {}".format(amplitude, output_file))

for process in processes:
    process.wait()

print("All processes completed.") 
