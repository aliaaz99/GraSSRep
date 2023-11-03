import subprocess

# vary the length of repeats and save them in the folder 'length':
repeat_length = [150, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
C_const = 25
for L_i in repeat_length:
    print("Data generating with  C=25 and L=" + str(L_i) + "")
    name_i = 'C25_L' + str(L_i)
    cmd_repeat_generate = 'python' + ' insertRepeat.py' + ' --name ' + name_i + ' --L ' + str(L_i) + ' --C ' + str(C_const) 
    p_repeat_generate = subprocess.run(cmd_repeat_generate, shell=True)

# vary the length of repeats and save them in the folder 'copyNum':
copy_number = [10, 20, 30 ,40 ,50 ,60, 70, 80, 90, 100, 125, 150]
L_const = 400
for C_i in copy_number:
    print("Data generating with  C=" + str(C_i) + " and L=400")
    name_i = 'L400_CN' + str(C_i)
    cmd_repeat_generate = 'python' + ' insertRepeat.py' + ' --name ' + name_i + ' --L ' + str(L_const) + ' --C ' + str(C_i) 
    p_repeat_generate = subprocess.run(cmd_repeat_generate, shell=True)


# Fix the length and copy number to vary the number of reads:
C_const = 25
L_const = 400
print("Data generating with  C=" + str(C_const) + " and L=" + str(L_const) + "")
name_i = 'C25_L400_coverage'
cmd_repeat_generate = 'python' + ' insertRepeat.py' + ' --name ' + name_i + ' --L ' + str(L_const) + ' --C ' + str(C_const) 
p_repeat_generate = subprocess.run(cmd_repeat_generate, shell=True)