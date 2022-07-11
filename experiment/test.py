import numpy as np

f = open('result.txt', mode = 'r')
# out = open('out.txt', mode = "w")

for i in range(9):
    data = f.readline()
    list = data.split(" ")[:-2]
    sub_core = [0] * 64
    for j in range(1, len(list)):
        sub_core[j % 64] += int(list[j]) - int(list[j-1])
    print("Matrix " + str(i) + ": ")
    # print(sub_core)
    print("[*]mean of num calculated by a subcore: " + str(np.mean(sub_core)))
    print("[*]std of num calculated by a subcore: " + str(np.std(sub_core)) + "\n")
f.close()