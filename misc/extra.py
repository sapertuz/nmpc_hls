out_min = [-15,  -3,  -3,  -3]
out_max = [15,   3,   3,   3]
in_max = [100, 100, 100, 100]
in_min = [-100, -100, -100, -100]

factor = []
for i in range(4):
    factor = (out_max[i] - out_min[i]) / (in_max[i] - in_min[i])
    print(factor)