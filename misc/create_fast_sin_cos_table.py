import matplotlib.pylab as plt
import numpy as np
import sys

# total arguments
n = len(sys.argv)

if n > 1:
    n_sine = int(sys.argv[1]) - 1
    x = np.linspace(0, 2*np.pi, n_sine)
    sinx = np.sin(x)
    plt.plot(x, sinx)
    plt.xlabel('Angle [rad]')
    plt.ylabel('sin(x)')
    plt.axis('tight')
    plt.show()
    
    converted_sinx = [str(element) for element in sinx]
    joined_string = ",".join(converted_sinx)

    with open("file.txt", "w") as output:
        output.write(str(joined_string))

else:
    print("ERROR: Enter number of points")