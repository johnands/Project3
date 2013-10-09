from scitools.std import *

infile1 = open('positions.dat', 'r')

# number of objects
n = int(infile1.readline())

# read fil and append to lists
# x30 is the third x-value of object zero:
# positions.dat is made in the following way: x00 y00 x01 y01 x02 y02 ... x0n-1 y0n-1
#                                             x10 y10 x11 y11 x12 y12 ... x1n-1 y1n-1
#                                             ..   ..  ..  ..  ..  .. ...  ..     ..
#                                             xN0 yN0 xN1 yN1 xN2 yN2 ... xNn-1 Nn-1
# I will make two lists x and y: 
# x = [[x00, x01, x02, ... , x0n-1], [x10, x11, x12, ... , x1n-1], ... , [xN0, xN1, xN2, ... , xNn-1]] 
                     
x = []
y = []
for line in infile1:
    words = line.split()
    xx = []
    yy = []
    for i in range(n):
        xx.append(float(words[2*i]))
        yy.append(float(words[2*i+1]))
    x.append(xx)
    y.append(yy)

infile1.close()

x = array(x) 	# [[x00 x01 x02], [x10 x11 x12] ...
y = array(y)

# the x-values of the Sun is then the first column of x, which is indexed as x[::,0]

for i in range(n):
    plot(x[::,i], y[::,i])
    hold('on')
hold('off')

if n == 2:
    legend('Sun', 'Earth')
elif n == 3:
    legend('Sun', 'Earth', 'Jupiter')
elif n == 10:
    legend('Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Neptune', 'Uranus', 'Pluto')

axis('equal')
#axis([-5, 3, -5, 3])

raw_input('Press Return key to quit: ')


infile2 = open('energymom.dat', 'r')

K = []
U = []
Lz = []
for line in infile2:
    words = line.split()
    K.append(float(words[0]))
    U.append(float(words[1]))
    Lz.append(float(words[2]))

infile2.close()

K = array(K)
U = array(U)
Lz = array(Lz)

plot(K)
hold('on')
plot(U)
hold('on')
plot(Lz)
hold('on')
plot(K+U)
legend('K', 'U', 'Lz', 'K+U')
raw_input('Press Return key to quit: ')



