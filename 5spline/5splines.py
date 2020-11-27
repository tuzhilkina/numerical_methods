import matplotlib.pyplot as plt
import numpy as np


def F0(t):
  if t < -1:
    f = 0
  elif t < 0:
    f = ((t + 1) ** 2) / 2
  elif t < 1:
    f = -(t + 1) ** 2 + 3 * (t + 1) - 3 / 2
  elif t < 2:
    f = ((t + 1) ** 2) / 2 - 3 * (t + 1) + 9 / 2
  elif t >= 2:
    f = 0
  return f

def P(t):
  p = [0.0, 0.0]
  for i in range(len(points)):
    p = [x + y for x, y in zip(p, [F0(t - i) * x for x in points[i]])]
  return p

points = [[0.285,0.088],[0.559,0.701],
            [0.112,1.490],[-1.044,1.815],
            [-2.422,1.172],[-3.252,-0.484],
            [-2.854,-2.638],[-1.008,-4.371],
            [1.845,-4.737],[4.687,-3.213]]
fig = plt.subplots(1)
t = np.arange(1.0, 9.0, 0.001)
f = [P(i) for i in t]
plt.plot([i[0] for i in points], [i[1] for i in points], 'ro')
plt.plot([i[0] for i in f], [i[1] for i in f])
plt.show()

for i in range(1, 9):
  temp = [0.5 * (x + y) for x, y in zip(points[i-1], points[i])]
  print('P({0}) = ({1:.4f}, {2:.4f}),    0.5 * (P_{3} + P_{4}) = ({5:.4f}, {6:.4f})'.format(i, P(i)[0], P(i)[1], i-1, i, temp[0], temp[1]))

points2 = [[-3.252,-0.484],[-2.854,-2.638],[-1.008,-4.371],
           [1.845,-4.737],[4.687,-3.213]]

x = [i[0] for i in points2]
y = [i[1] for i in points2]

def calc_c(i):
  left_side = (2 * x[i-1] * x[i] * (y[i] - c[i-1]) - x[i] ** 2 * (y[i-1] - c[i-1]) - x[i-1] ** 2 * (y[i] - c[i-1])) / (x[i-1] * (x[i] - x[i-1]))
  right_side = left_side - (x[i] ** 2 * y[i+1] - 2 * x[i] * x[i+1] * y[i] + x[i+1] ** 2 * y[i]) / (x[i+1] * (x[i+1] - x[i]))
  result = right_side * x[i+1] / (x[i] - x[i+1])
  return result

def calc_a(i):
  return (x[i] * (y[i+1] - c[i]) - x[i+1] * (y[i] - c[i])) / (x[i] * x[i+1] * (x[i+1] - x[i]))

def calc_b(i):
  return (x[i+1] ** 2 * (y[i] - c[i]) - x[i] ** 2 * (y[i+1] - c[i])) / (x[i] * x[i+1] * (x[i+1] - x[i])) 


c = [1]
a = [calc_a(0)]
b = [calc_b(0)]

for i in range(1, 4):
  c.append(calc_c(i))
  a.append(calc_a(i))
  b.append(calc_b(i))

for i in range(0, 4):
  print('{0} < x < {1}: {2}*x^2 + {3}*x + {4}'.format(x[i], x[i+1], a[i], b[i], c[i]))


def S(n):
  if (n < x[1]):
    s = a[0] * n ** 2 + b[0] * n + c[0]
  elif (n < x[2]):
    s = a[1] * n ** 2 + b[1] * n + c[1]
  elif (n < x[3]):
    s = a[2] * n ** 2 + b[2] * n + c[2]
  else:
    s = a[3] * n ** 2 + b[3] * n + c[3]
  return s

t = np.linspace(x[0], x[4], 100)
f = [S(i) for i in t]
plt.plot(x, y, 'ro')
plt.plot(t, f)
plt.show()