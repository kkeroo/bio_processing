state = [1, 0, 0, 0, 0, 0, 0, 0]

iters = 10

for _ in range(iters):
    y2 = state[0]
    y3 = state[1]
    y4 = state[2]
    y5 = state[3]
    y6 = state[4]
    y7 = state[5]
    y8 = state[6]
    y1 = state[7] ^ state[5] ^ state[3]

    state = [y1, y2, y3, y4, y5, y6, y7, y8]
    print(state)