import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import plotly.graph_objects as go

fig = go.Figure(go.Carpet(
    a = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
    b = [4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6],
    x = [2, 3, 4, 5, 2.2, 3.1, 4.1, 5.1, 2.7, 3.6, 4.6, 5.6, 1.5, 2.5, 3.5, 4.5],
    y = [1, 1.4, 1.6, 1.75, 2, 2.5, 2.7, 2.75,4, 4.5, 4.7, 4.75, 3, 3.5, 3.7, 3.75],
    aaxis = dict(
        tickprefix = 'a = ',
        smoothing = 0,
        minorgridcount = 9,
        type = 'linear'
    ),
    baxis = dict(
        tickprefix = 'b = ',
        smoothing = 0,
        minorgridcount = 9,
        type = 'linear'
    )
))

fig.show()

fig = go.Figure(go.Carpet(
    a = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
    b = [4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6],
    x = [2, 3, 4, 5, 2.2, 3.1, 4.1, 5.1, 1.5, 2.5, 3.5, 4.5],
    y = [1, 1.4, 1.6, 1.75, 2, 2.5, 2.7, 2.75, 3, 3.5, 3.7, 3.75],
    aaxis = dict(
        tickprefix = 'a = ',
        smoothing = 0,
        minorgridcount = 9,
        type = 'linear'
    ),
    baxis = dict(
        tickprefix = 'b = ',
        smoothing = 0,
        minorgridcount = 9,
        type = 'linear'
    )
))

fig.show()

# df1 = pd.read_excel('Aircraft Data.xlsx')
# print(df1)
# for column in df1:
#     if pd.isna(df1[column][1]) == False:
#         print(df1[column])
#         Variables = np.linspace((df1[column][1]), (df1[column][2]), int(df1[column][3]))
#         VariableName = column
#         print(VariableName)
# name = 'Takeoff Weight'
# for line in name.split():
#     print(line)
# array2d = [[1,2],
#            [1,3],
#            [1,4],
#            [1,5],
#            [1,6],
#            [1,7]]
# array1d = np.array([1,1,1,2,2,2])
# plt.contour(array2d, array1d)
# x = np.expand_dims(np.arange(1,11,1), axis=1)
# y = np.expand_dims(np.arange(2,21,2), axis=0)
# z = y * x

# print(x.shape)
# print(y.shape)
# print(z.shape)

# plt.figure()
# plt.contour(z)

x = np.linspace(750,1500,100)
y = np.linspace(0,1,101)
xx,yy = np.meshgrid(x,y)
row_vec = y[None]
col_vec = x[:, None]
print((np.ones([100,100])).shape)
print(yy)
print(xx)
a = xx.shape[0]
print(a)
#plt.plot(xx, yy, marker='.', color='k', linestyle='none')
contours = plt.contour(xx,yy,xx*yy, levels = 10)
plt.clabel(contours, inline=True, fontsize=8)
#plt.colorbar()


# plt.show()
