# 单像素互相关二阶量，近壁，原点783像素点，单像素瞬时速度场,分段
# -*- coding:utf-8 -*-
import skimage.io as io
import numpy as np
import time
from PIL import Image
from scipy import interpolate
np.set_printoptions(suppress=True)
start = time.perf_counter()
#####################################################################################################
# 基本信息
disk = 'J:'
loc0 = '单像素'
loc1 = 'p'
cal = 1.052/50
zero = 783
t = 1 / 400
nn0 = 4108 * 9
nn1 = 4108 * 10
nn = nn1 - nn0
n = 158
xp = 1280
yp = 800
yp1 = 790   # 上下边界
yp2 = 700
dx = 32     # 粒子位移
dy = 10
y1 = 710    # 相关计算边界
y2 = 780
delx = 64  # 左右边界
#####################################################################################################
# 文件位置
path1 = '%s\\%s\\%s' % (disk, loc0, loc1)
path2 = '%s\\单像素' % path1
path3 = '%s\\速度' % path2
#####################################################################################################
# 读数据
u = np.zeros((int(nn / n), y2 - y1, xp - 2 * delx))
v = np.zeros((int(nn / n), y2 - y1, xp - 2 * delx))
us = np.zeros((y2 - y1, xp - 2 * delx))
vs = np.zeros((y2 - y1, xp - 2 * delx))
for ii in range(nn0, nn1, n):
    print(ii + 1, ii + n)
    sumI = np.zeros((yp1 - yp2, xp))
    Imall = np.zeros((n, yp1 - yp2, xp))
    for i in range(n):
        im = np.array(Image.open('%s\\%d.jpg' % (path1, ii + i + 1)).convert('L'), 'f')
        im = np.array(im)[yp2 - 1: yp1 - 1, :]
        Imall[i, :, :] = im
        sumI += im
#####################################################################################################
    # 平均和均方根
    I1a = (sumI-Imall[0, :, :])/(n-1)
    I2a = (sumI-Imall[n - 1, :, :])/(n-1)
    dI1 = np.zeros((yp1 - yp2, xp))
    dI2 = np.zeros((yp1 - yp2, xp))
    for i in range(n - 1):
        dI1 += (Imall[i, :, :] - I1a) ** 2
        dI2 += (Imall[i + 1, :, :] - I2a) ** 2
    dI1 /= (n - 2)
    dI2 /= (n - 2)
    sigmaI1 = np.sqrt(dI1)
    sigmaI2 = np.sqrt(dI2)
#####################################################################################################
    # 相关计算
    # 倒着走
    for i in range(y1 - yp2, y2 - yp2):
        for j in range(delx, xp - delx):
            cor = np.zeros((dy * 2 + 1, dx + 1))
            for k in range(n - 1):
                im1 = Imall[k, i, j]
                im2 = Imall[k + 1, i - dy: i + dy + 1, j - dx: j + 1]
                cor += (im1 - I1a[i, j])*(im2 - I2a[i - dy: i + dy + 1, j - dx: j + 1])
            R = (cor / (sigmaI1[i, j] * sigmaI2[i - dy: i + dy + 1, j - dx: j + 1])) / (n - 2)
            loc = np.argwhere(R.max() == R)
            yy = loc[0][0]
            xx = loc[0][1]
            if 2 < yy < 2 * dy - 2 and 2 < xx < dx - 2:  # 不要超出边界
                # 二维插值
                xi = np.linspace(-2, 2, 5)
                yi = np.linspace(-2, 2, 5)
                xi, yi = np.meshgrid(xi, yi)
                zi = R[yy - 2: yy + 3, xx - 2: xx + 3]
                interfunc = interpolate.interp2d(xi, yi, zi, kind='cubic')
                xn = np.linspace(-1, 1, 500)
                yn = np.linspace(-1, 1, 500)
                Rn = interfunc(xn, yn)
                loc2 = np.argwhere(Rn.max() == Rn)
                rxx = xx + xn[loc2[0][1]]
                ryy = yy + yn[loc2[0][0]]
                xxx = dx - rxx
                yyy = ryy - dy
            u[int((ii - nn0) / n), i - (y1 - yp2), j - delx] = (xxx * cal / t) / 1000
            v[int((ii - nn0) / n), i - (y1 - yp2), j - delx] = (yyy * cal / t) / 1000
#####################################################################################################
for iii in range(nn0, nn1, n):
    u[int((iii - nn0) / n)].reshape(-1).tofile('%s\\u%d.bin' % (path3, int(iii/n) + 1))
    v[int((iii - nn0) / n)].reshape(-1).tofile('%s\\v%d.bin' % (path3, int(iii/n) + 1))
#####################################################################################################
end = time.perf_counter()
print('Running time: %s Seconds' % (end - start))
