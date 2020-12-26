# 单像素互相关二阶量，近壁，原点783像素点
# -*- coding:utf-8 -*-
import skimage.io as io
import numpy as np
import time
import csv
from PIL import Image
from scipy import interpolate
np.set_printoptions(suppress=True)
start = time.perf_counter()
#####################################################################################################
# 基本信息
cal = 1.052/50
zero = 783
t = 1 / 400
nn = 41080
n = 4108
num = 9
xp = 1280   # 图片信息
yp = 800
yp1 = 790   # 上下边界
yp2 = 700
y1 = 710    # 相关计算边界
y2 = 780
dx = 32     # 粒子位移
dy = 8
lx = 128     # 查询窗口
bc = 32
delx = 128  # 左右边界
#####################################################################################################
# 文件位置
disk = 'J:'
loc0 = '单像素'
loc1 = 'p'
path1 = '%s\\%s\\%s' % (disk, loc0, loc1)
path2 = '%s\\单像素' % path1
path3 = '%s\\41080-1220' % path2
path4 = '%s\\二阶量' % path3
#####################################################################################################
# 读数据
print('开始读取')
with open('%s\\单像素.csv' % path3, 'r') as f:
    data_a = csv.DictReader(f)
    ua = [row[' "u(m/s)"'] for row in data_a][1:]
    ua = [float(m) for m in ua]
    ua = np.array(ua)
with open('%s\\单像素.csv' % path3, 'r') as f:
    data_a = csv.DictReader(f)
    va = [row[' "v(m/s)" '] for row in data_a][1:]
    va = [float(m) for m in va]
    va = np.array(va)
Imall = np.zeros((n, yp1 - yp2, xp))
for i in range(n):
    im = np.array(Image.open('%s\\%d.jpg' % (path1, i + 1 + n * num)).convert('L'), 'f')
    im = np.array(im)[yp2 - 1: yp1 - 1, :]
    Imall[i, :, :] = im
print('组合完成')
#####################################################################################################
# 倒着走
#####################################################################################################
# 快照
uu = np.zeros(y2 - y1)
vv = np.zeros(y2 - y1)
uv = np.zeros(y2 - y1)
c = np.zeros(y2 - y1)
for k in range(n - 1):
#####################################################################################################
    # y
    for i in range(y1 - yp2, y2 - yp2):
        cor = np.zeros((dy * 2 + 1, dx + 1))
#####################################################################################################
        # x
        for j in range(delx + lx, xp - delx, bc):
#####################################################################################################
            # 平均值和均方根
            I1a = np.mean(Imall[k, i, j - lx + 1: j + 1])
            dI1 = np.mean((Imall[k, i, j - lx + 1: j + 1] - I1a) ** 2)
            dI1 /= lx
            sigmaI1 = np.sqrt(dI1)
            I2a = np.zeros((dy * 2 + 1, dx + 1))
            for jj in range(j - lx + 1, j + 1):
                I2a += Imall[k + 1, i - dy: i + dy + 1, jj - dx: jj + 1]
            I2a /= lx
            dI2 = np.zeros((dy * 2 + 1, dx + 1))
            for jj in range(j - lx + 1, j + 1):
                dI2 += (Imall[k + 1, i - dy: i + dy + 1, jj - dx: jj + 1] - I2a) ** 2
            dI2 /= lx
            sigmaI2 = np.sqrt(dI2)
#####################################################################################################
            # 相关计算
            for jj in range(j - lx + 1, j + 1):
                im1 = Imall[k, i, jj]
                im2 = Imall[k + 1, i - dy: i + dy + 1, jj - dx: jj + 1]
                cor += (im1 - I1a)*(im2 - I2a)
            R = (cor / (sigmaI1 * sigmaI2)) / lx
            loc = np.argwhere(R.max() == R)
            yy = loc[0][0]
            xx = loc[0][1]
#####################################################################################################
            # 二维差值
            if 2 < yy < 2 * dy - 2 and 2 < xx < dx - 2:  # 不要超出边界
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
                if xxx != 0:
                    fu = (xxx * cal / t) / 1000 - ua[[i - (y1 - yp2)]]
                    fv = (yyy * cal / t) / 1000 - va[[i - (y1 - yp2)]]
                    uu[i - (y1 - yp2)] += fu ** 2
                    vv[i - (y1 - yp2)] += fv ** 2
                    uv[i - (y1 - yp2)] += fu * fv
                    c[i - (y1 - yp2)] += 1
uu.tofile("%s\\uu%d.bin" % (path4, num))
vv.tofile("%s\\vv%d.bin" % (path4, num))
uv.tofile("%s\\uv%d.bin" % (path4, num))
c.tofile("%s\\c%d.bin" % (path4, num))
#####################################################################################################
end = time.perf_counter()
print('Running time: %s Seconds' % (end - start))
