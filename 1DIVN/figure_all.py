import os
import re
import matplotlib.pyplot as plt

# 点大小
pointSize = 4

# 是否创建单独的图像
singleFig = 0

# 创建文件夹
if not os.path.exists('figure'):
    os.mkdir('figure')

f = open("Result_utf8-0.05.1D", "r", encoding="utf-8")
fileContent = f.readlines()

reading = 0
q = []
p = []
n = []
index = 1
fn = 0
fq = 0
pl = []
nl = []
labels = []
for line in fileContent:
    if line.startswith('     characteristic axi. stage'):
        reading = 1
        # 保存流量和转速数据
        match = re.search(r'.*N =([\d\.]*).*QL =([\d\.]*)', line)
        fn = float(match.group(1))
        fq = float(match.group(2))
        continue
    if line.startswith('   NU TEK') and reading == 1:
        reading = 0
        # 创建标签
        labels.append('Q='+str(fq)+',N='+str(fn))
        # 创建图像
        plt.figure(1)
        pl.append(plt.scatter(q, p, s=pointSize))
        if singleFig > 0:
            plt.figure(3)
            plt.scatter(q, p, s=pointSize)
            plt.savefig("figure/P" + str(index) + ".png")
            plt.clf()
            plt.close()
        
        plt.figure(2)
        nl.append(plt.scatter(q, n, s=pointSize))
        if singleFig > 0:
            plt.figure(3)
            plt.scatter(q, n, s=pointSize)
            plt.savefig("figure/N" + str(index) + ".png")
            plt.clf()
            plt.close()

        # 保存
        # 清除临时数据
        q = []
        p = []
        n = []
        # 递增
        index = index + 1
        continue
    if reading > 0:
        if line.startswith('     Q(LQ'):
            continue
        nums = re.findall(r'[\d\.]+', line)
        q.append(float(nums[0]))
        p.append(float(nums[1]))
        n.append(float(nums[2]))



plt.figure(1)
plt.xlabel('Q')
plt.ylabel('P')
plt.legend(handles=pl, labels=labels, loc='best')
plt.savefig("figure/P.png")
plt.close()
plt.figure(2)
plt.xlabel('Q')
plt.ylabel('N')
plt.legend(handles=nl, labels=labels, loc='best')
plt.savefig("figure/N.png")
plt.close()

print(str(index-1) + ' figures successfully')
