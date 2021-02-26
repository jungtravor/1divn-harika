# -*- coding: utf-8 -*-

import os
import re
import matplotlib.pyplot as plt


def readContent(fileContent):
	reading = 0
	q = []
	p = []
	n = []
	for line in fileContent:
	    if line.startswith('     characteristic axi. stage'):
	    	# 开始读取当前流量转速的数据
	        reading = 1
	        # 保存当前流量和转速数据
	        match = re.search(r'.*N =([\d\.]*).*QL =([\d\.]*)', line)
	        fn = float(match.group(1))
	        fq = float(match.group(2))
	        continue
	    if line.startswith('   NU TEK') and reading == 1:
	    	# 该流量转速下数据获取完成
	        continue
	    if reading > 0:
	    	# 读取该行数据，依次为流量、压力、转速
	        if line.startswith('     Q(LQ'):
	            continue
	        nums = re.findall(r'[\d\.]+', line)
	        q.append(float(nums[0]))
	        p.append(float(nums[1]))
	        n.append(float(nums[2]))
	return q, p, n

# 点大小
pointSize = 2

# 创建文件夹
dirName = 'figure_allinone'
if not os.path.exists(dirName):
    os.mkdir(dirName)

# 获取结果文件列表
filelist = [name for name in os.listdir('.')
            if name.startswith('Result_utf8')]

# 逐一读取文件信息，并在图像上画点
plt.figure(1, figsize=(19.20, 10.80))
plt.figure(2, figsize=(19.20, 10.80))
pl = []
nl = []
labels = []
for filename in filelist:
	q = []
	p = []
	n = []
	# 读取文件
	with open(filename, "r", encoding="utf-8") as f:
		detail = re.search(r'Result_utf8-(.*)\.1D', filename).group(1)
		fileContent = f.readlines()
		(q, p, n) = readContent(fileContent)

	plt.figure(1)
	pl.append(plt.scatter(q, p, s=pointSize))
	plt.figure(2)
	nl.append(plt.scatter(q, n, s=pointSize))
	labels.append(detail)

# 处理坐标轴、系列等数据
plt.figure(1)
plt.xlabel('Q')
plt.ylabel('P')
plt.legend(handles=pl, labels=labels, loc='best')
plt.savefig(dirName + "/P.png")
plt.close()
plt.figure(2)
plt.xlabel('Q')
plt.ylabel('N')
plt.legend(handles=nl, labels=labels, loc='best')
plt.savefig(dirName + "/N.png")
plt.close()

print('figures saved successfully')
