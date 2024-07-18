# GROMACS中利用增强采样(伞型)算法计算分子穿膜的平均力势(PMF)
作者: Yunsen Zhang  编辑: Yunsen Zhang

## 研究背景

分子穿越细胞膜的生物物理过程涉及多种机制和动力学，具体取决于分子的性质和细胞膜的特性。以下是几种主要的穿越方式：

1. 简单扩散 (Simple Diffusion)
简单扩散是指小分子（如氧气、二氧化碳）通过细胞膜的脂双层以随机运动的方式移动。这种过程不需要能量输入，分子从高浓度区域移动到低浓度区域，直到达到平衡。

2. 易化扩散 (Facilitated Diffusion)
易化扩散是指通过特定的膜蛋白（如通道蛋白或载体蛋白）帮助分子穿越细胞膜。这种方式适用于大分子或带电分子（如葡萄糖、离子）。易化扩散也不需要能量输入，依赖浓度梯度。

3. 主动运输 (Active Transport)
主动运输是指分子逆浓度梯度通过细胞膜的过程，需要能量输入（通常来自ATP分解）。这种过程由特定的载体蛋白（如泵）介导，常见的例子包括钠钾泵和质子泵。

4. 内吞作用 (Endocytosis)
内吞作用是指细胞膜包裹外界大分子或颗粒形成囊泡，将其带入细胞内部的过程。这是一种主动运输方式，需要能量。内吞作用包括吞噬作用（吞噬大颗粒或细菌）和吞饮作用（吞饮液体和溶解的分子）。

5. 出胞作用 (Exocytosis)
出胞作用是指细胞通过囊泡将大分子或废物排出细胞外的过程。出胞作用也是一种主动运输方式，需要能量，常见于分泌蛋白质或激素的细胞。

6. 跨细胞运输 (Transcytosis)
跨细胞运输是指分子通过内吞作用进入细胞，再通过出胞作用排出细胞，横跨整个细胞的过程。这种方式常见于某些上皮细胞和内皮细胞，适用于转运抗体、激素等。

理论上说, 上述的6种不同穿越方式, 都可以用分子动力学模拟的方法来研究. 但计算思路各不相同.

今天我们主要探讨的**小分子(如乙醇)跨越细胞膜的分子动力学模拟过程**. 小分子在穿越细胞膜的过程中往往要翻越一个较高的能垒(Energe barrier). 所以常规的分子动力学模拟(Unbiased molecular dynacmis simulation)很难在短时间内捕捉到完整的穿膜过程. 而增强采样技术中的伞型采样算法(Umbrella sampling)通过对分子在特定的反应坐标(在分子穿膜过程中, 反应坐标往往是分子和细胞膜的**质心距离**)上施加限制势, 从而帮助分子动力学模拟算法在一些高能垒的反应坐标处采到丰富的样本.

## 计算流程
1. 小分子跨越细胞膜的初始构象获取.
   
首先使用在线免费药物制剂分子模拟平台FormulationMM(<http://formulationmm.computpharm.org>)和MembPy.py(本人开发的一个专门生各种复杂细胞膜和分子复合物的工具, 即将发布. 关注本人的github以获取最新的发布消息(https://github.com/YYYYYunsen)生成小分子和细胞膜的初始构象.

![Watch the video](https://github.com/YYYYYunsen/MembPy/blob/main/Example/title.jpg)(https://github.com/YYYYYunsen/MembPy/blob/main/Example/drug_membrane_generation.mp4)


2. 伞型采样过程的mdp设置.

首先讲一下mdp的设置:
```
;for pull
pull                     = yes
pull-ngroups             = 2
pull-group1-name         = MOL
pull-group2-name         = memb
pull-ncoords             = 1
pull-coord1-type         = umbrella
pull-coord1-geometry     = distance
pull-coord1-dim          = N N Y
pull-coord1-groups       = 1 2
pull-coord1-k            = 5000.0
pull-coord1-rate         = 0.0
pull-coord1-init         = DIST
pull-coord1-start        = no
pull-pbc-ref-prev-step-com = yes
pull-group1-pbcatom      = 16021 ;随便找一个MOL上面的原子即可
pull-group2-pbcatom      = 1     ;随便找一个细胞膜上面的原子即可
```
我们将pull-coord1-init的值设置为DIST, 这里的DIST代表的就是药物分子和细胞膜质心间的初始距离.

1. 伞型采样过程的代码设置.
```
#!/bin/bash
set -e


for ((i=0;i<60;i++)); do

	d=$(echo "0.1*$(($i))" | bc);

	sed 's/DIST/'$d'/g' mdp/em.mdp > grompp.mdp
	gmx grompp -o em.$i -pp em.$i -po em.$i -n index.ndx -maxwarn 99
	gmx mdrun -deffnm em.$i -pf pullf-em.$i -px pullx-em.$i -ntmpi 1 -ntomp 16 -v

	sed 's/DIST/'$d'/g' mdp/eq.mdp > grompp.mdp
	gmx grompp -o eq.$i -c em.$i -t em.$i -pp eq.$i -po eq.$i -n index.ndx -maxwarn 99
	gmx mdrun -deffnm eq.$i -pf pullf-eq.$i -px pullx-eq.$i -ntmpi  1 -ntomp 16 -update gpu -v 

	sed 's/DIST/'$d'/g' mdp/prd.mdp > grompp.mdp
	gmx grompp -o prd.$i -c eq.$i -t eq.$i -pp prd.$i -po prd.$i -n index.ndx -maxwarn 99
	gmx mdrun -deffnm prd.$i -pf pullf-prd.$i -px pullx-prd.$i -ntmpi  1 -ntomp 16 -update gpu -v

done

ls prd.*.tpr > tpr.dat

ls pullf-prd.*.xvg > pullf.dat

gmx wham -it tpr.dat -if pullf.dat -bsres -bins 200 -unit kJ -nBootstrap 100
```
该代码中变量i代表窗口数目, d代表当前窗口下, 药物分子和细胞膜质心间的距离, d的计算公式中的0.1代表每个窗口之间的间距. 复制上面的代码到空白的文本文件中, 保存为run.bsh. 复制到你的体系所在的文件中然后执行
```
bash run.bsh
```
你可以通过修改变量i和d, 来设置你想要模拟的窗口数量以及窗口间的距离.

4. 结果展示

a.histo图
```
xmgrace -nxy histo.xvg
```
![Image text](https://github.com/YYYYYunsen/MembPy/blob/main/Example/image.png)
从上图可见每一个窗口都合理的分布并重叠, 覆盖整个采样过程.

b.PMF图
```
xmgrace profile.xvg
```
![Image text](https://github.com/YYYYYunsen/MembPy/blob/main/Example/image-2.png)
上图的deltaG就是deltaG(permeability), 用来衡量分子的穿膜难易程度.

## 小结
本人将会在此公众号上持续发布分子动力学模拟教程及本人开发的一些开源计算程序. 公众号MD Mimic Cosmos, 微信号: zys1220098477
