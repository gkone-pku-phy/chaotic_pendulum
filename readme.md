# 混沌摆动力学模拟

This is a simulation of chaotic pendulum kinetic system in JavaScript 

使用方式：直接打开网页。

Usage: open the webpage directly.

更改参数：打开控制台，调用函数`livelyPropertyListener(key,val)`改变对应的参数。参数变化后，会重新初始化系统。

Change the parameter: Open the console, call the function `livelyPropertyListener(key,val)` to change the corresponding coefficients. The system will reinitialize after these changes.

可供改变的参数有：
```
count:int           小球个数
backgroundColor:str 背景色
showOrbit:[T/F]     是否显示轨道
orbitLength:int     轨道长度
standard:[T/F]      是否为标准模式（小球等质量，杆等长）。
非标准模式下，二者均在一定范围内随机产生。

```