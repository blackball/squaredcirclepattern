Squared Circle Pattern for Camera Calbration
---

The traditional chessboard is not very accurate due to the *corner* concept, and the motion blur from camera, etc. 
This new pattern was designed to combine the detection speed from rectangle with accuracy from circle fitting. 
And this detector could easily integrated into OpenCV's calibration module (I had done this in less than ten minutes).

The default pattern is showed below, and you could use the Qt tool in the tool/ folder to generate a new pattern image 
with a different size.

![alt clusters](https://github.com/blackball/SquaredCirclePattern/raw/master/pattern.png)

If you use SCons as your build tool, and have opencv installed, then simply type this to build:
```python
scons
```
