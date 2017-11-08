# Extented Kalman filter in C++
*Self-driving Car Nanodegree at Udacity.*

**Goals:**

- Learn about the Kalman filter
- Implement the Kalman filter algorithm in C++
- Learn C++

**Tools:**
- [CarND-Extended-Kalman-Filter-Project](https://github.com/udacity/CarND-Extended-Kalman-Filter-Project)
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [Xcode](https://developer.apple.com/support/xcode/)


## Kalman Filter:

The Kalman filter estimates the state of a system based on a set of observed system measurements. It expects gaussian distributions and linear motion and observation model. Essentially the algorithm has two steps: 

**update:** This step updates the last predicted state of a system and it's uncertainty according to the Kalman gain.

**predict:** Predicts the next state of the observed system and it's uncertainty.

In the context of self-driving cars, the Kalman filter is a popular technique used to track pedestrians, other cars or any relevant object that can potentially interfere with the self-driving car itself. The algorithm expects a set of measurements captured by the car's sensors - Radar, Lidar, etc.. - and outputs a gaussian distribution of the observed system state. 


## References

- [Udacity](https://br.udacity.com)
- [Michel van Biezen - The Kalman Filter](https://www.youtube.com/channel/UCiGxYawhEp4QyFcX0R60YdQ)
