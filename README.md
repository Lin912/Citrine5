# Citrine4

Author's Information

github Personal homepage: **https://github.com/Lin912**

ORCid (0000-0003-3820-3199) **[link](https://orcid.org/)**



## Version(5.1.0)

**Extended content:
CleanVersion(VV)  <Top Vel to Bottom Vel>
CleanVersion(VG)  <Top Vel to Bottom Gravity>

## 1.0  Theory
We base our assumptions on the following:
(1)The Euler-Bernoulli beam theory is applied to each element of the selected elastic cable. 
(2)The cross-section of the cable element remains homogeneous and round; 
(3)The tension is a single-valued function of the strain.


![1](images/Diagram.jpg)


## 4.0   Fix external flow fields and add interfaces
We realize the effect of the external flow on the motion of a flexible body cell by adding the Vel<1,2,3>

```
Vel<1,2,3>
double V1 = brr[3];            
double V2 = brr[4];
double V3 = brr[5];
```

## 5.0 Results Showcase:
We did a set of related pendant calculations based on the above theory and code, and the specific distribution state and mechanical properties of the flexible body are shown below:

![2](images/2D%20slushing1.jpg)
![3](images/2D2%20slushing.jpg)
![5](images/2D%20slushing2.jpg)

