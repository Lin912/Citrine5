# Citrine

Author's Information

github Personal homepage: **https://github.com/Lin912**

ORCid (0000-0003-3820-3199) **[link](https://orcid.org/)**



## Version(5.1.0)

**Extended content:
CleanVersion(VV)  <Top Vel to Bottom Vel>
CleanVersion(VG)  <Top Vel to Bottom Gravity>

## 1.0   Fix external flow fields and add interfaces
We realize the effect of the external flow on the motion of a flexible body cell by adding the Vel<1,2,3>

```
Vel<1,2,3>
double V1 = brr[3];            
double V2 = brr[4];
double V3 = brr[5];
```

## 2.0 Results Showcase:
We did a set of related pendant calculations based on the above theory and code, and the specific distribution state and mechanical properties of the flexible body are shown below:


Change the file "setting.json" to rebuild the project.
The line of CMake path, the CMake source path and the CMake test path.



-------------------------------------------------------
# 3.0 Presumption

First of all, we need to point out that the problem studied in this
article is a structural dynamics problem of a continuum that satisfies
the basic compatibility relationship. In the follow-up of this article,
the finite difference method will be used to discretely solve the
dynamic equations, in which the continuous structure is discretized into
discrete points in time and space. The finite difference method only
divides the continuous structure into a certain number of discrete
points through numerical approximation. It is a discretization process
of continuous features in space. The continuous features of space are
strongly dependent on the number of discrete points in space. Also for
finite difference methods, structural continuity is usually maintained
by defining displacement and velocity variables at the nodes. Therefore,
our subsequent discretization process of the dynamic system is not a
truncation process on the structure, but a description of the dynamic
characteristics of selected points on the continuous structure.

So far, we have given three hypotheses for subordinates. They all belong
to the physical description of Euler-Bernoulli beams in structural
dynamics. These three hypotheses can be found in the following
references.

1.  Each selected cross-section of the cable element is a homogeneous,
    fine element with a uniform surface and a regular circular or
    toroidal cross-section shape

2.  The effect of bending on the cable element is fully represented by
    the Euler-Bernoulli beam equation

3.  The force and deformation characteristics of the element satisfy the
    simple linear theory of elasticity, and the deformation of the
    element is the loaded force at both ends satisfies the following
    linear functional relationship: $d_{p} = (1 + e) d_{s}$

Where $d_{s}$ and $d_{p}$ are the lengths of the element before and
after stretched, respectively, and $e$ is the magnitude of the axial
strain.

But for more complex structural dynamics problems, for example, when we
want to consider the dynamic effects of beam shear deformation on the
structure, more complex beam structure assumptions will be more
suitable. One of the classic beam theories that considers shear
deformation is Timoshenko Beam Theory, also known as shear deformation
beam theory. Unlike the Euler-Bernoulli beam theory, Timoshenko's theory
introduces shear deformation and is therefore suitable for relatively
short, thick, or stiff beams.

We also provide a unified arrangement and summary of the dynamic
variables of the cable here: The mass and weight in air of the
unstretched cable element, represented by $m$ and $w$, respectively.
Each unstretched cable element is identified by the following
homogeneous physical characteristics: cross-sectional diameter $d$,
cross-sectional area $A$, cross-sectional moment of inertia $I$, and
cross-sectional center distance $I_p$. This paper addresses the
mechanical characteristics of cable elements using the Frenet-Serret
frame for continuous cable elements in space, where a set of spatially
distributed orthogonal bases with symbols $\vec{T}$ (Tangent), $\vec{N}$
(Normal), and $\vec{B}$ (Bi-normal) are employed. The distribution of
external forces that act on the cable elements is represented by $F_t$,
$F_n$, and $F_b$ in the three directions. In contrast, the moments of
the external forces around the axes are denoted by $M_t$, $M_n$, and
$M_b$. The water density, cable density, and elastic modulus are denoted
by $\rho$, $\rho_c$, and $E$, respectively.

# The coordinate transformation process(Euler Angle Rotation)

***The Euler angle*** (*Leonhard Euler*) was first proposed by the Swiss
mathematician Leonhard Euler in the 18th century. Related work is
scattered in several papers, but one important paper is \"On the Free
Motion of Rigid Bodies\" (\"On the Motion of Rigid Bodies\"), published
in 1765. *Euler angle* were introduced to describe the rotation of rigid
bodies, describing the orientation of a rigid body with respect to a
fixed coordinate system, which can represent the orientation of a moving
reference system in physics, or the orientation of a general basis in
three-dimensional linear algebra. Classical Euler angles usually use
$0 \deg$ to represent the angle of inclination in the vertical
direction.

In three-dimensional space, the spatial rotation relationship between
multiple coordinates is usually characterized by ***Euler angle***, and
spatial calculations and mathematical transformation relationships are
often involved in the process of calculating Euler angles. There are
usually dozens of different ways to define Euler angles, and different
authors may have their own inherent conventions for their definitions,
which are usually defined in the antecedent content.

Specific definition of *Euler angle*:

1.  Rotation around the *B-axis*($\vec{Z}$) is defined by the $\phi$

2.  Rotation around the *N-axis*($\vec{Y}$) is defined by the $\theta$

3.  Rotation around the *T-axis*($\vec{X}$) is defined by the $\varphi$

where for sub-rotation proofs around one axis there is usually the
following mathematical expression (in terms of rotations around each of
the *three axes*):

$$R_{Z}\phi = 
    \begin{bmatrix}
        cos\phi& -sin\phi& 0\\
        sin\phi& cos\phi& 0\\
        0& 0& 1\\
    \end{bmatrix}$$

$$R_{Y}\theta = 
    \begin{bmatrix}
        cos\theta& 0& sin\theta\\
        0& 1& 0\\
        -sin\theta& 0& cos\theta\\
    \end{bmatrix}$$

$$R_{X}\varphi = 
    \begin{bmatrix}
        1& 0& 0\\
        0&cos\varphi&-sin\varphi\\
        0&sin\varphi&cos\varphi\\
    \end{bmatrix}$$

After determining the method of defining the Euler angles, the process
of transforming the rigid body between coordinate systems is realized by
the matrix product method , for example(*Z-X-Y*):

$$R = R_{Z}R_{X}R_{Y} = 
    \begin{bmatrix}
        cos\phi& -sin\phi& 0\\
        sin\phi& cos\phi& 0\\
        0& 0& 1\\
    \end{bmatrix}
    \cdot
    \begin{bmatrix}
        1& 0& 0\\
        0&cos\varphi&-sin\varphi\\
        0&sin\varphi&cos\varphi\\
    \end{bmatrix}
    \cdot
    \begin{bmatrix}
        cos\theta& 0& sin\theta\\
        0& 1& 0\\
        -sin\theta& 0& cos\theta\\
    \end{bmatrix}$$

At this point, we get the overall Rotation matrix as:

$$R = 
    \begin{bmatrix}
        cos\phi cos\theta - sin\phi sin\theta sin\varphi&-cos\varphi sin\phi&cos\phi sin\theta + cos\theta sin\phi sin\varphi\\
        cos\theta sin\phi + cos\phi sin\theta sin\varphi&cos\phi cos\varphi&sin\phi sin\theta - cos\phi cos\theta sin\varphi\\
        -cos\varphi sin\theta&sin\varphi&cos\theta cos\varphi\\
    \end{bmatrix}$$

## Classical and Tait-Bryan formats for Euler angles

*Peter Guthrie Tait* and *George H. Bryan* introduced the ***Tait-Bryan
form*** for use in aeronautical and marine engineering, where degree $0$
denotes horizontal position. There is no specific paper cited for this
method, and it is usually used as an agreed upon concept to set the
stage and unfold for subsequent work.

For the ***Tait-Bryan form***, the three Euler angles define the motion
around each of the three axes, with the first angle typically
representing rotation around the z-axis, the second angle typically
representing rotation around the y-axis, and the third angle typically
representing rotation around the x-axis. For the classical form, the
motion around the three axes is then represented by two angles. For
example, the first rotation around the axis is a rotation around the
z-axis, the second rotation around the axis is a rotation around the
y-axis, and the third rotation around the axis is still a rotation
around the z-axis.

The order in which the axes are wound directly affects the exact form of
the rotation matrix: for example, there is a big difference between the
***Z-X-Y*** and ***Z-Y-X*** forms.

and the rotation matrix for ***Z-X-Y*** is of the form: $$R = 
    \begin{bmatrix}
        cos\phi cos\theta - sin\phi sin\theta sin\varphi&-cos\varphi sin\phi&cos\phi sin\theta + cos\theta sin\phi sin\varphi\\
        cos\theta sin\phi + cos\phi sin\theta sin\varphi&cos\phi cos\varphi&sin\phi sin\theta - cos\phi cos\theta sin\varphi\\
        -cos\varphi sin\theta&sin\varphi&cos\theta cos\varphi\\
    \end{bmatrix}$$

In the subsequent work of this paper, we take the rotation matrix of the
form *Eq.(7)* (***Z-X-Y***) as the standard form of the rotation matrix
and unfold the later drag-to process accordingly.

## Internal and external rotation

Intrinsic versus extrinsic rotation:In the intrinsic system, each
element rotation is performed on the coordinate system rotated by the
previous operation. In the extrinsic system, each rotation is performed
around the axes of the world coordinate system, which does not move. As
an example, suppose that the three angles of the Eulerian triad specify
rotations around the z, y, and x axes, and in that order. The first
elemental rotation around the z-axis is the same for both intrinsic and
extrinsic conventions. However, for the intrinsic convention, the second
elemental rotation is around the y-axis of the new position resulting
from the first rotation, whereas in the extrinsic convention it is
around the original (unrotated) y-axis. Similarly, the final rotation
around the x-axis will be performed around the x-axis rotated by the
first two operations in the intrinsic system and the original
(unrotated) x-axis in the extrinsic system. This paper will follow the
intrinsic convention:that is, the axis moves with each rotation.

The internal and external rotations described above correspond to the
left and right multiplication of the Euler angle matrix, where the
internal rotation of the Euler angle matrix is generally
right-multiplied, which corresponds to a rotation around one's own
coordinate system, where the base of the coordinate system is being
transformed in real time. The external rotation is realized by
left-multiplying the Euler angle matrix, which corresponds to a rotation
around a fixed coordinate system, where the coordinate vectors are being
transformed.

However, we need to point out the following relationships specifically
and explain them further:

$$P_{b} = R_{ba} \cdot P_{a}$$

For two coordinate systems $a$ and $b$ in space, we have a point denoted
$P_{a}$ under the a coordinate system and $P_{b}$ under the $b$
coordinate system, then the coordinate transformation relationship
between the two coordinate systems can be realized in the form of a
product of Euler matrices. For example, the relationship in the above
equation is the coordinate transformation relationship between from
$P_{a}$ to $P_{b}$.

And, Euler angle matrices are orthogonal arrays which generally satisfy
the relation $R^{T}= R^{-1}$. At this point, we may be able to stop
obsessing about whether or not left or right multiplication is the case.
For again, the following relationship of equations exists:

$$R^{T}_{ba} \cdot P_{b} = R^{-1}_{ba} \cdot P_{b} = R^{-1}_{ba} \cdot R_{ba} \cdot P_{a}$$

As an example, for a force in *space a*(X-Y-Z), its components in the
three directions can be expressed as:

$$f = 
    \begin{bmatrix}
        f_{z}\\
        f_{x}\\
        f_{y}\\
    \end{bmatrix}$$

and the corresponding rotation matrix R,target space b. We follow the
following relation:

$$P_{b} = R_{ba} \cdot P_{a}$$

An expression for the force f in *a-space*, in *b-space*, can be
obtained.

$$\begin{bmatrix}
        f_{z} (cos\phi cos\theta - sin\phi sin\theta sin\varphi) + f_{y} (cos\phi sin\theta + cos\theta sin\phi sin\varphi) - f_{x} cos\varphi sin\phi\\
        f_{z} (cos\theta sin\phi + cos\phi sin\theta sin\varphi) + f_{y} (sin\phi sin\theta - cos\phi cos\theta sin\varphi) + f_{x} cos\phi cos\varphi\\
        f_{x} sin\varphi + f_{y} cos\theta cos\varphi - f_{z} cos\varphi sin\theta\\
    \end{bmatrix}$$

## Angular velocity representation of Euler's angle matrix

### General form

For ***Eulerian angle matrices*** in space, we need specific notation.
In the Eulerian angular representation, the angular velocity $\omega$ is
usually represented by a vector, the components of which correspond to
the rotational velocities around the individual axes. Through the
relationship between the angular velocity vector and the Euler angle
matrix, this can be expressed using the following equation:

$$\omega = R^{-1} \cdot \frac{d R}{d t}$$

In general, $R^{-1}$ = $R^{T}$:

$$R^{-1} = R^{T} = 
    \begin{bmatrix}
        cos\phi cos\theta - sin\phi sin\theta sin\varphi&cos\theta sin\phi + cos\phi sin\theta sin\varphi&-cos\varphi sin\theta\\
        -cos\varphi sin\phi&cos\phi cos\varphi&sin\varphi\\
        cos\phi sin\theta + cos\theta sin\phi sin\varphi&sin\phi sin\theta - cos\phi cos\theta sin\varphi&cos\theta cos\varphi\\
    \end{bmatrix}$$

### Derivation form

Of course, we can also follow the following procedure for derivation.
For the transformation matrix R we have:

$$R = R_{zxy} = R_{z} \cdot R_{x} \cdot R_{y}$$

and for the transformation matrix, the derivatives are:

$$\frac{dR}{dt} = \frac{d(R_{z}\cdot R_{x}\cdot R_{y})}{dt}\\
    = \frac{dR_{z}}{dt}R_{x} R_{y} + R_{z} \frac{d R_{x}}{dt} R_{y} + R_{z} R_{x} \frac{d R_{y}}{dt}$$

In the above equation:

$$\frac{d R_{x}}{dt} = 
    \begin{bmatrix}
        0&0&0\\
        0& -sin \varphi& -cos \varphi\\
        0& cos \varphi& -sin \varphi\\
    \end{bmatrix}
    \cdot \frac{d \varphi}{dt}$$

$$\frac{d R_{y}}{dt} = 
    \begin{bmatrix}
        -sin \theta& 0& cos \theta\\
        0& 0& 0\\
        -cos \theta& 0& -sin \theta\\
    \end{bmatrix}
    \cdot \frac{d \theta}{dt}$$

$$\frac{d R_{z}}{dt} = 
    \begin{bmatrix}
        -sin \phi& -cos \phi& 0\\
        cos \phi& -sin \phi& 0\\
        0& 0& 0\\
    \end{bmatrix}
    \cdot \frac{d \phi}{dt}$$

And for the space curvature matrix $\Omega$, we have the following
relation:

$$\Omega_{0} = \frac{d R}{dt} \cdot R^{T} = (\frac{d R_{z}}{dt} R_{x} R_{y} + R_{z} \frac{d R_{x}}{dt} R_{y} + R_{z} R_{x} \frac{d R_{y}}{dt}) \cdot R^{T}$$

$$\Omega_{0} = \frac{d R_{z}}{dt} R^{T}_{z} + R_{z} \frac{d R_{x}}{dt} R^{T}_{x} R^{T}_{y} + R_{z} R_{x} \frac{d R_{y}}{dt} R^{T}_{y} R^{T}_{x} R^{T}_{z}$$

Based on the following equation, we obtain an expression for ***the
angular velocity component***:

$$\Omega_{0} = 
    \begin{bmatrix}
     cos \theta&0&-cos \varphi sin \theta\\
     0&1&sin \varphi\\
     sin \theta& 0& cos\varphi cos\theta\\
    \end{bmatrix}
    \cdot
    \begin{bmatrix}
    \dot \varphi\\
    \dot \theta\\
    \dot \phi\\
    \end{bmatrix}$$

$$\omega_{x} = cos \theta \dot \varphi - cos \varphi sin \theta \dot \phi\\ $$

$$\omega_{y} = \dot \theta\\$$

$$\omega_{z} = sin \theta \dot \varphi + cos \varphi cos \theta \dot \phi\\$$

# Representation of variables

The following physical and mechanical properties should be specified:
$m_0$ and $w_0$ denote the mass and weight of per unit unstretched
length. Moreover, let $d_0$, $A$, $I$, and $I_p$ denote, respectively,
the diameter, the cross-sectional area, the second moment, and the polar
moment of the unstretched cable. The density of water, the density of
cable, the modulus of elasticity, and the shear modulus are denoted by
$\rho_w$, $\rho_c$, $E$, $G$, respectively.

All the variables are defined to the local Lagrangian coordinate system
with the symbols $\vec{t}$, $\vec{n}$ and $\vec{b}$ respect denotes the
tangent, normal, and binormal vector to the axis.Fig 1 show the
relationship between $\vec{T}$, $\vec{T_1}$, $\vec{T_2}$, $\vec{T_3}$
and $\vec{M}$, $\vec{M_1}$, $\vec{M_2}$, $\vec{M_3}$ which denote
respectively the vectors of forces and moments. Furthermore,$\vec{T_1}$,
$\vec{T_2}$, $\vec{T_3}$ and $\vec{M_1}$, $\vec{M_2}$, $\vec{M_3}$
denotes the tension, the tangential force, the vertical force, the
moment due to tension and two orthogonal bending moments apply on the
cable. Finally, vector $\vec{R}$ denotes the external force arising from
mass and fluid effects.

# Construction of kinetic equations

The process of modeling the dynamics of a cable consists of three basic
equations:

1.  Equation of inertia based on string vibration theory and Newton's
    laws\
    *Refer to the following website:
    https://zhuanlan.zhihu.com/p/139238713*

2.  Moment equations obtained from the moment equilibrium relationship

3.  the compatibility equation

## Derivatives in time and space

To obtain the exact material derivatives in the Langrangian coordinate
system of arbitrary vector $\vec{A}$, the variation of the dynamic
components in time and space are divided into Time-related differential
term $dt$ and Spatial-related differential terms $ds$. Given the
arbitrary vector $\vec{A}$ in terms of the Langrangian coordinate
system: $$\vec{A} = A_1\vec{t} + A_2\vec{n} + A_3\vec{b}$$

The time derivative of $\vec{A}$ is written as:
$$\frac{D\vec{A}}{Dt} = \frac{\partial \vec{A}}{\partial t} + A_1\frac{\partial \vec{t}}{\partial t} + A_2\frac{\partial \vec{n}}{\partial t} + A_3\frac{\partial \vec{b}}{\partial t}$$

The rotation velocity is written as a Darboux vector form[@Hover1]
[@Hover2]:
$$\vec{\omega} = \omega_1\vec{t} + \omega_2\vec{n} + \omega_3\vec{b}$$

Terms $\omega_1$, $\omega_2$ and $\omega_3$ is the component of the
rotation velocity along coordinate axes $\vec{t}$, $\vec{n}$ and
$\vec{b}$.

Then the time derivatives term in Eq.(2) are expanded to symmetrical
properties form:
$$\frac{\partial \vec{t}}{\partial t} = \vec{\omega} \times \vec{t}$$

$$\frac{\partial \vec{n}}{\partial t} = \vec{\omega} \times \vec{n}$$

$$\frac{\partial \vec{b}}{\partial t} = \vec{\omega} \times \vec{b}$$

After introducing Eqs.(4)-(6) into Eq.(2), the latter becomes:
$$\frac{D\vec{A}}{Dt} = \frac{\partial \vec{A}}{\partial t} + \vec{\omega} \times \vec{A}$$

The space derivatives term can be written as the same form:
$$\frac{D\vec{A}}{Ds} = \frac{\partial \vec{A}}{\partial s} + \vec{\Omega} \times \vec{A}$$

where $\vec{\Omega}$ denote the Darboux vector of the rotation for a
space curve.

## Dynamic equilibrium of the cable element

### Inertia equation

The Newton's second law is applied to each discrete cable segment with
the conversion of mass yields [@Hover3]:
$$m\frac{D\vec{V}}{Dt} = \frac{D\vec{T}}{Ds} + (1+e)\sum \vec{R}$$ where
$V$ and $T$ is the sum of three vector of partial velocities and force
along the $\vec{t}$, $\vec{n}$ and $\vec{b}$ direction, can be denoted
by : $$\vec{V} = u\vec{t} + v\vec{n} + w\vec{b}$$
$$\vec{T} = T_e\vec{t} + S_n\vec{n} + S_b\vec{b}$$ The material
derivative function is brought into the conversion of mass equation:
$$m(\frac{\partial \vec{V}}{\partial t} + \vec{\omega} \times \vec{V}) = (\frac{\partial \vec{V}}{\partial s} + \vec{\Omega} \times \vec{T}) + (1+e)(\sum \vec{R})$$

### Balance of moments equation

The balance of moments equation can be expressed by the following vector
equation [@tjavaras1998mechanics] [@zhu1999mechanics]:
$$\frac{1}{1+e}\frac{D[\rho_c I \vec{\omega}]}{Dt} = \frac{1}{(1+e)^2}(\frac{\partial \vec{M}}{\partial s}+\vec{\Omega} \times \vec{M})+\vec{t} \times (1+e)\vec{T}$$

Where $\vec{M}$ is the sum of moment components along the three axes of
$\vec{t}, \vec{n} and \vec{b}$. $\rho_c I$ is the diagonal $3 \times 3$
matrix with $diag[\rho_c I] = (\rho_c I_p, \rho_c I, \rho_c I)$.
Expanding the full derivative into the material derivative yields:
$$\frac{\rho_c I}{1+e}(\frac{\partial \vec{\omega}}{\partial t} + \vec{\omega}\times\vec{\omega}) = \frac{1}{(1+e)^2}(\frac{\partial \vec{M}}{\partial s} + \vec{\Omega}\times\vec{M}) + \vec{t}\times(1+e)\vec{T}$$

### Compatibility relation

The compatibility relation is derived with the assume which the cable
segment is free of discontinuities. For each independent segment, with
the vector $\vec{S}$ represent each beam element, the function of
space-time continuum relationship can be expressed by:
$$\frac{D}{Dt}[\frac{D\vec{S}}{Ds}] = \frac{D}{Ds}[\frac{D\vec{S}}{Dt}]$$
Taking the relation $\vec{V} = \frac{D\vec{S}}{Dt}$ and
$(1+e)\vec{t} = \frac{D\vec{S}}{Ds}$ into account. The linear
stress-strain relation $e=\frac{T}{EA}$ is also brought into the
equation:
$$\frac{1}{EA}\frac{\partial T}{\partial t} \vec{t} + (1+e)\vec{\omega} \times \vec{t} = \frac{\partial \vec{V}}{\partial s} + \vec{\Omega} \times \vec{V}$$

# Governing equation

The governing equation of the cable element can be described by the
three Equations(45), (47) and (49), which represent the mass conversion,
balance of moment and compatibility relationship. It's time to expand
the cross term and give the detailed expression for the mechanical
distribution.

## Some assumptions

For the problem we are studying, we are not concerned with the role of
torque on the dynamics of the cable. This is because the torque action
in general has a much smaller effect on the dynamical characteristics of
the cable compared to the action of bending moments in other planes.
Accordingly, we have the following assumptions for the angle of torsion
$\varphi$ around the t-axis:

$$\varphi = 0$$

At this point, our rotation matrix is:

$$R = 
    \begin{bmatrix}
        cos \phi cos \theta & - sin \phi& cos \phi sin\theta\\
        cos \theta sin \phi & cos \phi& sin \phi sin \theta\\
        -sin \theta& 0& cos \theta\\
    \end{bmatrix}$$

## Distributed forces

For external forces distributed on the cable, these include gravity,
inertial forces and damping forces (hydrodynamic forces). There are some
differences in the form and direction of action of these three forces,
for example, gravity is a mass force whose direction of action does not
change with the direction of the cable section. The inertial force, on
the other hand, is only considered to act in the n and b directions, and
the magnitude of its action is somewhat related to the local flow
velocity. For the damping force, we use Morrison's equation to solve for
it.

### Gravity

For the mass force, we define its magnitude and direction in the global
coordinate system:

$$G = 
    \begin{bmatrix}
        G_{z}\\
        G_{x}\\
        G_{y}\\
    \end{bmatrix}$$

At this point, we transform it to the spatial follower coordinate
system(*t-n-b*) by rotating the matrix:

$$= [G]^T \cdot [R]$$

$$= 
    \begin{bmatrix}
        G_z cos \phi cos \theta + G_x cos\theta sin\phi - G_y sin\theta\\
        -G_z sin\phi + G_x cos\phi\\
        G_z cos\phi sin\theta + G_x sin\phi sin\theta + G_y cos\theta\\
    \end{bmatrix}$$

According to us, we do not consider torsion and define the original
direction of gravity in the x-direction($G_{y} = G_{z} = 0$):

$$= 
    \begin{bmatrix}
        G_{x} sin\phi cos\theta\\
        G_{x} cos\phi\\
        G_{x} sin\phi sin\theta\\
    \end{bmatrix}$$

For the under water cable, the distributed forces can be divided into
field force and fluid force. The field force is generated by the mass of
the cable element under the action of the gravity field, which is
denoted by the symbol $\vec{R_g}$. And $G_{x}$ is $w_{0}$, $w_{0}$
represent the underwater weight of the structure, with negative
values$(-1000)$ representing the positive direction of gravity and
positive values$(1000)$ representing negative gravity.
$$(1+e)\vec{R_g} = [w_{0} sin\phi cos\theta]\vec{t} + [w_{0} cos\phi]\vec{n} + [w_{0} sin\phi sin\theta]\vec{b}$$

### Fluid force

And the fluid force, calculated by the Morrsion's formula, is divide
into two part.The inertial force term$\vec{R_i}$ and flow resistance
term $\vec{R_r}$. Based on the assumption of neglecting the torsion of
the cable element around the $\vec{t}$-axis, all the forces acting on
the cable element are given by:

$$(1+e)\vec{R_i} = m_a \frac{\partial v_{n}}{\partial t} \vec{n} + m_a \frac{\partial v_{b}}{\partial t} \vec{b}$$

$$(1+e)\vec{R_r} = R_r{}_t \vec{t} + R_r{}_n \vec{n} + R_r{}_b \vec{b}$$

$$R_r{}_t = \frac{1}{2} \pi \rho d_0 C_d{}_t v_{t} |v_{t}| \sqrt{1+e}$$

$$R_r{}_n = \frac{1}{2} \rho d_0 C_d{}_n v_{n} \sqrt{v_{n}^2 + v_{b}^2} \sqrt{1+e}$$

$$R_r{}_b = \frac{1}{2} \rho d_0 C_d{}_b v_{b} \sqrt{v_{n}^2 + v_{b}^2} \sqrt{1+e}$$

where $M$ and $m_a$ are the mass of the cable element and the add mass
under the fluid, $g$ denote the gravitational acceleration, $\rho$ is
the density of the surrounding water. $C_d{}_t$, $C_d{}_n$ and $C_d{}_b$
denote the corresponding coefficients given by the Morrion's formula,
for Morrion's formula and $m_a$ negative value$(-1)$ represent
resistance. Finally the symbols $v_{t}$, $v_{n}$ and $v_{b}$ are the
current-related velocities can be calculated by $V_x$, $V_y$ and $V_z$
[@zhao2021numerical], the velocity components of the fluid around the
cable element:
$$v_{t} = u - (V_{z} cos\phi cos \theta + V_{x} cos \theta sin \phi - V_{y} sin \theta)$$

$$v_{n} = v - (V_{x} cos \phi - V_{z} sin \phi)$$

$$v_{b} = w - (V_{z} cos \phi sin \theta + V_{x} sin \phi sin \theta + V_{y} cos \theta)$$

## Final control equations

After organizing, we obtain the final set of kinetic equations:

### Inertia equations

$$m \frac{\partial u}{\partial t} + m(\omega_2 w - \omega_3 v) - \frac{\partial T}{\partial s} - (S_b \Omega_2 - S_n \Omega_3) - (w_{0}sin\phi cos\theta +\frac{1}{2} \pi \rho d_0 C_d{}_t v_{t} |v_{t}| \sqrt{1+e}) = 0$$

$$m \frac{\partial v}{\partial t} + m(\omega_3 u - \omega_1 w) - \frac{\partial S_n}{\partial s} - (T \Omega_3 - S_b \Omega_1) - (w_{0}cos \phi +m_a \frac{\partial v_{n}}{\partial t} +\frac{1}{2} \rho d_0 C_d{}_n v_{n} \sqrt{v_{n}^2 + v_{b}^2} \sqrt{1+e}) = 0$$

$$m \frac{\partial w}{\partial t} + m(\omega_1 v - \omega_2 u) - \frac{\partial S_b}{\partial s} - (S_n \Omega_1 - T \Omega_2) - (w_{0}sin \phi sin \theta +m_a \frac{\partial v_{b}}{\partial t} +\frac{1}{2} \rho d_0 C_d{}_b v_{b} \sqrt{v_{n}^2 + v_{b}^2} \sqrt{1+e}) = 0$$

### Balance of moments equations

$$(1+e)\rho I_p \frac{\partial \omega_1}{\partial t} - G I_p \frac{\partial \Omega_1}{\partial s} = 0$$

$$(1+e)\rho I \frac{\partial \omega_2}{\partial t} - E I \frac{\partial \Omega_2}{\partial s} - (G I_p - E I) + S_b (1+e)^3 = 0$$

$$(1+e)\rho I \frac{\partial \omega_3}{\partial t} - E I \frac{\partial \Omega_3}{\partial s} - (E I - G I_p) + S_n (1+e)^3 = 0$$

### Compatibility relation of the cable

$$\frac{1}{EA} \frac{\partial T}{\partial t} - \frac{\partial u}{\partial s} - (\Omega_2 w - \Omega_3 v) = 0$$

$$(1+e)\omega_3 - \frac{\partial v}{\partial s} - (\Omega_3 u - \Omega_1 w) = 0$$

$$(1+e)\omega_2 + \frac{\partial w}{\partial s} + (\Omega_1 v - \Omega_2 u) = 0$$

## The final control equation

It can be seen that there are nine control equations, however, we have
12 unknowns to be solved. The dynamic components are: $T$, $S_n$, $S_b$,
$u$ , $v$, $w$, $\omega_1$, $\omega_2$, $\omega_3$, $\Omega_1$,
$\Omega_2$, $\Omega_3$. It's worthless to note in Eqs.(9)-(14), the
angle velocity can be expressed by the Euler angle
$\phi, \theta, \varphi$. With the assume, the torsion of cable element
is negated, the relation of $\Omega_1$ and $\Omega_3$ can be found:
$\Omega_1 = -\Omega_3 tan \theta$. In order to supply the completeness
of the equation, Eqs. (12)-(14) are added to the system of control
functions. Thus, the system of 10 partial differential equations with 10
unknowns is derived. Where $"w_0"$ is underwater weight.
$$m \frac{\partial u}{\partial t} + m (w \frac{\partial \theta}{\partial t} - v \frac{\partial \phi}{\partial t} cos \theta) - \frac{\partial T}{\partial s} - (S_b \Omega_2 - S_n \Omega_3) - w_{0}sin\phi cos\theta -\frac{1}{2} \pi \rho d_0 C_d{}_t v_{t} |v_{t}| \sqrt{1+e} = 0$$

::: small
$$m \frac{\partial v}{\partial t} + m \frac{\partial \phi}{\partial t}(u cos\theta + w sin\theta) - \frac{\partial S_n}{\partial s} - (T \Omega_3 + S_b \Omega_3 tan\theta) - w_{0}cos \phi - m_a \frac{\partial v_{n}}{\partial t} - \frac{1}{2} \rho d_0 C_d{}_n v_{n} \sqrt{v_{n}^2 + v_{b}^2} \sqrt{1+e}  = 0$$
:::

::: small
$$m \frac{\partial w}{\partial t} - m (v \frac{\partial \phi}{\partial t} sin\theta + u \frac{\partial \theta}{\partial t}) - \frac{\partial S_b}{\partial s} + (S_n \Omega_3 tan\theta + T \Omega_2) - w_{0}sin \phi sin \theta - m_a \frac{\partial v_{b}}{\partial t} - \frac{1}{2} \rho d_0 C_d{}_b v_{b} \sqrt{v_{n}^2 + v_{b}^2} \sqrt{1+e}  = 0$$
:::

$$E I \frac{\partial \Omega_2}{\partial s} + E I \Omega_3^2 tan \theta - S_b (1+e)^3 = 0$$

$$E I \frac{\partial \Omega_3}{\partial s} - E I \Omega_2 \Omega_3 tan \theta + S_n (1+e)^3 = 0$$

$$\frac{1}{EA} \frac{\partial T}{\partial t} - \frac{\partial u}{\partial s} - (\Omega_2 w - \Omega_3 v) = 0$$

$$(1+e)\frac{\partial \phi}{\partial t} cos\theta - \frac{\partial v}{\partial s} - \Omega_3(u + w tan\theta) = 0$$

$$(1+e)\frac{\partial \theta}{\partial t} + \frac{\partial w}{\partial s} - (v \Omega_3 tan \theta + \Omega_2 u) = 0$$

$$\Omega_2 - \frac{\partial \theta}{\partial s} = 0$$

$$\Omega_3 - \frac{\partial \phi}{\partial s} cos \theta = 0$$

# The numerical solution method

We currently have a total of ten equations as ten unknowns to be
calculated. It is known that the equations satisfy the closure and can
be solved by a numerical calculation method. The vector of the unknown
can be written as:
$$Y = [u, v, w, T, S_{n}, S_{b}, \theta, \phi, \Omega_{2}, \Omega_{3}]^T$$
The equation of motion written in the vectoral form is:
$$M Y^t + N Y^s + Q = 0$$ where $Y^t$ is the time derivative of the
vector $Y$ to be solved and $Y^s$ is the spatial derivative of the
vector $Y$. The detailed expressions can be found in Appendix.

## Discrete Methods

A numerical discretization method called Keller-box method
[@meek1984nonlinear] is applied to solve the above equation. This
differential format with second-order accuracy in spatial and temporal
dimensions $O(\Delta t^2)+O(\Delta s^2)$ and good stability of the
solution was widely used in the solution of underwater-cable problems.
The motion equation can be discretized into the following form:
$$\textbf{M}\frac{\partial Y}{\partial t} + \textbf{N}\frac{\partial Y}{\partial s} + \textbf{Q} = 0$$
Where the bolded characters **M** and **N** are the time discrete term
correlation matrix and spatial discrete term correlation matrix,
respectively. $10 \times 10$ square matrices **M**, **N** and
$10 \times 1$ column vector **Q** are given in the Appendix at the end
of the article. The above equation is a set of nonlinear
non-simultaneous first-order partial differential equations, that can be
solved by the Newton-Raphson method. Since the solution process is a
quasi-static process, the number of computational layers set after
dividing the computational domain is the ratio of the solution time to
the time step. The spatial time step is advanced in each computational
layer to obtain a steady-state solution that satisfies the computational
residuals in that layer and then advanced to the next computational
layer. By modifying the boundary conditions of different calculation
layers, the motion and force of the cable under different external loads
can be obtained.
