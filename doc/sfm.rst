How to use the toolbox for Structure from Motion.

Creating the data
#################

Usually, you are given 2D measurements and you want to recover the 3D structures. generateToyAnimation is a function that can generate certain kinds of data (mostly some typical NRSFM data but also random rigid data).

.. code-block:: matlab

  anim=generateToyAnimation( 0,'nPoint',50,'nFrame',10);

This will create random rigid 3D data as well as random views and the corresponding orthographic projections.

You can also add 5% of noise as follows:

.. code-block:: matlab

  anim=anim.addNoise('noiseS','5','doFillW',true);

3D reconstruction
#################

Algorithms implemented
**********************

I implemented and tested the following algorithms:
  * projective camera
  
    * projective reconstruction:
    
      * Essential matrix decomposition. Reference: HZ2, p259, Result 9.19, and p. 294
      * Normalized 8-point algorithm. Reference: HZ2, p279 and alg 11.1 p. 282
      * Projective Factorization. Reference: HZ2, p445, Algorithm 18.2, Sturm and Triggs ECCV 96, A Factorization Based Algorithm for Multi-Image Projective Structure and Motion
      * Projective Factorization. Reference: Oliensis, Hartley, PAMI 07, Iterative Extensions of the Sturm/Triggs Algorithm Convergence and Nonconvergence
      
    * affine and metric upgrade. Reference: Chandraker IJCV 09: Globally Optimal Algorithms for Stratified Autocalibration (with fixes to bugs I found in the paper (as well as home made optimizations)
    
  * orthographic camera
  
    * affine reconstruction for calibrated cameras. Reference: HZ2, p437, Algorithm 18.1, Tomasi Kanade
    * metric reconstruction for calibrated cameras. Reference: Tomasi Kanade with the metric constraint
    * metric reconstruction for uncalibrated cameras. Not sure there is any paper on that but it's a simple SDP

Examples
********

If we want to recreate the 3D data from the projection, we just need to call:

.. code-block:: matlab

  animReconstructed = computeSMFromW( false, anim.W, 'method', Inf, 'nItrSBA', 100 );

Just type 

.. code-block:: matlab

  help computeSMFromW

to figure out that the above uses standard SfM algorithms and bundle adjustment to solve for the optimal reconstruction with calibrated orthographic cameras.

To compare with the original, different kinds of errors can then be computed using computeError:

.. code-block:: matlab

  err = animReconstructed.computeError('animGT', anim)

(animGT stands for anim Ground Truth)

Viewing the data
################

An Animation object can be displayed by just launching:

.. code-block:: matlab

  playAnim(anim)

The following parameters can be used:

  * 'nCam' [0] number of cameras to show before and after the current one (if -1 => not even the current one is displayed)
  * 'fps' [20] maximum number of frames to display per second. Use fps==0 to introduce no pause and have the movie play as fast as possible
  * 'nLoop' [1] number of time to nLoop video (may be inf), if neg plays video forward then backward then forward etc.
  * 'animGT' [] ground truth animation
  * 'frame' [] only display the specified frames
  * 'camMode' [1] choose the original camera mode as camMode
  * 'showFirst' [false] display the first point differently (with a square)
  * 'showTitle'  [true] display the title (frame number)
  * 'showGT' [false] display the ground truth from the very beginning
  * 'alignGT' [false] align onto the ground truth animation from the very beginning
  * 'showConn' [true] show the connectivities in anim
  * 'showPrettyAxes' [false] show pretty axes (useful for publications)

Other SfM problems
##################

Absolute orientation
********************

The function computeOrientation can help figuring out the absolute orientation between two sets of 3D points.

.. code-block:: matlab

  [ R t s ]=computeOrientation(x,xp,'absolute')

gives the rotation, translation and scale such that x=s*R*xp+t as much as possible. It uses Horn's algorithm as it is closed-form.

If you want to make sure that R has a positive determinant, use 'absoluteHard' instead of 'absolute' which will use Gloptipoly.

Exterior orientation
********************

Same function but with different arguments:

.. code-block:: matlab

  [ R t s ]=computeOrientation(x,xp,'exterior')

gives the rotation, translation and scale such that x=projection*(s*R*xp+t) as much as possible (x is now [ 2 x nPoint ]. This only works for the Orthographic case.

You can also submit an initial value for the rotation 

.. code-block:: matlab

  [ R t s ]=computeOrientation(x,xp,'exterior','RIni',rotationMatrix(rand(1,3)))

If you are dealing with a video sequence and have to do sequential exterior orientation, use the following:

.. code-block:: matlab

  [ R t s ]=computeOrientation(x,xp,'exteriorOrientation')

Triangulation
*************

Similaryl to computeSMFromW, give it isProj, the measurements W, the projection matrices P, and possible intrinsic camera matrices K and method method and that's it:

.. code-block:: matlab

  S = computeSFromWM( isProj, W, P, varargin );
