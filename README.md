# Yarn-cloth-sim

### Current result

I have implemented all the forces mentioned in the paper, and got the following result.

<img src="ani.png" width="60%">

This animation seems right but not robust as time steps increase.

Also, as for the cases when pinned from two corners, it doesn't converge and sometimes crash :(

I think maybe because I have not added boundary condition and just do it through clamping the position and velocity. Also, I think the energy formula mentioned in the paper is kind of strange. I don't know how to calculate the stiffness matrix for the friction. 

It is damn hard!
<<<<<<< HEAD

### 7.30
I change the implementation of boundary condition. And the behavior of stretch is correct now. 
But the determinant of sitffness matrix is too large, I cannot get a appropriate simulation parameter.

I will compute the stiffness matrix of friction later. And then add each energy one by one to check the bug. 

It's really complicated because the paper didn't give enough information...
=======
>>>>>>> 3bca2530fd26bb507288da2952314431b48255a5
