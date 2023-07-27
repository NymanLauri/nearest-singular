Accompanying code for "A Riemannian optimization method to compute the nearest singular pencil", F. Dopico, V. Noferini and L. Nyman.

This code has been tested with Manopt 7.1.

Example:

    >> importmanopt
    Manopt was added to Matlab's path.
    >> rng(0); n = 4;
    >> A = randn(n) + randn(n)*1i; 
    >> B = randn(n) + randn(n)*1i;
    >> [S,T,distance,~,Q,~] = nearest_complex(A, B);

The distance between B x + A and the singular T x + S is

    >> norm([A-S B-T],'f')
    ans =
      1.8432

which is also stored as the last element of the distance array:

    >> distance(end)
    ans =
      1.8432

We can bring the pencil T x + S into Schur from via Q:
   
    >> Q(:,:,1)*S*Q(:,:,2), Q(:,:,1)*T*Q(:,:,2)
    
    ans =
       0.0594 - 2.1844i   1.1663 - 1.4270i   2.0362 - 0.3444i   0.4781 - 0.5793i
       0.0000 - 0.0000i   0.0000 + 0.0000i   5.4332 + 1.8117i   0.6542 + 0.5363i
      -0.0000 + 0.0000i  -0.0000 - 0.0000i  -2.5804 + 0.5907i   1.0023 - 0.1099i
       0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.0979 + 1.9399i
    
    ans =
       1.2238 - 1.3734i   0.0158 + 0.4057i   0.6777 + 1.0727i  -0.7111 + 0.2648i
       0.0000 + 0.0000i  -0.0000 + 0.0000i   0.2778 + 0.7166i  -0.1375 - 0.9701i
      -0.0000 + 0.0000i   0.0000 - 0.0000i  -1.0159 - 1.5481i   3.8319 + 1.1145i
       0.0000 + 0.0000i  -0.0000 - 0.0000i  -0.0000 + 0.0000i   1.0218 - 1.8142i

The element (2,2) of the Schur form is zero, and hence the pencil is indeed singular.


You can also specify the maximum amount of iterations, maximum amount of time and the initial point:

    >> maxiter = 100;
    >> timemax_seconds = 100;
    >> x0 = eye(n); x0(:,:,2) = eye(n);
    >> [S,T,distance,~,Q,~] = nearest_complex(A, B, maxiter, timemax_seconds, x0);


In order to specify the minimal index of the output T x + S, use the function `nearest_with_minimal_index`:

    >> min_index = 0;
    >> [S,T,distance,~,Q,~] = nearest_with_minimal_index(A, B, min_index);
    >> Q(:,:,1)*S*Q(:,:,2), Q(:,:,1)*T*Q(:,:,2)

    ans =
      -0.0000 - 0.0000i   2.6596 + 0.3850i  -1.0184 + 2.8725i  -1.8485 - 0.9521i
      -0.0000 - 0.0000i   1.2888 - 0.3109i  -1.9930 + 0.8365i   0.6803 - 1.0732i
       0.0000 - 0.0000i   0.0000 - 0.0000i  -3.7515 + 0.4585i   1.5996 - 1.7822i
       0.0000 - 0.0000i  -0.0000 - 0.0000i   0.0000 - 0.0000i   1.0967 + 2.8371i
    
    ans =
      -0.0000 - 0.0000i   2.8909 + 0.9939i  -0.6366 - 0.1109i   0.8995 - 1.1392i
       0.0000 + 0.0000i  -2.5817 + 2.9456i   0.3045 + 0.2884i  -0.2515 + 0.2478i
       0.0000 - 0.0000i   0.0000 + 0.0000i  -0.2716 + 0.0589i   0.0241 + 0.4902i
       0.0000 + 0.0000i  -0.0000 - 0.0000i   0.0000 - 0.0000i  -0.5334 + 1.6781i


The file `example6.1.mat` contains the computed minimum for Example 6.1. of the paper "A Riemannian optimization method to compute the nearest singular pencil". 
