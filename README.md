Accompanying code for "A Riemannian optimization method to compute the nearest singular pencil", F. Dopico, V. Noferini and L. Nyman.

This code has been tested with Manopt 7.1.

Example:

    >> importmanopt
    Manopt was added to Matlab's path.
    >> rng(0); n = 4;
    >> A = randn(n) + randn(n)*1i; 
    >> B = randn(n) + randn(n)*1i;
    >> [S,T,distance,time_seconds,Q,infotable] = nearest_complex(A, B);

The distance between B x + A and the singular T x + S is

    >> norm([A-S B-T],'f')
    ans =
      1.8432

which is also stored as the last element of the distance array:

    >> distance(end)
    ans =
      1.8432

We can bring the pencil T*x + S into Schur from via Q:
   
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
    >> [S,T,distance,time_seconds,Q,infotable] = nearest_complex(A, B, maxiter, timemax_seconds, x0);
