function D=MakeLaplacian1D (N)
%make the diffusion matrix.
%mixed: homogeneous neumann at 0, homogeneous dirichlet on 1
%note: N must be large in this case!
N=N-1;
D=-2*eye(N,N)+diag(diag(eye(N-1,N-1)),1)+diag(diag(eye(N-1,N-1)),-1);
D(1,1)=-1;
D(1,2)=1;
D(N,N)=-2;
%D(N,N-1)=0;
%D(N,N-2)=0;
D(N,N-1)=1;

return