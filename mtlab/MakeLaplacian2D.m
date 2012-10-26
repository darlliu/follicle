function nabla=MakeLaplacian2D(N,M)
%make 2d laplacian with periodic on y (M nodes)
%and homogeneous neumann on z (N nodes)

tdx=-2*eye(N*M,N*M)+diag(diag(eye(N*(M-1),(M-1)*N)),N)+diag(diag(eye(N*(M-1),(M-1)*N)),-N) ...
+ diag(diag(eye(N,N)),(M-1)*N)+ diag(diag(eye(N,N)),-(M-1)*N);

D=-2*eye(N,N)+diag(diag(eye(N-1,N-1)),1)+diag(diag(eye(N-1,N-1)),-1);
D(1,1)=-1;
D(1,2)=1;
D(N-1,N-1)=-1;
D(N,N)=0;
D(N,N-1)=-1;
D(N,N-2)=1;
tdy=D;
for i= 2:M,
    tdy=blkdiag(tdy,D);
end
size(tdx)
size(tdy)
nabla=sparse(tdx+tdy);

return
