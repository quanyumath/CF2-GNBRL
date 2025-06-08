function Z = delta_fL(X,lambda,beta,A,L,S,psi)

   y1 = tnn_psi(L,psi); z1 = lambda*sum(y1(:));  %  y1 = abs(L).^p;
   y2 = tprod(A,L)+S-X; z2 = beta*norm(y2(:))^2/2;
   Z = z1+z2;

end