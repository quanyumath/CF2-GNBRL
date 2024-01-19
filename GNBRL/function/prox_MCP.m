function  y  = prox_MCP(z, lambda, alpha)

    size_z = size(z); z = z(:); y = z;
    i1 = find(alpha*lambda<abs(z) & abs(z)<=alpha);
    y(i1) = (sign(z(i1)).*(abs(z(i1))-lambda*alpha))./(1-lambda);
    i2 = find(abs(z)<=alpha*lambda);
    y(i2) = 0;
 
    y = reshape(y,size_z);
end

