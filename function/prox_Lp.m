function  y  = prox_Lp(z, lambda, p)

    size_z = size(z); z = z(:);
    if p == 0
        y = z;
        i1 = find(abs(z)<=sqrt(2*lambda));
        y(i1) = 0;
    elseif 0<p && p<1
        max_iter = 10;
        ABSTOL = 1e-6;
        y = zeros(length(z),1);
        beta = (1/(lambda*2*(1-p)))^(1/(p-2));
        f1 = lambda*p*beta.^(p-1)+beta-abs(z);
        i0 = find(f1<0);
        if ~isempty(i0) 
            z_u = abs(z(i0));
            y_u = z_u;
            for i=1:max_iter              
                deta_y = (lambda*p*y_u.^(p-1)+y_u-z_u)./(lambda*p*(p-1)*y_u.^(p-2)+1);
                y_u = y_u-deta_y;
                if norm(deta_y) < sqrt(length(y_u))*ABSTOL
                    break;
                end
            end
            y_u = y_u.*sign(z(i0));
            i1 = find(1/2*z_u.^2-lambda*abs(y_u).^p-1/2*(y_u-z(i0)).^2<0);
            y_u(i1) = 0;
            y(i0) = y_u;
        end
    elseif p == 1
        y = max(0,z-lambda)+min(0,z+lambda);
    end
y = reshape(y,size_z);
end

