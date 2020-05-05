function dual = primal2dual(prob)
%PRIMAL2DUAL Summary of this function goes here
%   Detailed explanation goes here

if logical(sum(prob.blc==prob.buc))&&((isempty(prob.blx)))&&((isempty(prob.bux)))

    [m,n] = size(prob.a);
    nb_cones = size(prob.cones,1);
    cl = struct2cell(cell2mat(prob.cones));
    ic = cell2mat(cl(2,:))';
    nc = length(ic);
    sizec = arrayfun(@(x) length(x{1}),cl(2,:))';
    sizecs = [0;cumsum(sizec)];
    
    dual.cones = cell(nb_cones,1);
    for t=1:nb_cones
        dual.cones{t}.type = prob.cones{t}.type;
        dual.cones{t}.sub = m+sizecs(t)+(1:sizec(t));
    end

    dual.c = [prob.buc;sparse(nc,1)];
    dual.blc = prob.c;
    dual.buc = dual.blc;
    I = sparse(ic,(1:nc)',1,n,nc);
    dual.a = [prob.a' I];
    dual.blx = [];
    dual.bux = [];


end

if logical(sum(prob.blc~=prob.buc))
    ip = find(prob.blc == -inf);
    dual.bux = inf*ones(m+nc,1);
    dual.bux(ip) = 0;
end
    
if not(isempty(prob.blx))
    [m,n] = size(prob.a);
    nb_cones = size(prob.cones,1);
    nc = 0;
    for t=1:nb_cones
        nc = nc + length(prob.cones{t}.sub);
    end
    
    ip = find(prob.blx==0);
    nx = length(ip);
    Ix = speye(n);
    
    dual.c = [prob.buc;sparse(nx+nc,1)];
    dual.blc = prob.c;
    dual.buc = dual.blc;
    dual.a = [prob.a' Ix(:,ip) [sparse(n-nc,nc);speye(nc)]];
    dual.blx = [-inf*ones(m,1);sparse(nx,1);-inf*ones(nc,1)];
    dual.bux = [];
    
    dual.cones = cell(nb_cones,1);
    for t=1:nb_cones
        dual.cones{t}.type = prob.cones{t}.type;
        dual.cones{t}.sub = m+nx-(n-nc)+prob.cones{t}.sub;
    end
end
    
end

