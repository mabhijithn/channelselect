function results = gurobi_direct(A, b, k, Galconn_mat, noflags)
    noflags = noflags+1;
    M = 10;
    n = size(A,2)/noflags;
    y = b;
    Aty = (A'*y)./size(A,1);
    Q = (A'*A)./size(A,1);
    vtype = [repmat('C',[size(A,2) 1]);repmat('B',[n 1])];
    b = k-n;
    
    if(noflags == 1)
        Aineq = zeros(1+2*n,2*n);
        Aineq(1,:) = [zeros(1,n), (-1)*ones(1,n)];    

        for i = 1:n
            Aineq(1+i,i) = 1;
            Aineq(1+i,n+i) = 1;
        end

        for i = n+1:2*n
            Aineq(1+i,i-n) = -1;
            Aineq(1+i,i) = 1;
        end
        b = [b; (M+1)*ones(2*n,1)];
    else
%         Aineq = zeros(1+2*size(A,2), size(A,2)+n);
%         Aineq(1,:) = [zeros(1,size(A,2)), (-1)*ones(1,n)];
%         for i = 1:n
%             for j = 1:noflags
%                 strt = (i-1)*noflags+j;
%                 stp = i*noflags;
%                 Aineq(1+strt,strt) = 1;
%                 Aineq(1+strt,size(A,2)+i) = 1;
%             end
%         end
%         
%         for i = 1:n
%             for j = 1:noflags
%                 strt = size(A,2)+(i-1)*noflags+j;
%                 stp = i*noflags;
%                 Aineq(strt,strt-size(A,2)) = -1;
%                 Aineq(strt,size(A,2)+i) = 1;
%             end
%         end
%         b = [b; (M+1)*ones(size(Aineq,1)-1,1)];
        Aineq = [zeros(1,n*noflags),ones(1,n),zeros(1,n);zeros(n,n*noflags+2*n)];
        offset = n*noflags;
        for i = 1:n
             Aineq(1+i,offset+i) = 1;
             Aineq(1+i,offset+n+i) = 1;
        end
    end
    sense = ['='; repmat('<',[size(Aineq,1)-1,1])];

    m.A = sparse(Aineq);
    m.Q = sparse([Q/2,zeros(size(A,2),n); zeros(n,size(A,2)+n)]);
    m.obj = [-Aty; zeros(n,1)];
    m.rhs = b;
    m.sense = sense;
    m.vtype = vtype;
    sos = struct('type',[],'index',[],'weight',[]);
    kk = 1;
    for i = 1:n
        for j = 1:noflags
            strt = (i-1)*noflags+j;
            sos_tmp.type = 1;
            sos_tmp.index = [strt, size(A,2)+i];
            sos_tmp.weight = [1, 2];
            if(kk==1)
                sos = sos_tmp;
                kk = 0;
            else
                sos = [sos; sos_tmp];
            end
        end
    end
    
    % Galvanic connection constraint as SOS constraint
    for i = 1:size(Galconn_mat, 1)
        for j = 1:i-1
           if(Galconn_mat(i,j)==1)
                Galconn_mat(i,j) = 0;
           end               
        end      
    end
    [r,c] = find(Galconn_mat==1);
    for i = 1:size(r,1)
        sos_tmp.type = 1;
        sos_tmp.index = [r(i),c(i)];
        sos_tmp.weight = [1,2];
        sos = [sos; sos_tmp];
    end
    
    m.sos = sos;
    gurobi_write(m, 'qp.lp');
    results = gurobi(m);
end